"""Time-domain simulation of the SRM + asymmetric half-bridge converter.

Steady state at constant speed is computed for one phase over one electrical
period (one rotor pole pitch). Because the phases are identical and displaced
by one stroke angle, the total torque is the sum of the single-phase torque
shifted by multiples of the stroke. If the current does not extinguish within
the period (continuous conduction) further periods are simulated until the
flux at the period boundary repeats.

Converter states (per phase, asymmetric half bridge):

* both switches on  -> v = +Vdc - 2 R_ds i            (magnetise)
* one switch on     -> v = -(R_ds i + V_d)            (soft-chopping freewheel)
* both switches off -> v = -(Vdc + 2 V_d)             (demagnetise via diodes)

Current regulation is a sampled hysteresis controller (soft or hard chopping);
at high speed the back-EMF prevents the current from reaching the reference and
operation becomes single-pulse automatically.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .geometry import SRMDesign
from .losses import ConverterParams, IronLossCoefficients, core_loss, mechanical_loss
from .magnetics import MagneticModel

B_SAT_MAX = 2.0  # T, practical saturation limit of non-oriented Si steel


@dataclass
class SRMDrive:
    design: SRMDesign
    v_dc: float
    i_max: float  # converter / thermal current limit for chopping reference
    converter: ConverterParams = field(default_factory=ConverterParams)
    iron: IronLossCoefficients = field(default_factory=IronLossCoefficients)
    f_sample: float = 40e3  # controller sampling frequency
    hysteresis_band: float = 0.2  # A, minimum band
    hysteresis_band_rel: float = 0.02  # band as a fraction of the reference
    soft_chopping: bool = True
    min_steps: int = 1800

    def __post_init__(self) -> None:
        self.model = MagneticModel.from_design(self.design)

    # ------------------------------------------------------------------
    @property
    def phases(self) -> int:
        return self.design.phases

    def turn_on_advance(self, n_rpm: float, i_ref: float) -> float:
        """Angle needed to build up ``i_ref`` in the unaligned inductance."""
        omega = n_rpm * 2 * math.pi / 60
        v = self.v_dc - 2 * self.converter.r_ds_on * i_ref - self.design.phase_resistance * i_ref
        return omega * self.model.l_u * i_ref / max(v, 1e-3)

    def default_angles(self, n_rpm: float, i_ref: float, off_fraction: float = 0.6) -> tuple[float, float]:
        """Turn-on at overlap minus build-up advance; turn-off a fraction into the
        overlap-to-aligned interval (earlier at speed, see ``optimise``)."""
        mm = self.model
        theta_ov = mm.overlap_start - mm.fringe / 2
        theta_on = theta_ov - self.turn_on_advance(n_rpm, i_ref)
        theta_off = theta_ov + off_fraction * (mm.tau_r / 2 - theta_ov)
        return theta_on, theta_off

    # ------------------------------------------------------------------
    def simulate(
        self,
        n_rpm: float,
        theta_on: float,
        theta_off: float,
        i_ref: float,
        max_periods: int = 6,
    ) -> "SimResult":
        m = self.phases
        mm = self.model
        cv = self.converter
        r_ph = self.design.phase_resistance
        omega = n_rpm * 2 * math.pi / 60
        tau = mm.tau_r
        period = tau / omega
        n = max(self.min_steps, int(math.ceil(period * self.f_sample)))
        n = int(math.ceil(n / m) * m)
        dt = period / n
        dth = tau / n
        dwell = (theta_off - theta_on) % tau
        h = 0.5 * max(self.hysteresis_band, self.hysteresis_band_rel * i_ref)
        vdc = self.v_dc

        psi0 = 0.0
        for _ in range(max_periods):
            th = np.empty(n)
            cur = np.empty(n)
            flux = np.empty(n)
            torque = np.empty(n)
            p_conv = np.empty(n)
            idc = np.empty(n)
            psi = psi0
            i = mm.current(psi, theta_on)
            sw_on = True
            transitions = 0
            for k in range(n):
                theta = theta_on + k * dth
                rel = k * dth
                if rel < dwell:
                    if i >= i_ref + h:
                        if sw_on:
                            transitions += 1
                        sw_on = False
                    elif i <= i_ref - h:
                        if not sw_on:
                            transitions += 1
                        sw_on = True
                    if sw_on:
                        v = vdc - 2 * cv.r_ds_on * i
                        pc = 2 * cv.r_ds_on * i * i
                        idc_k = i
                    elif self.soft_chopping:
                        v = -(cv.r_ds_on * i + cv.v_diode) if i > 0 else 0.0
                        pc = cv.r_ds_on * i * i + cv.v_diode * i
                        idc_k = 0.0
                    else:
                        v = -(vdc + 2 * cv.v_diode) if i > 0 else 0.0
                        pc = 2 * cv.v_diode * i
                        idc_k = -i
                else:
                    if sw_on and k > 0:
                        transitions += 1
                        sw_on = False
                    if i > 0:
                        v = -(vdc + 2 * cv.v_diode)
                        pc = 2 * cv.v_diode * i
                        idc_k = -i
                    else:
                        v, pc, idc_k = 0.0, 0.0, 0.0
                th[k] = theta
                cur[k] = i
                flux[k] = psi
                torque[k] = mm.torque(i, theta)
                p_conv[k] = pc
                idc[k] = idc_k
                psi = psi + dt * (v - r_ph * i)
                if psi <= 0.0:
                    psi = 0.0
                i = mm.current(psi, theta + dth, i)
            converged = abs(psi - psi0) < 1e-4 * max(flux.max(), 1e-9)
            psi0 = psi
            if converged:
                break
        return SimResult(
            drive=self,
            n_rpm=n_rpm,
            theta_on=theta_on,
            theta_off=theta_off,
            i_ref=i_ref,
            theta=th,
            current=cur,
            flux=flux,
            torque_phase=torque,
            p_conv_inst=p_conv,
            i_dc_phase=idc,
            transitions=transitions,
            continuous=psi0 > 0.0,
        )


@dataclass
class SimResult:
    drive: SRMDrive
    n_rpm: float
    theta_on: float
    theta_off: float
    i_ref: float
    theta: np.ndarray
    current: np.ndarray
    flux: np.ndarray
    torque_phase: np.ndarray
    p_conv_inst: np.ndarray
    i_dc_phase: np.ndarray
    transitions: int
    continuous: bool

    # ---- waveforms --------------------------------------------------------
    @property
    def torque_total(self) -> np.ndarray:
        m = self.drive.phases
        n = len(self.theta)
        shift = n // m
        return sum(np.roll(self.torque_phase, k * shift) for k in range(m))

    @property
    def i_dc_total(self) -> np.ndarray:
        m = self.drive.phases
        shift = len(self.theta) // m
        return sum(np.roll(self.i_dc_phase, k * shift) for k in range(m))

    # ---- scalar results ------------------------------------------------------
    @property
    def omega(self) -> float:
        return self.n_rpm * 2 * math.pi / 60

    @property
    def torque_avg(self) -> float:
        return float(self.drive.phases * self.torque_phase.mean())

    @property
    def torque_ripple(self) -> float:
        t = self.torque_total
        avg = t.mean()
        return float((t.max() - t.min()) / avg) if avg > 1e-9 else float("nan")

    @property
    def i_rms(self) -> float:
        return float(np.sqrt(np.mean(self.current**2)))

    @property
    def i_peak(self) -> float:
        return float(self.current.max())

    @property
    def psi_peak(self) -> float:
        return float(self.flux.max())

    @property
    def b_pole_peak(self) -> float:
        """Peak stator-pole flux density. psi/(T A) also counts leakage flux
        once the pole saturates, so it is capped at the saturation level."""
        d = self.drive.design
        return min(self.psi_peak / (d.turns_per_phase * d.pole_face_area), B_SAT_MAX)

    @property
    def switching_frequency(self) -> float:
        """Average device switching frequency per phase [Hz]."""
        period = self.drive.model.tau_r / self.omega
        return self.transitions / period

    def losses(self) -> dict:
        d = self.drive.design
        cv = self.drive.converter
        m = d.phases
        p_cu = m * d.phase_resistance * self.i_rms**2
        p_cond = m * float(self.p_conv_inst.mean())
        # each transition dissipates ~0.5 V I t_sw in the switching device
        i_sw = self.drive.i_max if self.transitions else 0.0
        i_sw = min(i_sw, self.i_peak)
        p_sw = m * self.switching_frequency * 0.5 * self.drive.v_dc * i_sw * cv.t_switch
        p_fe = core_loss(d, self.drive.iron, self.n_rpm, self.b_pole_peak)["total"]
        p_mech = mechanical_loss(self.n_rpm)
        return {
            "copper": p_cu,
            "core": p_fe,
            "converter_conduction": p_cond,
            "converter_switching": p_sw,
            "mechanical": p_mech,
        }

    @property
    def p_electromagnetic(self) -> float:
        return self.torque_avg * self.omega

    @property
    def p_shaft(self) -> float:
        return self.p_electromagnetic - self.losses()["core"] - self.losses()["mechanical"]

    @property
    def p_dc_circuit(self) -> float:
        """DC input from the circuit model (no core/switching losses)."""
        return float(self.drive.v_dc * self.i_dc_total.mean())

    def efficiency(self) -> dict:
        ls = self.losses()
        p_out = self.p_shaft
        p_in = p_out + sum(ls.values())
        motor_loss = ls["copper"] + ls["core"] + ls["mechanical"]
        conv_loss = ls["converter_conduction"] + ls["converter_switching"]
        return {
            "p_out": p_out,
            "p_in": p_in,
            "eta_system": p_out / p_in if p_in > 0 else 0.0,
            "eta_motor": p_out / (p_out + motor_loss) if p_out > 0 else 0.0,
            "eta_converter": (p_in - conv_loss) / p_in if p_in > 0 else 0.0,
        }

    def summary(self) -> dict:
        eff = self.efficiency()
        return {
            "speed_rpm": round(self.n_rpm, 1),
            "theta_on_deg": round(math.degrees(self.theta_on), 2),
            "theta_off_deg": round(math.degrees(self.theta_off), 2),
            "i_ref_A": round(self.i_ref, 2),
            "torque_avg_Nm": round(self.torque_avg, 3),
            "torque_ripple_pct": round(100 * self.torque_ripple, 1),
            "i_peak_A": round(self.i_peak, 2),
            "i_rms_A": round(self.i_rms, 2),
            "psi_peak_mWb": round(self.psi_peak * 1e3, 2),
            "b_pole_peak_T": round(self.b_pole_peak, 3),
            "f_switch_kHz": round(self.switching_frequency / 1e3, 2),
            "p_out_W": round(eff["p_out"], 1),
            "eta_motor_pct": round(100 * eff["eta_motor"], 1),
            "eta_system_pct": round(100 * eff["eta_system"], 1),
            "continuous_conduction": self.continuous,
        }
