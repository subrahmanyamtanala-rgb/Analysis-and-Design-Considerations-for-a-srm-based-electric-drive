"""First-order analytic surface-PM (SPM) benchmark sized to the same mission.

Purpose: a like-for-like comparison with the SRM at the *same* fidelity as the
SRM's analytical model (same envelope D_o x L, gear, battery voltage, converter
devices, iron-loss coefficients, mechanical loss and efficiency boundary).
It is intentionally simple and is NOT an FE-verified PM design:

* 12-slot / 10-pole fractional-slot concentrated winding (k_w = 0.933),
  inner rotor with surface NdFeB magnets (N42SH, B_r = 1.17 T at 100 C);
* air-gap field from the magnetic-circuit equation with Carter factor,
  fundamental from the magnet arc; no saturation, no slotting harmonics;
* turns chosen with the same voltage logic as the SRM: no-load line-to-line
  back-EMF at the 25-km/h base speed equals the minimum battery voltage (30 V),
  so no field weakening is needed inside the assist range;
* sinusoidal FOC (I_d = 0), 6-switch inverter with the same 100-V MOSFETs,
  synchronous conduction and 16-kHz switching;
* stator iron loss with the SRM's Bertotti coefficients (k_ns = 1.0 for the
  near-sinusoidal flux) plus an assumed magnet eddy-current loss of 10 % of the
  stator iron loss (FSCW space harmonics).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geometry import ALPHA_CU, DENSITY_CU, DENSITY_FE, MU0, RHO_CU_20
from .losses import ConverterParams, IronLossCoefficients, mechanical_loss

DENSITY_NDFEB = 7500.0


@dataclass
class SPMInputs:
    d_outer: float
    length: float
    v_dc_min: float
    base_speed_rpm: float
    slots: int = 12
    pole_pairs: int = 5
    k_w: float = 0.933
    split_ratio: float = 0.60  # stator bore / outer diameter
    air_gap: float = 0.5e-3
    h_magnet: float = 3.0e-3
    magnet_arc: float = 0.83
    b_r: float = 1.17
    mu_r: float = 1.05
    carter: float = 1.05
    b_tooth: float = 1.6
    b_yoke: float = 1.4
    fill_factor: float = 0.45
    winding_temp: float = 100.0
    f_pwm: float = 16e3
    magnet_loss_fraction: float = 0.10


class SPMBenchmark:
    def __init__(self, inp: SPMInputs, converter: ConverterParams | None = None, iron: IronLossCoefficients | None = None):
        self.inp = p = inp
        self.cv = converter or ConverterParams()
        self.iron = iron or IronLossCoefficients(k_nonsine=1.0)
        d_s = p.split_ratio * p.d_outer
        self.d_bore = d_s
        g_eff = p.carter * p.air_gap
        hm = p.h_magnet / p.mu_r
        self.b_gap = p.b_r * hm / (hm + g_eff)
        self.b1 = 4 / math.pi * self.b_gap * math.sin(p.magnet_arc * math.pi / 2)
        tau_p = math.pi * d_s / (2 * p.pole_pairs)
        self.phi1 = 2 / math.pi * self.b1 * tau_p * p.length
        # turns: no-load line-to-line peak EMF at base speed = V_dc,min
        w_e_base = p.pole_pairs * p.base_speed_rpm * 2 * math.pi / 60
        psi_target = p.v_dc_min / (math.sqrt(3) * w_e_base)
        coils_per_phase = p.slots // 3
        n_c = max(1, round(psi_target / (p.k_w * self.phi1) / coils_per_phase))
        self.turns = n_c * coils_per_phase
        self.turns_per_coil = n_c
        self.psi_m = self.turns * p.k_w * self.phi1
        self.k_t = 1.5 * p.pole_pairs * self.psi_m  # N m per A (peak)
        # stator geometry
        tau_s = math.pi * d_s / p.slots
        self.w_tooth = self.b_gap * tau_s / p.b_tooth
        phi_pole = self.b_gap * p.magnet_arc * tau_p * p.length
        self.yoke = phi_pole / 2 / (p.b_yoke * p.length)
        r_o, r_s = p.d_outer / 2, d_s / 2
        self.slot_depth = r_o - self.yoke - r_s
        a_slots = math.pi * ((r_o - self.yoke) ** 2 - r_s**2) - p.slots * self.w_tooth * self.slot_depth
        self.slot_area = a_slots / p.slots
        self.a_cond = p.fill_factor * self.slot_area / (2 * n_c)
        t_c = self.slot_area / (2 * self.slot_depth)
        self.mlt = 2 * (p.length + self.w_tooth) + math.pi * t_c
        rho = RHO_CU_20 * (1 + ALPHA_CU * (p.winding_temp - 20))
        self.r_ph = rho * self.turns * self.mlt / self.a_cond
        # synchronous inductance: air-gap (effective gap incl. magnet) + 50 % slot/tip leakage
        g_tot = g_eff + hm
        l_ag = 3 / math.pi * MU0 * (self.turns * p.k_w / p.pole_pairs) ** 2 * p.length * r_s / g_tot
        self.l_s = 1.5 * l_ag
        # masses
        self.m_teeth = DENSITY_FE * p.length * p.slots * self.w_tooth * self.slot_depth
        self.m_yoke = DENSITY_FE * p.length * math.pi * (r_o**2 - (r_o - self.yoke) ** 2)
        r_rotor = r_s - p.air_gap
        r_core = r_rotor - p.h_magnet
        self.m_magnet = DENSITY_NDFEB * p.length * p.magnet_arc * math.pi * (r_rotor**2 - r_core**2)
        rotor_yoke = self.yoke * 1.1
        self.m_rotor_fe = DENSITY_FE * p.length * math.pi * (r_core**2 - max(r_core - rotor_yoke, 7.5e-3) ** 2)
        self.m_cu = DENSITY_CU * 3 * self.turns * self.mlt * self.a_cond

    @property
    def mass_active(self) -> float:
        return self.m_teeth + self.m_yoke + self.m_magnet + self.m_rotor_fe + self.m_cu

    def operating_point(self, n_rpm: float, torque: float, v_dc: float) -> dict | None:
        p, cv = self.inp, self.cv
        w_m = n_rpm * 2 * math.pi / 60
        w_e = p.pole_pairs * w_m
        i_pk = torque / self.k_t
        i_rms = i_pk / math.sqrt(2)
        # voltage check (I_d = 0): |V| = sqrt((E + R I)^2 + (w L I)^2) <= V_dc / sqrt(3)
        e = w_e * self.psi_m
        v_ph = math.hypot(e + self.r_ph * i_pk, w_e * self.l_s * i_pk)
        if v_ph > v_dc / math.sqrt(3) * 0.98:
            return None
        f_e = w_e / (2 * math.pi)
        p_cu = 3 * self.r_ph * i_rms**2
        p_fe_s = self.m_teeth * self.iron.specific_loss(f_e, p.b_tooth) + self.m_yoke * self.iron.specific_loss(f_e, p.b_yoke)
        p_mag = p.magnet_loss_fraction * p_fe_s
        p_mech = mechanical_loss(n_rpm)
        p_out = torque * w_m - p_fe_s - p_mag - p_mech
        p_cond = 3 * cv.r_ds_on * i_rms**2  # synchronous conduction, one device per leg
        i_avg_abs = 2 * math.sqrt(2) / math.pi * i_rms
        p_sw = 3 * p.f_pwm * (v_dc * i_avg_abs * cv.t_switch + cv.q_rr * v_dc)
        p_motor_loss = p_cu + p_fe_s + p_mag + p_mech
        p_in = p_out + p_motor_loss + p_cond + p_sw
        return {
            "T": torque,
            "i_rms": i_rms,
            "eta_motor": p_out / (p_out + p_motor_loss) if p_out > 0 else 0.0,
            "eta_system": p_out / p_in if p_out > 0 else 0.0,
            "losses": {"copper": p_cu, "core": p_fe_s, "magnet": p_mag, "mechanical": p_mech, "converter_conduction": p_cond, "converter_switching": p_sw},
        }

    def drag_open_circuit(self, n_rpm: float) -> float:
        """No-load loss when coasting with the inverter off (no freewheel clutch) [W]."""
        p = self.inp
        f_e = p.pole_pairs * n_rpm / 60
        p_fe = self.m_teeth * self.iron.specific_loss(f_e, p.b_tooth) + self.m_yoke * self.iron.specific_loss(f_e, p.b_yoke)
        return p_fe * (1 + p.magnet_loss_fraction)

    def summary(self) -> dict:
        return {
            "D_o_mm": self.inp.d_outer * 1e3,
            "L_mm": self.inp.length * 1e3,
            "bore_mm": self.d_bore * 1e3,
            "B_gap_T": self.b_gap,
            "B1_T": self.b1,
            "turns_per_phase": self.turns,
            "turns_per_coil": self.turns_per_coil,
            "k_t_Nm_per_Apk": self.k_t,
            "R_ph_mohm": self.r_ph * 1e3,
            "L_s_uH": self.l_s * 1e6,
            "mass_magnet_kg": self.m_magnet,
            "mass_copper_kg": self.m_cu,
            "mass_iron_kg": self.m_teeth + self.m_yoke + self.m_rotor_fe,
            "mass_active_kg": self.mass_active,
        }


def spm_lookup(spm: SPMBenchmark, speeds, torques, v_dc: float, i_max_pk: float):
    """EfficiencyLookup-compatible map for the route simulation."""
    from .drive_cycle import EfficiencyLookup

    t_env = []
    for n in speeds:
        # torque limit: current limit or voltage limit
        lo, hi = 0.0, spm.k_t * i_max_pk
        if spm.operating_point(n, hi, v_dc) is None:
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if spm.operating_point(n, mid, v_dc) else (lo, mid)
        t_env.append(hi if spm.operating_point(n, hi, v_dc) else lo)
    eta = np.full((len(speeds), len(torques)), np.nan)
    for a, n in enumerate(speeds):
        for b, t in enumerate(torques):
            if t <= 0.97 * t_env[a]:
                op = spm.operating_point(n, t, v_dc)
                if op:
                    eta[a, b] = op["eta_system"]
    return EfficiencyLookup(np.asarray(speeds, float), np.asarray(t_env), np.asarray(torques, float), eta)
