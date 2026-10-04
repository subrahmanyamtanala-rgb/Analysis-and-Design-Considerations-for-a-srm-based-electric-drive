"""Loss models: core (Steinmetz), converter and mechanical losses."""

from __future__ import annotations

from dataclasses import dataclass

from .geometry import DENSITY_FE, SRMDesign


@dataclass
class IronLossCoefficients:
    """Bertotti-type coefficients for non-oriented steel (M270-35A class).

    p [W/kg] = k_h f B^2 + k_e (f B)^2, calibrated to ~2.7 W/kg at 1.5 T / 50 Hz.
    ``k_nonsine`` accounts for the unipolar, non-sinusoidal flux pulses of an
    SRM (higher dB/dt in the poles, minor loops in the yoke).
    """

    k_h: float = 0.018
    k_e: float = 1.3e-4
    k_nonsine: float = 1.3

    def specific_loss(self, f: float, b: float) -> float:
        return self.k_nonsine * (self.k_h * f * b * b + self.k_e * (f * b) ** 2)


@dataclass
class ConverterParams:
    """Asymmetric half-bridge built from 100 V trench MOSFETs."""

    r_ds_on: float = 4e-3  # ohm, at 100 C
    v_diode: float = 0.8  # body/fast diode forward drop
    t_switch: float = 60e-9  # combined rise + fall time per transition
    q_rr: float = 40e-9  # diode reverse-recovery charge [C]
    v_rating: float = 100.0
    i_rating: float = 60.0


def core_loss(design: SRMDesign, coeff: IronLossCoefficients, n_rpm: float, b_pole_peak: float) -> dict:
    """Split core loss into stator poles, stator yoke and rotor [W]."""
    f_ph = design.nr * n_rpm / 60.0
    f_rotor = design.ns * n_rpm / 120.0
    length = design.stack_length
    m_sp = DENSITY_FE * length * design.ns * design.stator_pole_width * design.stator_pole_height
    m_sy = design.mass_iron_stator - m_sp
    m_r = design.mass_iron_rotor
    b_sy = b_pole_peak * design.stator_pole_width / (2 * design.stator_yoke)
    b_r = b_pole_peak * design.stator_pole_width / design.rotor_pole_width
    p_sp = m_sp * coeff.specific_loss(f_ph, b_pole_peak)
    p_sy = m_sy * coeff.specific_loss(f_ph, b_sy)
    p_r = m_r * coeff.specific_loss(f_rotor, 0.7 * b_r)  # rotor yoke runs lighter
    return {"stator_poles": p_sp, "stator_yoke": p_sy, "rotor": p_r, "total": p_sp + p_sy + p_r}


def mechanical_loss(n_rpm: float) -> float:
    """Bearings + windage of a small sealed hub motor [W]."""
    k = n_rpm / 1000.0
    return 0.8 * k + 0.4 * k * k
