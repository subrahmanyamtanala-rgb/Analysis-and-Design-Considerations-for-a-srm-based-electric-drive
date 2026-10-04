"""End-to-end design flow: vehicle -> requirements -> machine -> converter."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .drive import SimResult, SRMDrive
from .geometry import SizingInputs, SRMDesign, size_srm
from .performance import max_torque, operating_point, required_peak_current
from .vehicle import DriveRequirements, EBikeSpec, derive_requirements


@dataclass
class DriveDesign:
    spec: EBikeSpec
    req: DriveRequirements
    inputs: SizingInputs
    design: SRMDesign
    drive: SRMDrive
    i_peak_required: float


def design_drive(
    spec: EBikeSpec | None = None,
    current_margin: float = 1.10,
    fea_table: str | None = None,
    **sizing_overrides,
) -> DriveDesign:
    """Size the machine and set the converter current limit.

    1. vehicle duty -> torque/speed requirements (corner speed for turns)
    2. analytical sizing (output equation, pole arcs, turns, slot fill)
    3. simulation-based current requirement: the smallest chopping current
       that produces the launch torque with optimal angles; the converter limit
       is set ``current_margin`` above it.

    With ``fea_table`` (an .npz written by ``scripts/run_fea.py``) the drive
    uses the FEA flux-linkage characteristic instead of the analytical one.
    """
    spec = spec or EBikeSpec()
    req = derive_requirements(spec)
    inputs = SizingInputs(
        peak_torque=req.peak_torque,
        base_speed_rpm=req.corner_speed_rpm,
        dc_voltage=spec.battery_voltage_min,
        **sizing_overrides,
    )
    design = size_srm(inputs)
    model = None
    if fea_table is not None:
        from .magnetics_fea import load_fea_model

        model = load_fea_model(design, fea_table)
    drive = SRMDrive(design, v_dc=spec.battery_voltage, i_max=100.0, magnetic_model=model)
    i_req = required_peak_current(drive, 0.1 * req.corner_speed_rpm, req.peak_torque, i_hi=100.0)
    drive.i_max = current_margin * i_req
    return DriveDesign(spec, req, inputs, design, drive, i_req)


# ---------------------------------------------------------------------------
# Converter and thermal design checks
# ---------------------------------------------------------------------------
def converter_ratings(dd: DriveDesign, peak: SimResult, rated: SimResult) -> dict:
    """Device and DC-link ratings for the asymmetric half bridge."""
    spec = dd.spec
    v_dev = 1.5 * spec.battery_voltage_max  # regen / ringing headroom
    i_dc_peak = peak.i_dc_total
    i_dc_rated = rated.i_dc_total
    # the battery supplies the mean; the capacitor carries the AC component
    i_cap_rms_peak = float(np.sqrt(np.mean((i_dc_peak - i_dc_peak.mean()) ** 2)))
    i_cap_rms_rated = float(np.sqrt(np.mean((i_dc_rated - i_dc_rated.mean()) ** 2)))
    # capacitor sized for <= 5 % voltage ripple at rated: C >= dQ / dV, with
    # dQ approximated by the charge exchanged in one stroke above the mean
    period = dd.drive.model.tau_r / rated.omega
    dt = period / len(i_dc_rated)
    ac = i_dc_rated - i_dc_rated.mean()
    q = np.cumsum(ac) * dt
    dq = float(q.max() - q.min())
    c_min = dq / (0.05 * spec.battery_voltage)
    return {
        "topology": f"asymmetric half bridge, {dd.design.phases} legs, "
        f"{2 * dd.design.phases} MOSFETs + {2 * dd.design.phases} diodes",
        "device_voltage_rating_V": round(v_dev, 0),
        "selected_device_voltage_V": 100,
        "phase_current_limit_A": round(dd.drive.i_max, 1),
        "device_rms_current_peak_A": round(peak.i_rms, 2),
        "battery_current_peak_A": round(float(i_dc_peak.mean()), 2),
        "battery_current_rated_A": round(float(i_dc_rated.mean()), 2),
        "dc_link_ripple_current_rms_peak_A": round(i_cap_rms_peak, 2),
        "dc_link_ripple_current_rms_rated_A": round(i_cap_rms_rated, 2),
        "dc_link_capacitance_min_uF": round(c_min * 1e6, 0),
        "switching_frequency_rated_kHz": round(rated.switching_frequency / 1e3, 2),
    }


def thermal_estimate(dd: DriveDesign, rated: SimResult, h_conv: float = 25.0, t_amb: float = 30.0) -> dict:
    """Lumped steady-state temperature rise at continuous rating.

    Heat leaves through the hub shell; ``h_conv`` lumps natural + forced
    convection from riding (~15-35 W/m^2K at 10-25 km/h).
    """
    d = dd.design
    ls = rated.losses()
    p_motor = ls["copper"] + ls["core"] + ls["mechanical"]
    r_shell = d.outer_diameter / 2 + 0.01
    l_shell = d.stack_length + 0.04  # end windings + covers
    area = 2 * math.pi * r_shell * l_shell + 2 * math.pi * r_shell**2
    dtheta = p_motor / (h_conv * area)
    j_rms = rated.i_rms / d.conductor_area / 1e6
    return {
        "motor_loss_W": round(p_motor, 1),
        "cooling_area_m2": round(area, 4),
        "temperature_rise_K": round(dtheta, 1),
        "winding_temp_estimate_C": round(t_amb + 1.2 * dtheta, 1),  # 20 % hot-spot
        "current_density_rated_A_mm2": round(j_rms, 2),
    }


def rated_and_peak_points(dd: DriveDesign) -> tuple[SimResult, SimResult]:
    req = dd.req
    rated = operating_point(dd.drive, req.base_speed_rpm, req.rated_torque)
    peak = max_torque(dd.drive, 0.1 * req.corner_speed_rpm)
    return rated, peak
