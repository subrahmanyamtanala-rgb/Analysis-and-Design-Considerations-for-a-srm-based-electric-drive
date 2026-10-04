"""E-bike road-load model and drive requirements.

All quantities are SI unless the name says otherwise (``_kmh``, ``_rpm``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

G = 9.81  # m/s^2


@dataclass
class EBikeSpec:
    """Vehicle, rider and powertrain parameters for a pedelec-class e-bike."""

    mass_total: float = 100.0  # rider (80 kg) + bike incl. battery (20 kg)
    wheel_radius: float = 0.33  # 26" wheel with tyre
    c_rr: float = 0.006  # rolling resistance coefficient (city tyre on asphalt)
    cd_a: float = 0.50  # drag area, upright rider [m^2]
    rho_air: float = 1.2  # kg/m^3
    rotating_mass_factor: float = 1.05  # equivalent inertia of wheels/rotor
    gear_ratio: float = 8.0  # motor speed / wheel speed (single-stage planetary, ring/sun = 7)
    gear_efficiency: float = 0.95  # mesh (load-dependent) efficiency
    gear_no_load_loss: float = 0.0  # churning/seal loss at base speed [W], ~ speed
    v_max_kmh: float = 25.0  # assist cut-off (EN 15194 / India CMVR e-cycle)
    rated_power: float = 250.0  # continuous rated motor output power [W]
    peak_power: float = 500.0  # short-term (hill / launch) motor output power [W]
    battery_voltage: float = 36.0  # 10S Li-ion, nominal
    battery_voltage_min: float = 30.0  # 10S at cut-off
    battery_voltage_max: float = 42.0  # 10S fully charged
    battery_capacity_ah: float = 10.0

    # ----- derived -------------------------------------------------------
    @property
    def v_max(self) -> float:
        return self.v_max_kmh / 3.6

    @property
    def battery_energy_wh(self) -> float:
        return self.battery_voltage * self.battery_capacity_ah

    def wheel_speed(self, v: float) -> float:
        """Wheel angular speed [rad/s] for road speed ``v`` [m/s]."""
        return v / self.wheel_radius

    def motor_speed(self, v: float) -> float:
        """Motor angular speed [rad/s] for road speed ``v`` [m/s]."""
        return self.gear_ratio * self.wheel_speed(v)

    def motor_speed_rpm(self, v: float) -> float:
        return self.motor_speed(v) * 60.0 / (2.0 * math.pi)


def road_load_force(spec: EBikeSpec, v: float, grade: float = 0.0, accel: float = 0.0) -> float:
    """Tractive force at the tyre contact patch [N].

    ``grade`` is the road slope as a fraction (0.05 = 5 %), ``accel`` in m/s^2.
    """
    alpha = math.atan(grade)
    f_roll = spec.c_rr * spec.mass_total * G * math.cos(alpha)
    f_grade = spec.mass_total * G * math.sin(alpha)
    f_aero = 0.5 * spec.rho_air * spec.cd_a * v * abs(v)
    f_accel = spec.rotating_mass_factor * spec.mass_total * accel
    return f_roll + f_grade + f_aero + f_accel


def gear_drag_torque(spec: EBikeSpec) -> float:
    """Motor-side drag torque of the gear's load-independent loss [N m]."""
    if spec.gear_no_load_loss <= 0:
        return 0.0
    return spec.gear_no_load_loss / spec.motor_speed(spec.v_max)


def motor_torque_demand(spec: EBikeSpec, force: float) -> float:
    """Motor shaft torque [N m] required to deliver tractive ``force`` (motoring).

    Gear model: load-dependent mesh efficiency plus a load-independent drag
    (no-load loss proportional to speed, i.e. a constant drag torque).
    """
    t_wheel = force * spec.wheel_radius
    if t_wheel >= 0:
        return t_wheel / (spec.gear_ratio * spec.gear_efficiency) + gear_drag_torque(spec)
    return t_wheel * spec.gear_efficiency / spec.gear_ratio


@dataclass
class DriveRequirements:
    """Motor-side requirements derived from the vehicle mission."""

    base_speed_rpm: float  # motor speed at the 25 km/h assist cut-off
    max_speed_rpm: float
    corner_speed_rpm: float  # end of the constant-torque region (T_peak * w = P_peak)
    rated_torque: float
    peak_torque: float
    cruise_torque: float
    cruise_power: float
    hill_torque: float
    start_torque: float


def derive_requirements(
    spec: EBikeSpec,
    hill_grade: float = 0.06,
    hill_speed_kmh: float = 15.0,
    start_grade: float = 0.08,
    start_accel: float = 0.5,
    assist_share: float = 1.0,
) -> DriveRequirements:
    """Translate vehicle-level duty points into motor torque/speed requirements.

    ``assist_share`` is the fraction of the tractive effort supplied by the motor
    (1.0 = throttle-only / worst case; a pedelec typically needs 0.5-0.7).
    """
    v_max = spec.v_max
    n_max = spec.motor_speed_rpm(v_max)

    f_cruise = road_load_force(spec, v_max)
    t_cruise = assist_share * motor_torque_demand(spec, f_cruise)
    p_cruise = t_cruise * spec.motor_speed(v_max)

    v_hill = hill_speed_kmh / 3.6
    t_hill = assist_share * motor_torque_demand(spec, road_load_force(spec, v_hill, hill_grade))

    t_start = assist_share * motor_torque_demand(
        spec, road_load_force(spec, 0.5, start_grade, start_accel)
    )

    rated_torque = spec.rated_power / spec.motor_speed(v_max)
    peak_torque = max(t_hill, t_start, 2.0 * rated_torque)
    # peak torque is only needed at low road speed; above the corner speed the
    # drive runs at constant (peak) power. Sizing the turns for this corner
    # instead of for 25 km/h halves the current needed for the launch torque.
    corner_rpm = min(n_max, spec.peak_power / peak_torque * 60 / (2 * math.pi))
    return DriveRequirements(
        base_speed_rpm=n_max,
        max_speed_rpm=1.2 * n_max,  # allow overspeed downhill with the motor idling
        corner_speed_rpm=corner_rpm,
        rated_torque=rated_torque,
        peak_torque=peak_torque,
        cruise_torque=t_cruise,
        cruise_power=p_cruise,
        hill_torque=t_hill,
        start_torque=t_start,
    )
