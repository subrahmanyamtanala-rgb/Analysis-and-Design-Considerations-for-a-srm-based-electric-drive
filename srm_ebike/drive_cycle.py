"""Pedelec route simulation: motor duty, battery energy and range."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .vehicle import EBikeSpec, motor_torque_demand, road_load_force


@dataclass
class Segment:
    length: float  # m
    grade: float  # fraction
    speed_kmh: float  # target speed
    stop_at_end: bool = False


# Mixed urban/hilly commute, ~6.4 km: start-stop city streets, a 6 % climb,
# descent and a flat run at the assist limit.
DEFAULT_ROUTE = [
    Segment(400, 0.00, 20, True),
    Segment(600, 0.01, 22, True),
    Segment(800, 0.00, 25, False),
    Segment(700, 0.06, 15, False),
    Segment(300, 0.03, 18, False),
    Segment(900, -0.04, 25, True),
    Segment(1500, 0.00, 25, False),
    Segment(500, 0.02, 20, True),
    Segment(700, 0.00, 25, True),
]


@dataclass
class EfficiencyLookup:
    """Interpolates system efficiency over a regular (speed, torque) grid."""

    speeds: np.ndarray
    t_envelope: np.ndarray  # max torque at each speed
    torques: np.ndarray
    eta: np.ndarray  # [n_speed, n_torque], NaN above the envelope

    def __post_init__(self) -> None:
        eta = np.array(self.eta, dtype=float)
        # extend each speed row past the envelope with its last solved value so
        # interpolation near the limit stays well defined
        for a in range(eta.shape[0]):
            row = eta[a]
            ok = ~np.isnan(row)
            if ok.any():
                eta[a] = np.interp(self.torques, self.torques[ok], row[ok])
            else:
                eta[a] = 0.5
        self._interp = RegularGridInterpolator(
            (self.speeds, self.torques), eta, bounds_error=False, fill_value=None
        )

    def t_max(self, n_rpm: float) -> float:
        return float(np.interp(n_rpm, self.speeds, self.t_envelope))

    def eta_at(self, n_rpm: float, torque: float) -> float:
        n = float(np.clip(n_rpm, self.speeds[0], self.speeds[-1]))
        t = float(np.clip(torque, self.torques[0], self.torques[-1]))
        return float(np.clip(self._interp([[n, t]])[0], 0.05, 0.99))


def speed_profile(route, accel: float = 0.6, decel: float = 1.0, dt: float = 0.5):
    """Kinematic speed trace (time, position, speed, grade) for the route."""
    t, x, v = 0.0, 0.0, 0.0
    out = []
    seg_start = 0.0
    for i, seg in enumerate(route):
        v_tgt = seg.speed_kmh / 3.6
        seg_end = seg_start + seg.length
        next_v = 0.0 if seg.stop_at_end else (route[i + 1].speed_kmh / 3.6 if i + 1 < len(route) else 0.0)
        while x < seg_end - 1e-6:
            dist_left = seg_end - x
            # braking distance to the speed required at the end of the segment
            v_brake = math.sqrt(next_v**2 + 2 * decel * dist_left)
            v_cap = min(v_tgt, v_brake)
            if v < v_cap:
                v_new = min(v + accel * dt, v_cap)
            else:
                v_new = max(v - decel * dt, v_cap)
            v_new = max(v_new, 0.3 if dist_left > 0.5 else 0.0)  # creep to the line
            acc = (v_new - v) / dt
            x += 0.5 * (v + v_new) * dt
            v = v_new
            t += dt
            out.append((t, x, v, seg.grade, acc))
        seg_start = seg_end
        if seg.stop_at_end:
            for _ in range(int(10 / dt)):  # 10 s wait at the stop
                t += dt
                out.append((t, x, 0.0, seg.grade, 0.0))
            v = 0.0
    arr = np.array(out)
    return {"t": arr[:, 0], "x": arr[:, 1], "v": arr[:, 2], "grade": arr[:, 3], "a": arr[:, 4]}


def simulate_route(
    spec: EBikeSpec,
    lookup: EfficiencyLookup,
    route=None,
    rider_power: float = 100.0,
    rider_force_max: float = 60.0,
    dt: float = 0.5,
) -> dict:
    """Energy simulation of a pedelec ride.

    The rider supplies up to ``rider_power`` (force limited at low speed); the
    motor supplies the remainder up to its torque envelope and only below the
    assist cut-off speed. Braking is mechanical (no regeneration).
    """
    prof = speed_profile(route or DEFAULT_ROUTE, dt=dt)
    v, grade, acc = prof["v"], prof["grade"], prof["a"]
    n = len(v)
    t_motor = np.zeros(n)
    n_motor = np.zeros(n)
    p_batt = np.zeros(n)
    p_rider = np.zeros(n)
    deficit = np.zeros(n)
    for k in range(n):
        if v[k] <= 0.0:
            continue
        f_req = road_load_force(spec, v[k], grade[k], acc[k])
        f_rider = min(rider_power / max(v[k], 0.5), rider_force_max)
        f_rider = min(f_rider, max(f_req, 0.0))
        p_rider[k] = f_rider * v[k]
        f_motor = f_req - f_rider
        if f_motor <= 0 or v[k] > spec.v_max + 1e-6:
            continue
        n_rpm = spec.motor_speed_rpm(v[k])
        t_dem = motor_torque_demand(spec, f_motor)
        t_avail = lookup.t_max(n_rpm)
        if t_dem > t_avail:
            deficit[k] = t_dem - t_avail
            t_dem = t_avail
        p_mech = t_dem * n_rpm * 2 * math.pi / 60
        t_motor[k] = t_dem
        n_motor[k] = n_rpm
        p_batt[k] = p_mech / lookup.eta_at(n_rpm, t_dem)
    dist_km = prof["x"][-1] / 1000
    e_batt = p_batt.sum() * dt / 3600
    e_mech = float(np.sum(t_motor * n_motor * 2 * math.pi / 60) * dt / 3600)
    usable = 0.9 * spec.battery_energy_wh
    return {
        **prof,
        "t_motor": t_motor,
        "n_motor": n_motor,
        "p_batt": p_batt,
        "p_rider": p_rider,
        "deficit": deficit,
        "distance_km": dist_km,
        "duration_min": prof["t"][-1] / 60,
        "energy_battery_Wh": e_batt,
        "energy_motor_Wh": e_mech,
        "energy_rider_Wh": p_rider.sum() * dt / 3600,
        "cycle_efficiency": e_mech / e_batt if e_batt > 0 else float("nan"),
        "Wh_per_km": e_batt / dist_km,
        "range_km": usable / (e_batt / dist_km) if e_batt > 0 else float("inf"),
        "torque_deficit_s": float((deficit > 1e-6).sum() * dt),
    }
