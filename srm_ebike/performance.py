"""Commutation-angle optimisation, torque-speed envelope and efficiency map."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .drive import SimResult, SRMDrive


@dataclass
class AngleGrid:
    """Search space for (theta_on, theta_off), relative to the overlap angle.

    theta_on  = theta_ov - advance, advance in [0, max_advance] (scaled with
                speed: the current needs time to build up in L_u)
    theta_off = theta_ov + f (tau_r/2 - theta_ov), f in off_fractions
    The dwell is capped at ``max_dwell`` of the rotor pole pitch so that the
    demagnetising interval fits in the period (no continuous conduction).
    """

    n_on: int = 5
    off_fractions: tuple[float, ...] = (0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
    max_dwell: float = 0.55


def _angle_candidates(drive: SRMDrive, n_rpm: float, i_ref: float, grid: AngleGrid):
    mm = drive.model
    theta_ov = mm.overlap_start - mm.fringe / 2
    span = mm.tau_r / 2 - theta_ov
    adv_nominal = drive.turn_on_advance(n_rpm, i_ref)
    adv_max = min(max(1.6 * adv_nominal, math.radians(2.0)), grid.max_dwell * mm.tau_r)
    for f in grid.off_fractions:
        theta_off = theta_ov + f * span
        for adv in np.linspace(0.0, adv_max, grid.n_on):
            theta_on = theta_ov - adv
            if theta_off - theta_on > grid.max_dwell * mm.tau_r:
                theta_on = theta_off - grid.max_dwell * mm.tau_r
            yield theta_on, theta_off


def max_torque(drive: SRMDrive, n_rpm: float, i_ref: float | None = None, grid: AngleGrid | None = None) -> SimResult:
    """Best (highest mean torque) operating point at the current limit."""
    grid = grid or AngleGrid()
    i_ref = drive.i_max if i_ref is None else i_ref
    best = None
    seen = set()
    for on, off in _angle_candidates(drive, n_rpm, i_ref, grid):
        key = (round(on, 6), round(off, 6))
        if key in seen:
            continue
        seen.add(key)
        res = drive.simulate(n_rpm, on, off, i_ref)
        if best is None or res.torque_avg > best.torque_avg:
            best = res
    return best


def operating_point(
    drive: SRMDrive, n_rpm: float, torque: float, grid: AngleGrid | None = None
) -> SimResult | None:
    """Minimum-input-power operating point delivering ``torque`` at ``n_rpm``.

    For each angle pair the chopping reference is solved with Brent's method;
    returns ``None`` if the torque is not reachable within ``drive.i_max``.
    """
    grid = grid or AngleGrid(n_on=3, off_fractions=(0.55, 0.7, 0.85, 0.95))
    best, best_pin = None, math.inf
    # the angle set depends on i_ref through the build-up advance; use the
    # advance for the current limit so the set is fixed during the root search
    for on, off in _angle_candidates(drive, n_rpm, drive.i_max, grid):
        def err(i):
            return drive.simulate(n_rpm, on, off, i).torque_avg - torque

        try:
            if err(drive.i_max) < 0:
                continue
            i_sol = brentq(err, 0.3, drive.i_max, xtol=0.02, rtol=1e-3, maxiter=30)
        except ValueError:
            continue
        res = drive.simulate(n_rpm, on, off, i_sol)
        p_in = res.efficiency()["p_in"]
        if p_in < best_pin:
            best, best_pin = res, p_in
    return best


def torque_speed_envelope(drive: SRMDrive, speeds_rpm, grid: AngleGrid | None = None) -> list[SimResult]:
    return [max_torque(drive, n, grid=grid) for n in speeds_rpm]


def efficiency_map(drive: SRMDrive, speeds_rpm, torques, envelope: list[SimResult] | None = None):
    """Efficiency over a regular (speed, torque) grid.

    Points above 97 % of the envelope torque at a speed are left as NaN.
    Returns (speeds, torques, eta_system[n_speed, n_torque], eta_motor, results).
    """
    envelope = envelope or torque_speed_envelope(drive, speeds_rpm)
    speeds = np.asarray(speeds_rpm, dtype=float)
    torques = np.asarray(torques, dtype=float)
    eta_s = np.full((len(speeds), len(torques)), np.nan)
    eta_m = np.full_like(eta_s, np.nan)
    results = {}
    for a, n in enumerate(speeds):
        t_lim = 0.97 * envelope[a].torque_avg
        for b, t in enumerate(torques):
            if t > t_lim:
                continue
            res = operating_point(drive, n, t)
            if res is None:
                continue
            eff = res.efficiency()
            eta_s[a, b] = eff["eta_system"]
            eta_m[a, b] = eff["eta_motor"]
            results[(a, b)] = res
    return speeds, torques, eta_s, eta_m, results


def required_peak_current(drive: SRMDrive, n_rpm: float, torque: float, i_hi: float = 80.0) -> float:
    """Smallest chopping current giving ``torque`` at ``n_rpm`` with optimal angles."""
    grid = AngleGrid(n_on=3, off_fractions=(0.75, 0.85, 0.95))

    def err(i):
        return max_torque(drive, n_rpm, i, grid).torque_avg - torque

    if err(i_hi) < 0:
        raise ValueError(f"{torque:.2f} N m not reachable at {n_rpm:.0f} rpm below {i_hi} A")
    return brentq(err, 1.0, i_hi, xtol=0.05)
