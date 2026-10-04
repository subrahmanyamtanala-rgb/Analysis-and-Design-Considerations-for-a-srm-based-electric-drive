"""Torque-sharing-function (TSF) control.

Each phase receives a torque reference T_k* = T* f_k(theta). The rising
function of the incoming phase and the falling function of the outgoing phase
overlap over theta_ov and always sum to one, so sum_k T_k* = T*. Phase
torque references are converted to current references through the inverse
static torque characteristic i(T, theta) and tracked by a three-level
hysteresis controller (+V, 0, -V) in the asymmetric half bridge.

Shapes (x in [0, 1] across the overlap): linear, cubic, cosine and exponential
(Xue et al., IEEE TPEL 2009; Vujicic, IEEE TPEL 2012).
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from .drive import SimResult, SRMDrive

SHAPES = {
    "linear": lambda x: x,
    "cubic": lambda x: x * x * (3 - 2 * x),
    "cosine": lambda x: 0.5 * (1 - math.cos(math.pi * x)),
    "exponential": lambda x: (1 - math.exp(-x * x / 0.15)) / (1 - math.exp(-1 / 0.15)),
}


class InverseTorque:
    """i(T, theta) from the static characteristic of the drive's magnetic model."""

    def __init__(self, drive: SRMDrive, i_max: float, n_theta: int = 451, n_i: int = 241):
        mm = drive.model
        self.tau = mm.tau_r
        self.i_max = i_max
        self.th = np.linspace(0.0, self.tau, n_theta)
        self.ii = np.linspace(0.0, i_max, n_i)
        self.t = np.array([[mm.torque(i, x) for i in self.ii] for x in self.th])

    def __call__(self, torque: float, theta: float) -> float:
        if torque <= 0.0:
            return 0.0
        x = theta % self.tau
        a = min(int(x / (self.th[1] - self.th[0])), len(self.th) - 1)
        row = np.maximum.accumulate(self.t[a])
        if torque >= row[-1]:
            return self.i_max
        return float(np.interp(torque, row, self.ii))


def tsf_profile(shape: str, eps: float, theta_ov: float):
    f = SHAPES[shape]

    def share(phi: float) -> float:
        if phi < 0:
            return 0.0
        if phi < theta_ov:
            return f(phi / theta_ov)
        if phi < eps:
            return 1.0
        if phi < eps + theta_ov:
            return 1.0 - f((phi - eps) / theta_ov)
        return 0.0

    return share


def simulate_tsf(drive: SRMDrive, inv: InverseTorque, n_rpm: float, t_ref: float, theta_on: float, theta_ov: float, shape: str) -> SimResult:
    eps = drive.design.stroke_angle
    share = tsf_profile(shape, eps, theta_ov)

    def ref(phi: float) -> float:
        return inv(t_ref * share(phi), theta_on + phi)

    return drive.simulate(n_rpm, theta_on, theta_on + eps + theta_ov, inv.i_max, control="tsf", ref_profile=ref)


def tsf_operating_point(
    drive: SRMDrive,
    inv: InverseTorque,
    n_rpm: float,
    torque: float,
    shape: str = "cubic",
    theta_on_deg=(3.0, 4.5, 6.0),
    theta_ov_deg=(2.5, 4.0, 5.5),
    objective: str = "ripple",
) -> SimResult | None:
    """Best TSF operating point meeting the mean torque.

    The torque command is trimmed with Brent's method so that the simulated mean
    torque equals ``torque``; among the (theta_on, theta_ov) grid the point
    with minimum ripple (or minimum input power) is returned.
    """
    best, best_val = None, math.inf
    for on in theta_on_deg:
        for ov in theta_ov_deg:
            th_on, th_ov = math.radians(on), math.radians(ov)

            def err(tr):
                return simulate_tsf(drive, inv, n_rpm, tr, th_on, th_ov, shape).torque_avg - torque

            try:
                hi = 2.5 * torque
                if err(hi) < 0:
                    continue
                t_cmd = brentq(err, 0.3 * torque, hi, xtol=1e-3 * torque, maxiter=25)
            except ValueError:
                continue
            res = simulate_tsf(drive, inv, n_rpm, t_cmd, th_on, th_ov, shape)
            val = res.torque_ripple if objective == "ripple" else res.efficiency()["p_in"]
            if val < best_val:
                best, best_val = res, val
                best.tsf = {"shape": shape, "theta_on_deg": on, "theta_ov_deg": ov, "t_cmd": t_cmd}
    return best
