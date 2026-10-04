"""Nonlinear flux-linkage model psi(i, theta) of one SRM phase.

The aligned magnetisation curve is an exponential knee model

    psi_a(i) = L_sat i + psi_k (1 - exp(-(L_a - L_sat) i / psi_k))

whose slope is the unsaturated aligned inductance L_a at i = 0 and L_sat in
deep saturation. The unaligned curve is linear, psi_u(i) = L_u i. The two are
blended by a position function w(theta) in [0, 1]:

    psi(i, theta) = psi_u(i) + w(theta) (psi_a(i) - psi_u(i))

w is 1 while the stator pole is fully covered by the (wider) rotor pole, 0 when
there is no overlap, and follows a C1 smooth-step across the overlap region,
which is widened by a fringing angle on either side. Because w multiplies an
i-only function, co-energy and torque have closed forms:

    W'(i, theta) = L_u i^2/2 + w(theta) (W'_a(i) - L_u i^2/2)
    T(i, theta)  = dw/dtheta (W'_a(i) - L_u i^2/2)

Phase angle convention: theta = 0 is the unaligned position, theta = tau_r/2
(half the rotor pole pitch) is aligned. Motoring torque is produced for
theta in (0, tau_r/2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geometry import SRMDesign


@dataclass
class MagneticModel:
    l_a: float
    l_u: float
    l_sat: float
    psi_k: float
    beta_s: float
    beta_r: float
    nr: int
    fringe: float = 0.0

    @classmethod
    def from_design(cls, d: SRMDesign) -> "MagneticModel":
        return cls(
            l_a=d.l_aligned,
            l_u=d.l_unaligned,
            l_sat=d.l_sat,
            psi_k=d.psi_knee,
            beta_s=d.beta_s,
            beta_r=d.beta_r,
            nr=d.nr,
            fringe=d.fringe_angle,
        )

    # ------------------------------------------------------------------
    @property
    def tau_r(self) -> float:
        return 2 * math.pi / self.nr

    @property
    def overlap_start(self) -> float:
        """Phase angle where the pole tips start to overlap (excl. fringing)."""
        return self.tau_r / 2 - (self.beta_s + self.beta_r) / 2

    @property
    def full_overlap(self) -> float:
        """Phase angle where the stator pole becomes fully covered."""
        return self.tau_r / 2 - abs(self.beta_r - self.beta_s) / 2

    def _x_limits(self) -> tuple[float, float]:
        # distance from the aligned position where w reaches 1 and 0
        x1 = abs(self.beta_r - self.beta_s) / 2
        x2 = (self.beta_s + self.beta_r) / 2
        # fringing widens the transition symmetrically but cannot push it past
        # the aligned or unaligned positions
        x1 = max(0.0, x1 - self.fringe / 2)
        x2 = min(self.tau_r / 2, x2 + self.fringe / 2)
        return x1, x2

    def w(self, theta: float) -> float:
        x = abs(theta % self.tau_r - self.tau_r / 2)
        x1, x2 = self._x_limits()
        if x <= x1:
            return 1.0
        if x >= x2:
            return 0.0
        s = (x2 - x) / (x2 - x1)
        return s * s * (3 - 2 * s)

    def dw_dtheta(self, theta: float) -> float:
        th = theta % self.tau_r
        x = abs(th - self.tau_r / 2)
        x1, x2 = self._x_limits()
        if x <= x1 or x >= x2:
            return 0.0
        s = (x2 - x) / (x2 - x1)
        dw_dx = -6 * s * (1 - s) / (x2 - x1)
        # dx/dtheta = -1 before alignment (th < tau/2), +1 after
        return -dw_dx if th < self.tau_r / 2 else dw_dx

    # ------------------------------------------------------------------
    def psi_aligned(self, i: float) -> float:
        k = (self.l_a - self.l_sat) / self.psi_k
        return self.l_sat * i + self.psi_k * (1 - math.exp(-k * i))

    def coenergy_aligned(self, i: float) -> float:
        k = (self.l_a - self.l_sat) / self.psi_k
        return 0.5 * self.l_sat * i * i + self.psi_k * i - self.psi_k / k * (1 - math.exp(-k * i))

    def psi(self, i: float, theta: float) -> float:
        i = abs(i)
        pu = self.l_u * i
        return pu + self.w(theta) * (self.psi_aligned(i) - pu)

    def torque(self, i: float, theta: float) -> float:
        i = abs(i)
        return self.dw_dtheta(theta) * (self.coenergy_aligned(i) - 0.5 * self.l_u * i * i)

    def current(self, psi: float, theta: float, i_guess: float | None = None) -> float:
        """Invert psi(i, theta) for i >= 0 (Newton, psi is monotonic in i)."""
        if psi <= 0.0:
            return 0.0
        w = self.w(theta)
        k = (self.l_a - self.l_sat) / self.psi_k
        lu, ls, pk = self.l_u, self.l_sat, self.psi_k
        i = i_guess if i_guess and i_guess > 0 else psi / (lu + w * (self.l_a - lu))
        for _ in range(30):
            e = math.exp(-k * i)
            f = (1 - w) * lu * i + w * (ls * i + pk * (1 - e)) - psi
            df = (1 - w) * lu + w * (ls + (self.l_a - ls) * e)
            step = f / df
            i -= step
            if i < 0:
                i = 0.5 * (i + step)  # back off towards the previous iterate
            if abs(step) < 1e-9 * (1 + i):
                break
        return i

    # ------------------------------------------------------------------
    def curves(self, i_max: float, thetas: np.ndarray, n: int = 100) -> tuple[np.ndarray, np.ndarray]:
        """Magnetisation curves psi(i) for each angle in ``thetas``."""
        ii = np.linspace(0, i_max, n)
        psi = np.array([[self.psi(i, th) for i in ii] for th in thetas])
        return ii, psi

    def static_torque(self, currents: np.ndarray, thetas: np.ndarray) -> np.ndarray:
        return np.array([[self.torque(i, th) for th in thetas] for i in currents])
