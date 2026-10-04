"""Flux-linkage model built from FEA tables psi(i, theta).

Exposes the same interface as :class:`magnetics.MagneticModel` (``psi``,
``torque``, ``current``, ``tau_r``, ``overlap_start``, ``fringe``, ``l_u``), so
it can replace the analytical model in :class:`drive.SRMDrive`.

Construction:
* the FEA half-period table (unaligned -> aligned) is mirrored to a full rotor
  pole pitch, interpolated with cubic splines in theta and monotone PCHIP in i,
  and extrapolated above the last FEA current with the final incremental
  inductance;
* the co-energy W'(i, theta) is integrated from that same psi table and torque
  is its theta-derivative, so the drive model is energy-consistent;
* dense lookup tables (psi, torque, and the inverse i(psi, theta)) give O(1)
  bilinear evaluation inside the time-stepping loop.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from .geometry import SRMDesign
from .magnetics import MagneticModel


class TabulatedMagneticModel:
    def __init__(
        self,
        design: SRMDesign,
        theta_deg: np.ndarray,
        currents: np.ndarray,
        psi_table: np.ndarray,
        i_max: float = 60.0,
        n_theta: int = 721,
        n_i: int = 601,
        n_psi: int = 601,
    ):
        ref = MagneticModel.from_design(design)
        self.nr = design.nr
        self.beta_s, self.beta_r = ref.beta_s, ref.beta_r
        self.fringe = ref.fringe
        self.tau_r = 2 * math.pi / self.nr
        self.overlap_start = ref.overlap_start
        self.full_overlap = ref.full_overlap

        th_half = np.radians(np.asarray(theta_deg, float))
        cur = np.concatenate([[0.0], np.asarray(currents, float)])
        psi_half = np.vstack([np.zeros(len(th_half)), np.asarray(psi_table, float).T]).T  # [theta, i]

        # mirror to a full pitch: psi(tau - th) = psi(th)
        th_full = np.concatenate([th_half, self.tau_r - th_half[-2::-1]])
        psi_full = np.vstack([psi_half, psi_half[-2::-1]])

        # dense current grid with linear extrapolation beyond the FEA range
        ii = np.linspace(0.0, i_max, n_i)
        psi_i = np.empty((len(th_full), n_i))
        for a in range(len(th_full)):
            f = PchipInterpolator(cur, psi_full[a])
            slope = (psi_full[a, -1] - psi_full[a, -2]) / (cur[-1] - cur[-2])
            psi_i[a] = np.where(ii <= cur[-1], f(np.minimum(ii, cur[-1])), psi_full[a, -1] + slope * (ii - cur[-1]))

        # dense theta grid (periodic cubic spline per current)
        tt = np.linspace(0.0, self.tau_r, n_theta)
        psi_t = np.empty((n_theta, n_i))
        for b in range(n_i):
            cs = CubicSpline(th_full, psi_i[:, b], bc_type="periodic")
            psi_t[:, b] = cs(tt)
        psi_t = np.maximum(psi_t, 0.0)
        # enforce monotonicity in i (numerical safety for the inverse table)
        psi_t = np.maximum.accumulate(psi_t, axis=1)

        di = ii[1] - ii[0]
        coenergy = np.concatenate([np.zeros((n_theta, 1)), np.cumsum(0.5 * (psi_t[:, 1:] + psi_t[:, :-1]) * di, axis=1)], axis=1)
        dth = tt[1] - tt[0]
        torque = (np.roll(coenergy, -1, axis=0) - np.roll(coenergy, 1, axis=0)) / (2 * dth)
        torque[0] = (coenergy[1] - coenergy[-2]) / (2 * dth)
        torque[-1] = torque[0]

        psi_max = float(psi_t.max())
        pp = np.linspace(0.0, psi_max, n_psi)
        inv = np.empty((n_theta, n_psi))
        for a in range(n_theta):
            row = psi_t[a]
            # strictly increasing for interpolation
            row = row + np.arange(n_i) * 1e-12
            inv[a] = np.interp(pp, row, ii, right=ii[-1] + (pp[-1] - row[-1]) / max(row[-1] - row[-2], 1e-9) * di)

        self._tt, self._ii, self._pp = tt, ii, pp
        self._psi, self._torque, self._inv = psi_t, torque, inv
        self._dth, self._di, self._dp = dth, di, pp[1] - pp[0]
        self.n_theta, self.n_i, self.n_psi = n_theta, n_i, n_psi

        # characteristic inductances for the controller heuristics
        a_u, a_a = 0, (n_theta - 1) // 2
        k = max(1, int(2.0 / di))
        self.l_u = float(psi_t[a_u, k] / ii[k])
        self.l_a = float(psi_t[a_a, k] / ii[k])
        self.l_sat = float((psi_t[a_a, -1] - psi_t[a_a, -k - 1]) / (ii[-1] - ii[-k - 1]))
        self.psi_k = float(psi_t[a_a, -1] - self.l_sat * ii[-1])

    # ------------------------------------------------------------------
    def _bilinear(self, table, x, y, dx, dy, nx, ny):
        fx = x / dx
        fy = y / dy
        ix = int(fx)
        iy = int(fy)
        if ix >= nx - 1:
            ix, fx = nx - 2, float(nx - 1)
        if iy >= ny - 1:
            iy = ny - 2
            # linear extrapolation along y beyond the table
        tx = fx - ix
        ty = fy - iy
        r0 = table[ix]
        r1 = table[ix + 1]
        a = r0[iy] + (r0[iy + 1] - r0[iy]) * ty
        b = r1[iy] + (r1[iy + 1] - r1[iy]) * ty
        return a + (b - a) * tx

    def psi(self, i: float, theta: float) -> float:
        return self._bilinear(self._psi, theta % self.tau_r, abs(i), self._dth, self._di, self.n_theta, self.n_i)

    def torque(self, i: float, theta: float) -> float:
        return self._bilinear(self._torque, theta % self.tau_r, abs(i), self._dth, self._di, self.n_theta, self.n_i)

    def current(self, psi: float, theta: float, i_guess: float | None = None) -> float:
        if psi <= 0.0:
            return 0.0
        return max(0.0, self._bilinear(self._inv, theta % self.tau_r, psi, self._dth, self._dp, self.n_theta, self.n_psi))

    def w(self, theta: float) -> float:  # for API compatibility in plots
        a = self.psi(2.0, theta)
        return (a - self.l_u * 2.0) / max((self.l_a - self.l_u) * 2.0, 1e-12)

    def curves(self, i_max: float, thetas: np.ndarray, n: int = 100):
        ii = np.linspace(0, i_max, n)
        return ii, np.array([[self.psi(i, th) for i in ii] for th in thetas])

    def static_torque(self, currents: np.ndarray, thetas: np.ndarray) -> np.ndarray:
        return np.array([[self.torque(i, th) for th in thetas] for i in currents])


def load_fea_model(design: SRMDesign, path) -> TabulatedMagneticModel:
    data = np.load(path)
    return TabulatedMagneticModel(design, data["theta_deg"], data["currents"], data["psi"])
