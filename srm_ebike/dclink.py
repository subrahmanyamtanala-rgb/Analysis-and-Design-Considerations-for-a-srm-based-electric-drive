"""DC-link analysis: ripple-current sharing between capacitor bank and battery.

The DC-side converter current i_dc(t) is periodic with the rotor pole pitch.
Its AC part divides between the capacitor branch Z_c = ESR + 1/(j w C) and the
battery branch Z_b = R_b + j w L_b (internal resistance + cable inductance),
harmonic by harmonic:

    I_c(k) = I_ac(k) Z_b / (Z_b + Z_c),     V_link(k) = I_c(k) Z_c

The simple charge-balance bound used for first sizing is
C >= dQ / dV with dQ = max - min of the integral of the AC current, which
assumes the battery supplies only the mean current (worst case).
"""

from __future__ import annotations

import numpy as np


def charge_bound(i_dc: np.ndarray, period: float, dv: float) -> dict:
    dt = period / len(i_dc)
    q = np.cumsum(i_dc - i_dc.mean()) * dt
    dq = float(q.max() - q.min())
    return {"dQ_uC": dq * 1e6, "C_min_uF": dq / dv * 1e6}


def sharing(i_dc: np.ndarray, period: float, c: float, esr: float, r_b: float, l_b: float) -> dict:
    n = len(i_dc)
    spec = np.fft.rfft(i_dc - i_dc.mean())
    k = np.arange(len(spec))
    w = 2 * np.pi * k / period
    w[0] = 1.0
    zc = esr + 1 / (1j * w * c)
    zb = r_b + 1j * w * l_b
    ic = spec * zb / (zb + zc)
    ib = spec * zc / (zb + zc)
    ic[0] = ib[0] = 0
    vc = ic * zc
    i_c = np.fft.irfft(ic, n)
    i_b = np.fft.irfft(ib, n)
    v = np.fft.irfft(vc, n)
    return {
        "C_uF": c * 1e6,
        "I_cap_rms_A": float(np.sqrt(np.mean(i_c**2))),
        "I_bat_ac_rms_A": float(np.sqrt(np.mean(i_b**2))),
        "V_ripple_pp_V": float(v.max() - v.min()),
        "P_esr_W": float(np.mean(i_c**2) * esr),
    }
