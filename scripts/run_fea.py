#!/usr/bin/env python3
"""FEA campaign for the baseline SRM (writes results/fea/*.npz|json).

1. psi(i, theta) and Arkkio torque tables for the baseline design
2. mutual-coupling check (two phases excited simultaneously)
3. air-gap tolerance (0.30 / 0.35 / 0.40 mm)
4. pole-arc variants
5. mesh convergence

Usage: python scripts/run_fea.py [--only table,mutual,gap,arcs,mesh]
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

from srm_ebike.design_flow import design_drive  # noqa: E402
from srm_ebike.fea import SRMFEA  # noqa: E402
from srm_ebike.geometry import size_srm  # noqa: E402

OUT = ROOT / "results" / "fea"
THETAS_DEG = np.arange(0.0, 22.5 + 1e-9, 1.5)
CURRENTS = np.array([1.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 23.0, 26.0, 30.0, 34.0, 38.0, 42.0])


def log(msg, t0):
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


def sweep(fea, thetas_deg, currents, phase_currents=lambda i: (i, 0.0, 0.0)):
    psi = np.zeros((len(thetas_deg), len(currents), 3))
    torque = np.zeros((len(thetas_deg), len(currents)))
    bmax = np.zeros_like(torque)
    for a, th in enumerate(thetas_deg):
        a_prev = None
        for b, i in enumerate(currents):
            r, a_prev = fea.solve(math.radians(th), phase_currents(i), a0=a_prev)
            psi[a, b] = r.psi
            torque[a, b] = r.torque
            bmax[a, b] = r.b_pole_max
    return psi, torque, bmax


def task_table(dd, t0):
    fea = SRMFEA(dd.design)
    psi, torque, bmax = sweep(fea, THETAS_DEG, CURRENTS)
    np.savez(
        OUT / "table.npz",
        theta_deg=THETAS_DEG,
        currents=CURRENTS,
        psi=psi[:, :, 0],
        psi_mutual=psi[:, :, 1:],
        torque=torque,
        bmax=bmax,
        l_end=fea.l_end,
    )
    log("table done", t0)


def task_mutual(dd, t0):
    """Phase A and B conducting together (overlap during commutation)."""
    fea = SRMFEA(dd.design)
    thetas = np.arange(6.0, 22.5 + 1e-9, 2.5)
    rows = []
    for th in thetas:
        for ia, ib in [(15.0, 15.0), (30.0, 30.0), (30.0, 10.0)]:
            ra, _ = fea.solve(math.radians(th), (ia, 0.0, 0.0))
            # phase B axis is +30 deg (one stator pole pitch); its phase angle
            # is theta - stroke when the rotor is at phase-A angle theta
            rb, _ = fea.solve(math.radians(th), (0.0, ib, 0.0))
            rab, _ = fea.solve(math.radians(th), (ia, ib, 0.0))
            rows.append(
                {
                    "theta_deg": float(th),
                    "iA": ia,
                    "iB": ib,
                    "psiA_alone_mWb": ra.psi[0] * 1e3,
                    "psiA_both_mWb": rab.psi[0] * 1e3,
                    "psiB_alone_mWb": rb.psi[1] * 1e3,
                    "psiB_both_mWb": rab.psi[1] * 1e3,
                    "T_sum_alone": ra.torque + rb.torque,
                    "T_both": rab.torque,
                }
            )
    (OUT / "mutual.json").write_text(json.dumps(rows, indent=1))
    log("mutual done", t0)


def variant(dd, **over):
    inp = copy.deepcopy(dd.inputs)
    for k, v in over.items():
        setattr(inp, k, v)
    d = size_srm(inp)
    # keep the baseline winding (same turns) so only the magnetics change
    d.turns_per_phase, d.turns_per_pole = dd.design.turns_per_phase, dd.design.turns_per_pole
    return d


def static_metrics(d, currents=(15.0, 30.0)):
    fea = SRMFEA(d)
    tau = 360.0 / d.nr
    th = np.linspace(0.0, tau / 2, 10)
    out = {}
    for i in currents:
        t = []
        a_prev = None
        for x in th:
            r, a_prev = fea.solve(math.radians(x), (i, 0.0, 0.0))
            t.append(r.torque)
        out[f"T_avg_motoring_{i:.0f}A"] = float(np.trapezoid(t, th) / (tau / 2))
        out[f"T_peak_{i:.0f}A"] = float(max(t))
    ra, _ = fea.solve(math.radians(tau / 2), (2.0, 0.0, 0.0))
    ru, _ = fea.solve(0.0, (2.0, 0.0, 0.0))
    out["L_a_mH"] = ra.psi[0] / 2.0 * 1e3
    out["L_u_mH"] = ru.psi[0] / 2.0 * 1e3
    out["ratio"] = out["L_a_mH"] / out["L_u_mH"]
    # average static torque over a full stroke with ideal flat-top current
    # = (strokes/rev)/(2 pi) * co-energy difference, approx from T_avg_motoring
    return out


def task_gap(dd, t0):
    rows = []
    for g in (0.30e-3, 0.35e-3, 0.40e-3):
        d = variant(dd, air_gap=g)
        rows.append({"air_gap_mm": g * 1e3, **static_metrics(d)})
        log(f"gap {g * 1e3:.2f} mm done", t0)
    (OUT / "airgap.json").write_text(json.dumps(rows, indent=1))


def task_arcs(dd, t0):
    rows = []
    for bs, br in [(15, 17), (15, 15), (15, 19), (16, 18), (17, 19), (15, 22)]:
        d = variant(dd, beta_s_deg=bs, beta_r_deg=br)
        rows.append({"beta_s": bs, "beta_r": br, **static_metrics(d)})
        log(f"arcs {bs}/{br} done", t0)
    (OUT / "arcs.json").write_text(json.dumps(rows, indent=1))


def task_mesh(dd, t0):
    rows = []
    for div, hi, ha in [(2, 2.0e-3, 2.5e-3), (3, 1.2e-3, 1.5e-3), (5, 0.8e-3, 1.0e-3)]:
        f = SRMFEA(dd.design, h_gap=dd.design.air_gap / div, h_iron=hi, h_air=ha)
        ra, _ = f.solve(math.radians(22.5), (15, 0, 0))
        rm, _ = f.solve(math.radians(11), (15, 0, 0))
        ru, _ = f.solve(0.0, (15, 0, 0))
        rows.append(
            {
                "h_gap_mm": dd.design.air_gap / div * 1e3,
                "nodes": ra.n_nodes,
                "psi_aligned_mWb": ra.psi_2d[0] * 1e3,
                "psi_11deg_mWb": rm.psi_2d[0] * 1e3,
                "psi_unaligned_mWb": ru.psi_2d[0] * 1e3,
                "T_11deg_Nm": rm.torque,
            }
        )
    (OUT / "mesh.json").write_text(json.dumps(rows, indent=1))
    log("mesh done", t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="table,mutual,gap,arcs,mesh")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    dd = design_drive()
    tasks = {"table": task_table, "mutual": task_mutual, "gap": task_gap, "arcs": task_arcs, "mesh": task_mesh}
    for name in args.only.split(","):
        tasks[name](dd, t0)


if __name__ == "__main__":
    main()
