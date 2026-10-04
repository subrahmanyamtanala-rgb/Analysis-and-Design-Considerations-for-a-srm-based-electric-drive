#!/usr/bin/env python3
"""Second-round review studies (needs results/fea/table.npz).

Tasks (each writes results/revision2/<task>.json):
  battery  - DC-link voltage 30/36/42 V: launch current, base-speed torque, speed limit,
             rated efficiency; FE efficiency maps per voltage; routes with a battery
             equivalent circuit (OCV(SOC), R_int, sag) and range to cut-off
  ironloss - iron-loss model x0.7 / x1.0 / x1.3: operating points, maps, routes
  coupled  - three-phase FE solves with the simulated phase currents (mutual coupling
             included) vs the uncoupled table model, for CCC and TSF points
  pm       - first-order analytic surface-PM benchmark at equal envelope and boundary
  corner   - quantitative test of P_c = k_c V I_pk: nine-design regression and a
             design-voltage sweep (24/36/48-V packs)

Usage: python scripts/run_revision2.py --only battery,ironloss,coupled,corner
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
import json
import math
import sys
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
warnings.filterwarnings("ignore")

from make_paper_figures import C, COL, INK2, LS, save  # noqa: E402

from srm_ebike.design_flow import design_drive  # noqa: E402
from srm_ebike.drive_cycle import (  # noqa: E402
    ROUTES,
    BatteryModel,
    EfficiencyLookup,
    VoltageLookup,
    simulate_route,
    simulate_route_battery,
)
from srm_ebike.performance import efficiency_map, max_torque, operating_point, required_peak_current, torque_speed_envelope  # noqa: E402
from srm_ebike.vehicle import EBikeSpec  # noqa: E402

TABLE = ROOT / "results" / "fea" / "table.npz"
OUT = ROOT / "results" / "revision2"
FIG = ROOT / "paper" / "figures"
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def eff(res):
    e = res.efficiency()
    return {"T": res.torque_avg, "i_rms": res.i_rms, "eta_motor": e["eta_motor"], "eta_system": e["eta_system"], "p_in": e["p_in"], "losses": res.losses()}


def build_lookup(dr, req, n_speed=10, n_torque=14):
    speeds = np.linspace(0.08, 1.2, n_speed) * req.base_speed_rpm
    env = torque_speed_envelope(dr, speeds)
    t_env = np.array([e.torque_avg for e in env])
    sp, tq, eta_s, eta_m, _ = efficiency_map(dr, speeds, np.linspace(0.25, t_env.max(), n_torque), env)
    return EfficiencyLookup(sp, t_env, tq, eta_s), float(np.nanmax(eta_s)), float(np.nanmax(eta_m))


def strip(d):
    return {k: v for k, v in d.items() if not isinstance(v, np.ndarray)}


# ---------------------------------------------------------------------------
def task_battery(dd):
    req = dd.req
    rows, lookups = [], {}
    for vdc in (30.0, 36.0, 42.0):
        dr = dataclasses.replace(dd.drive, v_dc=vdc)
        i_req = required_peak_current(dr, 0.1 * req.corner_speed_rpm, req.peak_torque)
        t_launch = max_torque(dr, 0.1 * req.corner_speed_rpm).torque_avg
        t_base = max_torque(dr, req.base_speed_rpm).torque_avg

        def short(n):
            return max_torque(dr, n).torque_avg - req.rated_torque

        n_hi = 4.0 * req.base_speed_rpm
        n_lim = brentq(short, req.base_speed_rpm, n_hi, xtol=10.0) if short(n_hi) < 0 else n_hi
        rated = operating_point(dr, req.base_speed_rpm, req.rated_torque)
        lk, es, em = build_lookup(dr, req)
        lookups[vdc] = lk
        rows.append(
            {
                "v_dc": vdc,
                "i_required_launch": i_req,
                "T_launch_at_Imax": t_launch,
                "T_max_base": t_base,
                "P_max_base": t_base * req.base_speed_rpm * 2 * math.pi / 60,
                "n_rated_torque_limit_rpm": n_lim,
                "v_rated_torque_limit_kmh": n_lim / req.base_speed_rpm * dd.spec.v_max_kmh,
                "rated": eff(rated) if rated else None,
                "eta_sys_max": es,
                "eta_mot_max": em,
            }
        )
        log(f"battery: {vdc} V done")
    vl = VoltageLookup(lookups)
    bat = BatteryModel()
    routes = []
    for name, route in ROUTES.items():
        fixed = simulate_route(dd.spec, lookups[36.0], route)
        row = {"route": name, "distance_km": fixed["distance_km"], "fixed36_Wh_per_km": fixed["Wh_per_km"], "fixed36_range_km": fixed["range_km"]}
        for soc0 in (1.0, 0.5, 0.2):
            r = simulate_route_battery(dd.spec, vl, bat, route, soc0=soc0)
            row[f"soc{int(100 * soc0)}"] = strip(r)
        full = simulate_route_battery(dd.spec, vl, bat, route, soc0=1.0, repeat_until_empty=True)
        row["range_to_cutoff_km"] = full["distance_km"]
        row["range_soc_end"] = full["soc_end"]
        row["range_deficit_s"] = full["torque_deficit_s"]
        routes.append(row)
        log(f"battery route {name} done")
    (OUT / "battery.json").write_text(json.dumps({"rows": rows, "routes": routes, "battery": dataclasses.asdict(bat)}, indent=1, default=float))


def task_ironloss(dd):
    req = dd.req
    out = []
    for f in (0.7, 1.0, 1.3):
        dr = copy.deepcopy(dd.drive)
        dr.iron = dataclasses.replace(dr.iron, k_nonsine=dr.iron.k_nonsine * f)
        pts = {
            "hill": operating_point(dr, req.base_speed_rpm * 15 / dd.spec.v_max_kmh, req.hill_torque),
            "rated": operating_point(dr, req.base_speed_rpm, req.rated_torque),
            "cruise": operating_point(dr, req.base_speed_rpm, req.cruise_torque),
        }
        lk, es, em = build_lookup(dr, req)
        routes = {name: simulate_route(dd.spec, lk, r)["Wh_per_km"] for name, r in ROUTES.items()}
        out.append({"factor": f, "points": {k: eff(v) for k, v in pts.items()}, "eta_sys_max": es, "eta_mot_max": em, "routes_Wh_per_km": routes})
        log(f"ironloss x{f} done")
    (OUT / "ironloss.json").write_text(json.dumps(out, indent=1, default=float))


def task_coupled(dd):
    """Apply the three simulated phase currents simultaneously in the FE model."""
    from srm_ebike.fea import SRMFEA
    from srm_ebike.tsf import InverseTorque, tsf_operating_point

    dr, req = dd.drive, dd.req
    mm = dr.model
    eps = dd.design.stroke_angle
    tau = mm.tau_r
    fea = SRMFEA(dd.design)
    inv = InverseTorque(dr, dr.i_max)
    cases = {
        "launch_ccc": max_torque(dr, 0.1 * req.corner_speed_rpm),
        "200rpm_3Nm_ccc": operating_point(dr, 200, 3.0),
        "200rpm_3Nm_tsf": tsf_operating_point(dr, inv, 200, 3.0, "cubic", theta_on_deg=(3.0, 4.5, 6.0, 7.5), theta_ov_deg=(1.5, 3.0, 4.5)),
        "rated_ccc": operating_point(dr, req.base_speed_rpm, req.rated_torque),
    }
    rows, waves = [], {}
    n_pos = 24
    for name, res in cases.items():
        ph = (res.theta - res.theta[0] + res.theta_on) % tau
        order = np.argsort(ph)
        ph_s, cur_s = ph[order], res.current[order]

        def i_at(x):
            return float(np.interp(x % tau, ph_s, cur_s, period=tau))

        # one stroke of rotor positions (phase-A angle), starting at phase-A turn-on
        thetas = res.theta_on + np.linspace(0, eps, n_pos, endpoint=False)
        t_c, t_u = [], []
        a_prev = None
        for th in thetas:
            # phase B (pole +30 deg) sees phase angle th + eps, phase C th - eps
            ia, ib, ic = i_at(th), i_at(th + eps), i_at(th - eps)
            r, a_prev = fea.solve(th % tau, (ia, ib, ic), a0=a_prev)
            t_c.append(r.torque)
            t_u.append(mm.torque(ia, th) + mm.torque(ib, th + eps) + mm.torque(ic, th - eps))
        t_c, t_u = np.array(t_c), np.array(t_u)

        def rip(t):
            return float((t.max() - t.min()) / t.mean() * 100)

        rows.append(
            {
                "case": name,
                "n_rpm": res.n_rpm,
                "T_mean_uncoupled": float(t_u.mean()),
                "T_mean_coupled": float(t_c.mean()),
                "mean_diff_pct": float(100 * (t_c.mean() - t_u.mean()) / t_u.mean()),
                "ripple_uncoupled_pct": rip(t_u),
                "ripple_coupled_pct": rip(t_c),
                "ripple_from_drive_sim_pct": 100 * res.torque_ripple,
                "mesh_retries_cumulative": getattr(fea, "n_mesh_retries", 0),
            }
        )
        waves[name] = (np.degrees(thetas - thetas[0]), t_u, t_c)
        log(f"coupled {name} done")
    (OUT / "coupled.json").write_text(json.dumps(rows, indent=1, default=float))
    fig, ax = plt.subplots(1, 2, figsize=(COL, 1.8), sharey=False)
    for k, (name, lab) in enumerate((("200rpm_3Nm_ccc", "CCC"), ("200rpm_3Nm_tsf", "Cubic TSF"))):
        x, tu, tc = waves[name]
        ax[k].plot(x, tu, color=C[0], ls=LS[0], lw=1.2, label="uncoupled (FE table)")
        ax[k].plot(x, tc, "o-", color=C[1], ms=2.5, lw=0.9, label="coupled 3-phase FE")
        ax[k].set_title(f"({'ab'[k]}) {lab}, 200 r/min, 3 N$\\cdot$m", loc="left", fontsize=7)
        ax[k].set_xlabel("Angle in stroke (deg)")
        ax[k].set_ylim(bottom=0)
    ax[0].set_ylabel("Torque (N$\\cdot$m)")
    ax[1].legend(fontsize=5.5, loc="lower center")
    fig.tight_layout(w_pad=0.6)
    save(fig, FIG, "fig_coupled")


def task_corner(dd):
    from run_analysis import trade_study

    trades = trade_study(False)
    mass = np.array([r["active_mass_kg"] for r in trades])
    i_req = np.array([r["I_peak_A"] for r in trades]) / 1.1
    p_c, v_min = dd.spec.peak_power, dd.spec.battery_voltage_min
    kc = p_c / (v_min * i_req)
    slope, icpt = np.polyfit(mass, i_req, 1)
    pred = icpt + slope * mass
    r2 = 1 - np.sum((i_req - pred) ** 2) / np.sum((i_req - i_req.mean()) ** 2)
    kc_mean = float(kc.mean())
    i_rule = p_c / (kc_mean * v_min)
    dev_rule = float(np.max(np.abs(i_req - i_rule) / i_rule) * 100)
    # design-voltage sweep (same G, A -> same frame), 24/36/48-V packs
    vrows = []
    for v_min_d, v_nom, v_max in ((20.0, 24.0, 28.0), (30.0, 36.0, 42.0), (40.0, 48.0, 54.6)):
        spec = EBikeSpec(battery_voltage_min=v_min_d, battery_voltage=v_nom, battery_voltage_max=v_max)
        d2 = design_drive(spec)
        vrows.append(
            {
                "v_min": v_min_d,
                "v_nom": v_nom,
                "turns_per_pole": d2.design.turns_per_pole,
                "mass_kg": d2.design.mass_active,
                "i_required": d2.i_peak_required,
                "k_c": p_c / (v_min_d * d2.i_peak_required),
                "i_times_v": d2.i_peak_required * v_min_d,
            }
        )
        log(f"corner: {v_nom} V pack done")
    lv, li = np.log([r["v_min"] for r in vrows]), np.log([r["i_required"] for r in vrows])
    exp_v = float(np.polyfit(lv, li, 1)[0])
    res = {
        "nine_designs": {
            "mass_kg": mass.tolist(),
            "i_required": i_req.tolist(),
            "k_c": kc.tolist(),
            "slope_A_per_kg": float(slope),
            "intercept_A": float(icpt),
            "r2_linear": float(r2),
            "rel_slope_pct_per_kg": float(100 * slope / i_req.mean()),
            "k_c_mean": kc_mean,
            "k_c_std": float(kc.std()),
            "i_rule_A": i_rule,
            "max_dev_from_rule_pct": dev_rule,
            "i_range_pct": float(100 * (i_req.max() - i_req.min()) / i_req.mean()),
            "mass_ratio": float(mass.max() / mass.min()),
        },
        "voltage_sweep": vrows,
        "voltage_exponent": exp_v,
    }
    (OUT / "corner.json").write_text(json.dumps(res, indent=1, default=float))
    fig, ax = plt.subplots(1, 2, figsize=(COL, 1.85))
    g = [r["gear_ratio"] for r in trades]
    for k, gg in enumerate(sorted(set(g))):
        sel = [j for j in range(len(g)) if g[j] == gg]
        ax[0].plot(mass[sel], i_req[sel], ["o", "s", "^"][k], color=C[k], ms=3.5, label=f"$G$={gg}")
    xx = np.linspace(mass.min() * 0.9, mass.max() * 1.05, 20)
    ax[0].plot(xx, icpt + slope * xx, color=INK2, lw=0.8, ls="--", label="linear fit")
    ax[0].axhline(i_rule, color=C[3], lw=1.0, label=r"$P_c/(\bar k_cV_{dc})$")
    ax[0].set_xlabel("Active mass (kg)")
    ax[0].set_ylabel("$I_{pk}$ (A)")
    ax[0].set_ylim(0, i_req.max() * 1.35)
    ax[0].legend(fontsize=5, loc="lower right", ncol=2)
    ax[0].set_title("(a)", loc="left")
    vv = np.array([r["v_min"] for r in vrows])
    ii = np.array([r["i_required"] for r in vrows])
    vx = np.linspace(18, 42, 50)
    ax[1].plot(vx, p_c / (kc_mean * vx), color=C[3], lw=1.0, label=r"$P_c/(\bar k_cV)$")
    ax[1].plot(vv, ii, "o", color=C[0], ms=4, label="redesigned")
    ax[1].set_xlabel("Design voltage $V_{dc,min}$ (V)")
    ax[1].set_ylabel("$I_{pk}$ (A)")
    ax[1].legend(fontsize=5.5)
    ax[1].set_title("(b)", loc="left")
    fig.tight_layout(w_pad=0.6)
    save(fig, FIG, "fig_corner_rule")


def task_pm(dd):
    """First-order analytic SPM benchmark vs the SRM at equal envelope and boundary."""
    from srm_ebike.pm_benchmark import SPMBenchmark, SPMInputs, spm_lookup

    req, d, spec = dd.req, dd.design, dd.spec
    pm = SPMBenchmark(SPMInputs(d_outer=d.outer_diameter, length=d.stack_length, v_dc_min=spec.battery_voltage_min, base_speed_rpm=req.base_speed_rpm))
    pts = {
        "launch": (0.1 * req.corner_speed_rpm, req.peak_torque),
        "hill": (req.base_speed_rpm * 15 / spec.v_max_kmh, req.hill_torque),
        "rated": (req.base_speed_rpm, req.rated_torque),
        "cruise": (req.base_speed_rpm, req.cruise_torque),
    }
    dd_a = design_drive()
    rows = {}
    for name, (n, t) in pts.items():
        srm_f = max_torque(dd.drive, n) if name == "launch" else operating_point(dd.drive, n, t)
        srm_a = max_torque(dd_a.drive, n) if name == "launch" else operating_point(dd_a.drive, n, t)
        rows[name] = {"pm": pm.operating_point(n, t, spec.battery_voltage), "srm_fe": eff(srm_f), "srm_analytic": eff(srm_a)}
    speeds = np.linspace(0.08, 1.2, 10) * req.base_speed_rpm
    i_lim = 1.1 * req.peak_torque / pm.k_t
    lk_pm = spm_lookup(pm, speeds, np.linspace(0.25, 6.8, 14), spec.battery_voltage, i_lim)
    lk_srm, _, _ = build_lookup(dd.drive, req)
    routes = {}
    for name, r in ROUTES.items():
        a, b = simulate_route(spec, lk_pm, r), simulate_route(spec, lk_srm, r)
        routes[name] = {"pm_Wh_per_km": a["Wh_per_km"], "srm_fe_Wh_per_km": b["Wh_per_km"], "pm_deficit_s": a["torque_deficit_s"], "srm_deficit_s": b["torque_deficit_s"]}
    out = {
        "pm_design": pm.summary(),
        "pm_i_peak_launch_A": req.peak_torque / pm.k_t,
        "pm_i_limit_A": i_lim,
        "pm_drag_open_circuit_W": {str(n): pm.drag_open_circuit(n) for n in (1608.0, 1930.0)},
        "srm_mass_active_kg": d.mass_active,
        "points": rows,
        "routes": routes,
    }
    (OUT / "pm.json").write_text(json.dumps(out, indent=1, default=float))
    log("pm done")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="battery,ironloss,coupled,corner")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    dd = design_drive(fea_table=str(TABLE))
    log(f"FE design ready, I_max = {dd.drive.i_max:.2f} A")
    for name in args.only.split(","):
        {"battery": task_battery, "ironloss": task_ironloss, "coupled": task_coupled, "corner": task_corner, "pm": task_pm}[name](dd)


if __name__ == "__main__":
    main()
