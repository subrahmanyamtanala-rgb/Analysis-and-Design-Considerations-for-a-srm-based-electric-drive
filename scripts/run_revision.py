#!/usr/bin/env python3
"""Studies added in response to the review (needs results/fea from run_fea.py).

Writes results/revision/revision.json and paper/figures/fig_*.pdf.
Usage: python scripts/run_revision.py [--only validation,drive,tsf,pwm,thermal,gear,dclink,dwell]
"""

from __future__ import annotations

import argparse
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
from scipy.interpolate import PchipInterpolator  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
warnings.filterwarnings("ignore")

from make_paper_figures import C, COL, DBL, INK, INK2, LS, save  # noqa: E402  (also sets rcParams)

from srm_ebike import dclink  # noqa: E402
from srm_ebike.design_flow import converter_ratings, design_drive  # noqa: E402
from srm_ebike.drive_cycle import DEFAULT_ROUTE, EfficiencyLookup, simulate_route  # noqa: E402
from srm_ebike.losses import core_loss, mechanical_loss  # noqa: E402
from srm_ebike.performance import AngleGrid, efficiency_map, max_torque, operating_point, torque_speed_envelope  # noqa: E402
from srm_ebike.thermal import ThermalModel, ThermalParams  # noqa: E402
from srm_ebike.tsf import InverseTorque, tsf_operating_point  # noqa: E402
from srm_ebike.vehicle import EBikeSpec  # noqa: E402

FEA = ROOT / "results" / "fea"
OUT = ROOT / "results" / "revision"
FIG = ROOT / "paper" / "figures"
TABLE = FEA / "table.npz"
R: dict = {}


def log(msg, t0=[time.time()]):
    print(f"[{time.time() - t0[0]:7.1f}s] {msg}", flush=True)


def summary(res):
    eff = res.efficiency()
    return {
        "n_rpm": res.n_rpm,
        "T": res.torque_avg,
        "theta_on_deg": math.degrees(res.theta_on),
        "theta_off_deg": math.degrees(res.theta_off),
        "i_ref": res.i_ref,
        "i_rms": res.i_rms,
        "i_peak": res.i_peak,
        "ripple_pct": 100 * res.torque_ripple,
        "f_sw_kHz": res.switching_frequency / 1e3,
        "eta_motor": eff["eta_motor"],
        "eta_system": eff["eta_system"],
        "p_out": eff["p_out"],
        "p_in": eff["p_in"],
        "losses": res.losses(),
    }


# ---------------------------------------------------------------------------
def study_validation(dd_a, dd_f):
    data = np.load(TABLE)
    th, cur, psi, tq = data["theta_deg"], data["currents"], data["psi"], data["torque"]
    mm = dd_a.drive.model
    cur0 = np.concatenate([[0.0], cur])
    rows = []
    for i in (8.0, 15.0, 23.0, 31.0):
        p_fea = np.array([PchipInterpolator(cur0, np.concatenate([[0.0], psi[a]]))(i) for a in range(len(th))])
        t_fea = np.array([PchipInterpolator(cur0, np.concatenate([[0.0], tq[a]]))(i) for a in range(len(th))])
        p_an = np.array([mm.psi(i, math.radians(x)) for x in th])
        t_an = np.array([mm.torque(i, math.radians(x)) for x in th])
        t_fea_tab = np.array([dd_f.drive.model.torque(i, math.radians(x)) for x in th])
        mean_f, mean_a = np.trapezoid(t_fea, th) / th[-1], np.trapezoid(t_an, th) / th[-1]
        rows.append(
            {
                "i": i,
                "psi_rms_err_pct": float(100 * np.sqrt(np.mean((p_an - p_fea) ** 2)) / p_fea.max()),
                "psi_max_err_pct": float(100 * np.max(np.abs(p_an - p_fea)) / p_fea.max()),
                "psi_aligned_err_pct": float(100 * (p_an[-1] - p_fea[-1]) / p_fea[-1]),
                "psi_unaligned_err_pct": float(100 * (p_an[0] - p_fea[0]) / p_fea[0]),
                "T_mean_fea": float(mean_f),
                "T_mean_an": float(mean_a),
                "T_mean_err_pct": float(100 * (mean_a - mean_f) / mean_f),
                "T_peak_fea": float(t_fea.max()),
                "T_peak_an": float(t_an.max()),
                "T_peak_err_pct": float(100 * (t_an.max() - t_fea.max()) / t_fea.max()),
                "T_coenergy_vs_arkkio_pct": float(100 * (np.trapezoid(t_fea_tab, th) / th[-1] - mean_f) / mean_f),
            }
        )
    la_f, lu_f = psi[-1, 1] / cur[1], psi[0, 1] / cur[1]
    R["validation"] = {
        "rows": rows,
        "L_a_fea_mH": la_f * 1e3,
        "L_u_fea_mH": lu_f * 1e3,
        "L_a_an_mH": mm.l_a * 1e3,
        "L_u_an_mH": mm.l_u * 1e3,
        "L_a_err_pct": 100 * (mm.l_a - la_f) / la_f,
        "L_u_err_pct": 100 * (mm.l_u - lu_f) / lu_f,
        "ratio_fea": la_f / lu_f,
        "l_end_mH": float(data["l_end"]) * 1e3,
        "bmax_T": float(data["bmax"].max()),
        "mesh": json.loads((FEA / "mesh.json").read_text()),
    }
    mut = json.loads((FEA / "mutual.json").read_text())
    R["mutual"] = {
        "rows": mut,
        "max_psiA_change_pct": max(100 * abs(r["psiA_both_mWb"] - r["psiA_alone_mWb"]) / r["psiA_alone_mWb"] for r in mut),
        "max_T_superposition_err_pct": max(100 * abs(r["T_both"] - r["T_sum_alone"]) / max(abs(r["T_both"]), 0.3) for r in mut),
        "max_T_superposition_err_Nm": max(abs(r["T_both"] - r["T_sum_alone"]) for r in mut),
    }
    for name in ("airgap", "arcs"):
        f = FEA / f"{name}.json"
        if f.exists():
            R[name] = json.loads(f.read_text())

    # figure: psi-i and torque, FEA markers vs analytic lines
    fig, ax = plt.subplots(1, 2, figsize=(COL, 2.0))
    sel = [0, 4, 7, 10, len(th) - 1]
    ii = np.linspace(0, cur[-1], 80)
    for k, a in enumerate(sel):
        x = math.radians(th[a])
        ax[0].plot(ii, [mm.psi(i, x) * 1e3 for i in ii], color=C[k], ls=LS[k], lw=1)
        ax[0].plot(cur, psi[a] * 1e3, "o", ms=2.2, color=C[k], label=f"{th[a]:.1f}$^\\circ$")
    ax[0].set_xlabel("Current (A)")
    ax[0].set_ylabel("Flux linkage (mWb)")
    ax[0].legend(fontsize=5.5, loc="upper left", handlelength=1)
    ax[0].set_title("(a)", loc="left")
    for k, i in enumerate((8.0, 15.0, 23.0, 31.0)):
        t_fea = [PchipInterpolator(cur0, np.concatenate([[0.0], tq[a]]))(i) for a in range(len(th))]
        xx = np.linspace(0, th[-1], 100)
        ax[1].plot(xx, [mm.torque(i, math.radians(x)) for x in xx], color=C[k], ls=LS[k], lw=1)
        ax[1].plot(th, t_fea, "o", ms=2.2, color=C[k], label=f"{i:.0f} A")
    ax[1].set_xlabel("Rotor angle (deg)")
    ax[1].set_ylabel("Torque (N$\\cdot$m)")
    ax[1].legend(fontsize=5.5, loc="upper left", handlelength=1)
    ax[1].set_title("(b)", loc="left")
    fig.tight_layout(w_pad=0.5)
    save(fig, FIG, "fig_fea_validation")

    # figure: field plot at mid-overlap, 30 A
    from srm_ebike.fea import SRMFEA, element_b

    fea = SRMFEA(dd_a.design)
    theta = math.radians(12.0)
    res, a = fea.solve(theta, (30.0, 0.0, 0.0))
    p, tri, phys = fea.build_mesh(theta)
    bx, by = element_b(p, tri, a)
    import matplotlib.tri as mtri

    fig, axf = plt.subplots(figsize=(COL * 0.62, COL * 0.62))
    iron = np.isin(phys, (1, 2))
    x, y = p[0] * 1e3, p[1] * 1e3
    axf.tripcolor(mtri.Triangulation(x, y, tri.T, mask=iron), facecolors=np.zeros(tri.shape[1]), cmap="Greys", vmin=0, vmax=12, rasterized=True)
    tp = axf.tripcolor(mtri.Triangulation(x, y, tri.T, mask=~iron), facecolors=np.hypot(bx, by), cmap="viridis", vmin=0, vmax=2.2, rasterized=True)
    axf.tricontour(x, y, tri.T, a, levels=18, colors="white", linewidths=0.35)
    axf.set_aspect("equal")
    axf.axis("off")
    cb = fig.colorbar(tp, ax=axf, fraction=0.045, pad=0.02)
    cb.set_label("|B| (T)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    save(fig, FIG, "fig_fea_field")
    R["validation"]["field_nodes"] = int(p.shape[1])
    log("validation done")


def study_drive(dd_a, dd_f):
    out = {}
    for tag, dd in (("analytic", dd_a), ("fea", dd_f)):
        req, dr = dd.req, dd.drive
        launch = max_torque(dr, 0.1 * req.corner_speed_rpm)
        pts = {
            "launch": launch,
            "hill": operating_point(dr, req.base_speed_rpm * 15 / dd.spec.v_max_kmh, req.hill_torque),
            "rated": operating_point(dr, req.base_speed_rpm, req.rated_torque),
            "cruise": operating_point(dr, req.base_speed_rpm, req.cruise_torque),
        }
        out[tag] = {
            "i_required": dd.i_peak_required,
            "i_max": dr.i_max,
            "points": {k: summary(v) for k, v in pts.items()},
            "launch_margin_pct": 100 * (launch.torque_avg / req.peak_torque - 1),
        }
        if tag == "fea":
            speeds = np.linspace(0.08, 1.2, 10) * req.base_speed_rpm
            env = torque_speed_envelope(dr, speeds)
            t_env = np.array([e.torque_avg for e in env])
            sp, tq, eta_s, eta_m, _ = efficiency_map(dr, speeds, np.linspace(0.25, t_env.max(), 14), env)
            lookup = EfficiencyLookup(sp, t_env, tq, eta_s)
            cyc = simulate_route(dd.spec, lookup, DEFAULT_ROUTE)
            out[tag]["envelope"] = {"speed": speeds.tolist(), "torque": t_env.tolist()}
            out[tag]["eta_sys_max"] = float(np.nanmax(eta_s))
            out[tag]["eta_mot_max"] = float(np.nanmax(eta_m))
            out[tag]["cycle"] = {k: v for k, v in cyc.items() if not isinstance(v, np.ndarray)}
            out[tag]["converter"] = converter_ratings(dd, launch, pts["rated"])
            # energy-balance check of the FEA-table drive model
            errs = []
            for n, i in [(80, 30), (400, 25), (800, 20), (1200, 18), (1600, 15), (1900, 12)]:
                on, off = dr.default_angles(n, i, 0.8)
                r = dr.simulate(n, on, off, i)
                ls = r.losses()
                errs.append(abs(r.p_dc_circuit - (r.p_electromagnetic + ls["copper"] + ls["converter_conduction"])) / r.p_dc_circuit)
            out[tag]["energy_balance_max_err_pct"] = 100 * max(errs)
            R["_cache"] = {"sp": sp, "tq": tq, "eta_s": eta_s, "env": env, "cyc": cyc, "lookup": lookup, "pts": pts}
            # paper figures regenerated with the FEA-based model
            from make_paper_figures import fig_cycle, fig_losses, fig_map, fig_torque_speed, fig_waveforms

            fig_waveforms([(launch, "Current chopping"), (pts["rated"], "Rated, single pulse")], FIG)
            fig_torque_speed(dd, env, FIG)
            fig_map(sp, tq, eta_s, env, cyc, FIG)
            fig_losses({"launch": launch, "hill 6%": pts["hill"], "rated": pts["rated"], "cruise": pts["cruise"]}, FIG)
            fig_cycle(cyc, FIG)
            R["_pts_fea"] = pts
        log(f"drive {tag} done")
    R["drive"] = out


def study_tsf(dd_f):
    dr = dd_f.drive
    inv = InverseTorque(dr, dr.i_max)
    rows = []
    wave = {}
    for n, t in [(200, 1.5), (200, 3.0), (500, 3.0), (800, 1.5), (800, 3.0), (1200, 1.5), (1608, 1.48)]:
        base = operating_point(dr, n, t)
        row = {"n_rpm": n, "T": t, "CCC": summary(base) if base else None}
        for shape in ("linear", "cubic", "exponential"):
            r = tsf_operating_point(dr, inv, n, t, shape, theta_on_deg=(3.0, 4.5, 6.0, 7.5), theta_ov_deg=(1.5, 3.0, 4.5))
            row[shape] = {**summary(r), **r.tsf} if r else None
            if r is not None and n == 200 and t == 3.0 and shape == "cubic":
                wave["tsf"] = r
        if n == 200 and t == 3.0:
            wave["ccc"] = base
        rows.append(row)
        log(f"tsf {n} rpm {t} Nm done")
    R["tsf"] = rows
    fig, ax = plt.subplots(2, 1, figsize=(COL, 2.6), sharex=True)
    for k, (lab, r) in enumerate((("Current-reference (CCC)", wave["ccc"]), ("Cubic TSF", wave["tsf"]))):
        x = np.degrees(r.theta - r.theta[0])
        ax[0].plot(x, r.current, color=C[k], ls=LS[k], lw=1, label=lab)
        ax[1].plot(x, r.torque_total, color=C[k], ls=LS[k], lw=1.2, label=f"{lab}, ripple {100 * r.torque_ripple:.0f}%")
    ax[0].set_ylabel("Phase current (A)")
    ax[1].set_ylabel("Total torque (N$\\cdot$m)")
    ax[1].set_xlabel("Rotor angle from turn-on (deg)")
    ax[1].legend(fontsize=6, loc="lower right")
    ax[1].set_ylim(bottom=0)
    fig.align_ylabels(ax)
    fig.tight_layout(h_pad=0.3)
    save(fig, FIG, "fig_tsf")


def study_pwm(dd_f):
    dr = dd_f.drive
    req = dd_f.req
    rows = []
    cases = {
        "launch": (0.1 * req.corner_speed_rpm, None),
        "hill": (req.base_speed_rpm * 15 / dd_f.spec.v_max_kmh, req.hill_torque),
        "rated": (req.base_speed_rpm, req.rated_torque),
    }
    for name, (n, t) in cases.items():
        ref = max_torque(dr, n) if t is None else operating_point(dr, n, t)
        t_target = ref.torque_avg if t is None else t
        row = {"point": name, "hysteresis": summary(ref)}
        for f in (4e3, 8e3, 16e3, 20e3):
            def err(i):
                return dr.simulate(n, ref.theta_on, ref.theta_off, i, control="pwm", f_pwm=f).torque_avg - t_target

            try:
                i_sol = brentq(err, 0.5, 1.3 * dr.i_max, xtol=0.02) if err(1.3 * dr.i_max) > 0 else 1.3 * dr.i_max
            except ValueError:
                i_sol = ref.i_ref
            r = dr.simulate(n, ref.theta_on, ref.theta_off, i_sol, control="pwm", f_pwm=f)
            row[f"{f / 1e3:.0f}kHz"] = summary(r)
        rows.append(row)
        log(f"pwm {name} done")
    R["pwm"] = rows


def losses_for_thermal(res, d, iron):
    ls = res.losses()
    cl = core_loss(d, iron, res.n_rpm, res.b_pole_peak)
    return (ls["copper"], cl["stator_poles"] + cl["stator_yoke"], cl["rotor"], mechanical_loss(res.n_rpm))


def study_thermal(dd_f):
    d, dr = dd_f.design, dd_f.drive
    pts = R.get("_pts_fea")
    if pts is None:
        req = dd_f.req
        pts = {
            "launch": max_torque(dr, 0.1 * req.corner_speed_rpm),
            "hill": operating_point(dr, req.base_speed_rpm * 15 / dd_f.spec.v_max_kmh, req.hill_torque),
            "rated": operating_point(dr, req.base_speed_rpm, req.rated_torque),
            "cruise": operating_point(dr, req.base_speed_rpm, req.cruise_torque),
        }
    L = {k: losses_for_thermal(v, d, dr.iron) for k, v in pts.items()}

    def rated(t):
        return L["rated"]

    def climb(t):  # 20-min 6 % climb, then cruise
        return L["hill"] if t < 1200 else L["cruise"]

    def launches(t):  # stop-and-go: 5 s launch every 60 s, cruise in between
        return L["launch"] if (t % 60.0) < 5.0 else L["cruise"]

    scen = {"rated_continuous": (rated, 3600), "climb_6pct_20min": (climb, 2400), "repeated_launch": (launches, 2400)}
    rows = []
    traces = {}
    for h in (15.0, 25.0, 40.0):
        for amb in (25.0, 35.0, 45.0):
            tm = ThermalModel(d, ThermalParams(h_conv=h, t_amb=amb))
            for name, (fn, t_end) in scen.items():
                res = tm.simulate(t_end, fn, dt=0.5)
                rows.append({"scenario": name, "h": h, "t_amb": amb, "hotspot_max_C": float(res["hotspot"].max()), "winding_max_C": float(res["winding"].max()), "shell_max_C": float(res["shell"].max())})
                if h == 25.0 and amb == 35.0:
                    traces[name] = res
    tm = ThermalModel(d)
    R["thermal"] = {
        "rows": rows,
        "losses_W": {k: list(v) for k, v in L.items()},
        "conductances_W_per_K": {"winding_stator": tm.g_ws, "stator_shell": tm.g_sh, "shell_ambient_h25": tm.g_amb, "rotor_shell": tm.p.g_rotor},
        "capacitances_J_per_K": tm.c.tolist(),
    }
    fig, ax = plt.subplots(figsize=(COL, 1.9))
    labels = {"rated_continuous": "250 W continuous", "climb_6pct_20min": "6% climb 20 min, then cruise", "repeated_launch": "launch 5 s / 60 s"}
    for k, (name, res) in enumerate(traces.items()):
        ax.plot(res["t"] / 60, res["hotspot"], color=C[k], ls=LS[k], label=labels[name])
    ax.axhline(155, color=INK2, lw=0.6, ls=":")
    ax.text(1, 150, "class F (155 $^\\circ$C)", fontsize=6, color=INK2, va="top")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Winding hot spot ($^\\circ$C)")
    ax.legend(fontsize=6, loc="center right")
    save(fig, FIG, "fig_thermal")
    log("thermal done")


def study_gear(dd_f):
    cache = R.get("_cache")
    lookup = cache["lookup"] if cache else None
    if lookup is None:
        dr, req = dd_f.drive, dd_f.req
        speeds = np.linspace(0.08, 1.2, 8) * req.base_speed_rpm
        env = torque_speed_envelope(dr, speeds)
        t_env = np.array([e.torque_avg for e in env])
        sp, tq, eta_s, _, _ = efficiency_map(dr, speeds, np.linspace(0.25, t_env.max(), 8), env)
        lookup = EfficiencyLookup(sp, t_env, tq, eta_s)
    from srm_ebike.vehicle import derive_requirements

    rows = []
    for label, eta, nl in [("0.90", 0.90, 0.0), ("0.93", 0.93, 0.0), ("0.95", 0.95, 0.0), ("0.97", 0.97, 0.0), ("0.97 mesh + 3 W drag", 0.97, 3.0), ("0.95 mesh + 5 W drag", 0.95, 5.0)]:
        spec = EBikeSpec(gear_efficiency=eta, gear_no_load_loss=nl)
        req = derive_requirements(spec)
        cyc = simulate_route(spec, lookup, DEFAULT_ROUTE)
        rows.append(
            {
                "gear": label,
                "launch_torque": req.peak_torque,
                "rated_torque_motor": req.rated_torque,
                "launch_within_envelope": bool(req.peak_torque <= lookup.t_max(0.1 * req.corner_speed_rpm)),
                "Wh_per_km": cyc["Wh_per_km"],
                "range_km": cyc["range_km"],
                "cycle_eff_motor_drive": cyc["cycle_efficiency"],
                "deficit_s": cyc["torque_deficit_s"],
            }
        )
    R["gear"] = rows
    log("gear done")


def study_dclink(dd_f):
    pts = R.get("_pts_fea")
    dr = dd_f.drive
    if pts is None:
        req = dd_f.req
        pts = {"launch": max_torque(dr, 0.1 * req.corner_speed_rpm), "rated": operating_point(dr, req.base_speed_rpm, req.rated_torque)}
    out = {}
    for name in ("rated", "launch"):
        r = pts[name]
        period = dr.model.tau_r / r.omega
        idc = r.i_dc_total
        b = {f"{int(round(100 * dv))}pct": dclink.charge_bound(idc, period, dv * dd_f.spec.battery_voltage) for dv in (0.02, 0.05)}
        share = [dclink.sharing(idc, period, c, esr=0.02, r_b=0.15, l_b=1.0e-6) for c in (470e-6, 1000e-6, 2200e-6, 3800e-6)]
        out[name] = {
            "i_dc_mean": float(idc.mean()),
            "i_dc_ac_rms": float(np.sqrt(np.mean((idc - idc.mean()) ** 2))),
            "f_stroke_Hz": dd_f.design.phases / period,
            "charge_bound": b,
            "sharing": share,
        }
    R["dclink"] = {"assumptions": {"ESR_ohm": 0.02, "R_battery_ohm": 0.15, "L_cable_H": 1.0e-6}, **out}
    log("dclink done")


def study_dwell(dd_f):
    dr, req = dd_f.drive, dd_f.req
    rows = []
    for md in (0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80):
        g = AngleGrid(max_dwell=md)
        r_max = max_torque(dr, req.base_speed_rpm, grid=g)
        r_op = operating_point(dr, req.base_speed_rpm, req.rated_torque, grid=AngleGrid(n_on=3, off_fractions=(0.55, 0.7, 0.85, 0.95), max_dwell=md))
        rows.append(
            {
                "max_dwell": md,
                "T_max_base": r_max.torque_avg,
                "continuous_at_Tmax": r_max.continuous,
                "eta_sys_rated": r_op.efficiency()["eta_system"] if r_op else None,
                "ripple_rated_pct": 100 * r_op.torque_ripple if r_op else None,
            }
        )
    R["dwell"] = rows
    log("dwell done")


def study_grid(dd_f):
    """Convergence of the angle search: default grid vs a dense grid."""
    dr, req = dd_f.drive, dd_f.req
    dense = AngleGrid(
        n_on=9,
        off_fractions=tuple(np.round(np.arange(0.40, 0.99, 0.04), 3)),
        dwell_fractions=tuple(np.round(np.arange(0.25, 0.551, 0.025), 3)),
    )
    rows = []
    for name, n, t in [("hill", req.base_speed_rpm * 15 / dd_f.spec.v_max_kmh, req.hill_torque), ("rated", req.base_speed_rpm, req.rated_torque), ("cruise", req.base_speed_rpm, req.cruise_torque)]:
        a = operating_point(dr, n, t)
        b = operating_point(dr, n, t, grid=dense)
        rows.append({"point": name, "eta_sys_default": a.efficiency()["eta_system"], "eta_sys_dense": b.efficiency()["eta_system"], "p_in_default": a.efficiency()["p_in"], "p_in_dense": b.efficiency()["p_in"]})
    a = max_torque(dr, 0.1 * req.corner_speed_rpm)
    b = max_torque(dr, 0.1 * req.corner_speed_rpm, grid=dense)
    rows.append({"point": "launch_Tmax", "T_default": a.torque_avg, "T_dense": b.torque_avg})
    a = max_torque(dr, req.base_speed_rpm)
    b = max_torque(dr, req.base_speed_rpm, grid=dense)
    rows.append({"point": "base_Tmax", "T_default": a.torque_avg, "T_dense": b.torque_avg})
    R["grid"] = rows
    log("grid done")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="validation,drive,tsf,pwm,thermal,gear,dclink,dwell,grid")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    dd_a = design_drive()
    dd_f = design_drive(fea_table=str(TABLE))
    log(f"designs ready: I_req analytic {dd_a.i_peak_required:.1f} A, FEA {dd_f.i_peak_required:.1f} A")
    todo = args.only.split(",")
    if "validation" in todo:
        study_validation(dd_a, dd_f)
    if "drive" in todo:
        study_drive(dd_a, dd_f)
    for name, fn in (("tsf", study_tsf), ("pwm", study_pwm), ("thermal", study_thermal), ("gear", study_gear), ("dclink", study_dclink), ("dwell", study_dwell), ("grid", study_grid)):
        if name in todo:
            fn(dd_f)
    path = OUT / f"revision_{'_'.join(todo) if len(todo) < 9 else 'all'}.json"
    path.write_text(json.dumps({k: v for k, v in R.items() if not k.startswith("_")}, indent=1, default=float))
    log(f"wrote {path}")


if __name__ == "__main__":
    main()
