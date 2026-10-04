#!/usr/bin/env python3
"""Column-width vector figures and key numbers for the IEEE paper (paper/).

Usage:  python scripts/make_paper_figures.py [--out paper/figures]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from run_analysis import trade_study  # noqa: E402

from srm_ebike.design_flow import (  # noqa: E402
    converter_ratings,
    design_drive,
    rated_and_peak_points,
    thermal_estimate,
)
from srm_ebike.drive_cycle import DEFAULT_ROUTE, EfficiencyLookup, simulate_route  # noqa: E402
from srm_ebike.performance import efficiency_map, operating_point, torque_speed_envelope  # noqa: E402
from srm_ebike.vehicle import motor_torque_demand, road_load_force  # noqa: E402

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
LS = ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d3"
COL = 3.45  # IEEE single-column width [in]
DBL = 7.16  # IEEE double-column width [in]

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["STIXGeneral"],
        "mathtext.fontset": "stix",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.edgecolor": INK2,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.3,
        "legend.frameon": False,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
    }
)


def save(fig, out: Path, name: str) -> None:
    fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


def fig_geometry(dd, out):
    d = dd.design
    fig, ax = plt.subplots(1, 2, figsize=(COL, 1.75), gridspec_kw={"width_ratios": [1, 1.15]})
    a = ax[0]
    th = np.linspace(0, 2 * np.pi, 721)
    r_o, r_b = d.outer_diameter / 2, d.bore_diameter / 2
    r_yi, r_r = r_b + d.stator_pole_height, r_b - d.air_gap
    r_ry, r_sh = r_r - d.rotor_pole_height, d.shaft_diameter / 2

    def ring(r1, r2, color):
        x = np.concatenate([r2 * np.cos(th), r1 * np.cos(th[::-1])])
        y = np.concatenate([r2 * np.sin(th), r1 * np.sin(th[::-1])])
        a.fill(x * 1e3, y * 1e3, color=color, lw=0)

    def poles(n, r_in, r_out, w, color):
        for k in range(n):
            ang = 2 * np.pi * k / n
            u, p = np.array([np.cos(ang), np.sin(ang)]), np.array([-np.sin(ang), np.cos(ang)])
            pts = np.array([r_in * u - w / 2 * p, r_out * u - w / 2 * p, r_out * u + w / 2 * p, r_in * u + w / 2 * p]) * 1e3
            a.fill(pts[:, 0], pts[:, 1], color=color, lw=0)

    ring(r_yi, r_o, "#7d7c77")
    poles(d.ns, r_b, r_yi + 5e-4, d.stator_pole_width, "#7d7c77")
    ring(r_sh, r_ry + 5e-4, "#b5b4ae")
    poles(d.nr, r_ry, r_r, d.rotor_pole_width, "#b5b4ae")
    for k in range(0, d.ns, d.phases):
        ang0 = 2 * np.pi * k / d.ns
        for s in (-1, 1):
            ang = ang0 + s * (d.beta_s / 2 + 0.6 * (np.pi / d.ns - d.beta_s / 2))
            rr = 0.5 * (r_b + r_yi) * 1e3
            a.plot(rr * np.cos(ang), rr * np.sin(ang), "s", ms=3.2, color=C[1])
    a.set_aspect("equal")
    a.axis("off")
    a.set_title("(a)", loc="left")

    b = ax[1]
    eps, tau = math.degrees(d.stroke_angle), math.degrees(d.rotor_pole_pitch)
    tri = np.array([[eps, eps], [eps, tau - eps], [tau / 2, tau / 2]])
    b.fill(tri[:, 0], tri[:, 1], color=C[0], alpha=0.18, lw=0)
    b.plot(*np.vstack([tri, tri[:1]]).T, color=C[0], lw=1)
    bs, br = math.degrees(d.beta_s), math.degrees(d.beta_r)
    b.plot(bs, br, "o", ms=5, color=C[1], mec="white", mew=1)
    b.annotate(rf"$({bs:.0f}^\circ,{br:.0f}^\circ)$", (bs, br), xytext=(4, 2), textcoords="offset points")
    b.set_xlabel(r"$\beta_s$ (deg)")
    b.set_ylabel(r"$\beta_r$ (deg)")
    b.set_title("(b)", loc="left")
    save(fig, out, "fig_geometry")


def fig_static(dd, out):
    mm = dd.drive.model
    fig, ax = plt.subplots(1, 2, figsize=(COL, 1.9))
    th_ov = mm.overlap_start - mm.fringe / 2
    fr = [0.0, 0.25, 0.5, 0.75, 1.0]
    ths = [th_ov + f * (mm.tau_r / 2 - th_ov) for f in fr]
    ii, psi = mm.curves(1.2 * dd.drive.i_max, np.array(ths))
    for k, t in enumerate(ths):
        ax[0].plot(ii, psi[k] * 1e3, color=C[k], ls=LS[k], label=rf"{math.degrees(t):.1f}$^\circ$")
    ax[0].set_xlabel("Current (A)")
    ax[0].set_ylabel("Flux linkage (mWb)")
    ax[0].legend(fontsize=5.5, handlelength=1.8, loc="upper left")
    th = np.linspace(0, mm.tau_r, 300)
    for k, f in enumerate([0.25, 0.5, 0.75, 1.0]):
        i = f * dd.drive.i_max
        ax[1].plot(np.degrees(th), [mm.torque(i, t) for t in th], color=C[k], ls=LS[k], label=f"{i:.0f} A")
    ax[1].set_xlabel(r"Rotor angle $\theta$ (deg)")
    ax[1].set_ylabel("Torque (N$\\cdot$m)")
    ax[1].legend(fontsize=5.5, handlelength=1.8, loc="upper right")
    ax[0].set_title("(a)", loc="left")
    ax[1].set_title("(b)", loc="left")
    fig.tight_layout(w_pad=0.6)
    save(fig, out, "fig_static")


def fig_waveforms(cases, out):
    fig, ax = plt.subplots(3, len(cases), figsize=(DBL, 3.4), sharex="col")
    for c, (res, lab) in enumerate(cases):
        th = np.degrees(res.theta)
        ax[0, c].plot(th, res.current, color=C[0])
        ax[0, c].axhline(res.i_ref, color=INK2, lw=0.6, ls="--")
        ax[1, c].plot(th, res.flux * 1e3, color=C[2])
        ax[2, c].plot(th, res.torque_phase, color=C[3], lw=0.9, ls="--", label="phase")
        ax[2, c].plot(th, res.torque_total, color=C[1], label="total")
        ax[2, c].axhline(res.torque_avg, color=INK2, lw=0.6, ls=":")
        for r in range(3):
            ax[r, c].axvline(math.degrees(res.theta_off), color=INK2, lw=0.6, ls=":")
        ax[0, c].set_title(
            f"{lab}: {res.n_rpm:.0f} r/min, "
            rf"$\theta_{{on}}$={math.degrees(res.theta_on):.1f}$^\circ$, "
            rf"$\theta_{{off}}$={math.degrees(res.theta_off):.1f}$^\circ$"
        )
        ax[2, c].set_xlabel("Rotor angle (deg)")
    ax[0, 0].set_ylabel("Current (A)")
    ax[1, 0].set_ylabel("Flux (mWb)")
    ax[2, 0].set_ylabel("Torque (N$\\cdot$m)")
    ax[2, 1].legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.12))
    fig.tight_layout(h_pad=0.3)
    save(fig, out, "fig_waveforms")


def fig_torque_speed(dd, env, out):
    spec, req = dd.spec, dd.req
    n = np.array([e.n_rpm for e in env])
    t = np.array([e.torque_avg for e in env])
    fig, ax = plt.subplots(figsize=(COL, 2.2))
    ax.plot(n, t, color=C[0], lw=1.8, label="SRM envelope ($I_{max}$)")
    v = np.linspace(1, spec.v_max_kmh, 60) / 3.6
    for k, g in enumerate([0.0, 0.03, 0.06, 0.08]):
        tq = [motor_torque_demand(spec, road_load_force(spec, x, g)) for x in v]
        ax.plot([spec.motor_speed_rpm(x) for x in v], tq, color=C[k + 1], ls=LS[k + 1], lw=1, label=f"steady, {g * 100:.0f}% grade")
    pts = [
        ("launch", 0.1 * req.corner_speed_rpm, req.peak_torque),
        ("hill", req.base_speed_rpm * 15 / spec.v_max_kmh, req.hill_torque),
        ("rated", req.base_speed_rpm, req.rated_torque),
    ]
    for lab, x, y in pts:
        ax.plot(x, y, "o", ms=4, color=INK, mec="white", mew=0.8)
        ax.annotate(lab, (x, y), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax2 = ax.secondary_xaxis("top", functions=(lambda r: r / spec.motor_speed_rpm(1 / 3.6), lambda k: k * spec.motor_speed_rpm(1 / 3.6)))
    ax2.set_xlabel("Road speed (km/h)")
    ax.set_xlabel("Motor speed (r/min)")
    ax.set_ylabel("Motor torque (N$\\cdot$m)")
    ax.legend(loc="upper right", fontsize=6)
    save(fig, out, "fig_torque_speed")
    return n, t


def fig_map(sp, tq, eta, env, cyc, out):
    fig, ax = plt.subplots(figsize=(COL, 2.3))
    lev = [0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8]
    cs = ax.contourf(sp, tq, eta.T, levels=lev, cmap="Blues", extend="both")
    cl = ax.contour(sp, tq, eta.T, levels=lev, colors="white", linewidths=0.5)
    ax.clabel(cl, fmt=lambda v: f"{v * 100:.0f}%", fontsize=6, colors=INK)
    ax.plot([e.n_rpm for e in env], [e.torque_avg for e in env], color=INK, lw=1.2, label="envelope")
    sel = cyc["t_motor"] > 0
    ax.plot(cyc["n_motor"][sel], cyc["t_motor"][sel], ".", ms=2, color=C[1], alpha=0.6, label="route points")
    cb = fig.colorbar(cs, ax=ax, pad=0.02)
    cb.set_label("System efficiency")
    cb.ax.tick_params(labelsize=6)
    ax.set_xlabel("Motor speed (r/min)")
    ax.set_ylabel("Torque (N$\\cdot$m)")
    ax.legend(loc="upper right", fontsize=6)
    ax.grid(False)
    save(fig, out, "fig_efficiency_map")


def fig_losses(points, out):
    cats = [("copper", "copper"), ("core", "core"), ("converter_conduction", "conv. cond."), ("converter_switching", "conv. sw."), ("mechanical", "mech.")]
    labels = list(points)
    fig, ax = plt.subplots(figsize=(COL, 1.7))
    left = np.zeros(len(labels))
    hatches = ["", "////", "....", "xxxx", "\\\\\\\\"]
    for k, (c, lab) in enumerate(cats):
        vals = np.array([points[l].losses()[c] for l in labels])
        ax.barh(labels, vals, left=left, color=C[k], edgecolor="white", linewidth=0.8, hatch=hatches[k], label=lab, height=0.6)
        left += vals
    for y, tot in enumerate(left):
        ax.text(tot + 3, y, f"{tot:.0f} W", va="center", fontsize=7)
    ax.set_xlim(0, left.max() * 1.18)
    ax.set_xlabel("Loss (W)")
    ax.grid(axis="y", visible=False)
    ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.45, -0.32), fontsize=6, handlelength=1.2, columnspacing=0.8)
    save(fig, out, "fig_losses")


def fig_cycle(cyc, out):
    t = cyc["t"] / 60
    fig, ax = plt.subplots(3, 1, figsize=(COL, 3.0), sharex=True)
    ax[0].plot(t, cyc["v"] * 3.6, color=C[0], lw=1)
    ax[0].set_ylabel("Speed\n(km/h)")
    ax[1].plot(t, cyc["t_motor"], color=C[1], lw=0.9)
    ax[1].set_ylabel("Torque\n(N$\\cdot$m)")
    ax[2].plot(t, cyc["p_batt"], color=C[2], lw=0.9, label="battery")
    ax[2].plot(t, cyc["p_rider"], color=C[6], lw=0.9, ls="--", label="rider")
    ax[2].set_ylabel("Power\n(W)")
    ax[2].set_xlabel("Time (min)")
    ax[2].legend(loc="upper right", ncol=2)
    fig.align_ylabels(ax)
    fig.tight_layout(h_pad=0.2)
    save(fig, out, "fig_cycle")


def fig_trade(rows, out):
    fig, ax = plt.subplots(figsize=(COL, 1.9))
    marks = ["o", "s", "^"]
    for k, g in enumerate(sorted({r["gear_ratio"] for r in rows})):
        rr = [r for r in rows if r["gear_ratio"] == g]
        ax.plot([r["active_mass_kg"] for r in rr], [r["eta_system_rated_pct"] for r in rr], marker=marks[k], ms=4, color=C[k], ls=LS[k], label=f"$G$ = {g}")
        for r in rr:
            ax.annotate(f"{r['A_peak_kA_per_m']:.0f}", (r["active_mass_kg"], r["eta_system_rated_pct"]), xytext=(3, -8), textcoords="offset points", fontsize=6, color=INK2)
    ax.set_xlabel("Active mass (kg)")
    ax.set_ylabel(r"$\eta_{sys}$ at rated (%)")
    ax.legend(loc="lower right")
    save(fig, out, "fig_trade")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "paper" / "figures"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    dd = design_drive()
    req, drive = dd.req, dd.drive
    rated, peak = rated_and_peak_points(dd)
    cruise = operating_point(drive, req.base_speed_rpm, req.cruise_torque)
    hill = operating_point(drive, req.base_speed_rpm * 15 / dd.spec.v_max_kmh, req.hill_torque)

    speeds = np.linspace(0.08, 1.2, 10) * req.base_speed_rpm
    env = torque_speed_envelope(drive, speeds)
    t_env = np.array([e.torque_avg for e in env])
    sp, tq, eta_s, eta_m, _ = efficiency_map(drive, speeds, np.linspace(0.25, t_env.max(), 14), env)
    cyc = simulate_route(dd.spec, EfficiencyLookup(sp, t_env, tq, eta_s), DEFAULT_ROUTE)
    trades = trade_study(False)
    print(f"[{time.time() - t0:5.1f}s] computations done")

    fig_geometry(dd, out)
    fig_static(dd, out)
    fig_waveforms([(peak, "Current chopping"), (rated, "Rated, single pulse")], out)
    fig_torque_speed(dd, env, out)
    fig_map(sp, tq, eta_s, env, cyc, out)
    fig_losses({"launch": peak, "hill 6%": hill, "rated": rated, "cruise": cruise}, out)
    fig_cycle(cyc, out)
    fig_trade(trades, out)

    k_c = [
        dd.spec.peak_power / (dd.spec.battery_voltage_min * r["I_peak_A"] / 1.1) for r in trades
    ]
    numbers = {
        "design": dd.design.summary(),
        "i_peak_required_A": dd.i_peak_required,
        "i_max_A": drive.i_max,
        "points": {k: v.summary() for k, v in {"launch": peak, "hill": hill, "rated": rated, "cruise": cruise}.items()},
        "losses": {k: v.losses() for k, v in {"launch": peak, "hill": hill, "rated": rated, "cruise": cruise}.items()},
        "envelope": {"speed_rpm": speeds.tolist(), "torque": t_env.tolist()},
        "eta_system_max": float(np.nanmax(eta_s)),
        "eta_motor_max": float(np.nanmax(eta_m)),
        "converter": converter_ratings(dd, peak, rated),
        "thermal": thermal_estimate(dd, rated),
        "cycle": {k: v for k, v in cyc.items() if not isinstance(v, np.ndarray)},
        "trade": trades,
        "k_c_range": [min(k_c), max(k_c)],
        "rotor_inertia": dd.design.rotor_inertia,
    }
    (out / "numbers.json").write_text(json.dumps(numbers, indent=2, default=float))
    print(f"[{time.time() - t0:5.1f}s] figures written to {out}")


if __name__ == "__main__":
    main()
