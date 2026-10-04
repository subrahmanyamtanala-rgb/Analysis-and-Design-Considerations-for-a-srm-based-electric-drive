#!/usr/bin/env python3
"""Run the complete SRM e-bike drive analysis and write figures + tables.

Usage:  python scripts/run_analysis.py [--quick] [--out results]
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from srm_ebike.design_flow import (  # noqa: E402
    converter_ratings,
    design_drive,
    rated_and_peak_points,
    thermal_estimate,
)
from srm_ebike.drive_cycle import DEFAULT_ROUTE, EfficiencyLookup, simulate_route  # noqa: E402
from srm_ebike.geometry import compare_pole_configs  # noqa: E402
from srm_ebike.performance import (  # noqa: E402
    efficiency_map,
    operating_point,
    torque_speed_envelope,
)
from srm_ebike.vehicle import EBikeSpec, motor_torque_demand, road_load_force  # noqa: E402

# categorical palette (fixed order) and recessive ink, see docs/figures note
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d3"

plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "font.size": 9,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 2.0,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "legend.frameon": False,
    }
)


def save(fig, out: Path, name: str) -> str:
    fig.tight_layout()
    path = out / name
    fig.savefig(path)
    plt.close(fig)
    return name


def md_table(rows: list[dict]) -> str:
    keys = list(rows[0].keys())
    lines = ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
    for r in rows:
        lines.append("| " + " | ".join(str(r[k]) for k in keys) + " |")
    return "\n".join(lines)


def kv_table(d: dict, head=("Quantity", "Value")) -> str:
    lines = [f"| {head[0]} | {head[1]} |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in d.items()]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_cross_section(dd, out):
    d = dd.design
    fig, ax = plt.subplots(figsize=(4.6, 4.6))
    r_o, r_b = d.outer_diameter / 2, d.bore_diameter / 2
    r_yi = r_b + d.stator_pole_height
    r_r = r_b - d.air_gap
    r_ry = r_r - d.rotor_pole_height
    r_sh = d.shaft_diameter / 2
    th = np.linspace(0, 2 * np.pi, 721)

    def ring(r1, r2, color, alpha=1.0):
        x = np.concatenate([r2 * np.cos(th), r1 * np.cos(th[::-1])])
        y = np.concatenate([r2 * np.sin(th), r1 * np.sin(th[::-1])])
        ax.fill(x * 1e3, y * 1e3, color=color, alpha=alpha, lw=0)

    def poles(n, r_in, r_out, width, color, offset=0.0):
        for k in range(n):
            a = offset + 2 * np.pi * k / n
            u = np.array([np.cos(a), np.sin(a)])
            p = np.array([-np.sin(a), np.cos(a)])
            pts = [r_in * u - width / 2 * p, r_out * u - width / 2 * p, r_out * u + width / 2 * p, r_in * u + width / 2 * p]
            pts = np.array(pts) * 1e3
            ax.fill(pts[:, 0], pts[:, 1], color=color, lw=0)

    steel_s, steel_r = "#8a8984", "#b5b4ae"
    ring(r_yi, r_o, steel_s)
    poles(d.ns, r_b, r_yi + 0.5e-3, d.stator_pole_width, steel_s)
    ring(r_sh, r_ry + 0.5e-3, steel_r)
    poles(d.nr, r_ry, r_r, d.rotor_pole_width, steel_r)
    ring(0, r_sh, "#ffffff")
    # phase-A coils on its stator poles
    for k in range(d.ns):
        if k % d.phases:
            continue
        a = 2 * np.pi * k / d.ns
        for side in (-1, 1):
            ang = a + side * (d.beta_s / 2 + 0.6 * (np.pi / d.ns - d.beta_s / 2))
            rr = 0.5 * (r_b + r_yi)
            ax.plot(rr * np.cos(ang) * 1e3, rr * np.sin(ang) * 1e3, "s", ms=7, color=C[1])
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel("mm")
    ax.set_title(f"{d.ns}/{d.nr} SRM lamination (phase-A coils marked)")
    return save(fig, out, "fig_cross_section.png")


def fig_feasible_triangle(dd, out):
    d = dd.design
    eps = math.degrees(d.stroke_angle)
    tau = math.degrees(d.rotor_pole_pitch)
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    # region: beta_s >= eps, beta_r >= beta_s, beta_s + beta_r <= tau
    tri = np.array([[eps, eps], [eps, tau - eps], [tau / 2, tau / 2]])
    ax.fill(tri[:, 0], tri[:, 1], color=C[0], alpha=0.18, lw=0)
    ax.plot(*np.vstack([tri, tri[:1]]).T, color=C[0], lw=1.5)
    bs, br = math.degrees(d.beta_s), math.degrees(d.beta_r)
    ax.plot(bs, br, "o", ms=8, color=C[1], mec="white", mew=2)
    ax.annotate(f"selected ({bs:.0f}°, {br:.0f}°)", (bs, br), xytext=(8, 6), textcoords="offset points", color=INK)
    ax.set_xlabel("stator pole arc βs [deg]")
    ax.set_ylabel("rotor pole arc βr [deg]")
    ax.set_title("Lawrenson feasible region for pole arcs")
    return save(fig, out, "fig_feasible_triangle.png")


def fig_magnetics(dd, out):
    mm = dd.drive.model
    i_max = 1.2 * dd.drive.i_max
    fracs = [0.0, 0.25, 0.5, 0.75, 1.0]
    theta_ov = mm.overlap_start - mm.fringe / 2
    thetas = [theta_ov + f * (mm.tau_r / 2 - theta_ov) for f in fracs]
    ii, psi = mm.curves(i_max, np.array(thetas))
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for k, th in enumerate(thetas):
        ax[0].plot(ii, psi[k] * 1e3, color=C[k], label=f"θ = {math.degrees(th):.1f}°")
    ax[0].set_xlabel("phase current [A]")
    ax[0].set_ylabel("flux linkage [mWb]")
    ax[0].set_title("Magnetisation curves ψ(i, θ)")
    ax[0].legend(fontsize=8)
    th = np.linspace(0, mm.tau_r, 400)
    currents = [0.25, 0.5, 0.75, 1.0]
    for k, f in enumerate(currents):
        i = f * dd.drive.i_max
        ax[1].plot(np.degrees(th), [mm.torque(i, t) for t in th], color=C[k], label=f"{i:.0f} A")
    ax[1].set_xlabel("rotor angle θ [deg] (0 = unaligned)")
    ax[1].set_ylabel("phase torque [N m]")
    ax[1].set_title("Static torque per phase")
    ax[1].legend(fontsize=8)
    return save(fig, out, "fig_magnetics.png")


def fig_waveforms(res, out, name, title):
    d = res.drive.design
    th = np.degrees(res.theta)
    fig, ax = plt.subplots(3, 1, figsize=(6.4, 6.2), sharex=True)
    ax[0].plot(th, res.current, color=C[0])
    ax[0].axhline(res.i_ref, color=INK2, lw=0.8, ls="--")
    ax[0].set_ylabel("phase current [A]")
    ax[1].plot(th, res.flux * 1e3, color=C[2])
    ax[1].set_ylabel("flux linkage [mWb]")
    ax[2].plot(th, res.torque_phase, color=C[3], lw=1.2, label="one phase")
    ax[2].plot(th, res.torque_total, color=C[1], label=f"total ({d.phases} phases)")
    ax[2].axhline(res.torque_avg, color=INK2, lw=0.8, ls="--")
    ax[2].set_ylabel("torque [N m]")
    ax[2].set_xlabel("rotor angle [deg]")
    ax[2].legend(fontsize=8, loc="upper right")
    for a in ax:
        a.axvline(math.degrees(res.theta_off), color=INK2, lw=0.8, ls=":")
    ax[0].set_title(
        f"{title}: {res.n_rpm:.0f} rpm, θon={math.degrees(res.theta_on):.1f}°, "
        f"θoff={math.degrees(res.theta_off):.1f}°, T={res.torque_avg:.2f} N m"
    )
    return save(fig, out, name)


def fig_torque_speed(dd, env, out):
    spec, req = dd.spec, dd.req
    n = np.array([e.n_rpm for e in env])
    t = np.array([e.torque_avg for e in env])
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.6))
    ax[0].plot(n, t, color=C[0], label="SRM envelope (I = I_max)")
    v = np.linspace(1, spec.v_max_kmh, 60) / 3.6
    for k, g in enumerate([0.0, 0.03, 0.06, 0.08]):
        tq = [motor_torque_demand(spec, road_load_force(spec, x, g)) for x in v]
        ax[0].plot([spec.motor_speed_rpm(x) for x in v], tq, color=C[k + 1], lw=1.2, label=f"steady, {g * 100:.0f} % grade")
    pts = {
        "launch 8 %": (0.1 * req.corner_speed_rpm, req.peak_torque),
        "hill 6 % @15 km/h": (req.base_speed_rpm * 15 / spec.v_max_kmh, req.hill_torque),
        "rated": (req.base_speed_rpm, req.rated_torque),
    }
    for k, (lab, (x, y)) in enumerate(pts.items()):
        ax[0].plot(x, y, "o", ms=8, color=INK, mec="white", mew=2)
        ax[0].annotate(lab, (x, y), xytext=(6, 4), textcoords="offset points", fontsize=8)
    ax[0].set_xlabel("motor speed [rpm]")
    ax[0].set_ylabel("motor torque [N m]")
    ax[0].set_title("Torque capability vs road-load demand")
    ax[0].legend(fontsize=7.5)
    ax[1].plot(n, t * n * 2 * np.pi / 60, color=C[0], label="max shaft-side power")
    ax[1].axhline(spec.rated_power, color=C[1], lw=1.2, ls="--", label="rated 250 W")
    ax[1].axhline(spec.peak_power, color=C[2], lw=1.2, ls="--", label="peak target 500 W")
    ax[1].axvline(req.base_speed_rpm, color=INK2, lw=0.8, ls=":")
    ax[1].set_xlabel("motor speed [rpm]")
    ax[1].set_ylabel("power [W]")
    ax[1].set_title("Power capability")
    ax[1].legend(fontsize=8)
    return save(fig, out, "fig_torque_speed.png")


def fig_eff_map(speeds, torques, eta, env, cycle, out):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    levels = [0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85]
    cs = ax.contourf(speeds, torques, eta.T, levels=levels, cmap="Blues", extend="both")
    cl = ax.contour(speeds, torques, eta.T, levels=levels, colors="white", linewidths=0.6)
    ax.clabel(cl, fmt=lambda v: f"{v * 100:.0f}%", fontsize=7, colors=INK)
    ax.plot([e.n_rpm for e in env], [e.torque_avg for e in env], color=INK, lw=1.5, label="torque envelope")
    if cycle is not None:
        sel = cycle["t_motor"] > 0
        ax.plot(cycle["n_motor"][sel], cycle["t_motor"][sel], ".", ms=3, color=C[1], alpha=0.5, label="drive-cycle points")
    fig.colorbar(cs, ax=ax, label="system efficiency (battery to shaft)")
    ax.set_xlabel("motor speed [rpm]")
    ax.set_ylabel("torque [N m]")
    ax.set_title("Drive-system efficiency map")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(False)
    return save(fig, out, "fig_efficiency_map.png")


def fig_losses(points: dict, out):
    cats = ["copper", "core", "converter_conduction", "converter_switching", "mechanical"]
    labels = list(points)
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    left = np.zeros(len(labels))
    for k, c in enumerate(cats):
        vals = np.array([points[l].losses()[c] for l in labels])
        ax.barh(labels, vals, left=left, color=C[k], edgecolor="white", linewidth=2, label=c.replace("_", " "), height=0.6)
        left += vals
    for y, tot in enumerate(left):
        ax.text(tot + 2, y, f"{tot:.0f} W", va="center", color=INK, fontsize=8)
    ax.set_xlabel("loss [W]")
    ax.set_title("Loss breakdown at key operating points")
    ax.legend(fontsize=7.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    ax.set_xlim(0, left.max() * 1.15)
    ax.grid(axis="y", visible=False)
    return save(fig, out, "fig_loss_breakdown.png")


def fig_cycle(cyc, out):
    t = cyc["t"] / 60
    fig, ax = plt.subplots(3, 1, figsize=(7.2, 6.0), sharex=True)
    ax[0].plot(t, cyc["v"] * 3.6, color=C[0], lw=1.2)
    ax[0].set_ylabel("speed [km/h]")
    ax[0].set_title(
        f"Route simulation: {cyc['distance_km']:.1f} km, {cyc['Wh_per_km']:.1f} Wh/km, "
        f"est. range {cyc['range_km']:.0f} km"
    )
    ax[1].plot(t, cyc["t_motor"], color=C[1], lw=1.0)
    ax[1].set_ylabel("motor torque [N m]")
    ax[2].plot(t, cyc["p_batt"], color=C[2], lw=1.0, label="battery (motor)")
    ax[2].plot(t, cyc["p_rider"], color=C[6], lw=1.0, label="rider")
    ax[2].set_ylabel("power [W]")
    ax[2].set_xlabel("time [min]")
    ax[2].legend(fontsize=8, loc="upper right")
    return save(fig, out, "fig_drive_cycle.png")


# ---------------------------------------------------------------------------
def trade_study(quick: bool) -> list[dict]:
    rows = []
    ratios = [5, 8] if quick else [5, 8, 11]
    loadings = [45e3] if quick else [60e3, 45e3, 35e3]
    for g in ratios:
        for a in loadings:
            spec = EBikeSpec(gear_ratio=g)
            dd = design_drive(spec, elec_loading_peak=a)
            op = operating_point(dd.drive, dd.req.base_speed_rpm, dd.req.rated_torque)
            s = dd.design.summary()
            rows.append(
                {
                    "gear_ratio": g,
                    "A_peak_kA_per_m": a / 1e3,
                    "D_out_mm": s["outer_diameter_mm"],
                    "L_stk_mm": s["stack_length_mm"],
                    "active_mass_kg": s["mass_active_kg"],
                    "I_peak_A": round(dd.drive.i_max, 1),
                    "eta_motor_rated_pct": round(100 * op.efficiency()["eta_motor"], 1) if op else "n/a",
                    "eta_system_rated_pct": round(100 * op.efficiency()["eta_system"], 1) if op else "n/a",
                }
            )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results")
    ap.add_argument("--quick", action="store_true", help="coarser grids for a fast smoke run")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    dd = design_drive()
    spec, req, d, drive = dd.spec, dd.req, dd.design, dd.drive
    print(f"[{time.time() - t0:5.1f}s] machine sized: {d.summary()['configuration']}, I_max={drive.i_max:.1f} A")

    rated, peak = rated_and_peak_points(dd)
    cruise = operating_point(drive, req.base_speed_rpm, req.cruise_torque)
    hill_n = req.base_speed_rpm * 15 / spec.v_max_kmh
    hill = operating_point(drive, hill_n, req.hill_torque)
    single_pulse = operating_point(drive, req.max_speed_rpm, 0.5 * req.rated_torque)

    n_speeds = 6 if args.quick else 10
    speeds = np.linspace(0.08, 1.2, n_speeds) * req.base_speed_rpm
    env = torque_speed_envelope(drive, speeds)
    print(f"[{time.time() - t0:5.1f}s] envelope done")
    t_env = np.array([e.torque_avg for e in env])
    torques = np.linspace(0.25, t_env.max(), 6 if args.quick else 14)
    sp, tq, eta_s, eta_m, _ = efficiency_map(drive, speeds, torques, env)
    print(f"[{time.time() - t0:5.1f}s] efficiency map done")
    lookup = EfficiencyLookup(sp, t_env, tq, eta_s)
    cyc = simulate_route(spec, lookup, DEFAULT_ROUTE)
    trades = trade_study(args.quick)
    print(f"[{time.time() - t0:5.1f}s] drive cycle + trade study done")

    figs = [
        fig_cross_section(dd, out),
        fig_feasible_triangle(dd, out),
        fig_magnetics(dd, out),
        fig_waveforms(peak, out, "fig_waveforms_low_speed.png", "Current chopping (launch)"),
        fig_waveforms(rated, out, "fig_waveforms_rated.png", "Rated point (25 km/h)"),
        fig_waveforms(single_pulse, out, "fig_waveforms_single_pulse.png", "Single-pulse (overspeed)"),
        fig_torque_speed(dd, env, out),
        fig_eff_map(sp, tq, eta_s, env, cyc, out),
        fig_losses({"launch (peak)": peak, "hill 6 %": hill, "rated 250 W": rated, "cruise 25 km/h": cruise}, out),
        fig_cycle(cyc, out),
    ]

    conv = converter_ratings(dd, peak, rated)
    therm = thermal_estimate(dd, rated)
    req_tbl = {
        "Assist cut-off speed": f"{spec.v_max_kmh:.0f} km/h",
        "Gear ratio (motor : wheel)": f"{spec.gear_ratio:.0f} : 1",
        "Motor speed at cut-off (base)": f"{req.base_speed_rpm:.0f} rpm",
        "Corner speed (end of constant torque)": f"{req.corner_speed_rpm:.0f} rpm",
        "Rated torque (250 W @ base)": f"{req.rated_torque:.2f} N m",
        "Cruise torque (flat, 25 km/h)": f"{req.cruise_torque:.2f} N m ({req.cruise_power:.0f} W)",
        "Hill torque (6 %, 15 km/h)": f"{req.hill_torque:.2f} N m",
        "Launch torque (8 %, 0.5 m/s^2)": f"{req.start_torque:.2f} N m",
        "Design peak torque": f"{req.peak_torque:.2f} N m",
    }
    op_rows = [
        {"point": k, **v.summary()}
        for k, v in {
            "launch (peak)": peak,
            "hill 6 %": hill,
            "rated": rated,
            "cruise": cruise,
            "overspeed": single_pulse,
        }.items()
    ]
    cyc_tbl = {
        "Distance": f"{cyc['distance_km']:.2f} km",
        "Duration": f"{cyc['duration_min']:.1f} min",
        "Rider energy": f"{cyc['energy_rider_Wh']:.1f} Wh",
        "Motor shaft energy": f"{cyc['energy_motor_Wh']:.1f} Wh",
        "Battery energy": f"{cyc['energy_battery_Wh']:.1f} Wh",
        "Cycle-average drive efficiency": f"{100 * cyc['cycle_efficiency']:.1f} %",
        "Consumption": f"{cyc['Wh_per_km']:.2f} Wh/km",
        "Estimated range (90 % of 360 Wh)": f"{cyc['range_km']:.0f} km",
        "Time with torque deficit": f"{cyc['torque_deficit_s']:.1f} s",
    }

    summary = {
        "requirements": req.__dict__,
        "design": d.summary(),
        "design_notes": d.notes,
        "i_peak_required_A": dd.i_peak_required,
        "operating_points": op_rows,
        "converter": conv,
        "thermal": therm,
        "drive_cycle": {k: v for k, v in cyc.items() if not isinstance(v, np.ndarray)},
        "pole_configs": compare_pole_configs(req.max_speed_rpm),
        "trade_study": trades,
        "figures": figs,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    md = [
        "# Generated results",
        "",
        "Produced by `python scripts/run_analysis.py`. Do not edit by hand.",
        "",
        "## Drive requirements",
        kv_table(req_tbl),
        "",
        "## Pole-configuration comparison (at max. speed)",
        md_table(summary["pole_configs"]),
        "",
        "## Machine design",
        kv_table(d.summary()),
        "",
        "## Operating points",
        md_table(op_rows),
        "",
        "## Converter ratings",
        kv_table(conv),
        "",
        "## Thermal estimate (continuous rated point)",
        kv_table(therm),
        "",
        "## Route simulation",
        kv_table(cyc_tbl),
        "",
        "## Trade study: gear ratio and electric loading",
        md_table(trades),
        "",
    ]
    (out / "results.md").write_text("\n".join(md))
    print(f"[{time.time() - t0:5.1f}s] wrote {out / 'results.md'} and {len(figs)} figures")


if __name__ == "__main__":
    main()
