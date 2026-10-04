"""Transient lumped-parameter thermal network of the hub SRM.

Nodes: winding (copper + insulation), stator core, rotor, and hub shell.

    winding --G_ws--> stator --G_sh--> shell --h A_shell--> ambient
                                  rotor --G_rs (air gap + bearings)--> shell

Copper loss is temperature dependent (R ~ 1 + alpha (T - 20)); core loss is
split between stator and rotor in the ratio given by the loss model.
All conductances are estimated from the geometry with stated material data;
they are design estimates, not calibrated values.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geometry import ALPHA_CU, SRMDesign

C_CU, C_FE, C_AL = 385.0, 460.0, 900.0  # J/(kg K)


@dataclass
class ThermalParams:
    h_conv: float = 25.0  # shell -> ambient, W/(m^2 K)
    t_amb: float = 30.0
    k_slot: float = 0.25  # W/(m K), equivalent slot liner + impregnation
    t_liner: float = 0.35e-3  # m
    h_contact: float = 1500.0  # stator -> shell interface, W/(m^2 K)
    g_rotor: float = 0.6  # rotor -> shell via air gap, end air and bearings, W/K
    m_shell: float = 1.1  # kg aluminium hub shell incl. covers
    hotspot_factor: float = 1.15  # hot-spot rise / mean winding rise


class ThermalModel:
    def __init__(self, d: SRMDesign, p: ThermalParams | None = None):
        self.d, self.p = d, p or ThermalParams()
        p = self.p
        r_o = d.outer_diameter / 2
        # winding <-> stator: slot perimeter (two pole sides + yoke) per slot
        slot_perim = 2 * d.stator_pole_height + 2 * math.pi * (d.bore_diameter / 2 + d.stator_pole_height) / d.ns
        a_ws = d.ns * slot_perim * d.stack_length
        self.g_ws = p.k_slot * a_ws / p.t_liner
        self.g_sh = p.h_contact * 2 * math.pi * r_o * d.stack_length
        r_shell = r_o + 0.01
        l_shell = d.stack_length + 0.04
        self.a_shell = 2 * math.pi * r_shell * l_shell + 2 * math.pi * r_shell**2
        self.g_amb = p.h_conv * self.a_shell
        self.c = np.array(
            [
                d.mass_copper * C_CU * 1.1,  # + insulation / varnish
                d.mass_iron_stator * C_FE,
                d.mass_iron_rotor * C_FE,
                p.m_shell * C_AL,
            ]
        )

    def simulate(self, t_end: float, loss_fn, dt: float = 0.5, t0=None) -> dict:
        """Integrate the network; ``loss_fn(t) -> (p_cu_at_100C, p_core_stator, p_core_rotor, p_mech)``."""
        p = self.p
        n = int(t_end / dt) + 1
        temp = np.full(4, p.t_amb) if t0 is None else np.array(t0, float)
        out = np.empty((n, 4))
        for k in range(n):
            out[k] = temp
            p_cu100, p_fs, p_fr, p_m = loss_fn(k * dt)
            p_cu = p_cu100 * (1 + ALPHA_CU * (temp[0] - 20)) / (1 + ALPHA_CU * 80)
            tw, ts, tr, tsh = temp
            q = np.array(
                [
                    p_cu - self.g_ws * (tw - ts),
                    p_fs + self.g_ws * (tw - ts) - self.g_sh * (ts - tsh),
                    p_fr + p_m - p.g_rotor * (tr - tsh),
                    self.g_sh * (ts - tsh) + p.g_rotor * (tr - tsh) - self.g_amb * (tsh - p.t_amb),
                ]
            )
            temp = temp + dt * q / self.c
        t = np.arange(n) * dt
        return {
            "t": t,
            "winding": out[:, 0],
            "stator": out[:, 1],
            "rotor": out[:, 2],
            "shell": out[:, 3],
            "hotspot": p.t_amb + p.hotspot_factor * (out[:, 0] - p.t_amb),
        }

    def steady_state(self, p_cu100, p_fs, p_fr, p_m) -> dict:
        res = self.simulate(6 * 3600, lambda t: (p_cu100, p_fs, p_fr, p_m), dt=2.0)
        return {k: float(v[-1]) for k, v in res.items() if k != "t"}
