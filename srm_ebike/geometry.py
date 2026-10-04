"""Analytical sizing of the switched reluctance machine.

The procedure follows the classical SRM design flow (Krishnan, *Switched
Reluctance Motor Drives*, ch. 3; Miller, *Switched Reluctance Motors and their
Control*):

1. pole-number selection,
2. bore diameter and stack length from the output equation,
3. pole arcs inside the Lawrenson feasible triangle,
4. pole heights / yoke thickness from flux continuity,
5. turns per phase from the volt-second balance at base speed,
6. slot fill -> conductor size -> phase resistance,
7. aligned / unaligned inductance from simple permeance models.

Angles are mechanical radians unless noted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

MU0 = 4e-7 * math.pi
RHO_CU_20 = 1.72e-8  # ohm m
ALPHA_CU = 0.00393  # 1/K
DENSITY_FE = 7650.0
DENSITY_CU = 8900.0


# ---------------------------------------------------------------------------
# Pole-number selection
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PoleConfig:
    ns: int
    nr: int

    @property
    def phases(self) -> int:
        # ns = 2*q*m with q pole pairs per phase; the common configurations
        # (6/4, 8/6, 12/8, 16/12 ...) all satisfy m = ns / gcd(ns, nr)
        return self.ns // math.gcd(self.ns, self.nr)

    @property
    def stroke_angle_deg(self) -> float:
        return 360.0 / (self.phases * self.nr)

    @property
    def strokes_per_rev(self) -> int:
        return self.phases * self.nr

    def phase_frequency(self, n_rpm: float) -> float:
        """Fundamental frequency of the phase current [Hz]."""
        return self.nr * n_rpm / 60.0


def compare_pole_configs(n_max_rpm: float) -> list[dict]:
    """Qualitative/quantitative comparison table for candidate topologies."""
    rows = []
    for ns, nr in [(6, 4), (8, 6), (12, 8), (16, 12), (24, 16)]:
        pc = PoleConfig(ns, nr)
        m = pc.phases
        rows.append(
            {
                "config": f"{ns}/{nr}",
                "phases": m,
                "stroke_deg": round(pc.stroke_angle_deg, 2),
                "strokes_per_rev": pc.strokes_per_rev,
                "f_phase_at_nmax_Hz": round(pc.phase_frequency(n_max_rpm), 1),
                "switches_asym_bridge": 2 * m,
                # m >= 3 always leaves a phase with positive dT/dtheta in either
                # direction (given Lawrenson pole arcs); m = 2 needs stepped gaps
                "self_start_both_dir": m >= 3,
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Design data
# ---------------------------------------------------------------------------
@dataclass
class SizingInputs:
    peak_torque: float  # N m, short-term torque to be met
    base_speed_rpm: float
    dc_voltage: float  # design (minimum) DC-link voltage
    ns: int = 12
    nr: int = 8
    # output-equation coefficients (Krishnan)
    b_pole: float = 1.7  # stator-pole flux density at aligned position, T
    elec_loading_peak: float = 45e3  # A/m, short-term (1-2 min) electric loading
    k_e: float = 0.85  # efficiency
    k_d: float = 0.90  # duty cycle (fraction of stroke the phase conducts)
    k_2: float = 0.70  # 1 - 1/lambda_u, aligned/unaligned inductance factor
    l_over_d: float = 0.5  # aspect ratio, pancake shape for a hub motor
    split_ratio: float = 0.53  # bore / outer diameter
    air_gap: float = 0.30e-3
    beta_s_deg: float | None = None  # stator pole arc (default: stroke angle)
    beta_r_deg: float | None = None  # rotor pole arc (default: beta_s + 2 deg)
    yoke_factor_stator: float = 0.75  # stator yoke thickness / stator pole width
    yoke_factor_rotor: float = 0.85  # rotor yoke thickness / stator pole width
    shaft_diameter: float = 15e-3
    dwell_strokes: float = 1.0  # conduction dwell at base speed in strokes
    fill_factor: float = 0.45  # net copper / slot area, bobbin-wound concentrated coils
    winding_temp: float = 100.0  # deg C used for resistance
    b_knee: float = 1.55  # flux density at the knee of the aligned curve, T
    iron_mmf_factor: float = 0.90  # reduction of aligned L by iron reluctance
    end_leakage_factor: float = 1.35  # end-winding + leakage on L_u


@dataclass
class SRMDesign:
    """Complete geometric/electrical description of the designed machine."""

    ns: int
    nr: int
    phases: int
    bore_diameter: float
    stack_length: float
    outer_diameter: float
    air_gap: float
    beta_s: float
    beta_r: float
    stator_pole_width: float
    rotor_pole_width: float
    stator_pole_height: float
    rotor_pole_height: float
    stator_yoke: float
    rotor_yoke: float
    shaft_diameter: float
    turns_per_phase: int
    turns_per_pole: int
    conductor_area: float
    wire_diameter: float
    mean_turn_length: float
    slot_area: float
    phase_resistance: float
    l_aligned: float  # unsaturated aligned inductance
    l_unaligned: float
    l_sat: float  # incremental inductance in deep saturation (aligned)
    psi_knee: float
    fringe_angle: float
    mass_iron_stator: float
    mass_iron_rotor: float
    mass_copper: float
    notes: list[str] = field(default_factory=list)

    @property
    def rotor_pole_pitch(self) -> float:
        return 2 * math.pi / self.nr

    @property
    def stroke_angle(self) -> float:
        return self.rotor_pole_pitch / self.phases

    @property
    def poles_per_phase(self) -> int:
        return self.ns // self.phases

    @property
    def pole_face_area(self) -> float:
        return self.stator_pole_width * self.stack_length

    @property
    def inductance_ratio(self) -> float:
        return self.l_aligned / self.l_unaligned

    @property
    def mass_active(self) -> float:
        return self.mass_iron_stator + self.mass_iron_rotor + self.mass_copper

    @property
    def rotor_inertia(self) -> float:
        """Approximate rotor inertia [kg m^2] (solid annulus + poles)."""
        r_o = self.bore_diameter / 2 - self.air_gap
        r_i = self.shaft_diameter / 2
        r_y = r_o - self.rotor_pole_height
        j_yoke = 0.5 * math.pi * DENSITY_FE * self.stack_length * (r_y**4 - r_i**4)
        m_poles = DENSITY_FE * self.stack_length * self.nr * self.rotor_pole_width * self.rotor_pole_height
        r_poles = 0.5 * (r_o + r_y)
        return j_yoke + m_poles * r_poles**2

    def summary(self) -> dict:
        mm = 1e3
        return {
            "configuration": f"{self.ns}/{self.nr}, {self.phases}-phase",
            "outer_diameter_mm": round(self.outer_diameter * mm, 1),
            "bore_diameter_mm": round(self.bore_diameter * mm, 1),
            "stack_length_mm": round(self.stack_length * mm, 1),
            "air_gap_mm": round(self.air_gap * mm, 2),
            "beta_s_deg": round(math.degrees(self.beta_s), 2),
            "beta_r_deg": round(math.degrees(self.beta_r), 2),
            "stator_pole_width_mm": round(self.stator_pole_width * mm, 2),
            "rotor_pole_width_mm": round(self.rotor_pole_width * mm, 2),
            "stator_pole_height_mm": round(self.stator_pole_height * mm, 2),
            "rotor_pole_height_mm": round(self.rotor_pole_height * mm, 2),
            "stator_yoke_mm": round(self.stator_yoke * mm, 2),
            "rotor_yoke_mm": round(self.rotor_yoke * mm, 2),
            "turns_per_phase": self.turns_per_phase,
            "turns_per_pole": self.turns_per_pole,
            "wire_diameter_mm": round(self.wire_diameter * mm, 3),
            "slot_area_mm2": round(self.slot_area * 1e6, 1),
            "mean_turn_length_mm": round(self.mean_turn_length * mm, 1),
            "phase_resistance_mohm": round(self.phase_resistance * 1e3, 1),
            "L_aligned_mH": round(self.l_aligned * 1e3, 3),
            "L_unaligned_mH": round(self.l_unaligned * 1e3, 3),
            "L_sat_mH": round(self.l_sat * 1e3, 3),
            "inductance_ratio": round(self.inductance_ratio, 2),
            "psi_knee_mWb": round(self.psi_knee * 1e3, 2),
            "mass_iron_kg": round(self.mass_iron_stator + self.mass_iron_rotor, 3),
            "mass_copper_kg": round(self.mass_copper, 3),
            "mass_active_kg": round(self.mass_active, 3),
            "rotor_inertia_kgm2": float(f"{self.rotor_inertia:.3e}"),
        }


# ---------------------------------------------------------------------------
# Sizing procedure
# ---------------------------------------------------------------------------
def output_equation_d2l(inp: SizingInputs) -> float:
    """Return D^2 L [m^3] from the SRM output equation.

    Krishnan's output equation P = k_e k_d k_1 k_2 B A D^2 L N_r with
    k_1 = pi^2/120 and N_r in rpm becomes, in torque form,
    T = (pi/4) k_e k_d k_2 B A D^2 L.
    """
    k = (math.pi / 4) * inp.k_e * inp.k_d * inp.k_2 * inp.b_pole * inp.elec_loading_peak
    return inp.peak_torque / k


def check_pole_arcs(beta_s: float, beta_r: float, ns: int, nr: int) -> list[str]:
    """Validate pole arcs against Lawrenson's feasible triangle."""
    m = PoleConfig(ns, nr).phases
    eps = 2 * math.pi / (m * nr)
    issues = []
    if min(beta_s, beta_r) < eps - 1e-12:
        issues.append("min(beta_s, beta_r) < stroke angle: no self-starting in all rotor positions")
    if beta_s > beta_r + 1e-12:
        issues.append("beta_s > beta_r: reduces slot area and unaligned clearance")
    if beta_s + beta_r >= 2 * math.pi / nr:
        issues.append("beta_s + beta_r >= rotor pole pitch: no zero-overlap (unaligned) region")
    return issues


def size_srm(inp: SizingInputs) -> SRMDesign:
    pc = PoleConfig(inp.ns, inp.nr)
    m = pc.phases
    notes: list[str] = []

    # --- 1. main dimensions -------------------------------------------------
    d2l = output_equation_d2l(inp)
    d = (d2l / inp.l_over_d) ** (1 / 3)
    length = inp.l_over_d * d
    d_out = d / inp.split_ratio
    notes.append(f"Output equation: D^2L = {d2l * 1e6:.1f} cm^3")

    # --- 2. pole arcs ---------------------------------------------------------
    eps = 2 * math.pi / (m * inp.nr)
    beta_s = math.radians(inp.beta_s_deg) if inp.beta_s_deg else eps
    beta_r = math.radians(inp.beta_r_deg) if inp.beta_r_deg else beta_s + math.radians(2.0)
    notes.extend(check_pole_arcs(beta_s, beta_r, inp.ns, inp.nr))

    r_bore = d / 2
    r_rotor = r_bore - inp.air_gap
    w_s = 2 * r_bore * math.sin(beta_s / 2)
    w_r = 2 * r_rotor * math.sin(beta_r / 2)

    # --- 3. yokes and pole heights -------------------------------------------
    # pole flux splits in two halves in the yoke -> y >= w/2; extra margin
    # lowers yoke saturation and acoustic noise (radial-force ovalisation)
    y_s = max(inp.yoke_factor_stator * w_s, 0.5 * w_s)
    y_r = max(inp.yoke_factor_rotor * w_s, 0.5 * w_s)
    r_out = d_out / 2
    h_s = r_out - y_s - r_bore
    r_shaft = inp.shaft_diameter / 2
    h_r = r_rotor - y_r - r_shaft
    if h_r < 0:
        raise ValueError("rotor does not fit: increase D or reduce shaft diameter")
    # deep rotor slots help L_u, but very tall poles are weak; Miller suggests
    # h_r ~ 20-30 x g. Cap at 30 g and give the rest to the rotor yoke.
    h_r_max = 30 * inp.air_gap
    if h_r > h_r_max:
        y_r += h_r - h_r_max
        h_r = h_r_max

    # --- 4. turns per phase from volt-seconds at base speed -------------------
    omega_b = inp.base_speed_rpm * 2 * math.pi / 60
    dwell = inp.dwell_strokes * eps
    psi_peak = inp.dc_voltage * dwell / omega_b
    phi_pole = inp.b_pole * w_s * length
    t_ph = psi_peak / phi_pole
    poles_ph = inp.ns // m
    n_c = max(1, round(t_ph / poles_ph))
    t_ph = n_c * poles_ph
    notes.append(f"Single-pulse flux at base speed: psi = {psi_peak * 1e3:.1f} mWb")

    # --- 5. slot, conductor and resistance -------------------------------------
    r_y_in = r_bore + h_s
    annulus = math.pi * (r_y_in**2 - r_bore**2)
    a_slot = (annulus - inp.ns * w_s * h_s) / inp.ns
    a_cond = inp.fill_factor * a_slot / (2 * n_c)  # two coil sides per slot
    d_wire = math.sqrt(4 * a_cond / math.pi)
    coil_thickness = a_slot / (2 * h_s)
    mlt = 2 * (length + w_s) + math.pi * coil_thickness
    rho = RHO_CU_20 * (1 + ALPHA_CU * (inp.winding_temp - 20))
    r_ph = rho * t_ph * mlt / a_cond

    # --- 6. inductances ------------------------------------------------------
    # aligned: two air gaps in series per pole pair, all poles of a phase in
    # series -> L = mu0 A T^2 / (p g) with fringing via (w + 2g)
    a_eff = (w_s + 2 * inp.air_gap) * length
    l_a = inp.iron_mmf_factor * MU0 * a_eff * t_ph**2 / (poles_ph * inp.air_gap)
    # unaligned: centre path into the rotor slot (g + h_r) in parallel with two
    # side paths to the neighbouring rotor-pole corners
    l_side = r_bore * (math.pi / inp.nr - (beta_s + beta_r) / 2) + inp.air_gap
    l_side = max(l_side, 2 * inp.air_gap)
    perm_ratio = inp.air_gap / (inp.air_gap + h_r) + inp.air_gap / l_side
    l_u = inp.end_leakage_factor * l_a / inp.iron_mmf_factor * perm_ratio
    psi_k = inp.b_knee * w_s * length * t_ph
    l_sat = 1.1 * l_u

    # --- 7. masses -------------------------------------------------------------
    m_s = DENSITY_FE * length * (math.pi * (r_out**2 - r_y_in**2) + inp.ns * w_s * h_s)
    r_ry = r_rotor - h_r
    m_r = DENSITY_FE * length * (math.pi * (r_ry**2 - r_shaft**2) + inp.nr * w_r * h_r)
    m_cu = DENSITY_CU * m * t_ph * mlt * a_cond

    return SRMDesign(
        ns=inp.ns,
        nr=inp.nr,
        phases=m,
        bore_diameter=d,
        stack_length=length,
        outer_diameter=d_out,
        air_gap=inp.air_gap,
        beta_s=beta_s,
        beta_r=beta_r,
        stator_pole_width=w_s,
        rotor_pole_width=w_r,
        stator_pole_height=h_s,
        rotor_pole_height=h_r,
        stator_yoke=y_s,
        rotor_yoke=y_r,
        shaft_diameter=inp.shaft_diameter,
        turns_per_phase=t_ph,
        turns_per_pole=n_c,
        conductor_area=a_cond,
        wire_diameter=d_wire,
        mean_turn_length=mlt,
        slot_area=a_slot,
        phase_resistance=r_ph,
        l_aligned=l_a,
        l_unaligned=l_u,
        l_sat=l_sat,
        psi_knee=psi_k,
        fringe_angle=2.0 * inp.air_gap / r_bore,
        mass_iron_stator=m_s,
        mass_iron_rotor=m_r,
        mass_copper=m_cu,
        notes=notes,
    )
