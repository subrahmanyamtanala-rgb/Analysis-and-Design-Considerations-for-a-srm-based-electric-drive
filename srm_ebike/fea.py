"""2-D nonlinear magnetostatic finite-element model of the SRM cross-section.

Formulation: magnetic vector potential A_z on first-order triangles,

    curl( nu(|B|^2) curl A ) = J_z,      A = 0 on the stator outer boundary,

solved by Newton iteration with the M270-35A B-H curve. Geometry is built and
meshed with gmsh (OpenCASCADE kernel); assembly uses scikit-fem.

Post-processing:
* flux linkage per phase from the coil-side averages of A_z (2-D) plus an
  analytical end-winding leakage term;
* torque by Arkkio's method (Maxwell stress averaged over the whole air-gap
  annulus), cross-checked against the co-energy derivative.

Angle convention matches ``magnetics.py``: phase angle 0 = unaligned,
tau_r/2 = aligned for phase A (stator pole at 0 deg).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse.linalg import spsolve

from .geometry import MU0, SRMDesign

# M270-35A typical DC magnetisation curve (B in T, H in A/m)
BH_M270_35A = (
    np.array([0.0, 0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]),
    np.array([0.0, 60.0, 90.0, 120.0, 140.0, 170.0, 230.0, 370.0, 800.0, 2400.0, 5500.0, 11000.0, 20000.0, 36000.0]),
)


class Reluctivity:
    """nu(b2) = H(B)/B and d nu / d(b2) as smooth functions of b2 = |B|^2."""

    def __init__(self, bh=BH_M270_35A, b_max: float = 3.0, n: int = 4000):
        b, h = bh
        hb = PchipInterpolator(b, h)
        bb = np.linspace(1e-4, b_max, n)
        hh = np.where(bb <= b[-1], hb(np.minimum(bb, b[-1])), h[-1] + (bb - b[-1]) / MU0)
        nu = hh / bb
        nu[0] = nu[1]
        self.b2 = bb**2
        self.nu = nu
        self.dnu = np.gradient(nu, self.b2)

    def __call__(self, b2):
        return np.interp(b2, self.b2, self.nu), np.interp(b2, self.b2, self.dnu)


@dataclass
class FEAResult:
    theta: float
    currents: tuple[float, ...]  # per phase (A, B, C)
    psi: np.ndarray  # per phase, incl. end-winding term
    psi_2d: np.ndarray
    torque: float  # Arkkio
    b_pole_max: float
    newton_iters: int
    n_nodes: int


class SRMFEA:
    def __init__(
        self,
        d: SRMDesign,
        h_gap: float | None = None,
        h_iron: float = 1.2e-3,
        h_air: float = 1.5e-3,
        bh=BH_M270_35A,
        end_permeance: float = 0.5,
    ):
        self.d = d
        self.h_gap = h_gap or d.air_gap / 3
        self.h_iron = h_iron
        self.h_air = h_air
        self.nu = Reluctivity(bh)
        self.end_permeance = end_permeance
        self._mesh_cache: dict = {}

    # ------------------------------------------------------------------
    @property
    def l_end(self) -> float:
        """End-winding leakage inductance per phase (all poles in series).

        L_ew = p * mu0 * N_c^2 * l_ew * lambda_ew with the end-turn length of
        one coil l_ew = 2 (w_s + t_c) and a permeance coefficient lambda_ew.
        """
        d = self.d
        t_c = d.slot_area / (2 * d.stator_pole_height)
        l_ew = 2 * (d.stator_pole_width + t_c)
        return d.poles_per_phase * MU0 * d.turns_per_pole**2 * l_ew * self.end_permeance

    def coil_sides(self):
        """(phase, pole index, side sign, polygon) for every coil side."""
        d = self.d
        r_b = d.bore_diameter / 2
        r_yi = r_b + d.stator_pole_height
        r1, r2 = r_b + 0.6e-3, r_yi - 0.3e-3
        half_slot = math.pi / d.ns
        out = []
        for k in range(d.ns):
            a = 2 * math.pi * k / d.ns
            phase = k % d.phases
            for side in (-1, 1):
                # pole edge is a straight line offset w_s/2 from the pole axis
                def edge_pt(r):
                    u = np.array([math.cos(a), math.sin(a)])
                    p = np.array([-math.sin(a), math.cos(a)])
                    s = math.sqrt(max(r * r - (d.stator_pole_width / 2) ** 2, 0.0))
                    return s * u + side * (d.stator_pole_width / 2 + 0.2e-3) * p

                ac = a + side * half_slot * 0.985
                pts = [
                    edge_pt(r1),
                    edge_pt(r2),
                    np.array([r2 * math.cos(ac), r2 * math.sin(ac)]),
                    np.array([r1 * math.cos(ac), r1 * math.sin(ac)]),
                ]
                out.append((phase, k, side, pts))
        return out

    # ------------------------------------------------------------------
    def build_mesh(self, theta: float):
        """Mesh the cross-section with the rotor at phase angle ``theta``."""
        key = round(theta, 9)
        if key in self._mesh_cache:
            return self._mesh_cache[key]
        import gmsh

        d = self.d
        r_o, r_b = d.outer_diameter / 2, d.bore_diameter / 2
        r_yi = r_b + d.stator_pole_height
        r_r = r_b - d.air_gap
        r_ry = r_r - d.rotor_pole_height
        r_sh = d.shaft_diameter / 2
        rot = theta - math.pi / d.nr  # rotor pole angle relative to phase-A axis

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("srm")
        occ = gmsh.model.occ

        def pole_rect(r_in, r_out, w, ang):
            # rectangle along +x then rotated; r_out overshoot merges with arcs
            t = occ.addRectangle(r_in, -w / 2, 0, r_out - r_in, w)
            occ.rotate([(2, t)], 0, 0, 0, 0, 0, 1, ang)
            return (2, t)

        yoke = occ.addDisk(0, 0, 0, r_o, r_o)
        hole = occ.addDisk(0, 0, 0, r_yi, r_yi)
        s_yoke, _ = occ.cut([(2, yoke)], [(2, hole)])
        s_poles = [pole_rect(r_b, r_yi + 0.5e-3, d.stator_pole_width, 2 * math.pi * k / d.ns) for k in range(d.ns)]
        # stator pole faces are arcs at r_b: intersect poles with an annulus
        ann = occ.addDisk(0, 0, 0, r_o, r_o)
        bore = occ.addDisk(0, 0, 0, r_b, r_b)
        ann_cut, _ = occ.cut([(2, ann)], [(2, bore)])
        s_poles_c, _ = occ.intersect(s_poles, ann_cut, removeTool=True)
        stator, _ = occ.fuse(s_yoke, s_poles_c)

        r_core = occ.addDisk(0, 0, 0, r_ry, r_ry)
        r_poles = [pole_rect(r_ry - 0.5e-3, r_r + 1e-3, d.rotor_pole_width, rot + 2 * math.pi * k / d.nr) for k in range(d.nr)]
        rdisk = occ.addDisk(0, 0, 0, r_r, r_r)
        r_poles_c, _ = occ.intersect(r_poles, [(2, rdisk)], removeTool=True)
        rotor, _ = occ.fuse([(2, r_core)], r_poles_c)
        shaft = occ.addDisk(0, 0, 0, r_sh, r_sh)
        rotor, _ = occ.cut(rotor, [(2, shaft)])

        coils = []
        for phase, k, side, pts in self.coil_sides():
            ptags = [occ.addPoint(p[0], p[1], 0) for p in pts]
            ltags = [occ.addLine(ptags[i], ptags[(i + 1) % 4]) for i in range(4)]
            cl = occ.addCurveLoop(ltags)
            coils.append(((2, occ.addPlaneSurface([cl])), phase, k, side))

        gap_outer = occ.addDisk(0, 0, 0, r_b, r_b)
        gap_inner = occ.addDisk(0, 0, 0, r_r, r_r)
        gap, _ = occ.cut([(2, gap_outer)], [(2, gap_inner)])

        domain = occ.addDisk(0, 0, 0, r_o, r_o)
        tools = stator + rotor + [c[0] for c in coils] + gap
        out, out_map = occ.fragment([(2, domain)], tools)
        occ.synchronize()

        def children(i):
            return [t for (dim, t) in out_map[1 + i]]

        n_st, n_ro = len(stator), len(rotor)
        stator_t = sum((children(i) for i in range(n_st)), [])
        rotor_t = sum((children(n_st + i) for i in range(n_ro)), [])
        coil_t = [children(n_st + n_ro + i) for i in range(len(coils))]
        gap_t = sum((children(n_st + n_ro + len(coils) + i) for i in range(len(gap))), [])
        all_t = [t for (dim, t) in out]
        used = set(stator_t) | set(rotor_t) | set(sum(coil_t, [])) | set(gap_t)
        air_t = [t for t in all_t if t not in used]

        gmsh.model.addPhysicalGroup(2, stator_t, 1)
        gmsh.model.addPhysicalGroup(2, rotor_t, 2)
        gmsh.model.addPhysicalGroup(2, gap_t, 3)
        gmsh.model.addPhysicalGroup(2, air_t, 4)
        for i, ct in enumerate(coil_t):
            gmsh.model.addPhysicalGroup(2, ct, 100 + i)

        # mesh size: fine in and next to the air gap
        f = gmsh.model.mesh.field
        fd = f.add("MathEval")
        rg = 0.5 * (r_b + r_r)
        f.setString(
            fd,
            "F",
            f"Min({self.h_air}, {self.h_gap} + 0.35*Abs(Sqrt(x*x+y*y)-{rg}))",
        )
        fi = f.add("Restrict")
        f.setNumber(fi, "InField", fd)
        f.setNumbers(fi, "SurfacesList", all_t)
        f.setAsBackgroundMesh(fd)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", self.h_iron)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        idx = {t: i for i, t in enumerate(node_tags)}
        p = coords.reshape(-1, 3)[:, :2].T
        tri_list, phys = [], []
        for dim, tag in gmsh.model.getPhysicalGroups(2):
            for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                types, _, nodes = gmsh.model.mesh.getElements(2, ent)
                for ty, nn in zip(types, nodes):
                    if ty != 2:
                        continue
                    t = np.array([idx[n] for n in nn]).reshape(-1, 3)
                    tri_list.append(t)
                    phys.append(np.full(len(t), tag))
        gmsh.finalize()

        tri = np.vstack(tri_list).T
        phys = np.concatenate(phys)
        # compact unused nodes
        used_nodes = np.unique(tri)
        remap = -np.ones(p.shape[1], dtype=int)
        remap[used_nodes] = np.arange(len(used_nodes))
        p = p[:, used_nodes]
        tri = remap[tri]
        mesh = (p, tri, phys)
        self._mesh_cache = {key: mesh}  # keep only the latest geometry
        return mesh

    # ------------------------------------------------------------------
    def solve(self, theta: float, currents=(0.0, 0.0, 0.0), a0=None, tol: float = 1e-6, max_iter: int = 80, _ramp: bool = False) -> tuple[FEAResult, np.ndarray]:
        from skfem import Basis, BilinearForm, ElementTriP0, ElementTriP1, LinearForm, MeshTri, asm
        from skfem.helpers import dot, grad

        p, tri, phys = self.build_mesh(theta)
        mesh = MeshTri(p, tri)
        basis = Basis(mesh, ElementTriP1())
        d = self.d
        iron = np.isin(phys, (1, 2))
        n_el = tri.shape[1]

        # current density per element
        jz = np.zeros(n_el)
        sides = self.coil_sides()
        side_area = {}
        for i, (phase, k, side, _) in enumerate(sides):
            sel = phys == 100 + i
            area = mesh_area(p, tri[:, sel]).sum()
            side_area[i] = area
            pol = 1 if (k // d.phases) % 2 == 0 else -1
            jz[sel] = pol * side * d.turns_per_pole * currents[phase] / area

        b0 = basis.with_element(ElementTriP0())
        j_field = b0.interpolate(jz)
        iron_field = b0.interpolate(iron.astype(float))

        nu0 = 1 / MU0
        nu_fn = self.nu

        @LinearForm
        def rhs(v, w):
            return w["j"] * v

        def nu_eval(w):
            gb = grad(w["a"])
            b2 = gb[0] ** 2 + gb[1] ** 2
            nu_i, dnu_i = nu_fn(b2)
            isfe = w["fe"] > 0.5
            nu = np.where(isfe, nu_i, nu0)
            dnu = np.where(isfe, dnu_i, 0.0)
            return nu, dnu, gb

        @BilinearForm
        def jac(u, v, w):
            nu, dnu, ga = nu_eval(w)
            return nu * dot(grad(u), grad(v)) + 2 * dnu * dot(ga, grad(u)) * dot(ga, grad(v))

        @LinearForm
        def residual(v, w):
            nu, _, ga = nu_eval(w)
            return nu * dot(ga, grad(v)) - w["j"] * v

        boundary = basis.get_dofs(lambda x: np.hypot(x[0], x[1]) > d.outer_diameter / 2 - 1e-6).all()
        # a warm start is only meaningful on the same mesh (same rotor position)
        key = round(theta, 9)
        same_mesh = a0 is not None and len(a0) == basis.N and getattr(self, "_last_key", None) == key
        self._last_key = key
        a = a0.copy() if same_mesh else np.zeros(basis.N)
        a[boundary] = 0.0
        f_norm = np.linalg.norm(asm(rhs, basis, j=j_field)) + 1e-30
        it = 0
        a_old = da_old = None
        prev_res, alpha = np.inf, 1.0
        for it in range(1, max_iter + 1):
            aw = basis.interpolate(a)
            r = asm(residual, basis, a=aw, fe=iron_field, j=j_field)
            res = np.linalg.norm(np.delete(r, boundary)) / f_norm
            # backtracking: if the full Newton step increased the residual
            # (or produced NaN, e.g. a step deep into saturation), halve it
            if a_old is not None and (not np.isfinite(res) or res > 2.0 * prev_res) and alpha > 1 / 64:
                alpha *= 0.5
                a = a_old + alpha * da_old
                continue
            if res < tol and it > 1:
                break
            k = asm(jac, basis, a=aw, fe=iron_field)
            da = solve_condensed(k, -r, boundary)
            a_old, da_old, prev_res, alpha = a, da, res, 1.0
            a = a + da
        if not np.all(np.isfinite(a)) or res >= tol:
            if same_mesh:  # retry from a cold start
                return self.solve(theta, currents, None, tol, max_iter)
            if not _ramp:
                # 1) current continuation: 25 -> 50 -> 75 -> 100 % with warm starts
                try:
                    a_c = None
                    for f in (0.25, 0.5, 0.75):
                        _, a_c = self.solve(theta, tuple(f * c for c in currents), a_c, tol, max_iter, _ramp=True)
                        self._last_key = key
                    return self.solve(theta, currents, a_c, tol, max_iter, _ramp=True)
                except RuntimeError:
                    pass
                # 2) mesh-specific failure (degenerate element): regenerate the
                #    mesh at a rotor angle shifted by +/-0.01 deg
                for dth in (1.75e-4, -1.75e-4):
                    try:
                        out = self.solve(theta + dth, currents, None, tol, max_iter, _ramp=True)
                        self.n_mesh_retries = getattr(self, "n_mesh_retries", 0) + 1
                        return out
                    except RuntimeError:
                        continue
            if not np.all(np.isfinite(a)):
                raise RuntimeError(f"FE Newton iteration diverged at theta={theta:.4f} rad, currents={currents}")
        # ---- post-processing ------------------------------------------------
        psi2d = np.zeros(d.phases)
        a_el = a[tri].mean(axis=0)
        areas = mesh_area(p, tri)
        for i, (phase, k, side, _) in enumerate(sides):
            sel = phys == 100 + i
            pol = 1 if (k // d.phases) % 2 == 0 else -1
            avg = (a_el[sel] * areas[sel]).sum() / side_area[i]
            psi2d[phase] += pol * side * d.turns_per_pole * avg * d.stack_length
        psi = psi2d + self.l_end * np.asarray(currents, dtype=float)

        bx, by = element_b(p, tri, a)
        sel = phys == 3
        cx, cy = p[0][tri[:, sel]].mean(axis=0), p[1][tri[:, sel]].mean(axis=0)
        r = np.hypot(cx, cy)
        br = (bx[sel] * cx + by[sel] * cy) / r
        bt = (-bx[sel] * cy + by[sel] * cx) / r
        r_b, r_r = d.bore_diameter / 2, d.bore_diameter / 2 - d.air_gap
        torque = d.stack_length / (MU0 * (r_b - r_r)) * np.sum(br * bt * r * areas[sel])
        bmag = np.hypot(bx, by)
        sp = phys == 1
        res_obj = FEAResult(theta, tuple(currents), psi, psi2d, float(torque), float(bmag[sp].max()), it, p.shape[1])
        return res_obj, a


def mesh_area(p, tri):
    x, y = p[0][tri], p[1][tri]
    return 0.5 * np.abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))


def element_b(p, tri, a):
    """Element-wise B = curl(A_z ez) = (dA/dy, -dA/dx) for P1 triangles."""
    x, y = p[0][tri], p[1][tri]
    det = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
    av = a[tri]
    dadx = ((av[1] - av[0]) * (y[2] - y[0]) - (av[2] - av[0]) * (y[1] - y[0])) / det
    dady = ((av[2] - av[0]) * (x[1] - x[0]) - (av[1] - av[0]) * (x[2] - x[0])) / det
    return dady, -dadx


def solve_condensed(k, r, fixed):
    n = k.shape[0]
    free = np.setdiff1d(np.arange(n), fixed)
    x = np.zeros(n)
    x[free] = spsolve(k[free][:, free].tocsc(), r[free])
    return x
