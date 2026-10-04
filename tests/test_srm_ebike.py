import math

import numpy as np
import pytest

from srm_ebike.drive import SRMDrive
from srm_ebike.drive_cycle import DEFAULT_ROUTE, speed_profile
from srm_ebike.geometry import PoleConfig, SizingInputs, check_pole_arcs, size_srm
from srm_ebike.magnetics import MagneticModel
from srm_ebike.performance import max_torque, operating_point
from srm_ebike.vehicle import EBikeSpec, derive_requirements, road_load_force


@pytest.fixture(scope="module")
def baseline():
    spec = EBikeSpec()
    req = derive_requirements(spec)
    design = size_srm(
        SizingInputs(peak_torque=req.peak_torque, base_speed_rpm=req.corner_speed_rpm, dc_voltage=spec.battery_voltage_min)
    )
    drive = SRMDrive(design, v_dc=spec.battery_voltage, i_max=30.0)
    return spec, req, design, drive


# --- vehicle -----------------------------------------------------------------
def test_road_load_flat_cruise():
    spec = EBikeSpec()
    f = road_load_force(spec, spec.v_max)
    # rolling 5.9 N + aero 14.5 N at 25 km/h
    assert f == pytest.approx(5.886 + 0.5 * 1.2 * 0.5 * spec.v_max**2, rel=1e-3)


def test_requirements_consistent():
    spec = EBikeSpec()
    req = derive_requirements(spec)
    assert req.corner_speed_rpm < req.base_speed_rpm
    assert req.peak_torque * req.corner_speed_rpm * 2 * math.pi / 60 == pytest.approx(spec.peak_power)
    assert req.rated_torque * spec.motor_speed(spec.v_max) == pytest.approx(spec.rated_power)


# --- geometry --------------------------------------------------------------------
def test_pole_config_phases():
    assert PoleConfig(12, 8).phases == 3
    assert PoleConfig(8, 6).phases == 4
    assert PoleConfig(6, 4).stroke_angle_deg == pytest.approx(30.0)


def test_pole_arc_rules():
    deg = math.radians
    assert check_pole_arcs(deg(15), deg(17), 12, 8) == []
    assert check_pole_arcs(deg(12), deg(17), 12, 8)  # below stroke angle
    assert check_pole_arcs(deg(20), deg(26), 12, 8)  # no unaligned region


def test_design_is_physical(baseline):
    _, _, d, _ = baseline
    assert check_pole_arcs(d.beta_s, d.beta_r, d.ns, d.nr) == []
    assert d.stator_pole_height > 0 and d.rotor_pole_height > 0
    assert 4 < d.inductance_ratio < 12
    assert d.l_sat < d.l_aligned
    # slot copper fits by construction
    assert 2 * d.turns_per_pole * d.conductor_area <= 0.46 * d.slot_area


# --- magnetics -------------------------------------------------------------------
def test_flux_inversion_roundtrip(baseline):
    mm = baseline[3].model
    for th in np.linspace(0, mm.tau_r, 17):
        for i in (0.5, 5.0, 20.0, 45.0):
            assert mm.current(mm.psi(i, th), th) == pytest.approx(i, rel=1e-6)


def test_torque_is_coenergy_derivative(baseline):
    mm = baseline[3].model
    i = 15.0
    for th in np.linspace(0.05, mm.tau_r - 0.05, 23):
        h = 1e-6
        w1 = mm.l_u * i * i / 2 + mm.w(th + h) * (mm.coenergy_aligned(i) - mm.l_u * i * i / 2)
        w0 = mm.l_u * i * i / 2 + mm.w(th - h) * (mm.coenergy_aligned(i) - mm.l_u * i * i / 2)
        assert mm.torque(i, th) == pytest.approx((w1 - w0) / (2 * h), abs=1e-4)


def test_motoring_region_sign(baseline):
    mm = baseline[3].model
    assert mm.torque(10, mm.tau_r / 4) > 0
    assert mm.torque(10, 3 * mm.tau_r / 4) < 0
    assert mm.torque(10, 0.0) == 0.0


def test_saturation(baseline):
    mm: MagneticModel = baseline[3].model
    th = mm.tau_r / 2
    assert mm.psi(1.0, th) == pytest.approx(mm.l_a * 1.0, rel=0.05)
    slope_high = mm.psi(60, th) - mm.psi(59, th)
    assert slope_high == pytest.approx(mm.l_sat, rel=0.05)


# --- drive simulation ------------------------------------------------------------
@pytest.mark.parametrize("rpm,i_ref", [(100, 25.0), (800, 20.0), (1600, 15.0)])
def test_energy_balance(baseline, rpm, i_ref):
    drive = baseline[3]
    on, off = drive.default_angles(rpm, i_ref, 0.8)
    res = drive.simulate(rpm, on, off, i_ref)
    ls = res.losses()
    p_circuit = res.p_electromagnetic + ls["copper"] + ls["converter_conduction"]
    assert res.p_dc_circuit == pytest.approx(p_circuit, rel=0.03)
    assert not res.continuous


def test_current_regulated(baseline):
    drive = baseline[3]
    on, off = drive.default_angles(100, 20.0)
    res = drive.simulate(100, on, off, 20.0)
    assert res.i_peak < 20.0 * 1.06 + 0.5
    assert res.transitions > 2


def test_peak_torque_reachable(baseline):
    _, req, _, drive = baseline
    res = max_torque(drive, 0.1 * req.corner_speed_rpm)
    assert res.torque_avg >= req.peak_torque


def test_operating_point_hits_target(baseline):
    _, req, _, drive = baseline
    res = operating_point(drive, req.base_speed_rpm, req.rated_torque)
    assert res is not None
    assert res.torque_avg == pytest.approx(req.rated_torque, rel=0.02)
    eff = res.efficiency()
    assert 0.6 < eff["eta_system"] < eff["eta_motor"] < 0.95


# --- drive cycle -------------------------------------------------------------------
def test_speed_profile_route_length():
    prof = speed_profile(DEFAULT_ROUTE)
    assert prof["x"][-1] == pytest.approx(sum(s.length for s in DEFAULT_ROUTE), rel=1e-3)
    assert prof["v"].max() <= 25 / 3.6 + 1e-9
    assert np.abs(prof["a"]).max() <= 1.0 + 1e-9
