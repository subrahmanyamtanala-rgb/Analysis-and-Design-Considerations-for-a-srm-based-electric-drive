import math

import numpy as np
import pytest

from srm_ebike import dclink
from srm_ebike.design_flow import design_drive
from srm_ebike.drive import SRMDrive
from srm_ebike.fea import Reluctivity
from srm_ebike.magnetics_fea import TabulatedMagneticModel
from srm_ebike.thermal import ThermalModel, ThermalParams
from srm_ebike.tsf import SHAPES, tsf_profile
from srm_ebike.vehicle import EBikeSpec, motor_torque_demand


@pytest.fixture(scope="module")
def dd():
    return design_drive()


def test_reluctivity_monotone_b_h():
    nu = Reluctivity()
    b = np.linspace(0.05, 2.6, 200)
    h = nu(b**2)[0] * b
    assert np.all(np.diff(h) > 0)
    # deep saturation approaches free space
    assert nu(2.6**2)[0] * 2.6 - nu(2.5**2)[0] * 2.5 == pytest.approx(0.1 / (4e-7 * math.pi), rel=0.05)


def test_tabulated_model_reproduces_source(dd):
    mm = dd.drive.model
    th = np.arange(0, 22.5 + 1e-9, 1.5)
    cur = np.array([1, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 23, 26, 30, 34, 38, 42.0])
    psi = np.array([[mm.psi(i, math.radians(x)) for i in cur] for x in th])
    tab = TabulatedMagneticModel(dd.design, th, cur, psi)
    for x in (3.0, 9.0, 14.0, 20.0):
        for i in (5.0, 15.0, 30.0):
            t = math.radians(x)
            assert tab.psi(i, t) == pytest.approx(mm.psi(i, t), rel=0.03, abs=1e-3)
            assert tab.current(tab.psi(i, t), t) == pytest.approx(i, rel=0.02)
    # torque from the tabulated co-energy matches the analytical torque on average
    xx = np.radians(np.linspace(2, 21, 40))
    t_tab = np.mean([tab.torque(20.0, x) for x in xx])
    t_an = np.mean([mm.torque(20.0, x) for x in xx])
    assert t_tab == pytest.approx(t_an, rel=0.05)


def test_tabulated_drive_energy_balance(dd):
    mm = dd.drive.model
    th = np.arange(0, 22.5 + 1e-9, 1.5)
    cur = np.array([1, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 23, 26, 30, 34, 38, 42.0])
    psi = np.array([[mm.psi(i, math.radians(x)) for i in cur] for x in th])
    tab = TabulatedMagneticModel(dd.design, th, cur, psi)
    dr = SRMDrive(dd.design, v_dc=36.0, i_max=30.0, magnetic_model=tab)
    on, off = dr.default_angles(800, 20.0, 0.8)
    r = dr.simulate(800, on, off, 20.0)
    ls = r.losses()
    assert r.p_dc_circuit == pytest.approx(r.p_electromagnetic + ls["copper"] + ls["converter_conduction"], rel=0.03)


@pytest.mark.parametrize("shape", list(SHAPES))
def test_tsf_shares_sum_to_one(shape):
    eps, ov = math.radians(15), math.radians(4)
    f = tsf_profile(shape, eps, ov)
    for phi in np.linspace(0, eps, 50):
        total = f(phi) + f(phi + eps) + f(phi + 2 * eps) + f(phi - eps)
        assert total == pytest.approx(1.0, abs=1e-9)


def test_pwm_switching_frequency(dd):
    dr = dd.drive
    r = dr.simulate(200, math.radians(4), math.radians(20), 15.0, control="pwm", f_pwm=16e3)
    # inside the chopping window the devices switch at most twice per PWM period
    assert r.switching_frequency <= 2 * 16e3 * 1.05
    assert r.i_peak < 15.0 * 1.15


def test_thermal_steady_state_balance(dd):
    tm = ThermalModel(dd.design, ThermalParams(h_conv=25, t_amb=30))
    ss = tm.steady_state(50.0, 15.0, 3.0, 1.0)
    p_total = 50.0 * (1 + 0.00393 * (ss["winding"] - 20)) / (1 + 0.00393 * 80) + 19.0
    assert tm.g_amb * (ss["shell"] - 30) == pytest.approx(p_total, rel=0.01)
    assert ss["winding"] > ss["stator"] > ss["shell"] > 30


def test_dclink_current_division():
    period = 1e-3
    t = np.linspace(0, period, 2000, endpoint=False)
    idc = 10 + 8 * np.sin(2 * np.pi * t / period * 3)
    big = dclink.sharing(idc, period, 1.0, esr=0.0, r_b=0.15, l_b=1e-6)
    assert big["I_cap_rms_A"] == pytest.approx(8 / math.sqrt(2), rel=0.01)
    small = dclink.sharing(idc, period, 1e-9, esr=0.0, r_b=0.15, l_b=1e-6)
    assert small["I_bat_ac_rms_A"] == pytest.approx(8 / math.sqrt(2), rel=0.01)


def test_gear_drag_adds_constant_torque():
    a = motor_torque_demand(EBikeSpec(), 50.0)
    b = motor_torque_demand(EBikeSpec(gear_no_load_loss=3.0), 50.0)
    spec = EBikeSpec()
    assert b - a == pytest.approx(3.0 / spec.motor_speed(spec.v_max))


def test_battery_terminal_voltage_and_sag():
    from srm_ebike.drive_cycle import BatteryModel

    b = BatteryModel()
    assert b.ocv(1.0) == pytest.approx(41.8)
    assert b.ocv(0.0) == pytest.approx(30.0)
    v, i = b.terminal(0.5, 300.0)
    assert v * i == pytest.approx(300.0, rel=1e-9)
    assert v == pytest.approx(b.ocv(0.5) - b.r_internal * i, rel=1e-9)


def test_voltage_lookup_interpolates_linearly():
    from srm_ebike.drive_cycle import EfficiencyLookup, VoltageLookup

    sp, tq = np.array([100.0, 2000.0]), np.array([0.5, 6.0])
    lo = EfficiencyLookup(sp, np.array([6.0, 4.0]), tq, np.full((2, 2), 0.6))
    hi = EfficiencyLookup(sp, np.array([7.0, 5.0]), tq, np.full((2, 2), 0.8))
    vl = VoltageLookup({30.0: lo, 42.0: hi})
    assert vl.eta_at(1000, 2.0, 36.0) == pytest.approx(0.7)
    assert vl.t_max(100, 36.0) == pytest.approx(6.5)
    assert vl.eta_at(1000, 2.0, 50.0) == pytest.approx(0.8)  # clamped


def test_battery_route_matches_fixed_voltage_when_flat():
    from srm_ebike.drive_cycle import BatteryModel, EfficiencyLookup, VoltageLookup, simulate_route, simulate_route_battery

    sp, tq = np.array([100.0, 2500.0]), np.array([0.2, 8.0])
    lk = EfficiencyLookup(sp, np.array([8.0, 8.0]), tq, np.full((2, 2), 0.75))
    vl = VoltageLookup({30.0: lk, 42.0: lk})
    spec = EBikeSpec()
    fixed = simulate_route(spec, lk)
    bat = simulate_route_battery(spec, vl, BatteryModel(), soc0=1.0)
    assert bat["energy_terminal_Wh"] == pytest.approx(fixed["energy_battery_Wh"], rel=1e-6)
    assert bat["energy_chemical_Wh"] > bat["energy_terminal_Wh"]


def test_spm_benchmark_consistency():
    from srm_ebike.pm_benchmark import SPMBenchmark, SPMInputs

    pm = SPMBenchmark(SPMInputs(d_outer=0.135, length=0.036, v_dc_min=30.0, base_speed_rpm=1608.0))
    # no-load line-to-line EMF at base speed close to the design voltage
    w_e = pm.inp.pole_pairs * 1608 * 2 * math.pi / 60
    assert math.sqrt(3) * w_e * pm.psi_m == pytest.approx(30.0, rel=0.08)
    op = pm.operating_point(1608, 1.5, 36.0)
    assert op is not None and 0.5 < op["eta_system"] < op["eta_motor"] < 1.0
    assert op["i_rms"] == pytest.approx(1.5 / pm.k_t / math.sqrt(2))
    assert pm.operating_point(5000, 1.5, 36.0) is None  # beyond the voltage limit
