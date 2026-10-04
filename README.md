# Analysis and Design Considerations for a Switched Reluctance Motor-based Drive System for an Electric Bike

This repository contains a design study and a reproducible Python toolchain for a
**250 W pedelec drive built around a 12/8, 3-phase switched reluctance motor
(SRM)**. The SRM drives the wheel through a planetary hub gear and is fed from a
36 V Li-ion battery through an asymmetric half-bridge converter.

The full write-up is in **[`docs/design_report.md`](docs/design_report.md)**.
The tables it quotes are regenerated into [`results/results.md`](results/results.md)
along with the figures.

An IEEE Transactions-format manuscript (IEEEtran, journal mode) is in
**[`paper/`](paper/)**: [`paper/main.pdf`](paper/main.pdf) is built from `main.tex`
with `make -C paper`, and its figures are regenerated with
`python scripts/make_paper_figures.py`.

![efficiency map](results/fig_efficiency_map.png)

## What the toolchain does

| Step | Module | Content |
|---|---|---|
| Vehicle mission | `srm_ebike/vehicle.py` | Road-load model (rolling, aero, grade, inertia) and motor torque/speed requirements, including the corner-speed choice |
| Topology & sizing | `srm_ebike/geometry.py` | Pole-number comparison, output equation, Lawrenson pole-arc rules, pole/yoke dimensions, turns from volt-seconds, slot fill, resistance, aligned/unaligned inductance |
| Magnetics | `srm_ebike/magnetics.py` | Nonlinear ψ(i, θ) with saturation, closed-form co-energy and torque |
| Drive simulation | `srm_ebike/drive.py` | Asymmetric half bridge (magnetise / freewheel / demagnetise) with device drops, hysteresis soft chopping, single-pulse operation, phase-shifted torque summation |
| Losses | `srm_ebike/losses.py` | Copper, iron (Bertotti, non-sinusoidal factor), MOSFET conduction/switching, mechanical |
| Control & performance | `srm_ebike/performance.py` | Turn-on/turn-off angle optimisation, torque–speed envelope, minimum-loss operating points, efficiency map |
| Design flow | `srm_ebike/design_flow.py` | End-to-end design, converter ratings (devices, DC-link capacitor), thermal check |
| Route | `srm_ebike/drive_cycle.py` | Pedelec ride over a hilly urban route: rider power + motor assist, battery Wh/km and range |

## Quick start

```bash
pip install -r requirements.txt
python -m pytest                       # 16 physics / consistency tests
python scripts/run_analysis.py         # full analysis -> results/  (~3-4 min)
python scripts/run_analysis.py --quick --out /tmp/srm   # coarse smoke run (~40 s)
```

Every design parameter is a dataclass field (`EBikeSpec`, `SizingInputs`,
`SRMDrive`, `ConverterParams`, `IronLossCoefficients`), so you can run trade
studies from a few lines of Python:

```python
from srm_ebike.design_flow import design_drive
from srm_ebike.vehicle import EBikeSpec

dd = design_drive(EBikeSpec(gear_ratio=11), elec_loading_peak=40e3, air_gap=0.25e-3)
print(dd.design.summary(), dd.drive.i_max)
```

## Baseline result at a glance

| | |
|---|---|
| Machine | 12/8 SRM, Ø135 mm × 36 mm stack, 0.3 mm air gap, 41 turns/pole, ≈3.4 kg active mass |
| Gear | 8 : 1 planetary; 25 km/h = 1608 rpm |
| Torque | 5.9 N·m launch requirement (8 % grade, 0.5 m/s²) met at ≈31 A phase current |
| Converter | 3 × asymmetric half bridge, 100 V MOSFETs, ≈31 A chopping limit |
| Efficiency | ≈ 76 % motor and ≈ 73 % battery-to-shaft at the 250 W rated point |
| Route | See `results/results.md` for Wh/km and estimated range |

These results come from analytical and lumped models. The modelling
assumptions and their limits are listed in section 9 of the report. Before
building hardware, check the magnetic design with 2-D FEA.
