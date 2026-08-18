# Solid-Rocker-Motor-Simulator

A MATLAB simulator for solid rocket motor (SRM) internal ballistics. It compares grain
geometries and propellant formulations, integrates the coupled chamber-pressure /
web-regression ODE with `ode45`, runs parametric trade studies, and validates output
against reference thrust-curve data.
![Internal ballistics dashboard — thrust, chamber pressure, web regression, and burn-rate phase portrait for the BATES/APCP reference design](docs/fig4.pdf)

> _Figure 4 — full internal ballistics dashboard. Regenerate with `SRM_SIMULATOR`._

---

## What it does

- **Three grain geometry engines** — BATES (cylindrical core), star/finocyl, and
  wagon-wheel — each computing burn area `Ab(w)`, propellant volume `Vp(w)`, and web
  regression.
- **Propellant thermochemistry database** — APCP (two blends), KNSB, KNSU, KNDX, with
  Saint-Robert's burn rate law `r = a·Pc^n` and temperature sensitivity bands.
- **Isentropic nozzle model** — exit Mach (bisection solver), thrust coefficient `Cf`,
  and characteristic velocity `c*`.
- **Coupled internal ballistics ODE** — chamber pressure and web regression integrated
  simultaneously with `ode45`, with automatic burnout termination and NAR/TRA motor
  classification.
- **Parametric trade studies** — grain-vs-grain, propellant-vs-propellant, and
  single-parameter sweeps (e.g. core diameter).
- **Validation module** — compares simulated thrust against `.eng` reference curves and
  reports R², RMSE, and total-impulse error.

The master script generates **8 figures** spanning all four phases in a single run.

---

## Quick start

Requires MATLAB (developed on R2023+). No toolboxes required — the nozzle exit-Mach
solve uses a built-in bisection loop rather than `fsolve`.

```matlab
% From the project folder:
addpath(pwd)
SRM_SIMULATOR      % runs all phases, produces all 8 figures + summary table
```

To run phases individually:

```matlab
phase1_run         % grain geometry engine only
phase234_run       % propellant model, ballistics ODE, trade studies
report_gen         % saves all figures + a text summary to a timestamped folder
```

---

## Repository structure

| File | Role |
|------|------|
| `SRM_SIMULATOR.m` | Master runner — executes all phases, produces all 8 figures |
| `BATES.m` | BATES cylindrical-core grain geometry |
| `star_grain.m` | Star / finocyl grain geometry (Shoelace port-area calc) |
| `wagon_wheel.m` | Wagon-wheel (hub + spokes) grain geometry |
| `grain_viz.m` | Grain cross-section and burn-profile plotting |
| `propellant_db.m` | Propellant thermochemistry lookup table |
| `burn_rate.m` | Saint-Robert's law + temperature sensitivity |
| `nozzle_model.m` | Isentropic nozzle flow, `Cf`, exit-Mach bisection solver |
| `ballistics_ode.m` | Core two-state internal ballistics ODE (`ode45`) |
| `trade_study.m` | Parametric sweeps and comparison studies |
| `validation.m` | `.eng` reference comparison and error metrics |
| `report_gen.m` | Batch figure export + performance summary |
| `phase1_run.m`, `phase234_run.m` | Per-phase runner scripts |

---

## The physics

The simulator integrates a two-state system — chamber pressure and web regression —
using a lumped-parameter mass balance on the chamber gas:

```
dPc/dt = (Rg·Tc / Vc) · (ρp·Ab·r − Pc·At / c*)
dw/dt  = r = a·Pc^n
```

where `Pc` is chamber pressure, `w` is web regression distance, `Ab(w)` is the
web-dependent burn area, `Vc(w)` is chamber free volume, `At` is throat area, `ρp` is
propellant density, `Rg·Tc` is the specific gas constant times flame temperature, and
`c*` is characteristic velocity.

The nonlinearity comes from the coupling: burn rate depends on pressure, and pressure
depends on burn area, which depends on the web the burn rate is advancing. The system
is mildly stiff at ignition (pressure rises on the order of tens of MPa/s in the first
few milliseconds), so `ode45` adapts its step through the transient and a terminal
event stops integration once the grain burns out and pressure decays to ambient.

---

## Output figures

1. Grain geometry comparison — cross-sections + burn area, volume, and profile classification
2. BATES core-diameter sensitivity sweep
3. Propellant characterization — burn rate (with temperature bands), log-log slope, nozzle `Cf`
4. Full ballistics dashboard — thrust, pressure, web, burn-rate phase portrait, summary
5. APCP vs KNSB head-to-head on a fixed grain
6. Grain geometry trade study — impulse, pressure, Isp, burn duration
7. Parametric sweep — how core diameter drives performance
8. Validation — simulated vs reference thrust curve with residuals and R²

---

## Validation

The validation module parses standard RASP `.eng` thrust-curve files (the format used
by AeroTech, Cesaroni, and the amateur rocketry community) and computes R², RMSE, peak-
thrust error, and total-impulse error against the simulation.

---
