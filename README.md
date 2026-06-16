---
editor_options: 
  markdown: 
    wrap: 72
---

# Model 1: ODE model of hospital transmission with testing and isolation

## Overview

This repository contains a deterministic compartmental ODE model for
in-hospital transmission of a pathogen, with explicit representation of:

- suspected case identification and isolation,
- test sensitivity and turnaround time,
- false negatives returning to general care,
- treated confirmed cases in isolation,
- discharge and mortality flows,
- cumulative missed cases and side-room occupancy.

The model is implemented in R with `deSolve` and run through a Quarto
workflow.

------------------------------------------------------------------------

## Model structure

### ![](Images/260422%20Model%201%20diagram/Slide1.PNG)

### State variables

- `S`: Susceptible (not colonised or infected)
- `C`: Colonised (carrying pathogen but not infected)
- `I`: Infected (symptomatic, on general ward, not yet identified)
- `Xtest`: Suspected cases either isolated (in a side room) or not
  isolated (on a general ward), awaiting test result
- `FN`: False negatives — symptomatic patients who received a false
  negative test result, returned to the general ward (not isolated),
  awaiting re-testing after `t_wait_retest` days
- `Xtreat`: Confirmed cases isolated in a side room, on treatment,
  awaiting recovery and return to a general ward or discharge
- `D`: Cumulative deaths (absorbing count state)

Total in-hospital living census is:

$$N = S + C + I + Xtest + FN + Xtreat$$

------------------------------------------------------------------------

## Key transitions

- `S -> C`: Force of colonisation (`lambda`)
- `C -> I`: Progression to infection (`alpha`)
- `C -> S`: Decolonisation (`mu`)
- `I -> Xtest`: Suspected-case identification (`gamma`)
- `I` not identified: Missed untested cases (`gamma0`)
- `Xtest -> Xtreat`: Test positive, case confirmed (`theta`)
- `Xtest -> FN`: False negative, return to general ward (`pi`)
- `FN -> Xtest`: False negatives to be re-tested (`sigma`)
- `Xtreat -> S`: Treatment and recovery to general ward (`delta`)
- Discharge and mortality from all relevant compartments
- Deaths accumulate in `D`

------------------------------------------------------------------------

## Case accounting outputs

The ODE function tracks the following diagnostic rates at each time
step, in addition to state variables:

- **Colonisation**: in-hospital colonisations, colonised admissions
- **Infection**: in-hospital infections, infected admissions
- **Case identification**: suspected cases (`gamma * I`), confirmed
  cases (`theta * Xtest`), false negatives (`pi * Xtest`), missed
  untested (`gamma0 * I`)
- **Outcomes**: recovered cases, deaths by state, side-room occupancy
- **Admissions and discharges** by state

Per-step, cumulative, and annual summary versions of all flow quantities
are computed in post-processing.

------------------------------------------------------------------------

## File structure

```         
R_Cdiff_Model/
├── Model1_ODE_Model_Run.qmd          # Main analysis document
├── render_run.R                       # Utility script for labelled runs
├── Scripts/
│   ├── Model1_ODE_Parameters.R       # Parameters and scenario builder
│   ├── Model1_ODE_Functions.R        # ODE, post-processing and plot functions
│   └── Model1_ODE_OAT_parameter_testing.R  # One-at-a-time sensitivity testing
└── Outputs/
    └── <run_label>/                   # One folder per labelled run
```

### `Scripts/Model1_ODE_Parameters.R`

Defines:

- fixed model parameters,
- helper `prob2rate()` conversion,
- validation helpers: `validate_fixed_parameters()`,
  `validate_scenario_inputs()`, `validate_rate_parameters()`,
- scenario builder `make_mod_parms(...)`.

### `Scripts/Model1_ODE_Functions.R`

Defines:

- `Model1A_ODE(time, state, parms)` for use with `deSolve::ode()`,
- `run_scenario(parms)`: convenience wrapper around `deSolve::ode()`,
- non-negative clamping of state variables,
- `add_flow_quantities()`: per-step and cumulative flow calculations,
- `make_annual_summary()`, `stabilisation_summary()`,
  `negative_value_rows()`, `make_final_outcomes()`,
  `make_final_outcomes_vs_base()`,
- long-format data prep functions for plotting,
- colour palette (`model_colour_palette`) and labels
  (`model_pretty_labels`),
- plot functions: `plot_states_line()`, `plot_states_area()`,
  `plot_flows_to_C_and_I()`, `plot_flows_cases()`,
  `plot_missed_cases_bars()`, `plot_side_room()`, `plot_final_bars()`,
  `plot_scenario_params()`.

### `Model1_ODE_Model_Run.qmd`

Main Quarto document. Runs the end-to-end workflow:

- setup and sourcing of scripts,
- scenario definition (edit values here for each run),
- model execution for all scenarios,
- stabilisation checks and annual summaries,
- negative-value diagnostics,
- final outcomes table (including differences vs base scenario),
- plots and CSV export to `Outputs/<run_label>/`.

Accepts a `run_label` parameter (see **Run instructions** below).

### `render_run.R`

Utility script for running the model with a labelled output folder.
Source this file by calling `source("render_run.R")` and call
`render_run("A")` from the R console.

------------------------------------------------------------------------

## Requirements

- R (\>= 4.1 recommended)
- Packages:
  - `tidyverse`
  - `deSolve`
  - `scales`
  - `knitr`
  - `quarto`

Install if needed:

``` r
install.packages(c("tidyverse", "deSolve", "scales", "knitr", "quarto"))
```

Quarto must also be installed and available — see <https://quarto.org>.

------------------------------------------------------------------------

## Run instructions

### Standard workflow (recommended)

1.  Open `Model1_ODE_Model_Run.qmd` and set the scenario parameter
    values in **Section 3**. Save the qmd (no need to change the
    filename).
2.  In the R console, from the project root:

``` r
source("render_run.R")
render_run("A")
```

This will:

- execute the QMD in a fresh R session,
- save all CSVs and plots to `Outputs/A/`,
- produce `Model1_ODE_Model_Run_A.html` in the project root.

To run a second set of scenarios, edit the scenario values again and
call `render_run("B")`, and so on. The `run_label` in the QMD YAML does
not need to be changed — `render_run()` overrides it.

### Interactive rendering (quick checks)

You can also click **Render** in RStudio, which uses the default
`run_label: "default"` from the YAML and saves outputs to
`Outputs/default/`. This is useful for quick checks but does not produce
a labelled output folder.

------------------------------------------------------------------------

## Outputs

All outputs are written to `Outputs/<run_label>/` (created
automatically), including:

- `model_parameters_all_scenarios.csv`
- `solution_S1.csv`, `solution_S2.csv`, etc. (full time-series by
  scenario)
- `annual_summary_all.csv`
- `stabilisation_summary_states_flows.csv`
- `negative_rows_*.csv`
- `final_outcomes_all.csv`
- `final_outcomes_vs_base.csv`
- `scenario_control_parameters.csv`
- Plot PNG files

------------------------------------------------------------------------

## Scenario definitions

Scenario parameter values are set in **Section 3** of
`Model1_ODE_Model_Run.qmd`. The four control variables are:

- `isolate_before_test`: whether suspected cases are isolated before the
  test result is available (`TRUE`/`FALSE`)
- `prop_I_suspected`: proportion of infected patients identified as
  suspected cases (0--1)
- `test_sens`: diagnostic test sensitivity (0--1)
- `t_test_turnaround`: days from test administration to result (and
  transfer to side room if `isolate_before_test = TRUE`)

------------------------------------------------------------------------

## Validation and safeguards

Implemented safeguards include:

- explicit `stop()` checks for invalid inputs and parameters,
- warnings for potentially unrealistic values,
- non-negative state clamping in the ODE function,
- `N <- max(N, 1e-12)` to prevent divide-by-zero,
- non-negative guard for `missed_untested` via `pmax`.

------------------------------------------------------------------------

## Model settings

- **Model run time** and **time step** are set at the top of **Section
  4** in the QMD. The default is 1 year (`365 * 1` days) at a half-day
  time step (`time_step = 0.5`).
- All model time intervals (e.g. `t_test_turnaround`, length-of-stay)
  and transition rates use the same day-based timescale.
- Annual summaries use 365-day intervals.
- Stabilisation is assessed automatically; plots show data up to the
  longest stabilisation period across scenarios.
