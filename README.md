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
- `I` not identified: Missed untested cases (rate governed by `disI` and
  `mortI`, as complement of `gamma`)
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
├── Model1_ODE_Model_Run.qmd           # Main analysis document
├── render_run.R                       # Utility script for labelled runs
├── Scripts/
│   ├── Model1_ODE_Parameters.R        # Parameters and scenario builder
│   ├── Model1_ODE_Functions.R         # ODE, post-processing and plot functions
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

- `validate_run_settings()`: checks model time and time step inputs,
- `Model1A_ODE(time, state, parms)`: ODE system for use with
  `deSolve::ode()`,
- `run_scenario(parms)`: convenience wrapper around `deSolve::ode()`,
- `add_quantities()`: per-step and cumulative flow calculations,
- `calculate_total_patient_days()`: integrates bed occupancy over the
  run,
- `negative_value_rows()`: flags any rows with negative state or flow
  values,
- `find_stabilisation_day()`, `stabilisation_summary()`: detect when
  model states reach steady state,
- `make_annual_summary()`: annual totals of all flow quantities,
- `make_final_outcomes()`: cumulative outcomes at end of run, including
  post-hoc attribution of hospital-acquired CDI and confirmed cases by
  test type,
- `make_research_outcomes_table()`: derives and formats all research
  question outcomes (counts, proportions, differences) into a single
  wide table,
- `make_final_outcomes_vs_base()`: absolute differences vs a chosen base
  scenario,
- `make_plot_long_states()`, `make_plot_long_bed_occupancy()`,
  `make_plot_long_flows()`: reshape outputs to long format for plotting,
- `model_colour_palette`, `model_pretty_labels`: named colour and label
  vectors used consistently across plots,
- plot functions: `plot_states_line()`, `plot_states_line_noS()`,
  `plot_states_area()`, `plot_states_area_noS()`,
  `plot_bed_occupancy()`, `plot_cumulative_bed_occupancy()`,
  `plot_flows_to_C_and_I()`, `plot_case_ascertainment_line()`,
  `plot_case_ascertainment_bar()`, `plot_testing_bar()`,
  `plot_case_confirmation_bar()`, `plot_missed_cases_bar()`,
  `plot_deaths_by_state_bar()`, `plot_deaths_by_state_area_noS()`,
  `plot_discharges_by_state_bar()`,
  `plot_discharges_by_state_area_noS()`, `plot_case_outcomes_bar()`,
  `plot_scenario_params()`.

### `Scripts/Model1_ODE_OAT_parameter_testing.R`

A standalone model verification script — **not for analysis**. Runs a
systematic one-at-a-time (OAT) sensitivity test to confirm that the ODE
model responds in the expected direction when each parameter is varied
individually, with all others held at their defaults.

**What it tests:**

For each of 22 parameters (fixed and scenario-level), the model is run
at two values under scenario S1:

- `low`: 0 for proportions and transmission parameters valid at 0; a
  small positive value for time and LOS parameters that must be \> 0;
  1.0 for `abx_mult` (no antibiotic effect).
- `high`: 1 for proportions, 5 × baseline for all other parameters.

A `baseline` column is included in the specification for documentation
only; model runs are performed at `low` and `high` only.

**Outputs assessed (high vs low):**

| Output               | Description                                 |
|----------------------|---------------------------------------------|
| `ss_I`               | Steady-state infected compartment occupancy |
| `yr_confirmed_cases` | Annual confirmed cases (last complete year) |
| `yr_missed_cases`    | Annual missed cases (last complete year)    |

Each output is given an expected direction of change (`"increases"`,
`"decreases"`, or `"complex"`). A parameter passes if the observed
direction matches the expected direction. Parameters marked `"complex"`
have theoretically ambiguous directions due to competing mechanisms and
are flagged for manual inspection rather than auto-assessed.

**Files produced:**

- `Outputs/oat_specs.csv` — test specification table (one row per
  parameter, with low/baseline/high values, expected directions, and
  rationale notes).
- `Outputs/combined_oat_plot.png` — four-panel figure comprising:
  tornado charts of percentage change (high vs low) for each of the
  three assessed outputs, and a pass/fail heatmap across all parameters
  and output checks.

### `Model1_ODE_Model_Run.qmd`

Main Quarto document. Runs the end-to-end workflow:

- setup and sourcing of scripts,
- scenario definition (edit values here for each run),
- model execution for all scenarios,
- stabilisation checks and annual summaries,
- negative-value diagnostics,
- final outcomes table (including differences vs base scenario),
- research question outcomes table (Section 11),
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
  - `cowplot`
  - `quarto`

Install if needed:

``` r
install.packages(c("tidyverse", "deSolve", "scales", "knitr", "cowplot", "quarto"))
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

## Model settings

- **Model run time** and **time step** are set at the top of **Section
  4** in the QMD. The default is 1 year (`365 * 1` days) at a half-day
  time step (`time_step = 0.5`).
- All model time intervals (e.g. `t_test_turnaround`, length-of-stay)
  and transition rates use the same day-based timescale.
- Annual summaries use 365-day intervals.
- Stabilisation is assessed automatically; plots show data up to the
  longest stabilisation period across scenarios.

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

## Outcomes of interest

All outcomes are derived from cumulative totals over the full model run
via `make_research_outcomes_table()` in
`Scripts/Model1_ODE_Functions.R`. The table below sets out how each
outcome is calculated.

### 1) Case ascertainment

|   | Outcome | Numerator | Denominator |
|------------------|------------------|------------------|-------------------|
| 1a | Total incident CDI cases | `cum_total_infected_incidence` | — |
| 1b | Confirmed | `cum_confirmed_cases` | `cum_total_infected_incidence` |
| 1c | Missed | `cum_total_missed_cases` | `cum_total_infected_incidence` |
| 1d | Confirmed via retest (delayed diagnosis)^i^ | `cum_confirmed_via_retest`\^ | `cum_total_infected_incidence` |

**Note**

The proportions for outcomes 1b and 1c do not sum to 1. The gap accounts
for two categories of patient that are neither confirmed nor missed
under the model definitions: (i) patients who die in the **Xtest** state
— they were tested but died before their result was acted on, so are
classified as neither confirmed nor missed; and (ii) patients still in
hospital at the end of the model run whose CDI journey is unresolved.

^i^`cum_confirmed_via_retest` is estimated post-hoc by apportioning
`cum_confirmed_cases` in proportion to the share of all tests that were
retests (`cum_retest_false_negatives / cum_total_test`). This assumes
the confirmation rate (`theta`) is equal for first tests and retests.

### 2) Transmission

|   | Outcome | Calculation |
|-------------------|----------------------|-------------------------------|
| 2a | Hospital-acquired CDI, n and % of all CDI | `cum_hospital_attributable_I`^ii^ |
| 2b | Hospital-acquired colonisations, n and % of all colonisations | `cum_colonised_in_hospital`; `prop_C_hospital_attributable` = `cum_colonised_in_hospital / cum_total_colonised_incidence` |
| 2c | Colonisation prevalence at admission, % | `cum_admitted_colonised / cum_total_admissions` |
| 2d | Colonisation prevalence at discharge, % | `cum_discharges_colonised / cum_total_discharges` |
| 2e | CDI prevalence at admission, % | `cum_admitted_infected / cum_total_admissions` |
| 2f | CDI prevalence at discharge, % | `cum_discharges_while_infected / cum_total_discharges` |
| 2g | Difference in CDI prevalence at discharge versus admission^i^, % | CDI prevalence at discharge - CDI prevalence at admission |

**Note**

^i^ A positive value indicates that the hospital is a net source of
colonisation or infection.

^ii^ `cum_hospital_attributable_I` is estimated post-hoc by attributing
`cum_infected_in_hospital` (C→I progressions within hospital) in
proportion to the share of all colonisation that is hospital-acquired
(`cum_colonised_in_hospital / cum_total_colonised_incidence`). See
`make_final_outcomes()` in `Scripts/Model1_ODE_Functions.R`.

### 3) Hospital resource use

|   | Outcome | Calculation |
|------------------------|------------------------|------------------------|
| 3a | Total hospital bed-days | `hospital_bed_days` from `calculate_total_patient_days()` |
| 3b | Average length of stay (days) | `hospital_bed_days / cum_total_admissions` |
| 3c | Side-room bed-days | `side_room_bed_days` from `calculate_total_patient_days()` |
| 3d | Side-room bed-days as % of total | `side_room_bed_days / hospital_bed_days` |
| 3e | Faecal specimens collected and tested (n) | `cum_total_test` |

### 4) Hospital patient mortality

|   | Outcome | Numerator | Denominator |
|------------------|------------------|------------------|------------------|
| 4a | Total hospital deaths (n) | `cum_total_deaths` | — |
| 4b | Deaths following CDI (n and %) | `cum_deaths_post_infection` | `cum_total_deaths` |
| 4c | Deaths among missed CDI cases^i^ (n and %) | `cum_deaths_infected + cum_deaths_FN` | `cum_total_deaths` |

**Note**

^i^Deaths among missed CDI cases covers patients who died without ever
receiving a CDI diagnosis: untested patients (died in **I**) and
false-negative patients (died in **FN**). Deaths in the **Xtest** state
are excluded as those patients had been tested, even though no diagnosis
was confirmed.

------------------------------------------------------------------------

## 11) Research question outcomes

Derive and display all outcomes of interest, as defined in the model
research questions. Proportions labelled `[1]` do not sum to 1 across
case ascertainment categories: the gap reflects deaths in Xtest
(patients tested but who died before their result was acted on) and
cases still unresolved in hospital at the end of the model run.
Differences labelled `[2]` are discharge minus admission prevalence in
percentage points (pp); a positive value indicates the hospital is a net
source.

```{r}
research_outcomes <- make_research_outcomes_table(final_outcomes, patient_days_summary)
write_csv(research_outcomes, file.path(output_dir, "research_outcomes.csv"))
knitr::kable(research_outcomes, caption = "Outcomes of interest by scenario")
```

------------------------------------------------------------------------

## 12) Plots

Create a series of plots to visualise the model states and flows over
time, as well as the final outcomes. The plotting functions take the
long-format data frames of states and flows, as well as the final
outcomes, and generate ggplot objects. Save each plot as a PNG file in
the output directory. All outputs are written to `Outputs/<run_label>/`
(created automatically), including:

- `model_parameters_all_scenarios.csv`
- `solution_S1.csv`, `solution_S2.csv`, etc. (full time-series by
  scenario)
- `annual_summary_all.csv`
- `stabilisation_summary_states_flows.csv`
- `negative_rows_*.csv`
- `final_outcomes_all.csv`
- `final_outcomes_vs_base.csv`
- `research_outcomes.csv`
- `scenario_control_parameters.csv`
- Plot PNG files
