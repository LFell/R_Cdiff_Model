---
editor_options: 
  markdown: 
    wrap: 72
---

# Model 1: ODE model of hospital transmission with testing and isolation

## Overview

This repository contains a deterministic compartmental ODE model for
in-hospital transmission of a pathogen, with explicit representation of:

- suspected case identification and optional isolation,
- test sensitivity and turnaround time,
- false negatives returning to general care,
- treated confirmed cases, optionally isolated,
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
- `Xtreat`: Confirmed cases on treatment, either isolated (in a side
  room) or not isolated (on a general ward), awaiting recovery and
  return to a general ward or discharge
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

## Model assumptions

Default parameter values are given in brackets; they are set in `Scripts/Model1_ODE_Parameters.R`.

*Pathogen and population*

- A single toxigenic strain is in circulation.
- The hospital is a single, well-mixed population: there are no separate wards or bays. The model is deterministic and its parameters don't change over time (no seasonality).
- Patients are not followed after discharge. Readmission is not modelled.

*Admissions and discharges*

- Admissions arrive at a constant rate (133/day). They do not replace discharges, so the hospital census can vary.
- Admissions are split between S, C and I in fixed proportions reflecting community prevalence (8.6% colonised, 0.2% infected, the rest susceptible). Community prevalence is fixed and isn't affected by discharges from the hospital.
- Nobody is admitted directly to Xtest or Xtreat.
- Each state has its own constant discharge rate, set by its length of stay:
  - S and C: 7 days;
  - I and FN: 14 days;
  - Xtreat: 20% of treated patients are discharged and 80% return to a general ward as susceptible, over 14 days (plus 0.5 days' transfer time when isolated).
- Patients awaiting a test result (Xtest) stay in hospital until the result arrives; they can't be discharged from Xtest.
- Death rates are constant:
  - S and C share a baseline rate (0.1% risk over a 7-day stay);
  - I, Xtest and FN die at 3× the baseline (`mortI_mult`);
  - Xtreat dies at 1.5× the baseline (`mortXtreat_mult`).

*Colonisation and progression*

- Susceptible patients become colonised through:
  - background acquisition from spores in the environment, at a constant rate that doesn't depend on the number of patients in each state (β₀ = 0.002);
  - frequency-dependent transmission from colonised and infected patients (λ = β₀ + Σβₖk/N).
- Colonised and infected patients all pose a transmission risk:
  - colonised (C): 0.008;
  - infected (I) and false negatives (FN): 0.04;
  - suspected cases (Xtest): 0.04 on a general ward;
  - confirmed cases on treatment (Xtreat): 0.02 on a general ward.
- Treatment is assumed to reduce shedding from 0.04 at the start to 0 at the end. The model uses the average over the treatment period (0.02) as a constant rate.
- Isolation in a side room reduces transmission to 5% of the general-ward value. It does not eliminate it.
- Progression from C to I is increased by antibiotic exposure among colonised patients (SCI structure, not PSCIX): the progression rate rises by a relative risk of 2 in the 37.3% of colonised patients receiving antibiotics. Without antibiotics, 50% of colonised patients progress, with a mean time to progression of 5 days.
- Colonised patients can lose colonisation without treatment (5% over 5 days).
- There is no immunity, and no recurrence or relapse. Patients who recover after treatment return to the susceptible state.
- Colonised patients are not screened, detected or decolonised.

*Case identification and testing*

- Only infected (symptomatic) patients are identified as suspected cases and tested. Patients in S and C don't develop diarrhoea that leads to testing. The test is assumed 100% specific, so there are no false positives. As a result, "faecal specimens tested" counts only tests of patients with CDI, and underestimates real testing volume.
- Testing coverage (the proportion of CDI cases tested for CDI) sets the identification rate, with a mean time to identification of 2 days from the onset of symptoms. Infected patients can be discharged or die before they're identified, so the share actually tested is lower than the nominal coverage.
- Test results arrive after a fixed turnaround time (2 days). Each test, including retests, returns a positive result with probability equal to the test sensitivity. Results of repeat tests on the same patient are independent.
- Patients with a false-negative result return to a general ward (FN), still symptomatic and transmitting. They are retested after a mean of 7 days unless they are discharged or die first.
- Confirmed cases are treated. Under P2 and P3 they are isolated; under P1 they stay on a general ward.
- Moving a patient between a general ward and a side room takes 0.5 days. This is added to the rate of any transition that involves such a move.
- Side-room capacity is unlimited.

*Outcome definitions*

- **Missed cases:** patients who leave I or FN by discharge or death without a confirmed diagnosis. Deaths while awaiting a result in Xtest count as neither missed nor confirmed.
- **Hospital-onset CDI:** all in-hospital C→I progressions, i.e. patients admitted colonised plus patients colonised in hospital who progress to CDI during their stay.
- **Hospital-acquired CDI:** in-hospital C→I progressions, multiplied by the share of all colonisations that were acquired in hospital (a subset of hospital-onset CDI). This assumes patients colonised in hospital and those admitted colonised progress at the same rate.

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
├── renv.lock, renv/, .Rprofile       # renv package library for this project
├── Scripts/
│   ├── Model1_ODE_Parameters.R        # Parameters and scenario builder
│   ├── Model1_ODE_Functions.R         # ODE, post-processing, scenario grid and plot functions
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

- `validate_run_settings()`: checks model time and time step inputs
  (both positive, time step smaller than model time),
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
  test type, and CDI-attributable deaths,
- `add_derived_outcomes()`: proportions confirmed, missed and confirmed
  via retest, admission and discharge prevalence and their differences,
  and specimens tested per confirmed case,
- `make_research_outcomes_table()`: derives and formats the research
  question outcomes (1a–4d) and ICER inputs into a single wide table,
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
  `plot_bed_days_bar()`, `plot_scenario_params()`.

### `Scripts/Model1_ODE_OAT_parameter_testing.R`

A standalone model verification script — **not for analysis**. Runs a
systematic one-at-a-time (OAT) sensitivity test to confirm that the ODE
model responds in the expected direction when each parameter is varied
individually, with all others held at their defaults.

**What it tests:**

For each of 23 fixed parameters (including `t_test_turnaround`), the
model is run at two values at a single baseline of test settings (50%
testing coverage, 50% test sensitivity, so the untested,
false-negative and retesting pathways are all active) under each of the
three isolation policies (`"none"`, `"confirmed"`,
`"suspected_and_confirmed"`). Testing coverage and test sensitivity
are not tested here, as the scenario grid in the QMD covers their full
range.

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
are flagged for manual inspection rather than auto-assessed. The
expectations in `oat_specs` are for `"suspected_and_confirmed"`. Where
the expected direction depends on the isolation policy (e.g.
`t_test_turnaround`, which has little effect when suspected cases are
isolated but increases transmission when they wait on the ward, or
`t_patient_transfer` and `beta_mult_isolated`, which have no effect
under `"none"`), `oat_specs_policy` overrides the expectation for the
named policy.

The % change (high vs low) in the ten cumulative outcomes compared across
the scenario grid (see section 13 of the QMD) is also reported for
inspection, without expected directions.

**Files produced:**

- `Outputs/oat_specs.csv` — test specification table (one row per
  parameter, with low/baseline/high values, expected directions, and
  rationale notes).
- `Outputs/oat_specs_policy_overrides.csv` — expected directions that
  differ under the `"confirmed"` or `"none"` policy.
- `Outputs/oat_results.csv` — low/high values, expected and actual
  directions and pass/fail for each parameter × isolation policy.
- `Outputs/oat_outcomes_pct_change.csv` — low/high values and % change
  in each cumulative outcome for each parameter × isolation policy.
- `Outputs/combined_oat_plot.png` — four-panel figure comprising:
  tornado charts of percentage change (high vs low) for each of the
  three assessed outputs, coloured by isolation policy, and a pass/fail
  heatmap across all parameters and output checks, faceted by isolation
  policy.
- `Outputs/oat_tornado_cumulative_outcomes.png` — tornado charts of
  percentage change in each of the ten cumulative outcomes, coloured by
  isolation policy.

### `Model1_ODE_Model_Run.qmd`

Main Quarto document. Runs the end-to-end workflow:

- setup and sourcing of scripts,
- scenario definition: vectors of control values, crossed into a
  scenario grid (edit values here for each run),
- model execution for every scenario in the grid (cached), plus full
  time series for a few example scenarios,
- negative-value and stabilisation checks for every scenario (Sections 7
  and 8),
- results for the example scenarios that help verify the model is
  working as expected (Sections 9–11):
  - time-series line and area plots, up to stabilisation (Section 9),
  - annual summary tables (Section 10),
  - stacked bar plots of annual summaries (Section 11),
- results for the write-up (Sections 12–14):
  - research question outcomes: the research outcomes table (annual
    totals, including rates per 10,000 patient bed-days) for every
    scenario, and differences between policies (P2 vs P1, P3 vs P1, P3 vs
    P2; Section 12),
  - outcomes by testing coverage and test sensitivity: plots comparing
    the isolation policies across the grid (Section 13),
  - change in outcomes due to increasing test sensitivity, e.g. from 50% to
    80% (Section 14),
- CSV and PNG export to `Outputs/<run_label>/`.

Accepts a `run_label` parameter (see **Run instructions** below).

### `render_run.R`

Utility script for running the model with a labelled output folder.
Source this file by calling `source("render_run.R")` and call
`render_run("A")` from the R console.

------------------------------------------------------------------------

## Requirements

- R 4.5.0 (the version recorded in `renv.lock`)
- Packages: `tidyverse`, `deSolve`, `scales`, `knitr`, `cowplot`,
  `quarto` (plus base R `parallel`), managed with **renv**. Opening the
  `.Rproj`, or starting R in the project root, activates the project
  library via `.Rprofile`. To install the recorded versions:

``` r
renv::restore()
```

  After adding a package, run `renv::snapshot()` to record it.

Quarto must also be installed and available — see <https://quarto.org>.
RStudio bundles a copy; outside RStudio, `quarto::quarto_render()` may
need the `QUARTO_PATH` environment variable set to that `quarto.exe`.

------------------------------------------------------------------------

## Run instructions

### Standard workflow (recommended)

1.  Open `Model1_ODE_Model_Run.qmd` and set the vectors of scenario
    control values in **Section 3**. Save the qmd (no need to change
    the filename).
2.  In the R console, from the project root:

``` r
source("render_run.R")
render_run("A")
```

This will:

- execute the QMD in a fresh R session,
- save all CSVs and plots to `Outputs/A/`,
- produce `YYMMDD_Model1_ODE_Model_Run_A.html` in the project root,
  prefixed with today's date (e.g. `261008_Model1_ODE_Model_Run_A.html`),
  so reports from earlier days are kept for comparison. Rendering the
  same label again on the same day overwrites that day's file.

The HTML is self-contained (`embed-resources: true` in the QMD YAML): the
plots are embedded in the file, so each dated report keeps its own plots
even after `Outputs/A/` is overwritten by a later run, and the HTML can be
moved or shared on its own (about 10 MB).

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

- **Model run time**, **time step** and **burn-in** are set at the top
  of **Section 4** in the QMD. The default is a 6-year run (`365 * 6`
  days) at a half-day time step (`time_step = 0.5`), with a 1-year
  burn-in (`burn_in = 365`), giving a 5-year analysis period.
- The burn-in is the period while the model moves from its initial
  conditions to steady state. It is the same for every scenario and is
  excluded from all cumulative totals, patient-days and annual summaries,
  so these cover the analysis period only. `add_quantities()` gives rows
  in the burn-in (including time 0) zero width and multiplies every other
  output row by the time step. Sections 10–14 of the QMD report annual
  totals (the analysis-period total divided by its length in years), for
  comparison with surveillance data.
- All model time intervals (e.g. `t_test_turnaround`, length-of-stay)
  and transition rates use the same day-based timescale.
- Annual summaries use 365-day intervals counted from the end of the
  burn-in (year 1 is the first year of the analysis period).
- Stabilisation is assessed automatically for every scenario in the grid
  (states changing by ≤ 1e-4 per time step for 10 consecutive time steps;
  the `consecutive_days` setting counts time steps, not days, so this is
  5 days at the half-day step). Section 8 checks that the burn-in is
  longer than every scenario's stabilisation period (with the defaults,
  all stabilise within 250 days) and stops the render with an error
  listing any that are not, so `burn_in` can be increased. Time-series plots show data up to
  the longest stabilisation period across the example scenarios, i.e.
  the approach to steady state during the burn-in.
- The OAT script uses the same 1-year burn-in on its 3-year runs.

------------------------------------------------------------------------

## Scenario definitions

Scenario control values are set as vectors in **Section 3** of
`Model1_ODE_Model_Run.qmd`, and every combination is run (the scenario
grid; by default 3 policies × 11 coverage values × 11 sensitivity values
= 363 scenarios, about 3–4 minutes on one core, about 30 seconds in
parallel). The three control variables are:

- `isolation_policy`: which patients are isolated in a side room,
  numbered from least to most isolation (set by `isolation_policy_levels`
  in `Scripts/Model1_ODE_Functions.R`):
  - P1 `"none"`: neither suspected nor confirmed cases (the reference
    policy)
  - P2 `"confirmed"`: confirmed cases (`Xtreat`) only
  - P3 `"suspected_and_confirmed"`: suspected cases awaiting a test
    result (`Xtest`) and confirmed cases (`Xtreat`), as recommended in
    best-practice guidance

  Differences between policies (policy minus comparator, negative =
  fewer) are reported for P2 vs P1, P3 vs P1 and P3 vs P2 (the
  incremental effect of also isolating suspected cases), set in
  `policy_comparisons`.
- `prop_I_suspected` (testing coverage): proportion of CDI cases tested
  for CDI (identified as suspected cases; 0--1). The model uses it to
  calculate the case identification rate `gamma = prop_I_suspected /
  t_identify_suspected` (plus `t_patient_transfer` under P3). As patients
  can be discharged or die from I first, the share of infected patients
  actually tested is lower (with the default parameters about 88% at 100%
  and 78% at 50%; 85% and 74% under P3). Plots and tables label it
  "Proportion of CDI cases tested for CDI".
- `test_sens`: diagnostic test sensitivity (0--1). The model uses it to
  calculate the confirmed case rate (`theta = test_sens /
  t_test_turnaround`) and false negative rate (`pi = (1 - test_sens) /
  t_test_turnaround`), adding `t_patient_transfer` where the patient moves
  between a general ward and an isolation bed (to `theta` under P2, to
  `pi` under P3). As transfer time is added to only one rate, the share
  of results that are positive differs slightly from `test_sens` under P2
  and P3 (at 50% sensitivity: 50% P1, 44% P2, 56% P3).

Test turnaround time (`t_test_turnaround`, days from test administration
to result) is a fixed parameter in `Scripts/Model1_ODE_Parameters.R`.

Scenario labels combine a policy code with the coverage and sensitivity
as percentages, e.g. `P2_confirmed_cov050_sens100` (`P1_none`,
`P2_confirmed`, `P3_susp_conf`). `example_settings` in Section 3 lists
the coverage and sensitivity settings for the example scenarios, in the
order shown in Sections 9 and 11 (by default 0%/0%, 50%/50%, 100% coverage/50%
sensitivity, 50% coverage/100% sensitivity and 100%/100%), each run
under all three policies (15 scenarios). They are used for time-series
plots, annual summaries and the outcome tables. At 0% coverage nobody is
tested, so the three policies give the same results (a no-testing
baseline).

In Section 9, the time-series line and area plots show the policies as columns and
the settings as rows; in Section 11, the stacked bar plots group the bars by setting, with P1,
P2 and P3 side by side and the setting label beneath each group.

Grid results are cached in `Outputs/<run_label>/grid_runs.rds` and
re-used on the next render if the grid, fixed parameters, run settings,
initial conditions and model/post-processing functions are unchanged.
Set `grid_clear_cache <- TRUE` in Section 6 to force a re-run.

The grid runs in parallel (base R `parallel`) on `grid_n_cores` cores,
set in Section 6 (default: all physical cores but one; about 30 seconds
on 15 cores). Set `grid_n_cores <- 1` to run scenarios one at a time.

`isolation_policy` determines:

- the transmission parameters for `Xtest` and `Xtreat`: on a general
  ward, suspected cases (`Xtest`) transmit at `betaI` and confirmed cases
  on treatment (`Xtreat`) at the fixed `betaXtreat` (default 0.02: the
  average across the treatment period, as treatment is assumed to reduce
  shedding from `betaI` (0.04) at the start to 0 at the end); patients
  in a side room transmit at their general-ward value ×
  `beta_mult_isolated` (fixed parameter, default 0.05);
- where patient transfer time (`t_patient_transfer`) is added. It is
  added to any transition that involves a move between a general ward
  and a side room: I → Xtest (`gamma`), Xtest → FN (`pi`) and
  FN → Xtest (`sigma`) under `"suspected_and_confirmed"`; Xtest → Xtreat
  (`theta`) under `"confirmed"`; and Xtreat → S (`delta`) under both
  `"confirmed"` and `"suspected_and_confirmed"`;
- side-room bed occupancy: none under `"none"`, `Xtreat` under
  `"confirmed"`, and `Xtest + Xtreat` under `"suspected_and_confirmed"`.

Section 13 of the QMD (outcomes by testing coverage and test sensitivity) compares the three isolation policies across
testing coverage and test sensitivity. The outcomes are the research
question outcomes defined in `grid_outcome_spec` in
`Scripts/Model1_ODE_Functions.R` (numbered as in the research questions
below):

| Outcome | Variable | Format |
|---|---|---|
| 1a) Incident CDI cases | `cum_total_infected_incidence` | count per year |
| 1b) Missed CDI cases | `cum_total_missed_cases` | count per year |
| 1c) Proportion of CDI cases confirmed | `prop_confirmed` | proportion |
| 1d) Proportion confirmed via retest (estimated) | `prop_confirmed_via_retest` | proportion |
| 2a) Hospital-acquired CDI cases (estimated) | `cum_hospital_attributable_I` | count per year |
| 2b) Discharge minus admission prevalence of colonisation | `diff_prev_C` | percentage points |
| 2b) Discharge minus admission prevalence of CDI | `diff_prev_CDI` | percentage points |
| 3a) Total patient bed-days | `hospital_bed_days` | count per year |
| 3c) Side-room patient bed-days | `side_room_bed_days` | count per year |
| 3e) Faecal specimens tested | `cum_total_test` | count per year |
| 3f) Faecal specimens tested per confirmed case | `tests_per_confirmed` | rate |

Counts are annual totals over the analysis period; proportions,
prevalence differences and specimens per confirmed case are over the
whole analysis period (derived in `add_derived_outcomes()`). Differences
in total CDI cases (1a) are driven entirely by differences in
hospital-acquired CDI cases (2a), so their difference plots are
identical: with fixed admissions, patients admitted with CDI and
progressions among patients admitted colonised are the same in every
scenario (testing and isolation act only on infected patients). The
rates per 10,000 patient bed-days, side-room and missed proportions and
missed cases by origin are only tabulated (`rate_spec`;
`grid_rates.csv`). Hospital-onset CDI is all in-hospital C→I
progressions (patients admitted colonised plus patients colonised in
hospital); hospital-acquired CDI is the part attributed to
hospital-acquired colonisation. Each outcome has its own tab with:

- a line plot of the outcome against test sensitivity, one line for each
  of 0%, 20%, 40%, 60%, 80% and 100% testing coverage
  (`grid_line_coverage` in `Scripts/Model1_ODE_Functions.R`; categorical
  legend, 0% yellow to 100% purple), faceted by policy with a shared y
  axis;
- the same plot for the differences from P1 (P2 vs P1, P3 vs P1) at the
  same coverage and sensitivity, placed under P2 and P3 (the space under
  P1 is left empty);
- heatmaps of test sensitivity (x) by testing coverage (y), for every
  grid value, in one figure: the outcome under each policy on the top row
  (shared colour scale, yellow lowest, purple highest), and the
  differences from P1 (P2 vs P1, P3 vs P1) on the bottom row under P2 and
  P3, with the space under P1 left empty (yellow largest negative
  difference, purple largest positive difference).

The incremental effect of also isolating suspected cases (P3 vs P2) is
saved in `grid_key_outcome_differences.csv`. Each panel has its own y axis so
the rows line up.

At very low test sensitivity (about 10% or below), testing can increase
CDI, hospital-acquired CDI, missed cases and deaths compared with no
testing. Most tested patients get a false negative result, return to
the ward (FN) still transmitting, and are retested sooner than they
would typically be discharged. There is no discharge from Xtest, so they
cycle between FN and Xtest, held back from treatment and isolation and
kept infectious in hospital for longer. This is intended, to show the
effect of a low-sensitivity assay.

Section 14 of the QMD (change in outcomes due to increasing test
sensitivity) shows, for the same outcomes as Section 13 (`gain = TRUE`
in `grid_outcome_spec`), the change from increasing test sensitivity
from `sens_from` (default 50%) to `sens_to` (default 80%), by isolation
policy (P1, P2 and P3 reported as sensitivities) and testing coverage.
How each outcome is reported is set by `direction` in
`sensitivity_gain_spec` in `Scripts/Model1_ODE_Functions.R`:

- `"avoided"` (most counts, 2b and 3f): value at `sens_from` minus value
  at `sens_to`, so a negative value is an increase;
- `"additional"` (side-room bed-days): value at `sens_to` minus value at
  `sens_from`, since they rise as more cases are confirmed and isolated;
- `"change"` (the proportions 1c and 1d): value at `sens_to` minus value
  at `sens_from`.

Counts are per year; proportions and rates are over the analysis
period. Each outcome has a tab
with a table (0%, 20%, ..., 100% coverage) and a chart (every coverage
value in the grid). Both values must be in the grid.

------------------------------------------------------------------------

## Validation and safeguards

Implemented safeguards include:

- explicit `stop()` checks for invalid inputs and parameters,
- warnings for potentially unrealistic values,
- non-negative state clamping in the ODE function,
- `N <- max(N, 1e-12)` to prevent divide-by-zero,
- non-negative guard for `missed_untested` via `pmax`.

------------------------------------------------------------------------

## Research questions

The model is designed to focus on diagnostic testing of hospital
patients with CDI. It does not include testing of symptomatic patients
without CDI, or screening asymptomatic patients.

It can be used to quantify the effect of changes in testing algorithm
sensitivity for different levels of testing coverage for symptomatic
hospital patients with CDI and different hospital policies on location
of suspected cases awaiting test results and/or confirmed cases while
receiving treatment (i.e. on a general ward versus isolated in a
side-room). The value of other parameters, such as test turnaround time
and time required for patient transfers between general ward and side
room beds can be varied for sensitivity analysis.

Differences in total and hospital-onset CDI cases between scenarios are
driven entirely by differences in hospital-acquired CDI because:

a)  the model uses fixed admission rates and fixed prevalence of
    colonisations and infections among admitted patients, and
b)  testing and isolation affect only infected patients, so progression
    to CDI among patients admitted colonised is the same in every
    scenario.

Where outcomes are described as estimates, they are calculated using
the shares of patient flows from one state to multiple states, or from
multiple states to one state. For example:

a)  cases confirmed through re-testing are estimated by apportioning the
    flow from Xtest to Xtreat by the share of flows into Xtest that come
    from FN (retests), rather than from I (first tests);
b)  hospital-acquired CDI is estimated by apportioning the flow from C
    to I by the share of all flows into C that come from S (colonisation
    in hospital), rather than from admissions to C.

### Outcomes of interest

1)  Case ascertainment:
    a)  the total number of incident CDI cases (to include patients that
        were colonised on admission that progress to CDI during their
        hospital stay, patients admitted with CDI and patients that
        become colonised and progress to CDI during their hospital
        stay);
    b)  the number of missed CDI cases, i.e. those that are discharged
        or die following progression to infection that either did not
        have a specimen tested for CDI or received a false negative test
        result;
    c)  the number and proportion of incident CDI cases amongst hospital
        patients that are confirmed through diagnostic testing;
    d)  the estimated number and proportion of incident CDI cases that
        are confirmed through re-testing following a false-negative test
        result (delayed diagnosis).
2)  Transmission of *C. difficile*:
    a)  the estimated number of hospital-acquired infections (patients
        that become colonised and progress to CDI during their hospital
        stay); and
    b)  the difference in prevalence of colonisation and infection
        amongst discharged versus admitted patients.
3)  Hospital resource use:
    a)  total hospital bed occupancy, reflecting any impacts on bed
        pressure and average patient length of stay;
    b)  general ward bed occupancy;
    c)  required side room bed occupancy;
    d)  the proportion of hospital bed occupancy that is required for
        patient isolation in a side room;
    e)  the number of faecal specimens collected and tested for CDI;
    f)  the number of faecal specimens collected and tested for CDI per
        CDI case confirmed;
    g)  the number of confirmed cases that are started on treatment for
        CDI (equal to 1c).
4)  CDI case outcomes:
    a)  the number of incident CDI cases that are treated and recover
        (then either returned to a general ward or discharged);
    b)  the number of incident CDI cases that are discharged untreated;
    c)  the number of incident CDI cases that die in hospital (with no
        specimen sent for testing, awaiting a test result, awaiting
        re-testing following a false negative test result, or during
        treatment);
    d)  number of CDI-attributable deaths in hospital (accounting for
        deaths that may have occurred due to other causes).

### Cost effectiveness analysis

This involves calculation of Incremental Cost Effectiveness Ratios
(ICERs) for specific changes in testing or isolation policy, that allows
comparison of multiple outcomes within a single metric (£ cost per
avoided QALY loss, or £ cost per avoided hospital-acquired CDI case). It
would require information on the hospital resource costs (specimen
collection and testing, side room and general ward patient care, CDI
treatment) and patient welfare values for untreated, treated/recovered
CDI cases and deaths from CDI (discounted QALY losses).

The analysis would not take account of outcomes for CDI cases that occur
in the community, as a result of:

a)  discharged colonisations that progress to infection;
b)  onwards transmission from discharged colonised or infected patients.

Model outcomes required for calculation of an ICER for a change in
testing or isolation are shown in the tables below.

| Model outcome (ICER numerator, change in costs) | ICER numerator |
|---|---|
| Number of patient days on a general ward | × general ward care, cost per day |
| Number of patient days in a side room | × side room care, cost per day |
| Number of confirmed CDI cases | × CDI antimicrobial therapy, cost per course of treatment |
| Number of faecal specimens tested for CDI | × specimen collection and testing, cost per specimen |

| Model outcome (ICER denominator, avoided QALY loss or avoided hospital-acquired CDI case) | ICER denominator |
|---|---|
| Number of treated, recovered CDI cases | × QALY loss for treated cases |
| Number of untreated, discharged cases | × QALY loss for untreated cases |
| Number of CDI-attributable deaths | × discounted QALY loss for deaths |
| Number of hospital-acquired CDI cases | |

Extending the model to include probabilistic analysis, with
distributions for uncertain parameters would allow value of information
analysis.

------------------------------------------------------------------------

## Research question outcomes table (QMD Section 12)

`make_research_outcomes_table()` in `Scripts/Model1_ODE_Functions.R`
derives every outcome of interest from the cumulative totals over the
analysis period after the burn-in (`make_final_outcomes()`), the
patient-days (`calculate_total_patient_days()`) and the derived
proportions (`add_derived_outcomes()`). With `per_year = TRUE` (as in the
QMD) counts and patient-days are divided by the length of the analysis
period. The table has one row per outcome, numbered as above, followed by
the ICER numerator and denominator inputs, and one column per scenario.
The QMD displays it for the example scenarios (`research_outcomes.csv`)
and saves it for every scenario in the grid (`research_outcomes_all.csv`).

| Outcome | Calculation |
|----|----|
| 1a) Incident CDI cases | `cum_total_infected_incidence` |
| 1b) Missed CDI cases | `cum_total_missed_cases` |
| 1c) CDI cases confirmed, n and % | `cum_confirmed_cases`; `prop_confirmed` = `cum_confirmed_cases / cum_total_infected_incidence` |
| 1d) Confirmed via retest (estimated)^i^, n and % | `cum_confirmed_via_retest`; `prop_confirmed_via_retest` = `cum_confirmed_via_retest / cum_total_infected_incidence` |
| 2a) Hospital-acquired CDI cases (estimated)^ii^ | `cum_hospital_attributable_I` |
| 2b) Colonisation prevalence at admission and discharge, and difference (pp)^iii^ | `cum_admitted_colonised / cum_total_admissions`; `cum_discharges_colonised / cum_total_discharges`; `diff_prev_C` = discharge − admission |
| 2b) CDI prevalence at admission and discharge, and difference (pp)^iii^ | `cum_admitted_infected / cum_total_admissions`; `cum_discharges_while_infected / cum_total_discharges`; `diff_prev_CDI` = discharge − admission |
| 3a) Patient bed-days and average length of stay | `hospital_bed_days`; `hospital_bed_days / cum_total_admissions` |
| 3b) General ward patient bed-days | `general_ward_bed_days` |
| 3c) Side-room patient bed-days | `side_room_bed_days` |
| 3d) Side-room patient bed-days, % of all | `side_room_bed_days / hospital_bed_days` |
| 3e) Faecal specimens tested | `cum_total_test` |
| 3f) Faecal specimens tested per confirmed case | `tests_per_confirmed` = `cum_total_test / cum_confirmed_cases` ("-" if none confirmed) |
| 3g) Confirmed cases started on treatment | `cum_confirmed_cases` (equal to 1c) |
| 4a) Treated and recovered | `cum_recovered_cases` |
| 4b) Discharged untreated | `cum_discharges_while_infected` |
| 4c) Died in hospital following CDI | `cum_deaths_post_infection` |
| 4d) CDI-attributable deaths in hospital^iv^ | `cum_deaths_CDI_attributable` |

ICER numerator (costs): general ward patient bed-days, side-room patient
bed-days, confirmed CDI cases (treatment courses), faecal specimens
tested. ICER denominator (outcomes): treated and recovered CDI cases,
untreated discharged CDI cases, CDI-attributable deaths,
hospital-acquired CDI cases (estimated).

**Notes**

The proportions confirmed (1c) and missed (1b) do not sum to 1: the gap
is mainly deaths while awaiting a test result (in **Xtest**: tested, but
neither confirmed nor missed).

^i^ `cum_confirmed_via_retest` is estimated post-hoc by apportioning
`cum_confirmed_cases` in proportion to the share of all tests that were
retests (`cum_retest_false_negatives / cum_total_test`). This assumes
the confirmation rate (`theta`) is equal for first tests and retests.

^ii^ `cum_hospital_attributable_I` is estimated post-hoc by attributing
`cum_infected_in_hospital` (C→I progressions within hospital) in
proportion to the share of all colonisation that is hospital-acquired
(`cum_colonised_in_hospital / cum_total_colonised_incidence`). See
`make_final_outcomes()`.

^iii^ A positive difference indicates that the hospital is a net source
of colonisation or infection.

^iv^ CDI-attributable deaths are deaths following CDI that would not
have occurred at the baseline death rate (`mortS`): deaths of untreated
CDI patients (I, Xtest, FN) × (1 − 1/`mortI_mult`) plus deaths of
treated patients (Xtreat) × (1 − 1/`mortXtreat_mult`). In the grid run
the multipliers are taken from each scenario's parameters
(`mortI / mortS`, `mortXtreat / mortS`).

`make_final_outcomes()` also returns outcomes not in the table (e.g.
in-hospital colonisations by route and source, deaths by state), which
are saved in `final_outcomes_all.csv`.

------------------------------------------------------------------------

## Outputs

All outputs are written to `Outputs/<run_label>/` (created
automatically), including:

All scenarios in the grid:

- `scenario_control_parameters.csv`
- `model_parameters_all_scenarios.csv`
- `grid_runs.rds` (cached grid results)
- `patient_days_summary.csv`
- `all_negative_rows_combined.csv`
- `stabilisation_summary_states_flows.csv`
- `final_outcomes_all.csv` (control values, cumulative outcomes,
  patient-days and derived proportions, one row per scenario)
- `grid_key_outcomes.csv` and `grid_key_outcome_differences.csv` (the
  Section 13 and 14 outcomes, counts per year, and their differences
  between policies including P3 vs P2, in long format)
- `grid_rates.csv` (rates per 10,000 patient bed-days, side-room and
  missed proportions and missed cases per year by origin, long format)
- `research_outcomes_all.csv` (per year)
- `grid_lines_<outcome>.png`, `grid_lines_diff_<outcome>.png` and
  `grid_heatmap_<outcome>.png` (heatmaps and difference heatmaps in one
  figure; one set per Section 13 outcome)
- `sensitivity_gain_<from>_to_<to>.csv` and
  `sensitivity_gain_<from>_to_<to>_<outcome>.png` (Section 14, e.g.
  `sensitivity_gain_050_to_080.csv`)

Example scenarios only:

- `solution_<scenario>.csv` (full time series)
- `annual_summary_all.csv` (every year of the analysis period; the QMD
  table shows year 1, as every year is the same at steady state)
- `research_outcomes.csv`
- time-series plots (Section 9): `plot_states_lines_facet.png`, `plot_states_lines_facet_noS.png`,
  `plot_states_area_facet.png`, `plot_states_area_facet_noS.png`,
  `plot_bed_occupancy.png` (area plot of general ward and side-room
  beds), `plot_flows_to_C_and_I.png`, `plot_case_ascertainment_line.png`,
  `plot_deaths_by_state_area_noS.png`, `plot_discharges_by_state_area_noS.png`
- stacked bar plots of annual totals (Section 11):
  `plot_colonisations_by_source.png` (background acquisition and
  patient-to-patient transmission from C, I, Xtest, FN and Xtreat),
  `plot_testing_bar.png` (first tests and retests),
  `plot_case_confirmation_bar.png` (confirmed via first test and retest),
  `plot_missed_cases_bar.png` (untested and false negatives) and
  `plot_case_outcomes_bar.png` (discharged untreated: untested from I and
  false negatives from FN; recovered after treatment; died following CDI)
  and `plot_bed_days_bar.png` (patient bed-days per year in general ward
  beds and in side rooms)

The rendered HTML (`YYMMDD_Model1_ODE_Model_Run_<run_label>.html`) embeds
these plots, so it can be moved or shared without the `Outputs/` folder.
