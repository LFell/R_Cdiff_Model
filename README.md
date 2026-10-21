---
editor_options: 
  markdown: 
    wrap: 72
---

# Hospital CDI case ascertainment and transmission model

## Overview

This repository contains a deterministic compartmental ODE model for
in-hospital transmission of *C. difficile*, with explicit representation
of:

- suspected case identification and optional isolation,
- test sensitivity and turnaround time,
- false negatives returning to general care,
- treated confirmed cases, optionally isolated,
- patient transfers between general ward and side-room beds,
- discharge and mortality flows,
- cumulative missed cases and side-room occupancy.

The model is implemented in R with `deSolve` and run through a Quarto
workflow. This README describes the current version, **Model 4**.

------------------------------------------------------------------------

## Model versions

Major model development steps are numbered and marked with git tags. The
current model's files are named with its number (`Model4_Run.qmd`,
`Scripts/Model4_Functions.R`, `Scripts/Model4_Parameters.R`). Earlier
versions are not kept on the main branch: check out their tag to run
them, e.g. `git worktree add ../R_Cdiff_Model_model2 model-2`. At each
tag the files still have the names used at the time (`Model1_ODE_...`).

| Model | Git tag (commit) | What changed |
|---|---|---|
| Model 1 | `model-1` (`5056b0e`) | First version on git: isolation before testing as a yes/no switch (`isolate_before_test`). |
| Model 2 | `model-2` (`9601ff6`) | Three isolation policies (P1 none, P2 confirmed cases, P3 suspected and confirmed cases); patient transfer time added to the rate of each transition involving a move; scenario grid of policy × testing coverage × test sensitivity; one-at-a-time (OAT) parameter test report (`Model1_ODE_OAT_Run.qmd`). |
| Model 3 | `model-3` (`52ebb3a`) | Patient transfer states (T1–T5) instead of adding transfer time to rates, so the share of positive results equals the test sensitivity under every policy; comparison report against Model 2 (`Model1_ODE_Transfer_Compare.qmd`). |
| Model 4 | `model-4` (`0dc2598`) | Case identification rate `gamma` set from the proportion of CDI cases sent for testing (`prop_I_sent_for_testing`), so the control variable is the share actually tested; new outcome 1c (CDI cases ever tested, point prevalence) for comparison with point-prevalence studies; under P3, tests counted on entry to Xtest and deaths in T3 counted as missed cases; renamed from `prop_I_suspected` / "testing coverage". Files renamed to `Model4_*` after the tag. |

Renamed copies of the Model 2 and Model 3 documents and scripts (e.g.
`Model2_OAT_Test.qmd`, `Model3_vs_Model2_Comparison.qmd`) are filed
locally in `Not_in_repo/Old_scripts/Model2/` and `Model3/` for
reference, and their rendered reports in `Not_in_repo/Old_outputs/`.
These copies do not run from there; use the git tags. The OAT test has
not yet been ported to Model 4.

------------------------------------------------------------------------

## Model structure

### ![](Images/261021%20P1-P3%20combined.png)

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
- `T1`–`T5`: Patient transfer states (see below)
- `D`: Cumulative deaths (absorbing count state)

Total in-hospital living census is:

$$N = S + C + I + Xtest + FN + Xtreat + T1 + T2 + T3 + T4 + T5$$

### Patient transfer states

Each move between a general ward and a side room passes through a
transfer state, so transfer time is not added to the rate of any other
transition. Patients in a transfer state are still in their original
bed, and leave it at rate `tau = 1 / t_patient_transfer`:

- P2: `T1` (Xtest → T1 by `theta`, T1 → Xtreat by `tau`) and `T2`
  (Xtreat → T2 by `delta`, T2 → S by `tau`)
- P3: `T3` (I → T3 by `gamma` and FN → T3 by `sigma`, T3 → Xtest by
  `tau`), `T4` (Xtest → T4 by `pi`, T4 → FN by `tau`) and `T5`
  (Xtreat → T5 by `delta`, T5 → S by `tau`)

| State | Location | Infectiousness | Discharge | Mortality | Force of colonisation (λ) |
|----|----|----|----|----|----|
| `T1` | General ward | `betaI` | None | `mortT1 = mortI` | Included |
| `T2` | Isolation | 0 | None | `mortT2 = mortS` | Census N only |
| `T3` | General ward | `betaI` | None | `mortT3 = mortI` | Included |
| `T4` | Isolation | `betaI * beta_mult_isolated` | `1 / losFN` | `mortT4 = mortI` | Included |
| `T5` | Isolation | 0 | None | `mortT5 = mortS` | Census N only |

All transfer states count in the census N (the denominator of λ). T2 and T5 are non-infectious, so they have no transmission parameter and add nothing to the numerator. States
not used by a policy stay empty, and P1 has no transfers.

------------------------------------------------------------------------

## Key transitions

- `S -> C`: Force of colonisation (`lambda`)
- `C -> I`: Progression to infection (`alpha`)
- `C -> S`: Decolonisation (`mu`)
- `I -> Xtest` (via `T3` under P3): Suspected-case identification
  (`gamma`)
- `I` not identified: Missed untested cases (discharge `disI` and death
  `mortI`, competing with `gamma`)
- `Xtest -> Xtreat` (via `T1` under P2): Test positive, case confirmed
  (`theta`)
- `Xtest -> FN` (via `T4` under P3): False negative, return to general
  ward (`pi`)
- `FN -> Xtest` (via `T3` under P3): False negatives to be re-tested
  (`sigma`)
- `Xtreat -> S` (via `T2` under P2, `T5` under P3): Treatment and
  recovery to general ward (`delta`)
- Transfer states to their destination (`tau`)
- Discharge and mortality from all relevant compartments
- Deaths accumulate in `D`

------------------------------------------------------------------------

## Model assumptions

Default parameter values are given in brackets; they are set in `Scripts/Model4_Parameters.R`.

**Pathogen and population**

- A single toxigenic strain is in circulation.
- The hospital is a single, well-mixed population: there are no separate wards or bays. The model is deterministic and its parameters don't change over time (no seasonality).
- Patients are not followed after discharge. Readmission is not modelled.

**Admissions and discharges**

- Admissions arrive at a constant rate (133/day). They do not replace discharges, so the hospital census can vary.
- Admissions are split between S, C and I in fixed proportions reflecting community prevalence (8.6% colonised, 0.2% infected, the rest susceptible). Community prevalence is fixed and isn't affected by discharges from the hospital.
- Nobody is admitted directly to Xtest or Xtreat.
- Each state has its own constant discharge rate, set by its length of stay:
  - S and C: 7 days;
  - I and FN: 14 days;
  - Xtreat: 20% of treated patients are discharged and 80% return to a general ward as susceptible, over 14 days (when isolated, those returning to the ward then wait in T2 or T5 for transfer).
- Patients awaiting a test result (Xtest) stay in hospital until the result arrives; they can't be discharged from Xtest.
- Patients awaiting transfer can't be discharged, except false negatives awaiting return to the ward (T4), who are discharged at the FN rate.
- Death rates are constant:
  - S and C share a baseline rate (0.1% risk over a 7-day stay);
  - I, Xtest and FN, and patients with CDI awaiting transfer (T1, T3, T4), die at 3× the baseline (`mortI_mult`);
  - Xtreat dies at 1.5× the baseline (`mortXtreat_mult`);
  - recovered patients awaiting transfer (T2, T5) die at the baseline rate.

**Colonisation and progression**

- Susceptible patients become colonised through:
  - background acquisition from spores in the environment, at a constant rate that doesn't depend on the number of patients in each state (β₀ = 0.002);
  - frequency-dependent transmission from colonised and infected patients (λ = β₀ + Σβₖk/N).
- Colonised and infected patients all pose a transmission risk:
  - colonised (C): 0.008;
  - infected (I) and false negatives (FN): 0.04;
  - suspected cases (Xtest): 0.04 on a general ward;
  - confirmed cases on treatment (Xtreat): 0.02 on a general ward;
  - patients awaiting transfer: as untreated infected patients in their current bed (T1 and T3 on a general ward: 0.04; T4 in a side room: 0.04 × 5%); recovered patients (T2, T5) don't transmit.
- Treatment is assumed to reduce shedding from 0.04 at the start to 0 at the end. The model uses the average over the treatment period (0.02) as a constant rate.
- Isolation in a side room reduces transmission to 5% of the general-ward value. It does not eliminate it.
- Progression from C to I is increased by antibiotic exposure among colonised patients (SCI structure, not PSCIX): the progression rate rises by a relative risk of 2 in the 37.3% of colonised patients receiving antibiotics. Without antibiotics, 50% of colonised patients progress, with a mean time to progression of 5 days.
- Colonised patients can lose colonisation without treatment (5% over 5 days).
- There is no immunity, and no recurrence or relapse. Patients who recover after treatment return to the susceptible state.
- Colonised patients are not screened, detected or decolonised.

**Case identification and testing**

- Only infected (symptomatic) patients are identified as suspected cases and tested. Patients in S and C don't develop diarrhoea that leads to testing. The test is assumed 100% specific, so there are no false positives. As a result, "faecal specimens tested" counts only tests of patients with CDI, and underestimates real testing volume.
- The proportion of CDI cases sent for testing sets the identification rate. It is the share of infected patients identified before they are discharged or die, so the mean time to identification is not set separately: it is shorter when more cases are sent for testing (see `prop_I_sent_for_testing` under **Scenario definitions**).
- Test results arrive after a fixed turnaround time (2 days). Each test, including retests, returns a positive result with probability equal to the test sensitivity. Results of repeat tests on the same patient are independent.
- Patients with a false-negative result return to a general ward (FN), still symptomatic and transmitting. They are retested after a mean of 7 days unless they are discharged or die first.
- Confirmed cases are treated. Under P2 and P3 they are isolated; under P1 they stay on a general ward.
- Moving a patient between a general ward and a side room takes a mean of 0.5 days. Patients wait for the move in a transfer state (T1-T5) in their original bed, and leave it at rate `tau = 1 / t_patient_transfer`. Transfer time is not added to the rates of any other transition.
- Side-room capacity is unlimited.

**Outcome definitions**

- **Missed cases:** patients who leave I, FN, T3 (awaiting a side room before the specimen is taken) or T4 (false negatives awaiting return to the ward) by discharge or death without a confirmed diagnosis. Deaths while awaiting a result (in Xtest, after the specimen is taken) count as neither missed nor confirmed.
- **Hospital-onset CDI:** all in-hospital C→I progressions, i.e. patients admitted colonised plus patients colonised in hospital who progress to CDI during their stay.
- **Hospital-acquired CDI:** in-hospital C→I progressions, multiplied by the share of all colonisations that were acquired in hospital (a subset of hospital-onset CDI). This assumes patients colonised in hospital and those admitted colonised progress at the same rate.

------------------------------------------------------------------------

## Case accounting outputs

The ODE function tracks the following diagnostic rates at each time
step, in addition to state variables:

- **Colonisation**: in-hospital colonisations, colonised admissions
- **Infection**: in-hospital infections, infected admissions
- **Case identification**: suspected cases (`gamma * I`), confirmed
  cases (`theta * Xtest`), false negatives (`pi * Xtest`)
- **Testing**: first tests and retests, counted on entry to Xtest (under
  P3, on leaving T3, split by the share of T3 that came from I or FN)
- **Missed cases**: untested (discharges and deaths from I, and deaths
  in T3 of patients from I) and false negatives (discharges and deaths
  from FN and T4, and deaths in T3 of patients from FN)
- **Current CDI cases and those ever tested** (`CDI_census`,
  `CDI_tested_census`), for the point-prevalence outcome 1c
- **Outcomes**: recovered cases, deaths by state, side-room occupancy
- **Admissions and discharges** by state

Per-step, cumulative, and annual summary versions of all flow quantities
are computed in post-processing.

------------------------------------------------------------------------

## File structure

```         
R_Cdiff_Model/
├── Model4_Run.qmd                     # Main analysis document (Model 4)
├── render_run.R                       # render_run() for labelled runs
├── renv.lock, renv/, .Rprofile        # renv package library for this project
├── Images/                            # Model diagrams used in the report
├── Scripts/
│   ├── Model4_Parameters.R            # Parameters and scenario builder
│   └── Model4_Functions.R             # ODE, post-processing, scenario grid and plot functions
├── Outputs/                           # (git-ignored)
│   └── Model4/
│       └── <run_label>/               # One folder per labelled model run
└── Not_in_repo/                       # (git-ignored) filed copies of earlier models, old outputs
```

### `Scripts/Model4_Parameters.R`

Defines:

- fixed model parameters,
- helper `prob2rate()` conversion,
- validation helpers: `validate_fixed_parameters()`,
  `validate_scenario_inputs()`, `validate_rate_parameters()`,
- scenario builder `make_mod_parms(...)`.

### `Scripts/Model4_Functions.R`

Defines:

- `validate_run_settings()`: checks model time and time step inputs
  (both positive, time step smaller than model time),
- `model_ode(time, state, parms)`: ODE system for use with
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
- `add_derived_outcomes()`: proportions sent for testing, ever tested
  (point prevalence, 1c), confirmed, missed and confirmed
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

### One-at-a-time (OAT) parameter test

Model 2 had a standalone verification document (`Model1_ODE_OAT_Run.qmd`
at tag `model-2`; filed locally as `Model2_OAT_Test.qmd`) that varied
each fixed parameter one at a time and checked that the outcomes moved
in the theoretically expected direction under each isolation policy. It
uses Model 2's parameters, so it does not run against Model 4; porting
it to Model 4 is a planned next step.

### `Model4_Run.qmd`

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
  - a table comparing the proportion of CDI cases sent for testing with
    the point-prevalence share ever tested (Section 12),
  - outcomes by proportion of CDI cases sent for testing and test
    sensitivity: plots comparing
    the isolation policies across the grid (Section 13),
  - change in outcomes due to increasing test sensitivity, e.g. from 50% to
    80% (Section 14),
- CSV and PNG export to `Outputs/Model4/<run_label>/`.

Accepts a `run_label` parameter (see **Run instructions** below).

### `render_run.R`

Utility script for rendering with a labelled output folder. Source this
file by calling `source("render_run.R")`, then call `render_run("A")`
from the R console. It renders `Model4_Run.qmd` only; earlier models are
run from their git tags.

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

1.  Open `Model4_Run.qmd` and set the vectors of scenario
    control values in **Section 3**. Save the qmd (no need to change
    the filename).
2.  In the R console, from the project root:

``` r
source("render_run.R")
render_run("A")
```

This will:

- execute the QMD in a fresh R session,
- save all CSVs and plots to `Outputs/Model4/A/`,
- produce `A_Model4_Run.html` in the project root. Rendering the
  same label again overwrites it (and `Outputs/Model4/A/`), so to keep an
  earlier report for comparison use a new label, e.g. a date
  (`render_run("20261014")`).

The HTML is self-contained (`embed-resources: true` in the QMD YAML): the
plots are embedded in the file, so each report keeps its own plots
even after `Outputs/Model4/A/` is overwritten by a later run, and the HTML can be
moved or shared on its own (about 10 MB).

To run a second set of scenarios, edit the scenario values again and
call `render_run("B")`, and so on. The `run_label` in the QMD YAML does
not need to be changed — `render_run()` overrides it.

### Interactive rendering (quick checks)

You can also click **Render** in RStudio, which uses the default
`run_label: "default"` from the YAML and saves outputs to
`Outputs/Model4/default/`. This is useful for quick checks but does not produce
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

------------------------------------------------------------------------

## Scenario definitions

Scenario control values are set as vectors in **Section 3** of
`Model4_Run.qmd`, and every combination is run (the scenario
grid; by default 3 policies × 11 values of the proportion sent for
testing (0–90% and 99%) × 11 sensitivity values
= 363 scenarios, about 3–4 minutes on one core, about 30 seconds in
parallel). The three control variables are:

- `isolation_policy`: which patients are isolated in a side room,
  numbered from least to most isolation (set by `isolation_policy_levels`
  in `Scripts/Model4_Functions.R`):
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
- `prop_I_sent_for_testing`: proportion of CDI cases sent for testing
  (0 to below 1): the share of patients entering I who are identified as
  suspected cases before they are discharged or die. Identification
  competes with discharge and death from I, so
  `prop_I_sent_for_testing = gamma / (gamma + disI + mortI)` and the
  model sets
  `gamma = prop_I_sent_for_testing * (disI + mortI) / (1 - prop_I_sent_for_testing)`.
  `gamma` is infinite at 1, so the grid's highest value is 99%. The mean
  time in I before identification, `1 / (gamma + disI + mortI)`, is
  longer when fewer cases are sent for testing (with the default
  parameters about 9 days at 33%, 7 days at 50% and 1.4 days at 90%).
  The baseline value is set from the share of hospital-onset diarrhoea
  patients ever tested for CDI during their stay in a point-prevalence
  study (Mawer et al. 2019, doi 10.1016/j.jhin.2019.05.001); outcome 1c
  is the model's equivalent, for comparison. Plots and tables label it
  "Proportion of CDI cases sent for testing".
- `test_sens`: diagnostic test sensitivity (0--1). The model uses it to
  calculate the confirmed case rate (`theta = test_sens /
  t_test_turnaround`) and false negative rate (`pi = (1 - test_sens) /
  t_test_turnaround`). Neither includes transfer time, so the share of
  results that are positive equals `test_sens` under every policy.

Test turnaround time (`t_test_turnaround`, days from test administration
to result) is a fixed parameter in `Scripts/Model4_Parameters.R`.

Scenario labels combine a policy code with the proportion sent for
testing and the sensitivity as percentages, e.g.
`P2_confirmed_test050_sens100` (`P1_none`, `P2_confirmed`,
`P3_susp_conf`). `example_settings` in Section 3 lists the settings for
the example scenarios, in the order shown in Sections 9 and 11 (by
default, sent for testing/sensitivity: 0%/0%, 50%/50%, 99%/50%, 50%/100%
and 99%/100%), each run under all three policies (15 scenarios). They
are used for time-series plots, annual summaries and the outcome tables.
When 0% are sent for testing nobody is tested, so the three policies
give the same results (a no-testing baseline).

In Section 9, the time-series line and area plots show the policies as columns and
the settings as rows; in Section 11, the stacked bar plots group the bars by setting, with P1,
P2 and P3 side by side and the setting label beneath each group.

Grid results are cached in `Outputs/Model4/<run_label>/grid_runs.rds` and
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
- the transfer states used (see **Patient transfer states**): none under
  `"none"`, T1 and T2 under `"confirmed"`, and T3, T4 and T5 under
  `"suspected_and_confirmed"`;
- side-room bed occupancy: none under `"none"`; `Xtreat` and `T2` under
  `"confirmed"`; and `Xtest`, `Xtreat`, `T4` and `T5` under
  `"suspected_and_confirmed"`.

Section 13 of the QMD (outcomes by proportion of CDI cases sent for testing and test sensitivity) compares the three isolation policies across
the proportion of CDI cases sent for testing and test sensitivity. The outcomes are the research
question outcomes defined in `grid_outcome_spec` in
`Scripts/Model4_Functions.R` (numbered as in the research questions
below):

| Outcome | Variable | Format |
|---|---|---|
| 1a) Incident CDI cases | `cum_total_infected_incidence` | count per year |
| 1b) Missed CDI cases | `cum_total_missed_cases` | count per year |
| 1c) CDI cases ever tested, point prevalence at model stabilisation | `prop_current_CDI_tested` | proportion |
| 1d) Proportion of CDI cases confirmed | `prop_confirmed` | proportion |
| 1e) Proportion confirmed via retest (estimated) | `prop_confirmed_via_retest` | proportion |
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
  of 0%, 20%, 40%, 60%, 80% and 99% of CDI cases sent for testing
  (`grid_line_sent` in `Scripts/Model4_Functions.R`; categorical
  legend, 0% yellow to 99% purple), faceted by policy with a shared y
  axis;
- the same plot for the differences from P1 (P2 vs P1, P3 vs P1) at the
  same proportion sent for testing and sensitivity, placed under P2 and P3 (the space under
  P1 is left empty);
- heatmaps of test sensitivity (x) by proportion of CDI cases sent for testing (y), for every
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
sensitivity) shows, for the outcomes in Section 13 except 1c (`gain = TRUE`
in `grid_outcome_spec`), the change from increasing test sensitivity
from `sens_from` (default 50%) to `sens_to` (default 80%), by isolation
policy (P1, P2 and P3 reported separately) and proportion sent for testing.
How each outcome is reported is set by `direction` in
`sensitivity_gain_spec` in `Scripts/Model4_Functions.R`:

- `"avoided"` (most counts, 2b and 3f): value at `sens_from` minus value
  at `sens_to`, so a negative value is an increase;
- `"additional"` (side-room bed-days): value at `sens_to` minus value at
  `sens_from`, since they rise as more cases are confirmed and isolated;
- `"change"` (the proportions 1d and 1e): value at `sens_to` minus value
  at `sens_from`.

Counts are per year; proportions and rates are over the analysis
period. Each outcome has a tab
with a table (0%, 20%, ..., 99% sent for testing) and a chart (every
value of the proportion sent for testing in the grid). `sens_from` and
`sens_to` must both be values in the grid.

------------------------------------------------------------------------

## Validation and safeguards

Implemented safeguards include:

- explicit `stop()` checks for invalid inputs and parameters,
- warnings for potentially unrealistic values,
- non-negative state clamping in the ODE function,
- `N <- max(N, 1e-12)` to prevent divide-by-zero,

------------------------------------------------------------------------

## Research questions

The model is designed to focus on diagnostic testing of hospital
patients with CDI. It does not include testing of symptomatic patients
without CDI, or screening asymptomatic patients.

It can be used to quantify the effect of changes in testing algorithm
sensitivity for different proportions of symptomatic hospital patients
with CDI sent for testing and different hospital policies on location
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
    c)  the proportion of CDI cases in hospital that have ever been
        tested for CDI, as a point prevalence at model stabilisation
        (for comparison with point-prevalence studies such as Mawer et
        al. 2019);
    d)  the number and proportion of incident CDI cases amongst hospital
        patients that are confirmed through diagnostic testing;
    e)  the estimated number and proportion of incident CDI cases that
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
        CDI (equal to 1d).
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

`make_research_outcomes_table()` in `Scripts/Model4_Functions.R`
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
| 1c) CDI cases ever tested, point prevalence (%) at model stabilisation^v^ | `prop_current_CDI_tested` = `CDI_tested_patient_days / CDI_patient_days` |
| CDI cases sent for testing (% of incident CDI cases; unnumbered) | `prop_sent_for_testing` = `cum_suspected_cases / cum_total_infected_incidence` (equals `prop_I_sent_for_testing` at steady state) |
| 1d) CDI cases confirmed, n and % | `cum_confirmed_cases`; `prop_confirmed` = `cum_confirmed_cases / cum_total_infected_incidence` |
| 1e) Confirmed via retest (estimated)^i^, n and % | `cum_confirmed_via_retest`; `prop_confirmed_via_retest` = `cum_confirmed_via_retest / cum_total_infected_incidence` |
| 2a) Hospital-acquired CDI cases (estimated)^ii^ | `cum_hospital_attributable_I` |
| 2b) Colonisation prevalence at admission and discharge, and difference (pp)^iii^ | `cum_admitted_colonised / cum_total_admissions`; `cum_discharges_colonised / cum_total_discharges`; `diff_prev_C` = discharge − admission |
| 2b) CDI prevalence at admission and discharge, and difference (pp)^iii^ | `cum_admitted_infected / cum_total_admissions`; `cum_discharges_while_infected / cum_total_discharges`; `diff_prev_CDI` = discharge − admission |
| 3a) Patient bed-days and average length of stay | `hospital_bed_days`; `hospital_bed_days / cum_total_admissions` |
| 3b) General ward patient bed-days | `general_ward_bed_days` |
| 3c) Side-room patient bed-days | `side_room_bed_days` |
| 3d) Side-room patient bed-days, % of all | `side_room_bed_days / hospital_bed_days` |
| 3e) Faecal specimens tested | `cum_total_test` |
| 3f) Faecal specimens tested per confirmed case | `tests_per_confirmed` = `cum_total_test / cum_confirmed_cases` ("-" if none confirmed) |
| 3g) Confirmed cases started on treatment | `cum_confirmed_cases` (equal to 1d) |
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

The proportions confirmed (1d) and missed (1b) do not sum to 1: the gap
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
CDI patients (I, Xtest, FN, and T1, T3 and T4 awaiting transfer) ×
(1 − 1/`mortI_mult`) plus deaths of
treated patients (Xtreat) × (1 − 1/`mortXtreat_mult`). In the grid run
the multipliers are taken from each scenario's parameters
(`mortI / mortS`, `mortXtreat / mortS`).

^v^ `prop_current_CDI_tested` is the share of current CDI cases who have
ever been tested, the model equivalent of a point-prevalence survey
such as Mawer et al. (2019): CDI patient-days in a tested state (Xtest,
FN, Xtreat, T1, T4, and the share of T3 that came from FN) divided by
all CDI patient-days (the same states plus I and the rest of T3),
summed over the analysis period. The specimen is taken on entry to
Xtest, so patients in T3 who came from I have not yet been tested. T3 is
not split by origin, so its share from FN is estimated from the inflows
to T3 (exact at steady state).

`make_final_outcomes()` also returns outcomes not in the table (e.g.
in-hospital colonisations by route and source, deaths by state), which
are saved in `final_outcomes_all.csv`.

------------------------------------------------------------------------

## Outputs

All outputs are written to `Outputs/Model4/<run_label>/` (created
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
- `sent_for_testing_point_prevalence.csv` (Section 12: proportion of CDI
  cases sent for testing and the point-prevalence share ever tested, by
  policy, at one test sensitivity)
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

The rendered HTML (`<run_label>_Model4_Run.html`) embeds
these plots, so it can be moved or shared without the `Outputs/` folder.
