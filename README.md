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
|----|----|----|
| Model 1 | `model-1` (`5056b0e`) | First version on git: isolation before testing as a yes/no switch (`isolate_before_test`). |
| Model 2 | `model-2` (`9601ff6`) | Three isolation policies (P1 none, P2 confirmed cases, P3 suspected and confirmed cases); patient transfer time added to the rate of each transition involving a move; scenario grid of policy × testing coverage × test sensitivity; one-at-a-time (OAT) parameter test report (`Model1_ODE_OAT_Run.qmd`). |
| Model 3 | `model-3` (`52ebb3a`) | Patient transfer states (T1–T5) instead of adding transfer time to rates, so the share of positive results equals the test sensitivity under every policy; comparison report against Model 2 (`Model1_ODE_Transfer_Compare.qmd`). |
| Model 4 | `model-4` (`0dc2598`) | Case identification rate `gamma` set from the proportion of CDI cases sent for testing (`prop_I_sent_for_testing`), so the control variable is the share actually tested; new outcome 1c (CDI cases ever tested, point prevalence) for comparison with point-prevalence studies; under P3, tests counted on entry to Xtest and deaths in T3 counted as missed cases; renamed from `prop_I_suspected` / "testing coverage". Files renamed to `Model4_*` after the tag. Changed after the tag (November 2026): `gamma = -log(1 - prob_test_daily)`, from the daily probability that an untested CDI case is sent for testing, so the mean time to testing, the share of cases sent for testing and the point-prevalence share ever tested (1c) are outputs (briefly `gamma = 1 / t_wait_test` in commit `fcc6ab4`); mean times to first test and to treatment added as outcomes 4e and 4f, and hospital-acquired CDI to the annual summaries. |

Renamed copies of the Model 2 and Model 3 documents and scripts (e.g.
`Model2_OAT_Test.qmd`, `Model3_vs_Model2_Comparison.qmd`) are filed
locally in `Not_in_repo/Old_scripts/Model2/` and `Model3/` for
reference, and their rendered reports in `Not_in_repo/Old_outputs/`.
These copies do not run from there; use the git tags. The OAT test has
been ported to Model 4 (`Model4_OAT_Test.qmd`).

------------------------------------------------------------------------

## Model structure

![](Images/261102%20P1%20P2%20P3%20combined.PNG)

This is the earlier hand-drawn diagram. The report now draws a flow diagram for each isolation policy
from the model code when it renders (Section 1a; `Scripts/Model4_Diagrams.R`), saved as
`flow_P1.png` to `flow_P3.png` in `Outputs/Model4/<run_label>/`.

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
  `tau`), `T4` (Xtest → T4 by `pi`, T4 → FN by `tau`) and `T5` (Xtreat →
  T5 by `delta`, T5 → S by `tau`)

| State | Location | Infectiousness | Discharge | Mortality | Force of colonisation (λ) |
|----|----|----|----|----|----|
| `T1` | General ward | `betaI` | None | `mortT1 = mortI` | Included |
| `T2` | Isolation | 0 | None | `mortT2 = mortS` | Census N only |
| `T3` | General ward | `betaI` | None | `mortT3 = mortI` | Included |
| `T4` | Isolation | `betaI * beta_mult_isolated` | `1 / losFN` | `mortT4 = mortI` | Included |
| `T5` | Isolation | 0 | None | `mortT5 = mortS` | Census N only |

All transfer states count in the census N (the denominator of λ). T2 and
T5 are non-infectious, so they have no transmission parameter and add
nothing to the numerator. States not used by a policy stay empty, and P1
has no transfers.

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

Default parameter values are given in brackets; they are set in
`Scripts/Model4_Parameters.R`.

**Pathogen and population**

- A single toxigenic strain is in circulation.
- The hospital is a single, well-mixed population: there are no separate
  wards or bays. The model is deterministic and its parameters don't
  change over time (no seasonality).
- Patients are not followed after discharge. Readmission is not
  modelled.

**Admissions and discharges**

- Admissions arrive at a constant rate (133/day). They do not replace
  discharges, so the hospital census can vary.
- Admissions are split between S, C and I in fixed proportions
  reflecting community prevalence (8.6% colonised, 0.2% infected, the
  rest susceptible). Community prevalence is fixed and isn't affected by
  discharges from the hospital.
- Nobody is admitted directly to Xtest or Xtreat.
- Each state has its own constant discharge rate, set by its length of
  stay:
  - S and C: 7 days;
  - I and FN: 14 days;
  - Xtreat: 20% of treated patients are discharged and 80% return to a
    general ward as susceptible, over 14 days (when isolated, those
    returning to the ward then wait in T2 or T5 for transfer).
- Patients awaiting a test result (Xtest) stay in hospital until the
  result arrives; they can't be discharged from Xtest.
- Patients awaiting transfer can't be discharged, except false negatives
  awaiting return to the ward (T4), who are discharged at the FN rate.
- Death rates are constant:
  - S and C share a baseline rate (a risk of 0.001 over a 7-day stay);
  - I, Xtest and FN, and patients with CDI awaiting transfer (T1, T3,
    T4), die at 3× the baseline (`mortI_mult`);
  - Xtreat dies at 1.5× the baseline (`mortXtreat_mult`);
  - recovered patients awaiting transfer (T2, T5) die at the baseline
    rate.

**Colonisation and progression**

- Susceptible patients become colonised through:
  - background acquisition from spores in the environment, at a constant
    rate that doesn't depend on the number of patients in each state (β₀
    = 0.002);
  - frequency-dependent transmission from colonised and infected
    patients (λ = β₀ + Σβₖk/N).
- Colonised and infected patients all pose a transmission risk:
  - colonised (C): 0.008;
  - infected (I) and false negatives (FN): 0.04;
  - suspected cases (Xtest): 0.04 on a general ward;
  - confirmed cases on treatment (Xtreat): 0.02 on a general ward;
  - patients awaiting transfer: as untreated infected patients in their
    current bed (T1 and T3 on a general ward: 0.04; T4 in a side room:
    0.04 × 5%); recovered patients (T2, T5) don't transmit.
- Treatment is assumed to reduce shedding from 0.04 at the start to 0 at
  the end. The model uses the average over the treatment period (0.02)
  as a constant rate.
- Isolation in a side room reduces transmission to 5% of the
  general-ward value. It does not eliminate it.
- Progression from C to I is increased by antibiotic exposure among
  colonised patients (SCI structure, not PSCIX): the progression rate
  rises by a relative risk of 2 in the 37.3% of colonised patients
  receiving antibiotics. Without antibiotics, 50% of colonised patients
  progress, with a mean time to progression of 5 days.
- Colonised patients can lose colonisation without treatment (5% over 5
  days).
- There is no immunity, and no recurrence or relapse. Patients who
  recover after treatment return to the susceptible state.
- There is no recovery from CDI without treatment (no flow from I back
  to C or S): patients with CDI leave the infected states only by
  treatment, discharge or death. Treatment duration and effectiveness
  are the same under every policy, so transfer delays under P2 (before
  treatment, T1) and P3 (before testing, T3, and after a false
  negative, T4) lengthen the duration of infection. The model does not
  capture any effect of these delays on severity or mortality, but
  quantifies them (Outcomes D5 and D6).
- Colonised patients are not screened, detected or decolonised.

**Case identification and testing**

- Only infected (symptomatic) patients are identified as suspected cases
  and tested. Patients in S and C don't develop diarrhoea that leads to
  testing. The test is assumed 100% specific, so there are no false
  positives. As a result, "faecal specimens tested" counts only tests of
  patients with CDI, and underestimates real testing volume.
- The daily probability that an untested CDI case is sent for testing
  (`prob_test_daily`) sets a constant identification rate,
  `-log(1 - prob_test_daily)`. Identification competes with discharge
  and death, so the share of cases ever tested and the actual mean time
  to testing among cases tested both follow from it (see
  `prob_test_daily` under
  **Scenario definitions**).
- The specimen is collected on entry to Xtest and treatment starts on
  entry to Xtreat, so patients waiting for transfer are not yet tested
  (T3, P3) or treated (T1, P2), and transfer time adds to the time to
  testing and/or treatment.
- Test results arrive after a mean turnaround time of 2 days
  (exponentially distributed, as for every transition). Each test,
  including retests, returns a positive result with probability equal to
  the test sensitivity. Results of repeat tests on the same patient are
  independent.
- Patients with a false-negative result return to a general ward (FN),
  still symptomatic and transmitting. They are retested after a mean of
  7 days unless they are discharged or die first.
- Confirmed cases are treated. Under P2 and P3 they are isolated; under
  P1 they stay on a general ward.
- Moving a patient between a general ward and a side room takes a mean
  of 0.5 days. Patients wait for the move in a transfer state (T1-T5) in
  their original bed, and leave it at rate
  `tau = 1 / t_patient_transfer`. Transfer time is not added to the
  rates of any other transition.
- Side-room capacity is unlimited.

**Outcome definitions**

- **Missed cases:** patients who leave I, FN, T3 (awaiting a side room
  before the specimen is taken) or T4 (false negatives awaiting return
  to the ward) by discharge or death without a confirmed diagnosis.
  Deaths while awaiting a result (in Xtest, after the specimen is taken)
  count as neither missed nor confirmed.
- **Hospital-onset CDI:** all in-hospital C→I progressions, i.e.
  patients admitted colonised plus patients colonised in hospital who
  progress to CDI during their stay.
- **Hospital-acquired CDI:** in-hospital C→I progressions, multiplied by
  the share of all colonisations that were acquired in hospital (a
  subset of hospital-onset CDI). This assumes patients colonised in
  hospital and those admitted colonised progress at the same rate.

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
  `CDI_tested_census`), for the point-prevalence Outcome A4
- **Outcomes**: recovered cases, deaths by state, side-room occupancy
- **Admissions and discharges** by state

Per-step, cumulative, and annual summary versions of all flow quantities
are computed in post-processing.

------------------------------------------------------------------------

## File structure

```         
R_Cdiff_Model/
├── Model4_Run.qmd                     # Main analysis document (Model 4)
├── Model4_OAT_Test.qmd                # One-at-a-time parameter test (verification)
├── render_run.R                       # render_run() / render_oat() for labelled runs
├── renv.lock, renv/, .Rprofile        # renv package library for this project
├── Images/                            # Model diagram used in the report (sources in Not_in_repo/Images/)
├── Scripts/
│   ├── Model4_Parameters.R            # Parameters and scenario builder
│   ├── Model4_Functions.R             # ODE, post-processing, scenario grid and plot functions
│   └── Model4_Diagrams.R              # Flow diagrams of each policy, drawn from a flow table
├── Outputs/                           # (git-ignored)
│   └── Model4/
│       ├── <run_label>/               # One folder per labelled model run
│       └── OAT/<run_label>/           # One folder per labelled OAT run
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
  (point prevalence, Outcome A4), confirmed, missed and confirmed via retest,
  admission and discharge prevalence and their differences, and
  specimens tested per confirmed case,
- `make_research_outcomes_table()`: derives and formats the research
  question outcomes (Outcome A1–Outcome D6) and ICER inputs into a single wide table,
- `make_final_outcomes_vs_base()`: absolute differences vs a chosen base
  scenario,
- `mean_wait_to_state(parms, target)`: mean time from CDI symptom onset
  (entry to I) to first test (`"Xtest"`) or to treatment (`"Xtreat"`),
  among the cases that get there, calculated exactly from the scenario's
  rates (the CDI pathway is a continuous-time Markov chain; includes
  transfer waits, test turnaround, false negatives and retests); added to
  the grid's final outcomes as `mean_wait_test` and
  `mean_wait_treatment` (Outcomes D5, D6),
- `wait_if_not_discharged()`, `mean_wait_to_testing()`,
  `make_testing_outputs()`, `format_testing_outputs()`: testing outputs
  that follow from `prob_test_daily` (mean time to testing if not
  discharged or dead first, actual mean time to testing among cases
  tested, shares of incident and prevalent cases ever tested), and the
  table of them in Sections 13 and 15,
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
  `plot_cdi_cases_by_origin_bar()`, `plot_bed_days_bar()`,
  `plot_scenario_params()`. All plots use one text size,
  `font_base_size` (9 pt, via `my_theme` and `grid_theme`); the QMD saves
  every plot 24 cm wide and shows it at the page width, so text appears
  the same size in every section.

### `Scripts/Model4_Diagrams.R`

Draws the flow diagram for each isolation policy (Graphviz, via `DiagrammeR`), sourced in Section 1a
of `Model4_Run.qmd`. Defines:

- `model_flows`: one row per flow in `model_ode()` (from, to, the routing flag that turns it on for
  a policy, its rate, or for a forked arrow the shared exit rate and the branch proportion, its label
  and its type: move, discharge or death); `model_states`: state labels, whether infectious, and
  location (general ward, side room or set by the policy),
- `scenario_flows()`: the flows used under a scenario, with their rates evaluated,
- `check_model_flows(parms)`: checks that the flow table gives the same derivative for every state
  as `model_ode()`, for a random state; stops if not, so the diagram cannot drift from the model,
- `diagram_pos`, `fixed_edge_labels`: fixed positions (inches) of the states, discharge/death nodes,
  junction dots, the formula key and the two labels between S and C, in three rows,
- `formula_key()`: the key box, with the formula for each rate on the diagram from the input
  parameters (lambda written out with the terms of the states the policy uses),
- `build_flow_dot(parms)`, `plot_flow_diagram(parms)`, `save_flow_diagram(parms, path)`: the
  Graphviz text, the diagram, and a PNG (2600 px wide). Only the policy in `parms` affects the
  diagram.

Forked arrows (one exit rate, split by proportions) are used for admissions (`adm_rate` split by
`prop_adm_C`, `prop_adm_I`), test results (`t_test_turnaround`, split by `test_sens`) and the end of
treatment (`t_treatment`, split by `prop_Xtreat_to_S`). When `model_ode()` changes, update
`model_flows` (and `diagram_pos` for any new state).

### `Model4_OAT_Test.qmd`

A standalone model verification document — **not for analysis**. Runs a
systematic one-at-a-time (OAT) sensitivity test to confirm that Model 4
responds in the expected direction when each fixed parameter is varied
individually, with all others held at their defaults. Render it with
`render_oat("<label>")` (see **Run instructions**); the HTML report
shows the pass/fail counts, any failed tests and the plots. It is ported
from the Model 2 test (`Model1_ODE_OAT_Run.qmd` at tag `model-2`); the
introduction of the report lists what changed.

**What it tests:** each of the 26 fixed parameters in
`Scripts/Model4_Parameters.R` (the render stops if any fixed parameter
has no specification) is run at a low and a high value, at a daily
probability of testing of 0.07 (`prob_test_daily`; about 50% of CDI cases sent
for testing) and 50% test sensitivity (so the untested,
false-negative and retesting pathways are all active), under each of the
three isolation policies.

- `low`: 0 for proportions and transmission parameters valid at 0; a
  small positive value for parameters that must be \> 0; 1.0 for
  `abx_mult` (no antibiotic effect); 0.1 × baseline for `adm_rate`.
- `high`: 1 for proportions, 5 × baseline for all other parameters.

**Outputs assessed (high vs low):** steady-state I (`ss_I`), annual
confirmed cases and annual missed cases (last complete year). Each has
an expected direction (`"increases"`, `"decreases"`, `"unchanged"` or
`"complex"`), with a 1% tolerance; `"complex"` marks theoretically
ambiguous directions, flagged for manual inspection rather than
auto-assessed. `oat_specs` gives the expectations for
`"suspected_and_confirmed"`, and `oat_specs_policy` overrides them where
they differ by policy. The % change in each Section 14 outcome
(`grid_outcome_spec`) is also reported, without expected directions.

**Files produced** (in `Outputs/Model4/OAT/<run_label>/`):
`oat_specs.csv`, `oat_specs_policy_overrides.csv`, `oat_results.csv`
(pass/fail per parameter × policy), `oat_outcomes_pct_change.csv`,
`combined_oat_plot.png` (tornado charts for the three checked outputs
and a pass/fail heatmap) and `oat_tornado_cumulative_outcomes_<n>.png`.

### `Model4_Run.qmd`

Main Quarto document. Runs the end-to-end workflow:

- setup and sourcing of scripts,
- scenario definition: vectors of control values, crossed into a
  scenario grid (edit values here for each run),
- model execution for every scenario in the grid (cached), plus full
  time series for a few example scenarios,
- negative-value and stabilisation checks for every scenario (Sections 8
  and 9),
- results for the example scenarios that help verify the model is
  working as expected (Sections 10–12):
  - time-series line and area plots, up to stabilisation (Section 10),
  - annual summary tables, including hospital-acquired CDI cases
    (Section 11),
  - stacked bar plots of annual summaries, including incident CDI cases
    by origin (Section 12),
- results for the write-up (Sections 13–17):
  - table of research question outcomes: the research outcomes table (annual
    totals, including rates per 10,000 patient bed-days) for every
    scenario, and differences between policies (P2 vs P1, P3 vs P1, P3
    vs P2; Section 13),
  - a table of testing outputs for each daily probability of testing
    (`prob_test_daily`): the mean time to testing if not discharged or
    dead first (days), the actual mean time to testing among cases
    tested (days), and the shares of incident and prevalent (Outcome A4, by
    policy) CDI cases ever tested (Section 13, repeated in Section 15),
  - outcomes by daily probability of testing and test sensitivity: plots
    comparing the isolation policies across the grid (Section 14),
  - change in outcomes due to increasing test sensitivity, e.g. from 50%
    to 80% (Section 15),
  - commentary on the effect of improving test sensitivity, by research
    question (1–4) and sub-outcome, with every figure calculated inline
    from the grid (Section 16),
  - responsiveness of outcomes to improving test sensitivity: a
    forest-style chart of the arc elasticity of every outcome by policy,
    at daily probabilities of testing of 0.03 and 0.40, and commentary
    grouping outcomes by how strongly they respond (Section 17),
- CSV and PNG export to `Outputs/Model4/<run_label>/`.

Accepts a `run_label` parameter (see **Run instructions** below).

### `render_run.R`

Utility script for rendering with a labelled output folder. Source this
file by calling `source("render_run.R")`, then call `render_run("A")`
(main model run) or `render_oat("A")` (OAT parameter test) from the R
console. Earlier models are run from their git tags.

------------------------------------------------------------------------

## Requirements

- R 4.5.0 (the version recorded in `renv.lock`)
- Packages: `tidyverse`, `deSolve`, `scales`, `knitr`, `cowplot`,
  `quarto`, and `DiagrammeR`, `DiagrammeRsvg`, `rsvg` for the flow
  diagrams (plus base R `parallel`), managed with **renv**. Opening the
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

1.  Open `Model4_Run.qmd` and set the vectors of scenario control values
    in **Section 4**. Save the qmd (no need to change the filename).
2.  In the R console, from the project root:

``` r
source("render_run.R")
render_run("A")
```

This will:

- execute the QMD in a fresh R session,
- save all CSVs and plots to `Outputs/Model4/A/`,
- produce `A_Model4_Run.html` in the project root. Rendering the same
  label again overwrites it (and `Outputs/Model4/A/`), so to keep an
  earlier report for comparison use a new label, e.g. a date
  (`render_run("20261014")`).

The HTML is self-contained (`embed-resources: true` in the QMD YAML):
the plots are embedded in the file, so each report keeps its own plots
even after `Outputs/Model4/A/` is overwritten by a later run, and the
HTML can be moved or shared on its own (about 10 MB).

To run a second set of scenarios, edit the scenario values again and
call `render_run("B")`, and so on. The `run_label` in the QMD YAML does
not need to be changed — `render_run()` overrides it.

### OAT parameter test

``` r
source("render_run.R")
render_oat("20261021")
```

This renders `Model4_OAT_Test.qmd`, saves its CSVs and plots to
`Outputs/Model4/OAT/20261021/` and produces
`20261021_Model4_OAT_Test.html` in the project root. As with
`render_run()`, re-using a label overwrites that run. Re-run it after
any change to the model structure or the fixed parameters.

### Interactive rendering (quick checks)

You can also click **Render** in RStudio, which uses the default
`run_label: "default"` from the YAML and saves outputs to
`Outputs/Model4/default/`. This is useful for quick checks but does not
produce a labelled output folder.

------------------------------------------------------------------------

## Model settings

- **Model run time**, **time step** and **burn-in** are set at the top
  of **Section 5** in the QMD. The default is a 6-year run (`365 * 6`
  days) at a half-day time step (`time_step = 0.5`), with a 1-year
  burn-in (`burn_in = 365`), giving a 5-year analysis period.
- The burn-in is the period while the model moves from its initial
  conditions to steady state. It is the same for every scenario and is
  excluded from all cumulative totals, patient-days and annual
  summaries, so these cover the analysis period only. `add_quantities()`
  gives rows in the burn-in (including time 0) zero width and multiplies
  every other output row by the time step. Sections 12–15 of the QMD
  report annual totals (the analysis-period total divided by its length
  in years), for comparison with surveillance data.
- All model time intervals (e.g. `t_test_turnaround`, length-of-stay)
  and transition rates use the same day-based timescale.
- Annual summaries use 365-day intervals counted from the end of the
  burn-in (year 1 is the first year of the analysis period).
- Stabilisation is assessed automatically for every scenario in the grid
  (states changing by ≤ 1e-4 per time step for 10 consecutive time
  steps; the `consecutive_days` setting counts time steps, not days, so
  this is 5 days at the half-day step). Section 9 checks that the
  burn-in is longer than every scenario's stabilisation period (with the
  defaults, all stabilise within 250 days) and stops the render with an
  error listing any that are not, so `burn_in` can be increased.
  Time-series plots show data up to the longest stabilisation period
  across the example scenarios, i.e. the approach to steady state during
  the burn-in.

------------------------------------------------------------------------

## Scenario definitions

Scenario control values are set as vectors in **Section 4** of
`Model4_Run.qmd`, and every combination is run (the scenario grid; by
default 3 policies × 12 values of the daily probability of testing (0,
1, 2, 3, 4, 5, 7, 10, 15, 20, 30 and 40%) × 11 sensitivity values = 396
scenarios, about 4 minutes on one core, about 30 seconds in
parallel). The three control variables
are:

- `isolation_policy`: which patients are isolated in a side room,
  numbered from least to most isolation (set by
  `isolation_policy_levels` in `Scripts/Model4_Functions.R`):

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

- `prob_test_daily`: daily probability that an untested CDI case (in I)
  is sent for testing (0 to below 1, 0 = no testing; multiples of 0.01,
  as scenario labels show it in hundredths); shown as a decimal (e.g. 0.03,
  `fmt_prob()`), not a percentage, so it is not confused with shares. It sets the case identification rate
  `gamma = prob2rate(prob_test_daily, 1) = -log(1 - prob_test_daily)`.
  Identification competes with discharge and death from I, so these are
  outputs (table in Sections 13 and 15, `testing_outputs.csv`):
  - the mean time to testing if not discharged or dead first,
    `1 / gamma` (about 2 days at 40%, 14 days at 7% and 100 days at 1%);
  - the actual mean time to testing among the cases tested,
    `1 / (gamma + disI + mortI)` (shorter than `1 / gamma`, as discharge
    and death cut off the longer waits);
  - the share of incident CDI cases ever sent for testing,
    `gamma / (gamma + disI + mortI)` (`prop_sent_for_testing`);
  - the share of prevalent CDI cases ever tested at model stabilisation
    (Outcome A4), by policy.

  The grid values are unevenly spaced, so the Section 14 plots show them
  as a categorical axis or legend. Mawer et al. (2019, doi
  10.1016/j.jhin.2019.05.001) found 75 of 230 hospital-onset diarrhoea
  patients in hospital (33%) had ever been tested for CDI. That is
  comparable in kind with Outcome A4, but its denominator is all diarrhoea, not
  CDI, so it is not used to calibrate `prob_test_daily`. Plots and
  tables label it "Daily probability of an untested CDI case being sent
  for testing".

- `test_sens`: diagnostic test sensitivity (0--1). The model uses it to
  calculate the confirmed case rate
  (`theta = test_sens / t_test_turnaround`) and false negative rate
  (`pi = (1 - test_sens) / t_test_turnaround`). Neither includes
  transfer time, so the share of results that are positive equals
  `test_sens` under every policy.

Test turnaround time (`t_test_turnaround`, days from test administration
to result) is a fixed parameter in `Scripts/Model4_Parameters.R`.

Scenario labels combine a policy code with the daily probability of
testing and the sensitivity in hundredths (whole percentages), e.g.
`P2_confirmed_test007_sens100` (`P1_none`, `P2_confirmed`,
`P3_susp_conf`; `test000` is no testing). `example_settings` in Section
3 lists the settings for the example scenarios, in the order shown in
Sections 10 and 12 (by default 0.03 and 0.40, about 33 and 2 days' mean time
to testing if not discharged or dead first, each at 50% and 100%
sensitivity), each run under all three policies (12
scenarios). They are used for time-series plots, annual summaries and
the outcome tables.

In Section 10, the time-series line and area plots show the policies as
columns and the settings as rows; in Section 12, the stacked bar plots
group the bars by setting, with P1, P2 and P3 side by side and the
setting label beneath each group.

Grid results are cached in `Outputs/Model4/<run_label>/grid_runs.rds`
and re-used on the next render if the grid, fixed parameters, run
settings, initial conditions and model/post-processing functions are
unchanged. Set `grid_clear_cache <- TRUE` in Section 7 to force a
re-run.

The grid runs in parallel (base R `parallel`) on `grid_n_cores` cores,
set in Section 7 (default: all physical cores but one; about 30 seconds
on 15 cores). Set `grid_n_cores <- 1` to run scenarios one at a time.

`isolation_policy` determines:

- the transmission parameters for `Xtest` and `Xtreat`: on a general
  ward, suspected cases (`Xtest`) transmit at `betaI` and confirmed
  cases on treatment (`Xtreat`) at the fixed `betaXtreat` (default 0.02:
  the average across the treatment period, as treatment is assumed to
  reduce shedding from `betaI` (0.04) at the start to 0 at the end);
  patients in a side room transmit at their general-ward value ×
  `beta_mult_isolated` (fixed parameter, default 0.05);
- the transfer states used (see **Patient transfer states**): none under
  `"none"`, T1 and T2 under `"confirmed"`, and T3, T4 and T5 under
  `"suspected_and_confirmed"`;
- side-room bed occupancy: none under `"none"`; `Xtreat` and `T2` under
  `"confirmed"`; and `Xtest`, `Xtreat`, `T4` and `T5` under
  `"suspected_and_confirmed"`.

Section 14 of the QMD (outcomes by daily probability of testing and test
sensitivity) compares the three isolation policies across the daily
probability of testing and test sensitivity. The outcomes are the research question outcomes defined in
`grid_outcome_spec` in `Scripts/Model4_Functions.R` (numbered as in the
research questions below):

| Outcome | Variable | Format |
|----|----|----|
| Outcome | Variable | Format | Section 15 `direction` |
|----|----|----|----|
| Outcome A1: Incident CDI cases | `cum_total_infected_incidence` | count per year | avoided |
| Outcome A2: Missed CDI cases | `cum_total_missed_cases` | count per year | avoided |
| Outcome A3: Incident CDI cases ever sent for testing | `prop_sent_for_testing` | proportion | change |
| Outcome A4: CDI cases ever tested, point prevalence at model stabilisation | `prop_current_CDI_tested` | proportion | change |
| Outcome A5: CDI cases confirmed | `cum_confirmed_cases` | count per year | additional |
| Outcome A5: Proportion of CDI cases confirmed | `prop_confirmed` | proportion | change |
| Outcome A6: CDI cases confirmed via retest (estimated) | `cum_confirmed_via_retest` | count per year | additional |
| Outcome A6: Proportion confirmed via retest (estimated) | `prop_confirmed_via_retest` | proportion | change |
| Outcome B1: Hospital-acquired colonisations | `cum_colonised_in_hospital` | count per year | avoided |
| Outcome B2: Hospital-acquired CDI cases (estimated) | `cum_hospital_attributable_I` | count per year | avoided |
| Outcome B3: Colonisation prevalence at discharge | `prop_discharged_C` | proportion | avoided |
| Outcome B3: Discharge minus admission prevalence of colonisation | `diff_prev_C` | percentage points | avoided |
| Outcome B3: CDI prevalence at discharge | `prop_discharged_CDI` | proportion | avoided |
| Outcome B3: Discharge minus admission prevalence of CDI | `diff_prev_CDI` | percentage points | avoided |
| Outcome C1: Total patient bed-days | `hospital_bed_days` | count per year | avoided |
| Outcome C1: Average length of stay | `avg_LOS` | days | avoided |
| Outcome C2: General ward patient bed-days | `general_ward_bed_days` | count per year | avoided |
| Outcome C3: Side-room patient bed-days | `side_room_bed_days` | count per year | additional |
| Outcome C4: Side-room share of patient bed-days | `prop_side_room_bed_days` | proportion | additional |
| Outcome C5: Faecal specimens tested | `cum_total_test` | count per year | avoided |
| Outcome C6: Faecal specimens tested per confirmed case | `tests_per_confirmed` | rate | avoided |
| Outcome D1: CDI cases treated and recovered | `cum_recovered_cases` | count per year | additional |
| Outcome D2: CDI cases discharged untreated | `cum_discharges_while_infected` | count per year | avoided |
| Outcome D3: CDI cases died in hospital | `cum_deaths_post_infection` | count per year | avoided |
| Outcome D4: CDI-attributable deaths in hospital | `cum_deaths_CDI_attributable` | count per year | avoided |
| Outcome D5: Mean time from CDI symptom onset to first test, among cases tested | `mean_wait_test` | days | avoided |
| Outcome D6: Mean time from CDI symptom onset to treatment, among cases treated | `mean_wait_treatment` | days | avoided |

Not plotted: colonisation and CDI prevalence at admission (Outcome B3; fixed
inputs, the same in every scenario) and confirmed cases started on
treatment (Outcome C7; equal to Outcome A5). Counts are annual totals over the analysis
period; proportions, prevalence differences, length of stay, specimens
per confirmed case and mean times are over the whole analysis period
(derived in `add_derived_outcomes()`). Differences
in total CDI cases (Outcome A1) are driven entirely by differences in
hospital-acquired CDI cases (Outcome B2), so their difference plots are
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
  daily probability of testing of 0.00 (no testing), 0.02, 0.05, 0.10, 0.20
  and 0.40 (`grid_line_sent` in `Scripts/Model4_Functions.R`; categorical
  legend, no testing yellow to most testing purple), faceted by policy
  with a shared y axis;
- the same plot for the differences from P1 (P2 vs P1, P3 vs P1) at the
  same daily probability of testing and sensitivity, placed under P2 and
  P3 (the space under P1 is left empty);
- heatmaps of test sensitivity (x) by daily probability of testing (y,
  one evenly spaced row per grid value, no testing at the bottom and the
  most testing at the top; tiles with no
  value, e.g. Outcome D5 and Outcome D6 with no testing, are grey), for every grid
  value, in one figure: the outcome under
  each policy on the top row (shared colour scale, yellow lowest, purple
  highest), and the differences from P1 (P2 vs P1, P3 vs P1) on the
  bottom row under P2 and P3, with the space under P1 left empty (yellow
  largest negative difference, purple largest positive difference).

The incremental effect of also isolating suspected cases (P3 vs P2) is
saved in `grid_key_outcome_differences.csv`. Each panel has its own y
axis so the rows line up.

At very low test sensitivity (about 10% or below), testing can increase
CDI, hospital-acquired CDI, missed cases and deaths compared with no
testing. Most tested patients get a false negative result, return to the
ward (FN) still transmitting, and are retested sooner than they would
typically be discharged. There is no discharge from Xtest, so they cycle
between FN and Xtest, held back from treatment and isolation and kept
infectious in hospital for longer. This is intended, to show the effect
of a low-sensitivity assay.

Section 15 of the QMD (change in outcomes due to increasing test
sensitivity) shows, for every outcome in Section 14 (`gain = TRUE` in
`grid_outcome_spec`), the change from increasing test sensitivity from
`sens_from` (default 50%) to `sens_to` (default 80%), by isolation
policy (P1, P2 and P3 reported separately) and daily probability of
testing. How each outcome is reported is set by the `direction` column
of `grid_outcome_spec` (table above), and the titles are built from it
in `sensitivity_gain_spec`:

- `"avoided"`: value at `sens_from` minus value at `sens_to`, so a
  negative value is an increase;
- `"additional"` (outcomes that rise as more cases are confirmed):
  value at `sens_to` minus value at `sens_from`;
- `"change"` (the shares of cases tested and confirmed): value at
  `sens_to` minus value at `sens_from`, in percentage points.

The arc elasticity of each outcome with respect to sensitivity uses the
midpoint method, `(ΔY / midpoint Y) / (ΔS / midpoint S)`: the % change in
the outcome for a 1% change in sensitivity, from the actual change in
the outcome (negative = falls as sensitivity rises). It is NA where the
outcome's midpoint is 0 and for the differences in prevalence
(Outcome B3), and is saved as `elasticity` in
`sensitivity_gain_<from>_to_<to>.csv` (`make_sensitivity_gain()`).
Section 17 compares the elasticities of all outcomes in one chart,
`plot_elasticities()`: outcomes in label order (Outcome A1 at the top),
one panel per daily probability of testing in `gain_table_sent`, a
point per policy (`policy_colours`), reference lines at 0 and ±1
(response in proportion to sensitivity) and solid lines between the
outcome groups A–D. Outcomes with no elasticity (the Outcome B3 differences in
prevalence) are left out.

Counts are per year; proportions and rates are over the analysis period.
Each outcome has a tab with a table of the values at `sens_from` and
`sens_to`, the change and the arc elasticity, by policy (rows ordered
by policy, then daily probability of testing of 0.03 and 0.40;
`gain_table_sent`, set in Section 14) and a chart,
`plot_sensitivity_levels()`, of the outcome itself over every daily
probability of testing in the grid: `sens_from` dashed and `sens_to`
solid, one panel per policy with a fixed y axis, the gap between them
shaded and a dotted line for no testing. The x axis is
`prob_test_daily` on a log scale from 0.01 (`sent_log_trans`; testing
increases left to right; no testing is the dotted horizontal line). The
testing outputs (shares of incident and prevalent cases ever tested,
mean times to first test and to treatment) are outputs, so they are
plotted as outcomes in their own Section 14 tabs rather than added to
this axis; the section also repeats the Section 13 table of testing
outputs. `sens_from` and `sens_to` must both be values
in the grid.

------------------------------------------------------------------------

## Validation and safeguards

Implemented safeguards include:

- explicit `stop()` checks for invalid inputs and parameters,
- warnings for potentially unrealistic values,
- non-negative state clamping in the ODE function,
- `N <- max(N, 1e-12)` to prevent divide-by-zero,

------------------------------------------------------------------------

## Potential model developments

Set out in Section 1 h) of `Model4_Run.qmd`:

- **Specimen collection and start of treatment while awaiting
  transfer (P2, P3):** count the test on entry to T3 and start
  treatment on entry to T1, with the mean time in Xtest shortened to
  `t_test_turnaround - t_patient_transfer` and in Xtreat to
  `t_treatment - t_patient_transfer`, so transfer no longer delays
  testing or treatment.
- **Environmental contamination that varies with prevalence:** a spore
  stock state `Env`, filled in proportion to bed-days in each infectious
  state (lower for C, Xtreat and isolated patients) and emptied at
  1 / persistence time (e.g. 180 days), with `beta0 = beta_env × Env`.
  Links environmental exposure to prevalence (with a lag), but hard to
  parameterise.
- Others: testing of non-CDI diarrhoea and colonised patients (with
  specificity), limited side-room capacity, recurrence and readmission,
  non-exponential times, and probabilistic analysis.

------------------------------------------------------------------------

## Research questions

The model is designed to focus on diagnostic testing of hospital
patients with CDI. It does not include testing of symptomatic patients
without CDI, or screening asymptomatic patients.

It can be used to quantify the effect of changes in testing algorithm
sensitivity for different proportions of symptomatic hospital patients
with CDI sent for testing and different hospital policies on location of
suspected cases awaiting test results and/or confirmed cases while
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

Where outcomes are described as estimates, they are calculated using the
shares of patient flows from one state to multiple states, or from
multiple states to one state. For example:

a)  cases confirmed through re-testing are estimated by apportioning the
    flow from Xtest to Xtreat by the share of flows into Xtest that come
    from FN (retests), rather than from I (first tests);
b)  hospital-acquired CDI is estimated by apportioning the flow from C
    to I by the share of all flows into C that come from S (colonisation
    in hospital), rather than from admissions to C.

### Outcomes of interest

Each outcome is evaluated across the full range of the control
variables (plots in QMD Section 14), and compared for an improvement in
test sensitivity from 50% to 80% for P1, P2 and P3 (plots across the
daily probability of testing, and tables at 0.03 and 0.40, in Section 15).

**A) Case ascertainment**

- Outcome A1: the total number of incident CDI cases (patients colonised on
  admission who progress to CDI during their stay, patients admitted
  with CDI and patients colonised in hospital who progress to CDI);
- Outcome A2: the number of missed CDI cases, untested or after a false negative
  result, discharged or dying;
- Outcome A3: the proportion of incident CDI cases that are ever sent for
  testing;
- Outcome A4: the proportion of CDI cases in hospital that have ever been tested
  for CDI, as a point prevalence at model stabilisation (for comparison
  with point-prevalence studies such as Mawer et al. 2019);
- Outcome A5: the number and proportion of incident CDI cases that are confirmed
  through diagnostic testing;
- Outcome A6: the estimated number and proportion of incident CDI cases that are
  confirmed through re-testing following a false-negative test result
  (delayed diagnosis).

**B) Transmission of *C. difficile***

- Outcome B1: the number of hospital-acquired colonisations;
- Outcome B2: the estimated number of hospital-acquired infections (patients
  colonised in hospital who progress to CDI during their stay);
- Outcome B3: the prevalence of colonisation and infection amongst discharged
  patients, and its difference from the prevalence amongst admitted
  patients.

**C) Hospital resource use**

- Outcome C1: total hospital bed occupancy (patient bed-days) and average
  patient length of stay;
- Outcome C2: general ward bed occupancy;
- Outcome C3: required side-room bed occupancy;
- Outcome C4: the proportion of hospital bed occupancy that is required for
  patient isolation in a side room;
- Outcome C5: the number of faecal specimens collected and tested for CDI;
- Outcome C6: the number of faecal specimens tested per CDI case confirmed;
- Outcome C7: the number of confirmed cases that are started on treatment (equal
  to Outcome A5).

**D) CDI case outcomes**

- Outcome D1: the number of incident CDI cases that are treated and recover
  (then either returned to a general ward or discharged);
- Outcome D2: the number of incident CDI cases that are discharged untreated;
- Outcome D3: the number of incident CDI cases that die in hospital (with no
  specimen sent for testing, awaiting a test result, awaiting re-testing
  following a false negative test result, or during treatment);
- Outcome D4: the number of CDI-attributable deaths in hospital (accounting for
  deaths that may have occurred due to other causes);
- Outcome D5: the mean time from CDI symptom onset to first test, among the
  cases tested;
- Outcome D6: the mean time from CDI symptom onset to the start of treatment,
  among the cases treated.

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
|----|----|
| Number of patient days on a general ward | × general ward care, cost per day |
| Number of patient days in a side room | × side room care, cost per day |
| Number of confirmed CDI cases | × CDI antimicrobial therapy, cost per course of treatment |
| Number of faecal specimens tested for CDI | × specimen collection and testing, cost per specimen |

| Model outcome (ICER denominator, avoided QALY loss or avoided hospital-acquired CDI case) | ICER denominator |
|----|----|
| Number of treated, recovered CDI cases | × QALY loss for treated cases |
| Number of untreated, discharged cases | × QALY loss for untreated cases |
| Number of CDI-attributable deaths | × discounted QALY loss for deaths |
| Number of hospital-acquired CDI cases |  |

Extending the model to include probabilistic analysis, with
distributions for uncertain parameters would allow value of information
analysis.

------------------------------------------------------------------------

## Research question outcomes table (QMD Section 13)

`make_research_outcomes_table()` in `Scripts/Model4_Functions.R` derives
every outcome of interest from the cumulative totals over the analysis
period after the burn-in (`make_final_outcomes()`), the patient-days
(`calculate_total_patient_days()`) and the derived proportions
(`add_derived_outcomes()`). With `per_year = TRUE` (as in the QMD)
counts and patient-days are divided by the length of the analysis
period. The table has one row per outcome, numbered as above, followed
by the ICER numerator and denominator inputs, and one column per
scenario. The QMD displays it for the example scenarios
(`research_outcomes.csv`) and saves it for every scenario in the grid
(`research_outcomes_all.csv`).

| Outcome | Calculation |
|----|----|
| Outcome A1: Incident CDI cases | `cum_total_infected_incidence` |
| Outcome A2: Missed CDI cases | `cum_total_missed_cases` |
| Outcome A3: Incident CDI cases ever sent for testing (%) | `prop_sent_for_testing` = `cum_suspected_cases / cum_total_infected_incidence` (equals `gamma / (gamma + disI + mortI)` at steady state) |
| Outcome A4: CDI cases ever tested, point prevalence (%) at model stabilisation^v^ | `prop_current_CDI_tested` = `CDI_tested_patient_days / CDI_patient_days` |
| Outcome A5: CDI cases confirmed, n and % | `cum_confirmed_cases`; `prop_confirmed` = `cum_confirmed_cases / cum_total_infected_incidence` |
| Outcome A6: Confirmed via retest (estimated)^i^, n and % | `cum_confirmed_via_retest`; `prop_confirmed_via_retest` = `cum_confirmed_via_retest / cum_total_infected_incidence` |
| Outcome B2: Hospital-acquired CDI cases (estimated)^ii^ | `cum_hospital_attributable_I` |
| Outcome B3: Colonisation prevalence at admission and discharge, and difference (pp)^iii^ | `cum_admitted_colonised / cum_total_admissions`; `cum_discharges_colonised / cum_total_discharges`; `diff_prev_C` = discharge − admission |
| Outcome B3: CDI prevalence at admission and discharge, and difference (pp)^iii^ | `cum_admitted_infected / cum_total_admissions`; `cum_discharges_while_infected / cum_total_discharges`; `diff_prev_CDI` = discharge − admission |
| Outcome C1: Patient bed-days and average length of stay | `hospital_bed_days`; `hospital_bed_days / cum_total_admissions` |
| Outcome C2: General ward patient bed-days | `general_ward_bed_days` |
| Outcome C3: Side-room patient bed-days | `side_room_bed_days` |
| Outcome C4: Side-room patient bed-days, % of all | `side_room_bed_days / hospital_bed_days` |
| Outcome C5: Faecal specimens tested | `cum_total_test` |
| Outcome C6: Faecal specimens tested per confirmed case | `tests_per_confirmed` = `cum_total_test / cum_confirmed_cases` ("-" if none confirmed) |
| Outcome C7: Confirmed cases started on treatment | `cum_confirmed_cases` (equal to Outcome A5) |
| Outcome D1: Treated and recovered | `cum_recovered_cases` |
| Outcome D2: Discharged untreated | `cum_discharges_while_infected` |
| Outcome D3: Died in hospital following CDI | `cum_deaths_post_infection` |
| Outcome D4: CDI-attributable deaths in hospital^iv^ | `cum_deaths_CDI_attributable` |
| Outcome D5: Mean time from CDI symptom onset to first test, among cases tested (days) | `mean_wait_test` (entry to I to first entry to Xtest; `mean_wait_to_state()`) |
| Outcome D6: Mean time from CDI symptom onset to treatment, among cases treated (days) | `mean_wait_treatment` (entry to I to entry to Xtreat, including false negatives and retests) |

ICER numerator (costs): general ward patient bed-days, side-room patient
bed-days, confirmed CDI cases (treatment courses), faecal specimens
tested. ICER denominator (outcomes): treated and recovered CDI cases,
untreated discharged CDI cases, CDI-attributable deaths,
hospital-acquired CDI cases (estimated).

**Notes**

The proportions confirmed (Outcome A5) and missed (Outcome A2) do not sum to 1: the gap
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
CDI patients (I, Xtest, FN, and T1, T3 and T4 awaiting transfer) × (1 −
1/`mortI_mult`) plus deaths of treated patients (Xtreat) × (1 −
1/`mortXtreat_mult`). In the grid run the multipliers are taken from
each scenario's parameters (`mortI / mortS`, `mortXtreat / mortS`).

^v^ `prop_current_CDI_tested` is the share of current CDI cases who have
ever been tested, the model equivalent of a point-prevalence survey such
as Mawer et al. (2019): CDI patient-days in a tested state (Xtest, FN,
Xtreat, T1, T4, and the share of T3 that came from FN) divided by all
CDI patient-days (the same states plus I and the rest of T3), summed
over the analysis period. The specimen is taken on entry to Xtest, so
patients in T3 who came from I have not yet been tested. T3 is not split
by origin, so its share from FN is estimated from the inflows to T3
(exact at steady state).

`make_final_outcomes()` also returns outcomes not in the table (e.g.
in-hospital colonisations by route and source, deaths by state), which
are saved in `final_outcomes_all.csv`.

------------------------------------------------------------------------

## Outputs

All outputs are written to `Outputs/Model4/<run_label>/` (created
automatically), including:

- `flow_P1.png`, `flow_P2.png`, `flow_P3.png` (flow diagram of each
  isolation policy; Section 1a)

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
  Section 14 and 15 outcomes, counts per year, and their differences
  between policies including P3 vs P2, in long format)
- `grid_rates.csv` (rates per 10,000 patient bed-days, side-room and
  missed proportions and missed cases per year by origin, long format)
- `research_outcomes_all.csv` (per year)
- `testing_outputs.csv` (Section 13: for each daily probability of
  testing, the mean time to testing if not discharged or dead first in
  days, the actual mean time to testing among cases tested in days, the
  share of incident CDI cases sent for testing and the point-prevalence share ever tested (Outcome A4) by policy, at
  one test sensitivity)
- `grid_lines_<outcome>.png`, `grid_lines_diff_<outcome>.png` and
  `grid_heatmap_<outcome>.png` (heatmaps and difference heatmaps in one
  figure; one set per Section 14 outcome)
- `sensitivity_gain_<from>_to_<to>.csv` and
  `sensitivity_levels_<from>_and_<to>_<outcome>.png` (values at both
  sensitivities) (Section 15, e.g. `sensitivity_gain_050_to_080.csv`)
- `sensitivity_elasticity_<from>_to_<to>.png` (arc elasticities of all
  outcomes; Section 17)

Example scenarios only:

- `solution_<scenario>.csv` (full time series)
- `annual_summary_all.csv` (every year of the analysis period; the QMD
  table shows year 1, as every year is the same at steady state)
- `research_outcomes.csv`
- time-series plots (Section 10): `plot_states_lines_facet.png`,
  `plot_states_lines_facet_noS.png`, `plot_states_area_facet.png`,
  `plot_states_area_facet_noS.png`, `plot_bed_occupancy.png` (area plot
  of general ward and side-room beds), `plot_flows_to_C_and_I.png`,
  `plot_case_ascertainment_line.png`,
  `plot_deaths_by_state_area_noS.png`,
  `plot_discharges_by_state_area_noS.png`
- stacked bar plots of annual totals (Section 12):
  `plot_colonisations_by_source.png` (background acquisition and
  patient-to-patient transmission from C, I, Xtest, FN and Xtreat),
  `plot_testing_bar.png` (first tests and retests),
  `plot_case_confirmation_bar.png` (confirmed via first test and
  retest), `plot_missed_cases_bar.png` (untested and false negatives)
  and `plot_case_outcomes_bar.png` (discharged untreated: untested from
  I and false negatives from FN; recovered after treatment; died
  following CDI) and `plot_bed_days_bar.png` (patient bed-days per year
  in general ward beds and in side rooms)

The rendered HTML (`<run_label>_Model4_Run.html`) embeds these plots, so
it can be moved or shared without the `Outputs/` folder.
