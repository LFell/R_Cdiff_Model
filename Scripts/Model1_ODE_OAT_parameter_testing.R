# Model 1 OAT (One-At-a-Time) Parameter Testing ####

# 1) Introduction ####

# This script runs a systematic one-at-a-time (OAT) sensitivity test to verify
# that the ODE model responds in the expected direction when each parameter is
# varied individually, with all others held at their scenario/fixed defaults.

# For each parameter, the model is run at two values (low and high) under each of
# the four scenarios defined in Model1_ODE_Model_Run.qmd (S1-S4). The test
# specifications record a third reference column (baseline) for documentation,
# but model runs are only performed at low and high.

#   - low:      a value below the baseline (0 where valid; small positive otherwise)
#   - baseline: the default value in the fixed parameter list (reference only)
#   - high:     1 for proportions (prop_*), 5 × baseline for all other parameters

# Key outputs are extracted from each run and the direction of change (high vs low)
# is compared against the expected direction. A PASS/FAIL flag is reported for each
# parameter × scenario. A parameter passes overall only if it passes in all four
# scenarios.

# IMPORTANT: This script is intended for model verification, not for analysis.
# Results here are not used to draw conclusions about the model's predictions.

# =========================================================== #

# 2) Setup ####

rm(list = ls())

library(tidyverse)
library(deSolve)
library(scales)
library(cowplot)

options(scipen = 1000)

source("Scripts/Model1_ODE_Functions.R")
source("Scripts/Model1_ODE_Parameters.R")

# Model run settings - use a shorter run than the main analysis to save time,
# but long enough for the model to stabilise.
model_time <- 365 * 3
time_step  <- 0.5
mod.t      <- seq(0, model_time, by = time_step)

validate_run_settings(model_time, time_step)

# Initial conditions (matching Model1_ODE_Model_Run.qmd)
mod.init <- c(S = 909, C = 86, I = 5, Xtest = 0, FN = 0, Xtreat = 0, D = 0)

# Scenarios: the four scenarios as defined in Model1_ODE_Model_Run.qmd
scenarios <- list(
  S1 = list(isolate_before_test = FALSE, prop_I_suspected = 1.00, test_sens = 1.00, t_test_turnaround = 2),
  S2 = list(isolate_before_test = FALSE, prop_I_suspected = 1.00, test_sens = 0.50, t_test_turnaround = 2),
  S3 = list(isolate_before_test = FALSE, prop_I_suspected = 0.50, test_sens = 1.00, t_test_turnaround = 2),
  S4 = list(isolate_before_test = FALSE, prop_I_suspected = 0.50, test_sens = 0.50, t_test_turnaround = 2)
)

# =========================================================== #

# 3) Wrapper function: run_with_overrides() ####

# Runs the model for a single set of parameter values under a given scenario.
#
# Arguments:
#   overrides - a named list of parameters to override. Any parameter accepted
#               by make_mod_parms() can be overridden (fixed or scenario-level).
#               All others default to the values in `scenario` plus the fixed list.
#   scenario  - a named list of scenario control parameters with elements:
#               isolate_before_test, prop_I_suspected, test_sens, t_test_turnaround.
#               Defaults to S1.
#
# Returns a one-row tibble of summary statistics.

run_with_overrides <- function(overrides = list(), scenario = scenarios$S1) {
  
  base_args <- c(scenario, list(scenario_label = "oat_test"))
  all_args  <- modifyList(base_args, overrides)
  
  parms  <- do.call(make_mod_parms, all_args)
  result <- run_scenario(parms) |> add_quantities(time_step)
  
  # Steady-state values: use the last time point
  ss <- tail(result, 1)
  
  # Annual flows: use the last complete year
  annual  <- make_annual_summary(result, "oat_test")
  last_yr <- tail(annual, 1)
  
  tibble(
    # Steady-state compartment occupancy
    ss_S      = ss$S,
    ss_C      = ss$C,
    ss_I      = ss$I,
    ss_Xtest  = ss$Xtest,
    ss_Xtreat = ss$Xtreat,
    ss_lambda = ss$lambda,
    # Last-year cumulative flows
    yr_infected_incidence = last_yr$annual_total_infected_incidence,
    yr_confirmed_cases    = last_yr$annual_confirmed_cases,
    yr_missed_cases       = last_yr$annual_total_missed_cases,
    yr_side_room_days     = last_yr$annual_side_room_bed_days,
    yr_total_deaths       = last_yr$annual_total_deaths
  )
}

# =========================================================== #

# 4) OAT test specifications ####

# Each row defines one parameter test.
#
# Columns:
#   param             - parameter name (as passed to make_mod_parms())
#   low / baseline / high - values to test
#   exp_ss_I          - expected direction of change in steady-state I, comparing
#                       HIGH vs LOW: "increases", "decreases", or "complex"
#   exp_yr_confirmed  - expected direction for annual confirmed cases (high vs low)
#   exp_yr_missed     - expected direction for annual missed cases (high vs low)
#   notes             - brief rationale

# Notes on "complex" entries:
#   "complex" is used where the direction of effect is theoretically ambiguous
#   (e.g. competing mechanisms). These cases should be inspected manually.

# Notes on low values:
#   - Proportions (prop_*) and rates/betas that are valid at 0 use low = 0.
#   - Parameters that must be strictly > 0 (time and LOS parameters, abx_mult) use a
#     small but valid positive value; 0 would fail validation in make_mod_parms().
#   - abx_mult low = 1.0 (no antibiotic-mediated increase in progression risk).
#
# Notes on high values:
#   - Proportions (prop_*): high = 1 (upper bound).
#   - All other parameters: high = 5 × baseline.

# Note on ~notes column width: some notes are long; the CSV output (Section 4b)
# is the intended human-readable reference; the tribble is the source of truth.

oat_specs <- tribble(
  ~param,                   ~low,    ~baseline,  ~high,   ~exp_ss_I,    ~exp_yr_confirmed,  ~exp_yr_missed,  ~notes,
  
  # ~exp_ss_I means the expected direction of change in steady-state I when comparing HIGH vs LOW parameter value.
  # ~exp_yr_confirmed means the expected direction of change in annual confirmed cases when comparing HIGH vs LOW parameter value.
  # ~exp_yr_missed means the expected direction of change in annual missed cases when comparing HIGH vs LOW parameter value.
  
  # --- Admission prevalences ---
  # Note: prop_adm_C + prop_adm_I must be <= 1. High values are capped so that each
  # parameter's high + the other's baseline remains <= 1.
  "prop_adm_C",              0,       0.086,      0.99,    "increases",  "increases",        "increases",     "More colonised admissions -> more C -> more progression to I. High capped at 0.99 (joint constraint with prop_adm_I baseline)",
  "prop_adm_I",              0,       0.002,      0.90,    "increases",  "increases",        "increases",     "More infected admissions -> directly increases I. High capped at 0.90 (joint constraint with prop_adm_C baseline)",
  
  # --- Transmission parameters ---
  "beta0",                   0,       0.005,      0.025,   "increases",  "increases",        "increases",     "Higher background transmission -> more S->C -> more I",
  "betaC",                   0,       0.008,      0.040,   "increases",  "increases",        "increases",     "More transmission from C -> more S->C -> more I",
  "betaI",                   0,       0.040,      0.200,   "increases",  "increases",        "increases",     "More transmission from I -> more S->C -> more I",
  "betaXtreat",              0,       0.004,      0.020,   "increases",  "increases",        "increases",     "More transmission from Xtreat -> more S->C -> more I",
  
  # --- Progression C -> I ---
  "prop_progress_C_to_I",    0,       0.50,       1,       "increases",  "increases",        "increases",     "Higher progression proportion -> higher alpha -> more I",
  "abx_prop",                0,       0.373,      1,       "increases",  "increases",        "increases",     "More C on antibiotics -> higher effective alpha -> more I",
  "abx_mult",                1.00,    2.00,       10.00,   "increases",  "increases",        "increases",     "Higher antibiotic RR -> higher effective alpha -> more I. Low = 1.0 (no abx effect; must be > 0)",
  "t_progress_C_to_I",       0.50,    5.00,       25.00,   "decreases",  "decreases",        "decreases",     "Slower progression (lower alpha) -> C discharged before progressing. Low > 0 required",
  
  # --- Decolonisation C -> S ---
  "prop_decolonise_C_to_S",  0,       0.05,       1,       "decreases",  "decreases",        "decreases",     "More decolonisation -> lower C -> lower alpha*C -> less I",
  "t_decolonise_C_to_S",     0.50,    5.00,       25.00,   "increases",  "increases",        "increases",     "Slower decolonisation (lower mu) -> C persists longer -> more I. Low > 0 required",
  
  # --- Length of stay ---
  "losS",                    1.00,    7.00,       35.00,   "complex",    "complex",          "complex",       "Longer S LOS: lower throughput vs longer exposure time; net effect unclear. Low > 0 required",
  "losC",                    1.00,    7.00,       35.00,   "increases",  "increases",        "increases",     "Longer C LOS -> more C in hospital at any time -> more I. Low > 0 required",
  "losI",                    1.00,    14.00,      70.00,   "increases",  "increases",        "decreases",     "Longer I LOS -> lower disI -> smaller missed fraction (disI+mortI)/(gamma+disI+mortI) -> fewer I exit before identification. ss_I and confirmed both increase. Low > 0 required",
  "losXtreat",               1.00,    14.00,      70.00,   "complex",    "complex",          "complex",       "Longer Xtreat LOS -> more Xtreat occupancy -> more betaXtreat transmission -> more I (detectable in S1/S3 where test_sens=1). In S2/S4 (test_sens=0.5) FN transmission dominates and dilutes betaXtreat effect below 1% threshold. Note: confirmed_cases = theta*Xtest (inflow), not outflow. Low > 0 required",
  
  # --- Transfer and recovery ---
  "t_patient_transfer",      0.10,    0.50,       2.50,    "increases",  "complex",          "complex",       "Longer transfer -> slower gamma (I->Xtest) -> I accumulates; also affects discharge and recovery rates. Low > 0 required",
  "prop_Xtreat_to_S",        0,       0.80,       1,       "complex",    "complex",          "complex",       "More return to ward vs discharge -> more S in hospital; net effect on I depends on transmission dynamics",
  
  # --- Mortality ---
  "prob_mort_S",             0,       0.001,      0.005,   "complex",    "complex",          "complex",       "Higher S mortality -> more exits from S -> more admissions needed; small indirect effect on C and I",
  "mortI_mult",              0,       2.00,       10.00,   "complex",    "complex",          "complex",       "Higher I mortality -> more I deaths counted as missed (missed = (disI+mortI)*I); detectable in S1/S3 (test_sens=1). In S2/S4 (test_sens=0.5) FN pathway dominates total missed cases, diluting the extra I deaths below 1% threshold. Effect on ss_I/confirmed also directionally negative but mortI <2% of disI, undetectable at 1%",
  "mortXtreat_mult",         0,       1.20,       6.00,    "complex",    "complex",          "complex",       "Higher Xtreat mortality multiplier -> lower Xtreat; small indirect effect on I",
  
  # --- Time to identification ---
  "t_identify_suspected",    0.50,    2.00,       10.00,   "increases",  "complex",          "complex",       "Slower identification -> I stays in I longer -> higher ss_I; effect on confirmed depends on test dynamics. Low > 0 required",
  
  # --- Scenario-specific parameters ---
  "prop_I_suspected",        0,       1.00,       1,       "decreases",  "increases",        "decreases",     "HIGH(=1) vs LOW(=0): at low=0 gamma=0 so confirmed=0 and all infected are missed; at high=1 gamma>0 so confirmed>0 and fewer missed. High = baseline so only low tested.",
  "test_sens",               0,       1.00,       1,       "decreases",  "increases",        "decreases",     "Lower sensitivity -> more false negatives return to I; fewer confirmed, more missed. High = baseline.",
  "t_test_turnaround",       1.00,    2.00,       10.00,   "complex",    "complex",          "complex",       "Longer turnaround -> lower theta but large Xtest accumulation; betaXtest=betaI (no pre-test isolation) so Xtest drives transmission -> more I -> more confirmed, offsetting slower processing. Net direction of confirmed depends on betaXtest magnitude. Low > 0 required"
)

## b) Save OAT test specifications as CSV ----

write_csv(
  oat_specs |>
    select(param, low, baseline, high,
           exp_ss_I, exp_yr_confirmed, exp_yr_missed, notes),
  "Outputs/oat_specs.csv"
)

message("OAT test specifications saved to Outputs/oat_specs.csv")

# =========================================================== #

# 5) Run OAT tests ####

# For each parameter, run the model at low and high values under scenario S1.
# All other parameters are held at S1 defaults plus the fixed list.
# Scenario-specific parameters (prop_I_suspected, test_sens, t_test_turnaround)
# are tested by overriding their S1 baseline values directly.
# Suppress validation messages during the loop for readability.

message("Running OAT tests under scenario S1... (validation messages suppressed)")

oat_results <- oat_specs |>
  rowwise() |>
  mutate(
    result_low = list(suppressMessages(
      run_with_overrides(
        overrides = setNames(list(low),  param),
        scenario  = scenarios$S1
      )
    )),
    result_high = list(suppressMessages(
      run_with_overrides(
        overrides = setNames(list(high), param),
        scenario  = scenarios$S1
      )
    ))
  ) |>
  ungroup()

# =========================================================== #

# 6) Assess direction of change and compare to expected ####

# Helper: classify direction of change as "increases", "decreases", or "unchanged",
# with a relative tolerance of 1%.
classify_direction <- function(val_high, val_low, tol = 0.01) {
  if (is.na(val_high) | is.na(val_low)) return(NA_character_)
  if (val_low == 0 & val_high == 0)     return("unchanged")
  rel_change <- (val_high - val_low) / (abs(val_low) + 1e-12)
  if      (rel_change >  tol) "increases"
  else if (rel_change < -tol) "decreases"
  else                        "unchanged"
}

# Assess direction and pass/fail per parameter × scenario
oat_assessed <- oat_results |>
  mutate(
    ss_I_low          = map_dbl(result_low,  "ss_I"),
    ss_I_high         = map_dbl(result_high, "ss_I"),
    yr_confirmed_low  = map_dbl(result_low,  "yr_confirmed_cases"),
    yr_confirmed_high = map_dbl(result_high, "yr_confirmed_cases"),
    yr_missed_low     = map_dbl(result_low,  "yr_missed_cases"),
    yr_missed_high    = map_dbl(result_high, "yr_missed_cases"),
    
    act_ss_I         = map2_chr(ss_I_high,         ss_I_low,         classify_direction),
    act_yr_confirmed = map2_chr(yr_confirmed_high,  yr_confirmed_low, classify_direction),
    act_yr_missed    = map2_chr(yr_missed_high,     yr_missed_low,    classify_direction),
    
    pass_ss_I         = (exp_ss_I         == "complex") | (act_ss_I         == exp_ss_I),
    pass_yr_confirmed = (exp_yr_confirmed  == "complex") | (act_yr_confirmed == exp_yr_confirmed),
    pass_yr_missed    = (exp_yr_missed     == "complex") | (act_yr_missed    == exp_yr_missed),
    
    overall_pass = pass_ss_I & pass_yr_confirmed & pass_yr_missed
  )

# =========================================================== #

# 7) Summary table ####

# Per parameter results (single scenario: S1)
oat_summary_long <- oat_assessed |>
  select(
    param, low, baseline, high,
    ss_I_low, ss_I_high, exp_ss_I, act_ss_I, pass_ss_I,
    yr_confirmed_low, yr_confirmed_high, exp_yr_confirmed, act_yr_confirmed, pass_yr_confirmed,
    yr_missed_low, yr_missed_high, exp_yr_missed, act_yr_missed, pass_yr_missed,
    overall_pass, notes
  )

oat_summary <- oat_assessed |>
  select(param, low, baseline, high,
         exp_ss_I, exp_yr_confirmed, exp_yr_missed,
         act_ss_I, act_yr_confirmed, act_yr_missed,
         overall_pass, notes)

n_params  <- nrow(oat_summary)
n_pass    <- sum(oat_summary$overall_pass)
n_fail    <- n_params - n_pass
n_complex <- sum(oat_specs$exp_ss_I == "complex" |
                   oat_specs$exp_yr_confirmed == "complex" |
                   oat_specs$exp_yr_missed == "complex")

message("\n===== OAT Test Results (S1) =====")
message("Parameters tested: ", n_params, "  |  Total model runs: ", nrow(oat_assessed) * 2)
message("PASS: ", n_pass, " / ", n_params)
message("FAIL: ", n_fail, " / ", n_params)
message("Parameters with at least one 'complex' expected direction: ", n_complex,
        " (manual inspection recommended)")

failures <- oat_summary |> filter(!overall_pass)
if (nrow(failures) > 0) {
  message("\nFailed parameters:")
  print(failures |> select(param, low, high, exp_ss_I, exp_yr_confirmed, exp_yr_missed,
                            act_ss_I, act_yr_confirmed, act_yr_missed))
} else {
  message("\nAll non-complex tests passed.")
}

# =========================================================== #

# 8) Visualise results ####

font_base_size <- 8

# Compute % change for all three outputs, order each plot independently
tornado_base <- oat_assessed |>
  mutate(
    pct_change_ss_I         = 100 * (ss_I_high         - ss_I_low)         / (ss_I_low         + 1e-12),
    pct_change_yr_confirmed = 100 * (yr_confirmed_high - yr_confirmed_low) / (yr_confirmed_low + 1e-12),
    pct_change_yr_missed    = 100 * (yr_missed_high    - yr_missed_low)    / (yr_missed_low    + 1e-12)
  )

## a) Tornado chart: % change in steady-state I (high vs low) ----

plot_tornado_ssI <- ggplot(tornado_base |>
         mutate(param = fct_reorder(param, abs(pct_change_ss_I))),
       aes(x = pct_change_ss_I, y = param)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  labs(
    title    = "OAT sensitivity: % change in steady-state I",
    subtitle = "Each parameter varied individually under scenario S1",
    x        = "% change (high vs low parameter value)",
    y        = NULL
  ) +
  theme_classic(base_size = font_base_size)

## b) Tornado chart: % change in annual confirmed cases ----

plot_tornado_yr_confirmed <- ggplot(tornado_base |>
         mutate(param = fct_reorder(param, abs(pct_change_yr_confirmed))),
       aes(x = pct_change_yr_confirmed, y = param)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  labs(
    title    = "OAT sensitivity: % chng. in annual confirmed cases",
    subtitle = "Each parameter varied individually under scenario S1",
    x        = "% change (high vs low parameter value)",
    y        = NULL
  ) +
  theme_classic(base_size = font_base_size)

## c) Tornado chart: % change in annual missed cases (high vs low) ----

plot_tornado_yr_missed <- ggplot(tornado_base |>
         mutate(param = fct_reorder(param, abs(pct_change_yr_missed))),
       aes(x = pct_change_yr_missed, y = param)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  labs(
    title    = "OAT sensitivity: % chng. in annual missed cases",
    subtitle = "Each parameter varied individually under scenario S1",
    x        = "% change (high vs low parameter value)",
    y        = NULL
  ) +
  theme_classic(base_size = font_base_size)

## d) Pass/fail heatmap: param × output check ----

pass_plot_data <- oat_assessed |>
  select(param, pass_ss_I, pass_yr_confirmed, pass_yr_missed) |>
  pivot_longer(
    cols      = starts_with("pass"),
    names_to  = "check",
    values_to = "passed"
  ) |>
  mutate(
    check  = recode(check,
                    pass_ss_I         = "Steady-state I",
                    pass_yr_confirmed = "Annual confirmed",
                    pass_yr_missed    = "Annual missed"),
    result = if_else(passed, "PASS", "FAIL")
  )

plot_heatmap_oat <- ggplot(pass_plot_data,
       aes(x = check, y = param, fill = result)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(PASS = "#4DAF4A", FAIL = "#E41A1C")) +
  labs(
    title    = "OAT test pass/fail (scenario S1)",
    subtitle = "Green = PASS (or 'complex' expected direction); Red = FAIL",
    x        = NULL,
    y        = NULL,
    fill     = NULL
  ) +
  theme_classic(base_size = font_base_size)

# Combine the 4 plots using cowplot
  
combined_oat_plot <- plot_grid(plot_tornado_ssI, 
                               plot_tornado_yr_confirmed,
                               plot_tornado_yr_missed,
                               plot_heatmap_oat,
                               ncol = 2)

combined_oat_plot

# save as png

ggsave(("Outputs/combined_oat_plot.png"),
             combined_oat_plot, width = 32, height = 18, units = "cm", dpi = 300)

# =========================================================== #
# END OF SCRIPT ####

