# Model 1 OAT (One-At-a-Time) Parameter Testing ####

# 1) Introduction ####

# This script runs a systematic one-at-a-time (OAT) sensitivity test to verify
# that the ODE model responds in the expected direction when each parameter is
# varied individually, with all others held at their scenario/fixed defaults.

# For each parameter, the model is run at two values (low and high) at a single baseline of test
# settings (50% testing coverage, 50% test sensitivity, so that untested, false-negative and
# retesting pathways are all active) under each of the three isolation policies:
#   "none", "confirmed" and "suspected_and_confirmed".
# Testing coverage (prop_I_suspected) and test sensitivity are not tested here, as the scenario
# grid in Model1_ODE_Model_Run.qmd covers their full range.
# The test specifications record a third reference column (baseline) for documentation,
# but model runs are only performed at low and high.

#   - low:      a value below the baseline (0 where valid; small positive otherwise)
#   - baseline: the default value in the fixed parameter list (reference only)
#   - high:     1 for proportions (prop_*), 5 Ã— baseline for all other parameters

# Key outputs are extracted from each run and the direction of change (high vs low)
# is compared against the expected direction. A PASS/FAIL flag is reported for each
# parameter Ã— isolation policy. The % change in the ten cumulative outcomes compared across the
# scenario grid (grid_outcome_spec) is also reported, for inspection only (no pass/fail).

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
# but long enough for the model to stabilise. As in the main analysis, the first year is a burn-in
# excluded from cumulative and annual outcomes (here a 2-year analysis period).
model_time <- 365 * 3
time_step  <- 0.5
burn_in    <- 365
mod.t      <- seq(0, model_time, by = time_step)

validate_run_settings(model_time, time_step, burn_in)

# Initial conditions (matching Model1_ODE_Model_Run.qmd)
mod.init <- c(S = 909, C = 86, I = 5, Xtest = 0, FN = 0, Xtreat = 0, D = 0)

# Scenarios: baseline test settings under each isolation policy, named by policy
oat_baseline <- list(prop_I_suspected = 0.50, test_sens = 0.50)

scenarios <- map(isolation_policy_levels,
                 \(pol) c(list(isolation_policy = pol), oat_baseline)) |>
  set_names(isolation_policy_levels)

# =========================================================== #

# 3) Wrapper function: run_with_overrides() ####

# Runs the model for a single set of parameter values under a given scenario.
#
# Arguments:
#   overrides - a named list of parameters to override. Any parameter accepted
#               by make_mod_parms() can be overridden (fixed or scenario-level).
#               All others default to the values in `scenario` plus the fixed list.
#   scenario  - a named list of scenario control parameters with elements:
#               isolation_policy, prop_I_suspected, test_sens.
#               Defaults to the suspected_and_confirmed policy at the baseline test settings.
#
# Returns a one-row tibble of summary statistics, plus the cumulative outcomes in
# grid_outcome_spec (totals over the whole run).

run_with_overrides <- function(overrides = list(), scenario = scenarios$suspected_and_confirmed) {

  base_args <- c(scenario, list(scenario_label = "oat_test"))
  all_args  <- modifyList(base_args, overrides)

  parms  <- do.call(make_mod_parms, all_args)
  result <- run_scenario(parms) |> add_quantities(time_step, burn_in = burn_in)

  # Steady-state values: use the last time point
  ss <- tail(result, 1)

  # Annual flows: use the last complete year
  annual  <- make_annual_summary(result, "oat_test")
  last_yr <- tail(annual, 1)

  # Cumulative outcomes over the run (as compared across the scenario grid)
  cum_outcomes <- make_final_outcomes(list(oat_test = result)) |>
    bind_cols(calculate_total_patient_days(result)) |>
    select(all_of(grid_outcome_spec$variable))

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
  ) |>
    bind_cols(cum_outcomes)
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
#   - All other parameters: high = 5 Ã— baseline.

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
  "beta0",                   0,       0.002,      0.010,   "increases",  "increases",        "increases",     "Higher background transmission -> more S->C -> more I",
  "betaC",                   0,       0.008,      0.040,   "increases",  "increases",        "increases",     "More transmission from C -> more S->C -> more I",
  "betaI",                   0,       0.040,      0.200,   "increases",  "increases",        "increases",     "More transmission from I -> more S->C -> more I",
  "betaXtreat",              0,       0.020,      0.100,   "increases",  "increases",        "increases",     "More transmission from treated confirmed cases (Xtreat) -> more S->C -> more I. Under the isolation policies Xtreat transmits at betaXtreat * beta_mult_isolated, so the effect is smaller",
  "beta_mult_isolated",      0,       0.05,       1,       "increases",  "increases",        "increases",     "Higher multiplier -> more transmission from isolated Xtest (betaI) and Xtreat (betaXtreat) -> more S->C -> more I. High = 1 (isolation has no effect on transmission)",
  
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
  "losXtreat",               1.00,    14.00,      70.00,   "complex",    "complex",          "complex",       "Longer Xtreat LOS -> more Xtreat occupancy -> more betaXtreat transmission -> more I (detectable when test_sens = 1). At the 50% sensitivity baseline FN transmission dominates and dilutes the betaXtreat effect below the 1% threshold. Note: confirmed_cases = theta*Xtest (inflow), not outflow. Low > 0 required",
  
  # --- Transfer and recovery ---
  "t_patient_transfer",      0.10,    0.50,       2.50,    "increases",  "complex",          "complex",       "Longer transfer -> slower gamma (I->Xtest) under suspected_and_confirmed -> I accumulates; also slows pi, sigma and delta (Xtreat->S). Discharge rates do not include transfer time. Low > 0 required",
  "prop_Xtreat_to_S",        0,       0.80,       1,       "complex",    "complex",          "complex",       "More return to ward vs discharge -> more S in hospital; net effect on I depends on transmission dynamics",
  
  # --- Mortality ---
  "prob_mort_S",             0,       0.001,      0.005,   "complex",    "complex",          "complex",       "Higher S mortality -> more exits from S -> more admissions needed; small indirect effect on C and I",
  "mortI_mult",              0,       2.00,       10.00,   "complex",    "complex",          "complex",       "Higher I mortality -> more I deaths counted as missed (missed = (disI+mortI)*I); detectable when test_sens = 1. At the 50% sensitivity baseline the FN pathway dominates total missed cases, diluting the extra I deaths below the 1% threshold. Effect on ss_I/confirmed also directionally negative but mortI <2% of disI, undetectable at 1%",
  "mortXtreat_mult",         0,       1.20,       6.00,    "complex",    "complex",          "complex",       "Higher Xtreat mortality multiplier -> lower Xtreat; small indirect effect on I",
  
  # --- Time to identification ---
  "t_identify_suspected",    0.50,    2.00,       10.00,   "increases",  "complex",          "complex",       "Slower identification -> I stays in I longer -> higher ss_I; effect on confirmed depends on test dynamics. Low > 0 required",
  
  # --- Test turnaround (fixed parameter; no longer a scenario control variable) ---
  "t_test_turnaround",       1.00,    2.00,       10.00,   "complex",    "complex",          "complex",       "suspected_and_confirmed: longer turnaround -> lower theta and pi, so Xtest accumulates, but isolated Xtest transmits at betaI * beta_mult_isolated, so the effect on I is small (<1% at baseline). Direction depends on beta_mult_isolated. See oat_specs_policy for the expectation under the confirmed and none policies. Low > 0 required"
)

## a2) Policy-specific expected directions ----

# Rows here override the expected directions (and notes) in oat_specs for runs under the named
# isolation policy. oat_specs gives the expectations for "suspected_and_confirmed"; parameters not
# listed here use the oat_specs expectations under every policy.

oat_specs_policy <- tribble(
  ~isolation_policy, ~param,                ~exp_ss_I,    ~exp_yr_confirmed,  ~exp_yr_missed,  ~notes,
  "confirmed",       "t_test_turnaround",   "increases",  "increases",        "complex",       "Confirmed policy: suspected cases wait on the ward and transmit at betaXtest = betaI. Longer turnaround -> more Xtest on the ward -> more S->C -> more I -> more confirmed cases. Missed cases are complex: more I increases untested missed cases, but when test_sens < 1 patients spend longer in Xtest (no discharge) and fewer are discharged as false negatives, so FN missed cases fall. Low > 0 required",
  "confirmed",       "t_patient_transfer",  "increases",  "complex",          "complex",       "Confirmed policy: gamma, pi and sigma do not include transfer time; longer transfer slows theta (Xtest stays on the ward longer, transmitting at betaI) and delta. Low > 0 required",
  "confirmed",       "beta_mult_isolated",  "increases",  "increases",        "increases",     "Confirmed policy: only Xtreat is isolated, so the multiplier applies to betaXtreat only -> more transmission from Xtreat -> more I",
  "none",            "t_test_turnaround",   "increases",  "increases",        "complex",       "No isolation: suspected cases wait on the ward and transmit at betaXtest = betaI, as under the confirmed policy. Longer turnaround -> more Xtest on the ward -> more S->C -> more I -> more confirmed cases. Missed cases complex (see confirmed policy). Low > 0 required",
  "none",            "t_patient_transfer",  "unchanged",  "unchanged",        "unchanged",     "No isolation: no patient moves between a general ward and a side room, so transfer time enters no rate",
  "none",            "beta_mult_isolated",  "unchanged",  "unchanged",        "unchanged",     "No isolation: no patients are in side rooms, so the isolation multiplier has no effect"
)

## b) Save OAT test specifications as CSV ----

write_csv(
  oat_specs |>
    select(param, low, baseline, high,
           exp_ss_I, exp_yr_confirmed, exp_yr_missed, notes),
  "Outputs/oat_specs.csv"
)
write_csv(oat_specs_policy, "Outputs/oat_specs_policy_overrides.csv")

message("OAT test specifications saved to Outputs/oat_specs.csv")

# =========================================================== #

# 5) Run OAT tests ####

# Each parameter is run at its low and high value at the baseline test settings under each
# isolation policy. All other parameters are held at the baseline test settings plus the fixed list
# (which includes t_test_turnaround). Expected directions come from oat_specs, with overrides from
# oat_specs_policy for the named policies.
# Suppress validation messages during the loop for readability.

oat_grid <- expand_grid(isolation_policy = names(scenarios), oat_specs) |>
  left_join(oat_specs_policy,
            by = c("isolation_policy", "param"), suffix = c("", "_override")) |>
  mutate(
    isolation_policy = factor(isolation_policy, levels = isolation_policy_levels),
    exp_ss_I         = coalesce(exp_ss_I_override,         exp_ss_I),
    exp_yr_confirmed = coalesce(exp_yr_confirmed_override, exp_yr_confirmed),
    exp_yr_missed    = coalesce(exp_yr_missed_override,    exp_yr_missed),
    notes            = coalesce(notes_override,            notes)
  ) |>
  select(-ends_with("_override"))

message("Running OAT tests under each isolation policy... (validation messages suppressed)")

oat_results <- oat_grid |>
  rowwise() |>
  mutate(
    result_low = list(suppressMessages(
      run_with_overrides(
        overrides = setNames(list(low),  param),
        scenario  = scenarios[[as.character(isolation_policy)]]
      )
    )),
    result_high = list(suppressMessages(
      run_with_overrides(
        overrides = setNames(list(high), param),
        scenario  = scenarios[[as.character(isolation_policy)]]
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

# Assess direction and pass/fail per parameter Ã— isolation policy
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

# Per parameter Ã— isolation policy results
oat_summary_long <- oat_assessed |>
  select(
    isolation_policy, param, low, baseline, high,
    ss_I_low, ss_I_high, exp_ss_I, act_ss_I, pass_ss_I,
    yr_confirmed_low, yr_confirmed_high, exp_yr_confirmed, act_yr_confirmed, pass_yr_confirmed,
    yr_missed_low, yr_missed_high, exp_yr_missed, act_yr_missed, pass_yr_missed,
    overall_pass, notes
  )

write_csv(oat_summary_long, "Outputs/oat_results.csv")

oat_summary <- oat_assessed |>
  select(isolation_policy, param, low, baseline, high,
         exp_ss_I, exp_yr_confirmed, exp_yr_missed,
         act_ss_I, act_yr_confirmed, act_yr_missed,
         overall_pass, notes)

oat_pass_counts <- oat_summary |>
  group_by(isolation_policy) |>
  summarise(n_params = n(), n_pass = sum(overall_pass), n_fail = n_params - n_pass,
            .groups = "drop")

n_params  <- nrow(oat_specs)
n_pass    <- sum(oat_summary$overall_pass)
n_fail    <- nrow(oat_summary) - n_pass
n_complex <- sum(oat_grid$exp_ss_I == "complex" |
                   oat_grid$exp_yr_confirmed == "complex" |
                   oat_grid$exp_yr_missed == "complex")

message("\n===== OAT Test Results (baseline coverage ", percent(oat_baseline$prop_I_suspected),
        ", sensitivity ", percent(oat_baseline$test_sens), "; all isolation policies) =====")
message("Parameters tested: ", n_params, "  |  Total model runs: ", nrow(oat_assessed) * 2)
message("PASS: ", n_pass, " / ", nrow(oat_summary), " parameter Ã— policy tests")
message("FAIL: ", n_fail, " / ", nrow(oat_summary))
message("Parameter Ã— policy tests with at least one 'complex' expected direction: ", n_complex,
        " (manual inspection recommended)")
print(oat_pass_counts)

failures <- oat_summary |> filter(!overall_pass)
if (nrow(failures) > 0) {
  message("\nFailed parameters:")
  print(failures |> select(isolation_policy, param, low, high, exp_ss_I, exp_yr_confirmed, exp_yr_missed,
                            act_ss_I, act_yr_confirmed, act_yr_missed))
} else {
  message("\nAll non-complex tests passed.")
}

# =========================================================== #

# 8) Visualise results ####

font_base_size <- 8

oat_subtitle <- paste0("Each parameter varied individually at ",
                       percent(oat_baseline$prop_I_suspected), " testing coverage and ",
                       percent(oat_baseline$test_sens), " test sensitivity")

# % change (high vs low). NA where the low value is 0 (% change undefined), e.g. side-room
# bed-days under the "none" policy.
pct_change <- function(val_high, val_low) {
  if_else(val_low == 0, NA_real_, 100 * (val_high - val_low) / abs(val_low))
}

# Compute % change for all three checked outputs
tornado_base <- oat_assessed |>
  mutate(
    pct_change_ss_I         = pct_change(ss_I_high,         ss_I_low),
    pct_change_yr_confirmed = pct_change(yr_confirmed_high, yr_confirmed_low),
    pct_change_yr_missed    = pct_change(yr_missed_high,    yr_missed_low)
  )

## a-c) Tornado charts: % change (high vs low), by isolation policy ----

# Parameters are ordered by their largest absolute % change across policies;
# bars coloured by isolation policy
plot_oat_tornado <- function(var, title) {
  ggplot(tornado_base |> mutate(param = fct_reorder(param, abs(.data[[var]]), .fun = max, .na_rm = TRUE)),
         aes(x = .data[[var]], y = param, fill = isolation_policy)) +
    geom_col(position = position_dodge()) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_discrete(labels = \(x) isolation_policy_labels[x]) +
    labs(
      title    = title,
      subtitle = oat_subtitle,
      x        = "% change (high vs low parameter value)",
      y        = NULL,
      fill     = "Isolation policy"
    ) +
    theme_classic(base_size = font_base_size) +
    theme(legend.position = "bottom")
}

plot_tornado_ssI          <- plot_oat_tornado("pct_change_ss_I",         "OAT sensitivity: % change in steady-state I")
plot_tornado_yr_confirmed <- plot_oat_tornado("pct_change_yr_confirmed", "OAT sensitivity: % chng. in annual confirmed cases")
plot_tornado_yr_missed    <- plot_oat_tornado("pct_change_yr_missed",    "OAT sensitivity: % chng. in annual missed cases")

## d) Pass/fail heatmap: param Ã— output check, by isolation policy ----

pass_plot_data <- oat_assessed |>
  select(isolation_policy, param, pass_ss_I, pass_yr_confirmed, pass_yr_missed) |>
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
  facet_wrap(~ isolation_policy, nrow = 1, labeller = as_labeller(isolation_policy_labels)) +
  labs(
    title    = "OAT test pass/fail, by isolation policy",
    subtitle = "Green = PASS (or 'complex' expected direction); Red = FAIL",
    x        = NULL,
    y        = NULL,
    fill     = NULL
  ) +
  theme_classic(base_size = font_base_size) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the 4 plots using cowplot

combined_oat_plot <- plot_grid(plot_tornado_ssI,
                               plot_tornado_yr_confirmed,
                               plot_tornado_yr_missed,
                               plot_heatmap_oat,
                               ncol = 2)

combined_oat_plot

# save as png

ggsave("Outputs/combined_oat_plot.png",
       combined_oat_plot, width = 48, height = 36, units = "cm", dpi = 300)

## e) Tornado chart: % change in the cumulative outcomes compared across the scenario grid ----

# For inspection only (no expected directions or pass/fail). One panel per outcome; parameters
# ordered by their largest absolute % change across outcomes and policies.

oat_outcomes_long <- map(grid_outcome_spec$variable, \(v)
  oat_assessed |>
    transmute(isolation_policy, param, low, high,
              outcome    = v,
              value_low  = map_dbl(result_low,  v),
              value_high = map_dbl(result_high, v))) |>
  bind_rows() |>
  mutate(
    pct_change = pct_change(value_high, value_low),
    outcome    = factor(outcome, levels = grid_outcome_spec$variable,
                        labels = grid_outcome_spec$label)
  )

write_csv(oat_outcomes_long, "Outputs/oat_outcomes_pct_change.csv")

plot_tornado_outcomes <- oat_outcomes_long |>
  mutate(param = fct_reorder(param, abs(pct_change), .fun = max, .na_rm = TRUE)) |>
  ggplot(aes(x = pct_change, y = param, fill = isolation_policy)) +
  geom_col(position = position_dodge()) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  facet_wrap(~ outcome, nrow = 2, scales = "free_x") +
  scale_fill_discrete(labels = \(x) isolation_policy_labels[x]) +
  labs(
    title    = "OAT sensitivity: % change in cumulative outcomes over the model run",
    subtitle = paste0(oat_subtitle, ". Bars missing where the low value is 0."),
    x        = "% change (high vs low parameter value)",
    y        = NULL,
    fill     = "Isolation policy"
  ) +
  theme_classic(base_size = font_base_size) +
  theme(legend.position = "bottom")

plot_tornado_outcomes

ggsave("Outputs/oat_tornado_cumulative_outcomes.png",
       plot_tornado_outcomes, width = 60, height = 30, units = "cm", dpi = 300)

# =========================================================== #
# END OF SCRIPT ####

