# Run CPE Compartmental Model ####

# 1) Initial setup ####

rm(list = ls())

library(tidyverse)
library(deSolve)
library(scales)

options(scipen = 1000)

# Ensure output folders exist
dir.create("Outputs", showWarnings = FALSE, recursive = TRUE)

# Source model files
source("Scripts/CPE_ODE_Functions.R")
source("Scripts/CPE_ODE_Parameters.R")

# ========================================================== #

# 2) Model run settings ####

model.time <- 365 * 6   # 5-year horizon + 1-year burn-in
time.step <- 1
burn_in_days <- 365
mod.t <- seq(from = 0, to = model.time, by = time.step)

# ========================================================== #

# 3) Initial conditions ####

Sg_0 <- 640
Cg_0 <- 0
Ig_0 <- 0
St_0 <- 0
Ct_0 <- 0
It_0 <- 0
Sx_0 <- 0
Cx_0 <- 0
Ix_0 <- 0
Sy_0 <- 0
Cy_0 <- 0
Iy_0 <- 0

Tot_hos_0 <- sum(Sg_0, Cg_0, Ig_0, St_0, Ct_0, It_0, Sx_0, Cx_0, Ix_0, Sy_0, Cy_0, Iy_0)

Sdxa_0 <- 0
Cdxa_0 <- 0
Sdxb_0 <- 0
Cdxb_0 <- 0
Sda_0 <- 0
Cda_0 <- 0
Sdb_0 <- 0
Cdb_0 <- 0

Tot_com_0 <- sum(Sdxa_0, Cdxa_0, Sdxb_0, Cdxb_0, Sda_0, Cda_0, Sdb_0, Cdb_0)

mod.init <- c(
  Sg = Sg_0, Cg = Cg_0, Ig = Ig_0,
  St = St_0, Ct = Ct_0, It = It_0,
  Sx = Sx_0, Cx = Cx_0, Ix = Ix_0,
  Sy = Sy_0, Cy = Cy_0, Iy = Iy_0,
  Tot_hos = Tot_hos_0,
  Sdxa = Sdxa_0, Cdxa = Cdxa_0, Sdxb = Sdxb_0, Cdxb = Cdxb_0,
  Sda = Sda_0, Cda = Cda_0, Sdb = Sdb_0, Cdb = Cdb_0,
  Tot_com = Tot_com_0
)

# ========================================================== #

# 4) Utility functions ####

run_scenario <- function(parms) {
  as.data.frame(ode(
    y = mod.init,
    times = mod.t,
    func = CPE.ODE.model,
    parms = parms
  ))
}

apply_burn_in <- function(df, burn_days = 365) {
  df %>%
    filter(time > burn_days) %>%
    mutate(time = time - burn_days)
}

add_derived_outputs <- function(df) {
  df %>%
    mutate(year = if_else(time == 0, 1, ((time - 1) %/% 365) + 1)) %>%
    relocate(year, .before = time) %>%
    mutate(
      new_case_ci = new_case_c + new_case_i,
      new_dis_total = new_dis_s + new_dis_c + new_dis_i,
      new_mort_hosp_total = new_mort_hosp_s + new_mort_hosp_c + new_mort_hosp_i,
      new_readmit30_total = new_readmit30_s_flag + new_readmit30_c_flag + new_readmit30_s + new_readmit30_c,
      new_readmit365_total = new_readmit365_s_flag + new_readmit365_c_flag + new_readmit365_s + new_readmit365_c,
      new_readmit_total = new_readmit30_total + new_readmit365_total,
      new_readmit_flag_total = new_readmit30_s_flag + new_readmit30_c_flag + new_readmit365_s_flag + new_readmit365_c_flag,
      new_admit_com_total = new_admit_com_s + new_admit_com_c + new_admit_com_i,
      hos_daily_costs = bed_c + screen_test_c + screen_staff_c + inf_trt_c + inf_trt_meet_c +
        tox_initial_test_c + tox_ongoing_test_c + ppe_equip_c + ipc_staff_c + room_clean_c +
        stock_disp_c + infect_waste_c + mort_rev_c + outbreak_c
    ) %>%
    group_by(year) %>%
    mutate(
      cum_case_s = cumsum(new_case_s),
      cum_case_c = cumsum(new_case_c),
      cum_new_case_known_c = cumsum(new_case_known_c),
      cum_case_i = cumsum(new_case_i),
      cum_new_case_known_i = cumsum(new_case_known_i),
      cum_case_ci = cumsum(new_case_ci),
      cum_dis_s = cumsum(new_dis_s),
      cum_dis_c = cumsum(new_dis_c),
      cum_dis_i = cumsum(new_dis_i),
      cum_dis_total = cumsum(new_dis_total),
      cum_mort_hosp_s = cumsum(new_mort_hosp_s),
      cum_mort_hosp_c = cumsum(new_mort_hosp_c),
      cum_mort_hosp_i = cumsum(new_mort_hosp_i),
      cum_mort_hosp_total = cumsum(new_mort_hosp_total),
      cum_readmit30_s_flag = cumsum(new_readmit30_s_flag),
      cum_readmit30_c_flag = cumsum(new_readmit30_c_flag),
      cum_readmit30_s = cumsum(new_readmit30_s),
      cum_readmit30_c = cumsum(new_readmit30_c),
      cum_readmit30_total = cumsum(new_readmit30_total),
      cum_readmit365_s_flag = cumsum(new_readmit365_s_flag),
      cum_readmit365_c_flag = cumsum(new_readmit365_c_flag),
      cum_readmit365_s = cumsum(new_readmit365_s),
      cum_readmit365_c = cumsum(new_readmit365_c),
      cum_readmit365_total = cumsum(new_readmit365_total),
      cum_readmit_total = cumsum(new_readmit_total),
      cum_readmit_flag_total = cumsum(new_readmit_flag_total),
      cum_admit_com_s = cumsum(new_admit_com_s),
      cum_admit_com_c = cumsum(new_admit_com_c),
      cum_admit_com_i = cumsum(new_admit_com_i),
      cum_admit_com_total = cumsum(new_admit_com_total),
      cum_hos_daily_costs = cumsum(hos_daily_costs)
    ) %>%
    ungroup()
}

negative_value_rows <- function(df) {
  df %>% filter(if_any(where(is.numeric), ~ .x < 0))
}

negative_value_summary <- function(df) {
  df %>%
    summarise(across(where(is.numeric), ~ sum(.x < 0, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "num_negative_values") %>%
    filter(num_negative_values > 0)
}

summarise_outcomes <- function(df) {
  list(
    infect.days.total = sum(df$Ig, df$It, df$Ix, df$Iy),
    new.colon.import.total = sum(df$new_case_import_c),
    new.colon.acquired.total = sum(df$new_case_acquired_c),
    new.colon.total = sum(df$new_case_c),
    new.infect.import.total = sum(df$new_case_import_i),
    new.infect.acquired.total = sum(df$new_case_acquired_i),
    new.infect.total = sum(df$new_case_i),
    susc.deaths.total = sum(df$new_mort_hosp_s),
    colon.deaths.total = sum(df$new_mort_hosp_c),
    infect.deaths.total = sum(df$new_mort_hosp_i),
    deaths.total = sum(df$new_mort_hosp_s, df$new_mort_hosp_c, df$new_mort_hosp_i),
    discharged.total = sum(df$D),
    readmit.total = sum(df$R),
    admit.total = sum(df$H),
    screen.tests.total = sum(df$screen_tests),
    clinical.tests.total = sum(df$clinical_tests),
    hos.costs.total = sum(df$hos_daily_costs)
  )
}

plot_hospital_compartments <- function(df, title_text, out_file) {
  dat <- df %>%
    select(time, Sg, Cg, Ig, St, Ct, It, Sx, Cx, Ix, Sy, Cy, Iy) %>%
    pivot_longer(cols = -time, names_to = "compartment", values_to = "count") %>%
    mutate(compartment = factor(compartment, levels = c("Sg","Cg","Ig","St","Ct","It","Sx","Cx","Ix","Sy","Cy","Iy")))
  
  p <- ggplot(dat, aes(x = time, y = count, fill = compartment)) +
    geom_area() +
    labs(x = "Time (days)", y = "Number of patients", fill = "Compartment", title = title_text) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  ggsave(out_file, plot = p, width = 10, height = 6, dpi = 300)
  p
}

plot_community_compartments <- function(df, title_text, out_file) {
  dat <- df %>%
    select(time, Sdxa, Cdxa, Sdxb, Cdxb, Sda, Cda, Sdb, Cdb) %>%
    pivot_longer(cols = -time, names_to = "compartment", values_to = "count") %>%
    mutate(compartment = factor(compartment, levels = c("Sdxa","Cdxa","Sdxb","Cdxb","Sda","Cda","Sdb","Cdb")))
  
  p <- ggplot(dat, aes(x = time, y = count, fill = compartment)) +
    geom_area() +
    labs(x = "Time (days)", y = "Number of individuals", fill = "Compartment", title = title_text) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  ggsave(out_file, plot = p, width = 10, height = 6, dpi = 300)
  p
}

# ========================================================== #

# 5) Run scenarios (i0, i1, i2) ####

scenario_parms <- list(
  i0 = mod.parms.i0,
  i1 = mod.parms.i1,
  i2 = mod.parms.i2
)

scenario_raw <- lapply(scenario_parms, run_scenario)
scenario_abso <- lapply(scenario_raw, function(x) x %>% apply_burn_in(burn_in_days) %>% add_derived_outputs())

# Keep original object names for compatibility
CPE.ODE.solution.i0.abso <- scenario_abso$i0
CPE.ODE.solution.i1.abso <- scenario_abso$i1
CPE.ODE.solution.i2.abso <- scenario_abso$i2

# Save full outputs
write_csv(CPE.ODE.solution.i0.abso, "Outputs/CPE_ODE_solution_i0_abso.csv")
write_csv(CPE.ODE.solution.i1.abso, "Outputs/CPE_ODE_solution_i1_abso.csv")
write_csv(CPE.ODE.solution.i2.abso, "Outputs/CPE_ODE_solution_i2_abso.csv")

# ========================================================== #

# 6) Annual summary tables (after burn-in) ####

# Purpose:
# Create compact annual cumulative metrics for quick plausibility checks and
# cross-scenario comparison (i0 baseline, i1, i2).

make_annual_summary <- function(df, scenario_label) {
  df %>%
    filter(time %% 365 == 0 & time != 0) %>%
    mutate(year_after_burnin = time / 365) %>%   # row label for year (1..5)
    transmute(
      scenario = scenario_label,                  # row label for scenario
      year_after_burnin,
      cum_case_s,
      cum_case_c,
      cum_new_case_known_c,
      prop_known_c = if_else(cum_case_c > 0, cum_new_case_known_c / cum_case_c, NA_real_),
      cum_case_i,
      cum_new_case_known_i,
      prop_known_i = if_else(cum_case_i > 0, cum_new_case_known_i / cum_case_i, NA_real_),
      cum_case_ci,
      cum_dis_total,
      cum_mort_hosp_c,
      cum_mort_hosp_i,
      cum_mort_hosp_total,
      cum_readmit30_total,
      cum_readmit365_total,
      cum_readmit_total,
      cum_admit_com_total
    )
}

annual_summary_i0 <- make_annual_summary(CPE.ODE.solution.i0.abso, "i0_baseline")
annual_summary_i1 <- make_annual_summary(CPE.ODE.solution.i1.abso, "i1")
annual_summary_i2 <- make_annual_summary(CPE.ODE.solution.i2.abso, "i2")

# Combined table (useful for side-by-side checks)
annual_summary_all <- bind_rows(annual_summary_i0, annual_summary_i1, annual_summary_i2) %>%
  arrange(scenario, year_after_burnin)

# Combined table (useful for side-by-side checks)
annual_summary_all_pretty <- annual_summary_all %>%
  mutate(
    across(starts_with("cum_"), ~ round(.x, 2)),
    prop_known_c = round(prop_known_c, 4),
    prop_known_i = round(prop_known_i, 4)
  )

if (interactive()) {
  View(annual_summary_all_pretty)
}

write_csv(annual_summary_all_pretty, "Outputs/CPE_ODE_solution_annual_summary_all_pretty.csv")

# ========================================================== #
# 7) Stabilisation summary (when key states stop changing) ####

# Detect first day where all selected states change less than tolerance
find_stabilisation_day <- function(df, states, tol = 1e-4, consecutive_days = 30) {
  stopifnot(all(states %in% names(df)))
  
  x <- df %>%
    arrange(time) %>%
    select(time, all_of(states))
  
  # absolute day-to-day changes
  d <- x %>%
    mutate(across(all_of(states), ~ abs(.x - lag(.x))))
  
  # true if all states changed by <= tol on that day
  stable_flag <- d %>%
    mutate(stable = if_all(all_of(states), ~ !is.na(.x) & .x <= tol)) %>%
    pull(stable)
  
  # run-length logic: first index where we have 'consecutive_days' TRUEs
  r <- rle(stable_flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  idx <- which(r$values & r$lengths >= consecutive_days)[1]
  
  if (is.na(idx)) {
    return(NA_real_)
  } else {
    # first day in the first qualifying stable run
    first_stable_row <- starts[idx]
    return(x$time[first_stable_row])
  }
}

# Summarise state values at stabilisation day (or NA if never stabilised)
stabilisation_summary <- function(df, scenario_label,
                                  states = c("Sg","Cg","Ig"),
                                  tol = 1e-4,
                                  consecutive_days = 30) {
  
  t_star <- find_stabilisation_day(df, states, tol = tol, consecutive_days = consecutive_days)
  
  if (is.na(t_star)) {
    return(tibble(
      scenario = scenario_label,
      stabilisation_day = NA_real_,
      tolerance = tol,
      consecutive_days = consecutive_days,
      stabilised = FALSE
    ))
  }
  
  row_star <- df %>%
    filter(time == t_star) %>%
    slice(1) %>%
    select(all_of(states))
  
  tibble(
    scenario = scenario_label,
    stabilisation_day = t_star,
    tolerance = tol,
    consecutive_days = consecutive_days,
    stabilised = TRUE
  ) %>%
    bind_cols(row_star)
}

# Build summary for all scenarios
stability_i0 <- stabilisation_summary(CPE.ODE.solution.i0.abso, "i0")
stability_i1 <- stabilisation_summary(CPE.ODE.solution.i1.abso, "i1")
stability_i2 <- stabilisation_summary(CPE.ODE.solution.i2.abso, "i2")

stability_all <- bind_rows(stability_i0, stability_i1, stability_i2)

# Optional rounded version
stability_all_pretty <- stability_all %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))


if (interactive()) {
  View(stability_all_pretty)
}

write_csv(stability_all, "Outputs/CPE_ODE_stabilisation_summary.csv")
write_csv(stability_all_pretty, "Outputs/CPE_ODE_stabilisation_summary_pretty.csv")

# ========================================================== #

# 8) Negative values checks ####

# If any negative values are present, these functions will identify which variables 
# and how many instances of negative values there are, to help with debugging and model validation. 
# The results are saved to CSV files for review.

# If there are negative values for any of the model states, this could indicate an issue with 
# the model equations, parameter values, or numerical solver.

neg_i0 <- negative_value_rows(CPE.ODE.solution.i0.abso)
neg_i1 <- negative_value_rows(CPE.ODE.solution.i1.abso)
neg_i2 <- negative_value_rows(CPE.ODE.solution.i2.abso)

write_csv(neg_i0, "Outputs/negative_values_i0.csv")
write_csv(neg_i1, "Outputs/negative_values_i1.csv")
write_csv(neg_i2, "Outputs/negative_values_i2.csv")

neg_sum_i0 <- negative_value_summary(CPE.ODE.solution.i0.abso)
neg_sum_i1 <- negative_value_summary(CPE.ODE.solution.i1.abso)
neg_sum_i2 <- negative_value_summary(CPE.ODE.solution.i2.abso)

write_csv(neg_sum_i0, "Outputs/negative_values_i0_summary.csv")
write_csv(neg_sum_i1, "Outputs/negative_values_i1_summary.csv")
write_csv(neg_sum_i2, "Outputs/negative_values_i2_summary.csv")

# ========================================================== #

# 9) Plots ####

plot_hospital_compartments(
  CPE.ODE.solution.i0.abso,
  "Number of patients in each hospital compartment over time - Base case",
  "Outputs/hospital_compartments_counts_plot_i0.png"
)

plot_hospital_compartments(
  CPE.ODE.solution.i1.abso,
  "Number of patients in each hospital compartment over time - Intervention i1",
  "Outputs/hospital_compartments_counts_plot_i1.png"
)

plot_hospital_compartments(
  CPE.ODE.solution.i2.abso,
  "Number of patients in each hospital compartment over time - Intervention i2",
  "Outputs/hospital_compartments_counts_plot_i2.png"
)

plot_community_compartments(
  CPE.ODE.solution.i0.abso,
  "Number of individuals in each community compartment over time - Base case",
  "Outputs/community_compartments_counts_plot_i0.png"
)

plot_community_compartments(
  CPE.ODE.solution.i1.abso,
  "Number of individuals in each community compartment over time - Intervention i1",
  "Outputs/community_compartments_counts_plot_i1.png"
)

plot_community_compartments(
  CPE.ODE.solution.i2.abso,
  "Number of individuals in each community compartment over time - Intervention i2",
  "Outputs/community_compartments_counts_plot_i2.png"
)

# ========================================================== #

# 10) Summarise outcomes for CEA ####

i0 <- summarise_outcomes(CPE.ODE.solution.i0.abso)
i1 <- summarise_outcomes(CPE.ODE.solution.i1.abso)
i2 <- summarise_outcomes(CPE.ODE.solution.i2.abso)

# Incremental analysis vs baseline (i0)
incr.costs.i1 <- i1$hos.costs.total - i0$hos.costs.total
infect.avoid.i1 <- i0$new.infect.total - i1$new.infect.total
icer.i1 <- incr.costs.i1 / infect.avoid.i1

incr.costs.i2 <- i2$hos.costs.total - i0$hos.costs.total
infect.avoid.i2 <- i0$new.infect.total - i2$new.infect.total
icer.i2 <- incr.costs.i2 / infect.avoid.i2

# ========================================================== #

# 11) Create CEA results table ####

cea_results <- tibble(
  intervention = c("Universal - culture (i0)", "Universal - PCR (i1)", "Intervention i2"),
  run.days = c(nrow(CPE.ODE.solution.i0.abso), nrow(CPE.ODE.solution.i1.abso), nrow(CPE.ODE.solution.i2.abso)),
  colonisations.imported = comma(round(c(i0$new.colon.import.total, i1$new.colon.import.total, i2$new.colon.import.total), 4), accuracy = 0.0001),
  colonisations.acquired = comma(round(c(i0$new.colon.acquired.total, i1$new.colon.acquired.total, i2$new.colon.acquired.total), 4), accuracy = 0.0001),
  colonisations.total = comma(round(c(i0$new.colon.total, i1$new.colon.total, i2$new.colon.total), 4), accuracy = 0.0001),
  infections.imported = comma(round(c(i0$new.infect.import.total, i1$new.infect.import.total, i2$new.infect.import.total), 4), accuracy = 0.0001),
  infections.acquired = comma(round(c(i0$new.infect.acquired.total, i1$new.infect.acquired.total, i2$new.infect.acquired.total), 4), accuracy = 0.0001),
  infections.total = comma(round(c(i0$new.infect.total, i1$new.infect.total, i2$new.infect.total), 4), accuracy = 0.0001),
  deaths.s = comma(round(c(i0$susc.deaths.total, i1$susc.deaths.total, i2$susc.deaths.total), 4), accuracy = 0.0001),
  deaths.c = comma(round(c(i0$colon.deaths.total, i1$colon.deaths.total, i2$colon.deaths.total), 4), accuracy = 0.0001),
  deaths.i = comma(round(c(i0$infect.deaths.total, i1$infect.deaths.total, i2$infect.deaths.total), 4), accuracy = 0.0001),
  deaths.total = comma(round(c(i0$deaths.total, i1$deaths.total, i2$deaths.total), 4), accuracy = 0.0001),
  discharged = comma(round(c(i0$discharged.total, i1$discharged.total, i2$discharged.total), 4), accuracy = 0.0001),
  readmitted = comma(round(c(i0$readmit.total, i1$readmit.total, i2$readmit.total), 4), accuracy = 0.0001),
  admitted = comma(round(c(i0$admit.total, i1$admit.total, i2$admit.total), 4), accuracy = 0.0001),
  screen.tests = comma(round(c(i0$screen.tests.total, i1$screen.tests.total, i2$screen.tests.total), 4), accuracy = 0.0001),
  clinical.tests = comma(round(c(i0$clinical.tests.total, i1$clinical.tests.total, i2$clinical.tests.total), 4), accuracy = 0.0001),
  hos.costs = comma(round(c(i0$hos.costs.total, i1$hos.costs.total, i2$hos.costs.total), 4), accuracy = 0.0001),
  infect.avoid.vs.i0 = c(
    "",
    comma(round(infect.avoid.i1, 4), accuracy = 0.0001),
    comma(round(infect.avoid.i2, 4), accuracy = 0.0001)
  ),
  incr.costs.vs.i0 = c(
    "",
    comma(round(incr.costs.i1, 4), accuracy = 0.0001),
    comma(round(incr.costs.i2, 4), accuracy = 0.0001)
  ),
  icer.vs.i0 = c(
    "",
    comma(round(icer.i1, 4), accuracy = 0.0001),
    comma(round(icer.i2, 4), accuracy = 0.0001)
  )
)

write_csv(cea_results, "Outputs/cea_results.csv")

# Optional interactive viewing
if (interactive()) {
  View(cea_results)
}

# ========================================================== #

# 12) Clear workspace ####

rm(list = ls())

# Clear console

cat("\014")

# Clear plots

dev.off()

# END OF SCRIPT ####