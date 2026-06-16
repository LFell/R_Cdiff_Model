# Model 1A ODE Functions ####

# 1) Introduction ####

# This script defines the ODE function for Model 1A, which is a compartmental model of hospital
# transmission of a pathogen with the following states:

# The model compartments are:
# S: Susceptible, not colonised or infected
# C: Colonised, carrying the pathogen but not infected
# I: Infected, with symptoms caused by the pathogen, not yet suspected cases or false negatives
# Xtest: Suspected cases, isolated or on a general ward, awaiting test result
# Xtreat: Confirmed cases, isolated and receiving treatment, awaiting recovery & discharge or return
# to a general ward
# D: Dead (cumulative count of deaths, absorbing state)

# The parameters are defined in the accompanying script "Model1_ODE_Parameters.R". Together, these
# scripts are used to run the model and generate outputs in the main quarto document
# (Model1_ODE_Model_Run.qmd).

# ================================================================================================ #

# 2) Validate model settings ####

# Ensure that the model settings are valid before running the ODE function.

validate_run_settings <- function(model_time, time_step) {
  if (!is.numeric(model_time) ||
      length(model_time) != 1 ||
      is.na(model_time) || model_time <= 0) {
    stop("model_time must be a single numeric value > 0.")
  }
  if (!is.numeric(time_step) ||
      length(time_step) != 1 ||
      is.na(time_step) || time_step <= 0) {
    stop("time_step must be a single numeric value > 0.")
  }
}

# ================================================================================================ #

# 3) ODE function ####

Model1A_ODE <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    # a) Ensure non-negative state variables ####
    
    # Avoid negative values for S, C and I as these are used in the force of colonisation
    # calculation and negative values can lead to negative force of colonisation and thus
    # negative flow from S to C.
    
    S <- max(S, 0)
    C <- max(C, 0)
    I <- max(I, 0)
    
    Xtest <- max(Xtest, 0)
    Xtreat <- max(Xtreat, 0)
    
    D <- max(D, 0)
    
    # b) Calculate total in-hospital population (N) ####
    
    # Total in-hospital census - the number of patients in the hospital at a given time (N)
    # is the sum of those in the S, C, I, Xtest and Xtreat states.
    
    N <- S + C + I + Xtest + Xtreat
    
    # c) Ensure non-negative, non-zero N ####
    
    # Ensure the number of patients in the hospital (N) is not negative or zero to avoid
    # issues with division in the force of colonisation calculation.
    
    N <- max(N, 1e-12)
    
    # d) Calculate the force of colonisation (lambda) ####
    
    # The force of colonisation (lambda) is the rate at which susceptible patients become colonised.
    
    # It is calculated as the sum of colonisation due to environmental exposure (beta0) and
    # contributions from colonised, infected, and isolated patients, weighted by their respective
    # transmission rates (betaC, betaI, betaXtest, betaXtreat) and normalized by the total hospital
    # census (N).
    
    lambda <- beta0 + ((betaC * C + betaI * I + betaXtest * Xtest + betaXtreat * Xtreat) / N)
    
    # e) Calculate total exits ####
    
    # Total exits (discharges + deaths) from the hospital are calculated as the sum of the exit
    # rates (discharge + mortality) for each state multiplied by the number of patients in that
    # state.
    
    # It is assumed there are no discharges from Xtest as patients are awaiting test results.
    
    exits_total <- (disS + mortS) * S +
      (disC + mortC) * C +
      (disI + mortI) * I +
      (mortXtest) * Xtest +
      (disXtreat + mortXtreat) * Xtreat
    
    # f) Assign admissions to S, C and I ####
    
    # The total number of admissions to the hospital is determined by the total exits (discharges +
    # deaths) to maintain a fixed census.
    
    # The admissions are split between the S, C and I states, based on the a fixed import prevalence
    # of colonised and infected patients, with the remaining admissions to the S state.
    
    # Assume there are no admissions direct to Xtest or Xtreat as patients are only moved there
    # after admission and subsequent identification as suspected/confirmed cases.
    
    admI <- exits_total * prop_adm_I
    admC <- exits_total * prop_adm_C
    admS <- exits_total - admC - admI
    
    # g) Ensure non-negative admissions to S ####
    
    # Ensure that admissions to S are not negative, which can occur if the import prevalence of
    # colonised and infected patients is very high and/or if discharge and death rates are low,
    # leading to a situation where the calculated admissions to C and I exceed the total exits.
    
    admS <- max(admS, 0)
    
    # h) ODEs ####
    
    # Formulae for the derivatives of each state variable (dS, dC, dI, dXtest, dXtreat, dD),
    # based on the flows into and out of each state.
    
    dS <- admS + mu * C + delta * Xtreat - lambda * S - (disS + mortS) * S
    dC <- admC + lambda * S - (alpha + mu + disC + mortC) * C
    dI <- admI + alpha * C + pi * Xtest - (gamma + disI + mortI) * I
    dXtest <- gamma * I - (theta + pi + mortXtest) * Xtest
    dXtreat <- theta * Xtest - (delta + disXtreat + mortXtreat) * Xtreat
    dD <- mortS * S + mortC * C + mortI * I + mortXtest * Xtest + mortXtreat * Xtreat
    
    # i) Flows (instantaneous rates per day) ####
    
    # These are rates per day, evaluated at each ODE time step (half-day).
    
    # Admissions
    admitted_susceptible_rate <- admS
    admitted_colonised_rate <- admC
    admitted_infected_rate <- admI
    total_admissions_rate <- admS + admC + admI
    
    # Colonisation and decolonisation
    colonised_in_hospital_rate <- lambda * S
    decolonised_in_hospital_rate <- mu * C
    total_colonised_incidence_rate <- colonised_in_hospital_rate + admitted_colonised_rate
    
    # Infection
    infected_in_hospital_rate <- alpha * C
    total_infected_incidence_rate <- infected_in_hospital_rate + admitted_infected_rate
    
    # Case identification and recovery
    suspected_cases_rate <- gamma * I
    confirmed_cases_rate <- theta * Xtest
    recovered_cases_rate <- (delta * Xtreat) + (disXtreat * Xtreat)
    
    # Missed cases
    missed_untested_cases_rate <- gamma0 * I
    missed_false_negative_cases_rate <- pi * Xtest
    total_missed_cases_rate <- missed_untested_cases_rate + missed_false_negative_cases_rate
    
    # Discharges
    discharges_susceptible_rate <- disS * S
    discharges_colonised_rate <- disC * C
    discharges_infected_rate <- disI * I
    # no discharges from Xtest as patients are awaiting test results
    discharges_recovered_rate <- disXtreat * Xtreat
    total_discharges_rate <- discharges_susceptible_rate + discharges_colonised_rate +
      discharges_infected_rate + discharges_recovered_rate
    
    # Deaths
    deaths_susceptible_rate <- mortS * S
    deaths_colonised_rate <- mortC * C
    deaths_infected_rate <- mortI * I
    deaths_suspected_cases_rate <- mortXtest * Xtest
    deaths_confirmed_cases_rate <- mortXtreat * Xtreat
    deaths_post_infection_rate <- deaths_infected_rate + deaths_suspected_cases_rate + deaths_confirmed_cases_rate
    total_deaths_rate <- (mortS * S) + (mortC * C) + (mortI * I) + (mortXtest * Xtest) + (mortXtreat * Xtreat)
    
    total_exits_rate <- total_discharges_rate + total_deaths_rate
    
    # j) Side-room occupancy ####
    
    # Number of patients in Xtest and/or Xtreat at each time step (half-day)
    # Include Xtest and Xtreat if isolate_before_test = TRUE
    # Include only Xtreat when isolate_before_test = FALSE
    
    if (isolate_before_test) {
      side_room_patient_days <- Xtest + Xtreat
    } else {
      side_room_patient_days <- Xtreat
    }
    
    # k) Return derivatives and outputs as a list ####
    list(
      c(dS, dC, dI, dXtest, dXtreat, dD),
      c(
        N = N,
        lambda = lambda,
        
        # Admissions
        admitted_susceptible_rate = admitted_susceptible_rate,
        admitted_colonised_rate = admitted_colonised_rate,
        admitted_infected_rate = admitted_infected_rate,
        total_admissions_rate = total_admissions_rate,
        
        # Colonisation and decolonisation
        colonised_in_hospital_rate = colonised_in_hospital_rate,
        decolonised_in_hospital_rate = decolonised_in_hospital_rate,
        total_colonised_incidence_rate = total_colonised_incidence_rate,
        
        # Infection
        infected_in_hospital_rate = infected_in_hospital_rate,
        total_infected_incidence_rate = total_infected_incidence_rate,
        
        # Case identification and recovery
        suspected_cases_rate = suspected_cases_rate,
        confirmed_cases_rate = confirmed_cases_rate,
        recovered_cases_rate = recovered_cases_rate,
        
        # Missed cases
        missed_untested_cases_rate = missed_untested_cases_rate,
        missed_false_negative_cases_rate = missed_false_negative_cases_rate,
        total_missed_cases_rate = total_missed_cases_rate,
       
        # Discharges
        discharges_susceptible_rate = discharges_susceptible_rate,
        discharges_colonised_rate = discharges_colonised_rate,
        discharges_infected_rate = discharges_infected_rate,
        discharges_recovered_rate = discharges_recovered_rate,
        total_discharges_rate = total_discharges_rate,
        
        # Deaths
        deaths_susceptible_rate = deaths_susceptible_rate,
        deaths_colonised_rate = deaths_colonised_rate,
        deaths_infected_rate = deaths_infected_rate,
        deaths_suspected_cases_rate = deaths_suspected_cases_rate,
        deaths_confirmed_cases_rate = deaths_confirmed_cases_rate,
        deaths_post_infection_rate = deaths_post_infection_rate,
        total_deaths_rate = total_deaths_rate,
        
        total_exits_rate = total_exits_rate,
        
        side_room_patient_days = side_room_patient_days
    )
    )
  })
}

# ================================================================================================ #

# 4) Run the ODE model ####

# For a given set of parameters and return the output as a data frame.

run_scenario <- function(parms) {
  as.data.frame(ode(
    y = mod.init,
    times = mod.t,
    func = Model1A_ODE,
    parms = parms
  ))
}

# ================================================================================================ #

# 5) Post model run utility functions ####

## a) Add per time_step and cumulative flow quantities ####

add_flow_quantities <- function(df, time_step) {
  df %>%
    mutate(
      # If time is 0, assign year 1; otherwise calculate year based on time/365
      # using ceiling() to bin days into years: (0-365 -> Year 1, 365.5-730 -> Year 2)
      year = if_else(time == 0, 1, ceiling(time / 365)),
      
      # Per-step quantities based directly on the ODE function outputs
      admitted_susceptible_step = admitted_susceptible_rate * time_step,
      admitted_colonised_step = admitted_colonised_rate * time_step,
      admitted_infected_step = admitted_infected_rate * time_step,
      total_admissions_step = total_admissions_rate * time_step,
      
      colonised_in_hospital_step = colonised_in_hospital_rate * time_step,
      decolonised_in_hospital_step = decolonised_in_hospital_rate * time_step,
      total_colonised_incidence_step = total_colonised_incidence_rate * time_step,
      
      infected_in_hospital_step = infected_in_hospital_rate * time_step,
      total_infected_incidence_step = total_infected_incidence_rate * time_step,
      
      suspected_cases_step = suspected_cases_rate * time_step,
      confirmed_cases_step = confirmed_cases_rate * time_step,
      recovered_cases_step = recovered_cases_rate * time_step,
      
      missed_untested_cases_step = missed_untested_cases_rate * time_step,
      missed_false_negative_cases_step = missed_false_negative_cases_rate * time_step,
      total_missed_cases_step = total_missed_cases_rate * time_step,
      
      discharges_susceptible_step = discharges_susceptible_rate * time_step,
      discharges_colonised_step = discharges_colonised_rate * time_step,
      discharges_infected_step = discharges_infected_rate * time_step,
      discharges_recovered_step = discharges_recovered_rate * time_step,
      total_discharges_step = total_discharges_rate * time_step,
      
      deaths_susceptible_step = deaths_susceptible_rate * time_step,
      deaths_colonised_step = deaths_colonised_rate * time_step,
      deaths_infected_step = deaths_infected_rate * time_step,
      deaths_suspected_cases_step = deaths_suspected_cases_rate * time_step,
      deaths_confirmed_cases_step = deaths_confirmed_cases_rate * time_step,
      deaths_post_infection_step = deaths_post_infection_rate * time_step,
      total_deaths_step = total_deaths_rate * time_step,
      
      total_exits_step = total_exits_rate * time_step,
      
      side_room_patient_days_step = side_room_patient_days * time_step,
      
      # Cumulative quantities
      cum_admitted_susceptible = cumsum(admitted_susceptible_step),
      cum_admitted_colonised = cumsum(admitted_colonised_step),
      cum_admitted_infected = cumsum(admitted_infected_step),
      cum_total_admissions = cumsum(total_admissions_step),
      
      cum_colonised_in_hospital = cumsum(colonised_in_hospital_step),
      cum_decolonised_in_hospital = cumsum(decolonised_in_hospital_step),
      cum_total_colonised_incidence = cumsum(total_colonised_incidence_step),
      
      cum_infected_in_hospital = cumsum(infected_in_hospital_step),
      cum_total_infected_incidence = cumsum(total_infected_incidence_step),
      
      cum_suspected_cases = cumsum(suspected_cases_step),
      cum_confirmed_cases = cumsum(confirmed_cases_step),
      cum_recovered_cases = cumsum(recovered_cases_step),
      
      cum_missed_untested_cases = cumsum(missed_untested_cases_step),
      cum_missed_false_negative_cases = cumsum(missed_false_negative_cases_step),
      cum_total_missed_cases = cumsum(total_missed_cases_step),
      
      cum_discharges_susceptible = cumsum(discharges_susceptible_step),
      cum_discharges_colonised = cumsum(discharges_colonised_step),
      cum_discharges_infected = cumsum(discharges_infected_step),
      cum_discharges_recovered = cumsum(discharges_recovered_step),
      cum_total_discharges = cumsum(total_discharges_step),
      
      cum_deaths_susceptible = cumsum(deaths_susceptible_step),
      cum_deaths_colonised = cumsum(deaths_colonised_step),
      cum_deaths_infected = cumsum(deaths_infected_step),
      cum_deaths_suspected_cases = cumsum(deaths_suspected_cases_step),
      cum_deaths_confirmed_cases = cumsum(deaths_confirmed_cases_step),
      cum_deaths_post_infection = cumsum(deaths_post_infection_step),
      cum_total_deaths = cumsum(total_deaths_step),
      
      cum_total_exits = cumsum(total_exits_step),
      
      cum_side_room_patient_days = cumsum(side_room_patient_days_step)
    ) %>%
    relocate(year, .before = time)
}

## b) Identify any rows with negative values ####

# Return any rows with negative values for inspection, and print a warning message if any are found.
negative_value_rows <- function(df) {
  # 1. Save the filtered data to a new variable (e.g., neg_df)
  neg_df <- df %>% filter(if_any(where(is.numeric), ~ .x < 0))
  
  # 2. Check the row count of the FILTERED data (neg_df), not the original df
  if (nrow(neg_df) == 0) {
    message("No negative values found in the model output.")
  } else {
    message(
      "Negative values found in the model output - check Outputs/negative_rows.csv for details."
    )
    # 3. Write the FILTERED data to the CSV
    write_csv(neg_df, "Outputs/negative_rows.csv")
  }
  
  # 4. Explicitly return the filtered dataframe so lapply() captures it
  return(neg_df)
}


## c) Find the day when the model states stabilise ####

# Defined as the first day when all specified states change by less than a given tolerance for a
# specified number of consecutive days.
find_stabilisation_day <- function(df,
                                   states,
                                   tol = 1e-4,
                                   consecutive_days = 30) {
  stopifnot(all(states %in% names(df)))
  
  x <- df %>% arrange(time) %>% select(time, all_of(states))
  d <- x %>% mutate(across(all_of(states), ~ abs(.x - lag(.x))))
  stable_flag <- d %>%
    mutate(stable = if_all(all_of(states), ~ !is.na(.x) &
                             .x <= tol)) %>%
    pull(stable)
  
  r <- rle(stable_flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  idx <- which(r$values & r$lengths >= consecutive_days)[1]
  
  # If no stabilisation period is found, print message and return NA
  if (is.na(idx)) {
    message("No stabilisation period of ", consecutive_days, 
            " days found within the model run time frame of ", max(df$time), " days.")
    return(NA_real_)
  }
  
  # If stabilisation is found, save the time, print message, and explicitly RETURN it
  stab_time <- x$time[starts[idx]]
  
  message("Stabilisation achieved in ", consecutive_days, 
          " consecutive steps starting at day ", stab_time, 
          " (tolerance = ", tol, ").")
  
  return(stab_time)
}

## d) Summary table of model states and flows at stabilisation ####

# Along with the stabilisation day and parameters used for the stabilisation assessment.

stabilisation_summary <- function(df,
                                  scenario_label,
                                  state_vars = c("S", "C", "I", "Xtest", "Xtreat", "D"),
                                  flow_vars = c(
                                    "admitted_susceptible_rate",
                                    "admitted_colonised_rate",
                                    "admitted_infected_rate",
                                    "total_admissions_rate",
                                    
                                    "colonised_in_hospital_rate",
                                    "decolonised_in_hospital_rate",
                                    "total_colonised_incidence_rate",
                                    
                                    "infected_in_hospital_rate",
                                    "total_infected_incidence_rate",
                                    
                                    "suspected_cases_rate",
                                    "confirmed_cases_rate",
                                    "recovered_cases_rate",
                                    
                                    "missed_untested_cases_rate",
                                    "missed_false_negative_cases_rate",
                                    "total_missed_cases_rate",
                                    
                                    "discharges_susceptible_rate",
                                    "discharges_colonised_rate",
                                    "discharges_infected_rate",
                                    "discharges_recovered_rate",
                                    "total_discharges_rate",
                                    
                                    "deaths_susceptible_rate",
                                    "deaths_colonised_rate",
                                    "deaths_infected_rate",
                                    "deaths_suspected_cases_rate",
                                    "deaths_confirmed_cases_rate",
                                    "deaths_post_infection_rate",
                                    "total_deaths_rate",
                                    
                                    "total_exits_rate",
                                    
                                    "side_room_patient_days"
                                  ),
                                  extra_vars = c("N", "lambda"),
                                  tol = 1e-4,
                                  consecutive_days = 10) {
  
  # Stabilisation assessed on living/transmission states (i.e. not "Dead")
  stab_states <- c("S", "C", "I", "Xtest", "Xtreat")
  
  stopifnot(all(stab_states %in% names(df)))
  
  t_star <- find_stabilisation_day(df,
                                   stab_states,
                                   tol = tol,
                                   consecutive_days = consecutive_days)
  
  if (is.na(t_star)) {
    return(
      tibble(
        scenario = scenario_label,
        stabilisation_day = NA_real_,
        tolerance = tol,
        consecutive_days = consecutive_days,
        stabilised = FALSE
      )
    )
  }
  
  keep_vars <- unique(c(state_vars, flow_vars, extra_vars))
  keep_vars <- keep_vars[keep_vars %in% names(df)]
  
  row_star <- df %>%
    filter(time == t_star) %>%
    slice(1) %>%
    select(all_of(keep_vars))
  
  tibble(
    scenario = scenario_label,
    stabilisation_day = t_star,
    tolerance = tol,
    consecutive_days = consecutive_days,
    stabilised = TRUE
  ) %>%
    bind_cols(row_star)
}

## e) Create annual summary table ####

make_annual_summary <- function(df, scenario_label) {
  df %>%
    group_by(year) %>%
    summarise(
      annual_admitted_susceptible = sum(admitted_susceptible_step),
      annual_admitted_colonised = sum(admitted_colonised_step),
      annual_admitted_infected = sum(admitted_infected_step),
      annual_total_admissions = sum(total_admissions_step),
      
      annual_colonised_in_hospital = sum(colonised_in_hospital_step),
      annual_decolonised_in_hospital = sum(decolonised_in_hospital_step),
      annual_total_colonised_incidence = sum(total_colonised_incidence_step),
      
      annual_infected_in_hospital = sum(infected_in_hospital_step),
      annual_total_infected_incidence = sum(total_infected_incidence_step),
      
      annual_suspected_cases = sum(suspected_cases_step),
      annual_confirmed_cases = sum(confirmed_cases_step),
      annual_recovered_cases = sum(recovered_cases_step),
      
      annual_missed_untested_cases = sum(missed_untested_cases_step),
      annual_missed_false_negative_cases = sum(missed_false_negative_cases_step),
      annual_total_missed_cases = sum(total_missed_cases_step),
      
      annual_discharges_susceptible = sum(discharges_susceptible_step),
      annual_discharges_colonised = sum(discharges_colonised_step),
      annual_discharges_infected = sum(discharges_infected_step),
      annual_discharges_recovered = sum(discharges_recovered_step),
      annual_total_discharges = sum(total_discharges_step),
      
      annual_deaths_susceptible = sum(deaths_susceptible_step),
      annual_deaths_colonised = sum(deaths_colonised_step),
      annual_deaths_infected = sum(deaths_infected_step),
      annual_deaths_suspected_cases = sum(deaths_suspected_cases_step),
      annual_deaths_confirmed_cases = sum(deaths_confirmed_cases_step),
      annual_deaths_post_infection = sum(deaths_post_infection_step),
      annual_total_deaths = sum(total_deaths_step),
      
      annual_total_exits = sum(total_exits_step),
      
      annual_side_room_patient_days = sum(side_room_patient_days_step),
      .groups = "drop"
    ) %>%
    mutate(scenario = scenario_label) %>%
    select(scenario, year, everything())
}


# End of script ####