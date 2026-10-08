# Model 1A ODE Functions ####

# 1) Introduction ####

# This script defines the ODE function for Model 1A, which is a compartmental model of hospital
# transmission of a pathogen with the following states:

# The model compartments are:
# S: Susceptible, not colonised or infected
# C: Colonised, carrying the pathogen but not infected
# I: Infected, with symptoms caused by the pathogen, not yet suspected cases or false negatives
# Xtest: Suspected cases, isolated or on a general ward, awaiting test result
# FN: False negative cases, still with symptoms but not isolated, awaiting re-testing
# Xtreat: Confirmed cases, isolated or on a general ward and receiving treatment, awaiting recovery
# & discharge or return to a general ward
# Whether Xtest and/or Xtreat are in side rooms is set by the scenario control isolation_policy
# ("none", "confirmed" or "suspected_and_confirmed").
# D: Dead (cumulative count of deaths, absorbing state)

# The parameters are defined in the accompanying script "Model1_ODE_Parameters.R". Together, these
# scripts are used to run the model and generate outputs in the main quarto document
# (Model1_ODE_Model_Run.qmd).

# ================================================================================================ #

# 2) Validate model settings ####

# Ensure that the model settings are valid before running the ODE function.

validate_run_settings <- function(model_time, time_step, burn_in = 0) {
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
  if (time_step >= model_time) {
    stop("time_step must be smaller than model_time.")
  }
  if (!is.numeric(burn_in) ||
      length(burn_in) != 1 ||
      is.na(burn_in) || burn_in < 0 || burn_in >= model_time) {
    stop("burn_in must be a single numeric value >= 0 and smaller than model_time.")
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
    FN <- max(FN, 0)
    Xtreat <- max(Xtreat, 0)
    
    D <- max(D, 0)
    
    # b) Calculate total in-hospital population (N) ####
    
    # Total in-hospital census - the number of patients in the hospital at a given time (N)
    # is the sum of those in the S, C, I, Xtest, FN and Xtreat states.
    
    N <- S + C + I + Xtest + FN + Xtreat
    
    # c) Ensure non-negative, non-zero N ####
    
    # Ensure the number of patients in the hospital (N) is not negative or zero to avoid
    # issues with division in the force of colonisation calculation.
    
    if (N <= 0) {
      stop("N <= 0, total hospital census must be positive.")
    }
    
    # d) Calculate the force of colonisation (lambda) ####
    
    # The force of colonisation (lambda) is the rate at which susceptible patients become colonised.
    
    # It is calculated as the sum of colonisation due to environmental exposure (beta0) and
    # contributions from colonised, infected, and isolated patients, weighted by their respective
    # transmission rates (betaC, betaI, betaXtest, betaFN, betaXtreat) and normalized by the total hospital
    # census (N).
    
    lambda <- beta0 + ((betaC * C + betaI * I + betaXtest * Xtest + betaFN * FN + betaXtreat * Xtreat) / N)
    
    # e) Assign admissions to S, C and I ####
    
    # The total number of admissions to the hospital is set at a fixed rate (adm_rate) per day
    
    # The admissions are split between the S, C and I states, based on the a fixed import prevalence
    # of colonised and infected patients, with the remaining admissions to the S state.
    
    # Assume there are no admissions direct to Xtest or Xtreat as patients are only moved there
    # after admission and subsequent identification as suspected/confirmed cases.
    
    admI <- adm_rate * prop_adm_I
    admC <- adm_rate * prop_adm_C
    admS <- adm_rate - admC - admI
    
    # f) Differential equations ####
    
    # Formulae for updating the numbers of patients in each state variable (dS, dC, dI, dXtest, 
    # dXtreat, dD), based on the flows into and out of each state per time step.
    
    dS      <- admS + mu * C + delta * Xtreat - lambda * S - (disS + mortS) * S
    dC      <- admC + lambda * S - (alpha + mu + disC + mortC) * C
    dI      <- admI + alpha * C - (gamma + disI + mortI) * I
    dXtest  <- gamma * I + sigma * FN - (theta + pi + mortXtest) * Xtest
    dFN     <- pi * Xtest - (sigma + disFN + mortFN) * FN
    dXtreat <- theta * Xtest - (delta + disXtreat + mortXtreat) * Xtreat
    dD      <- mortS * S + mortC * C + mortI * I + mortXtest * Xtest + mortFN * FN + mortXtreat * Xtreat
    
    # g) Bed and side-room bed occupancy ####
    
    # Bed occupancy is calculated from the number of patients in each state at each time step
    
    # Total bed occupancy
    
    # Total bed occupancy at each time step is the sum of patients in S, C, I, Xtest, FN 
    # and Xtreat.
    hospital_bed_occupancy <- N
    
    # Side-room bed days
    
    # Number of patients in side rooms at each time step, depending on isolation_policy:
    #   "none"                    - no side rooms used
    #   "confirmed"               - Xtreat only
    #   "suspected_and_confirmed" - Xtest and Xtreat
    
    side_room_Xtest_occupancy  <- if (isolation_policy == "suspected_and_confirmed") Xtest else 0
    side_room_Xtreat_occupancy <- if (isolation_policy %in% c("confirmed", "suspected_and_confirmed")) Xtreat else 0
    side_room_bed_occupancy    <- side_room_Xtest_occupancy + side_room_Xtreat_occupancy
    
    general_ward_bed_occupancy <- hospital_bed_occupancy - side_room_bed_occupancy
    
    # h) Flows (instantaneous rates per day) ####
    
    # These are rates per day, evaluated at each ODE time step (half-day).
    
    # Admissions
    admitted_susceptible_rate <- admS
    admitted_colonised_rate <- admC
    admitted_infected_rate <- admI
    total_admissions_rate <- admS + admC + admI
    
    # Colonisation and decolonisation
    colonised_in_hospital_rate <- lambda * S
    
    # Colonisation in hospital by source: background acquisition (beta0) and contact with each
    # patient group (patient-to-patient transmission). The sources sum to colonised_in_hospital_rate.
    colonised_background_rate  <- beta0 * S
    colonised_from_C_rate      <- (betaC * C / N) * S
    colonised_from_I_rate      <- (betaI * I / N) * S
    colonised_from_Xtest_rate  <- (betaXtest * Xtest / N) * S
    colonised_from_FN_rate     <- (betaFN * FN / N) * S
    colonised_from_Xtreat_rate <- (betaXtreat * Xtreat / N) * S
    colonised_patient_to_patient_rate <- colonised_from_C_rate + colonised_from_I_rate +
      colonised_from_Xtest_rate + colonised_from_FN_rate + colonised_from_Xtreat_rate
    decolonised_in_hospital_rate <- mu * C
    total_colonised_incidence_rate <- colonised_in_hospital_rate + admitted_colonised_rate
    
    # Infection
    infected_in_hospital_rate <- alpha * C
    total_infected_incidence_rate <- infected_in_hospital_rate + admitted_infected_rate
    
    # Case identification and recovery
    suspected_cases_rate <- gamma * I
    confirmed_cases_rate <- theta * Xtest
    recovered_cases_rate <- (delta * Xtreat) + (disXtreat * Xtreat)
    
    # Testing
    first_test_rate <- gamma * I
    retest_false_negatives_rate <- sigma * FN
    total_test_rate <- first_test_rate + retest_false_negatives_rate
    false_negatives_rate <- pi * Xtest
    
    # Missed cases
    # Missed untested cases are those that leave the I state without being tested, either through 
    # discharge or death
    missed_untested_cases_rate <- (disI + mortI) * I
    # Missed false negative cases are those that leave the FN state without being retested, either 
    # through discharge or death
    missed_false_negative_cases_rate <- (disFN + mortFN) * FN
    total_missed_cases_rate <- missed_untested_cases_rate + missed_false_negative_cases_rate
    
    # Discharges
    discharges_susceptible_rate <- disS * S
    discharges_colonised_rate <- disC * C
    discharges_infected_rate <- disI * I
    # no discharges from Xtest as patients are awaiting test results
    discharges_FN_rate <- disFN * FN
    discharges_recovered_rate <- disXtreat * Xtreat
    discharges_while_infected_rate <-  discharges_infected_rate + discharges_FN_rate
    total_discharges_rate <- discharges_susceptible_rate + discharges_colonised_rate +
                          discharges_infected_rate + discharges_FN_rate + discharges_recovered_rate
    
    # Deaths
    deaths_susceptible_rate <- mortS * S
    deaths_colonised_rate <- mortC * C
    deaths_infected_rate <- mortI * I
    deaths_suspected_cases_rate <- mortXtest * Xtest
    deaths_FN_rate <- mortFN * FN
    deaths_confirmed_cases_rate <- mortXtreat * Xtreat
    deaths_post_infection_rate <- deaths_infected_rate + deaths_suspected_cases_rate +
      deaths_FN_rate + deaths_confirmed_cases_rate
    total_deaths_rate <- deaths_susceptible_rate + deaths_colonised_rate + deaths_infected_rate +
      deaths_suspected_cases_rate + deaths_FN_rate + deaths_confirmed_cases_rate
    
    total_exits_rate <- total_discharges_rate + total_deaths_rate
    
    # i) Return state levels and flows for each time step as a list ####
    
    list(
      c(dS, dC, dI, dXtest, dFN, dXtreat, dD),
      c(
        N = N,
        
        # Force of colonisation
        lambda = lambda,
        
        # Bed occupancy
        hospital_bed_occupancy = hospital_bed_occupancy,
        side_room_bed_occupancy = side_room_bed_occupancy,
        side_room_Xtest_occupancy = side_room_Xtest_occupancy,
        side_room_Xtreat_occupancy = side_room_Xtreat_occupancy,
        general_ward_bed_occupancy = general_ward_bed_occupancy,
        
        # Admissions
        admitted_susceptible_rate = admitted_susceptible_rate,
        admitted_colonised_rate = admitted_colonised_rate,
        admitted_infected_rate = admitted_infected_rate,
        total_admissions_rate = total_admissions_rate,
        
        # Colonisation and decolonisation
        colonised_in_hospital_rate = colonised_in_hospital_rate,
        colonised_background_rate = colonised_background_rate,
        colonised_from_C_rate = colonised_from_C_rate,
        colonised_from_I_rate = colonised_from_I_rate,
        colonised_from_Xtest_rate = colonised_from_Xtest_rate,
        colonised_from_FN_rate = colonised_from_FN_rate,
        colonised_from_Xtreat_rate = colonised_from_Xtreat_rate,
        colonised_patient_to_patient_rate = colonised_patient_to_patient_rate,
        decolonised_in_hospital_rate = decolonised_in_hospital_rate,
        total_colonised_incidence_rate = total_colonised_incidence_rate,
        
        # Infection
        infected_in_hospital_rate = infected_in_hospital_rate,
        total_infected_incidence_rate = total_infected_incidence_rate,
        
        # Case identification and recovery
        suspected_cases_rate = suspected_cases_rate,
        confirmed_cases_rate = confirmed_cases_rate,
        recovered_cases_rate = recovered_cases_rate,
        
        # Testing
        first_test_rate = first_test_rate,
        retest_false_negatives_rate = retest_false_negatives_rate,
        total_test_rate = total_test_rate,
        false_negatives_rate = false_negatives_rate,
        
        # Missed cases
        missed_untested_cases_rate = missed_untested_cases_rate,
        missed_false_negative_cases_rate = missed_false_negative_cases_rate,
        total_missed_cases_rate = total_missed_cases_rate,
        
        # Discharges
        discharges_susceptible_rate = discharges_susceptible_rate,
        discharges_colonised_rate = discharges_colonised_rate,
        discharges_infected_rate = discharges_infected_rate,
        discharges_FN_rate = discharges_FN_rate,
        discharges_recovered_rate = discharges_recovered_rate,
        discharges_while_infected_rate = discharges_while_infected_rate,
        total_discharges_rate = total_discharges_rate,
        
        # Deaths
        deaths_susceptible_rate = deaths_susceptible_rate,
        deaths_colonised_rate = deaths_colonised_rate,
        deaths_infected_rate = deaths_infected_rate,
        deaths_suspected_cases_rate = deaths_suspected_cases_rate,
        deaths_FN_rate = deaths_FN_rate,
        deaths_confirmed_cases_rate = deaths_confirmed_cases_rate,
        deaths_post_infection_rate = deaths_post_infection_rate,
        total_deaths_rate = total_deaths_rate,
        
        # Exits
        total_exits_rate = total_exits_rate
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

## a) Add per time step and cumulative quantities ####

# burn_in (days): the period at the start of the run while the model moves from its initial
# conditions to steady state. It is excluded from all per-step flows, patient-days and cumulative
# totals, which therefore cover the analysis period (burn_in, model_time].
add_quantities <- function(df, time_step, burn_in = 0) {
  df %>%
    mutate(
      # Year of the analysis period, using ceiling() to bin days after the burn-in into years:
      # (burn_in, burn_in + 365] -> Year 1, and so on. Rows in the burn-in are year 0.
      year = if_else(time <= burn_in, 0, ceiling((time - burn_in) / 365)),

      # Width of the interval each row represents. Each row value is multiplied by the step width
      # (a right-endpoint sum: the row at time t covers the interval (t - time_step, t]). Rows up to
      # and including the end of the burn-in (always including time 0) have width 0, so cumulative
      # totals and patient-days cover exactly the analysis period.
      step_width = if_else(time <= burn_in, 0, time_step),
      
      # Patient-days per state
      S_days     = S * step_width,
      C_days     = C * step_width,
      I_days     = I * step_width,
      Xtest_days = Xtest * step_width,
      FN_days    = FN * step_width,
      Xtreat_days= Xtreat * step_width,
      general_ward_bed_days = general_ward_bed_occupancy * step_width,
      side_room_bed_days    = side_room_bed_occupancy * step_width,
      side_room_Xtest_bed_days  = side_room_Xtest_occupancy * step_width,
      side_room_Xtreat_bed_days = side_room_Xtreat_occupancy * step_width,
      hospital_bed_days     = hospital_bed_occupancy * step_width,

      # Per-step flows
      
      # Admissions
      admitted_susceptible_step = admitted_susceptible_rate * step_width,
      admitted_colonised_step = admitted_colonised_rate * step_width,
      admitted_infected_step = admitted_infected_rate * step_width,
      total_admissions_step = total_admissions_rate * step_width,
      
      # Colonisation and decolonisation
      colonised_in_hospital_step = colonised_in_hospital_rate * step_width,
      colonised_background_step = colonised_background_rate * step_width,
      colonised_from_C_step = colonised_from_C_rate * step_width,
      colonised_from_I_step = colonised_from_I_rate * step_width,
      colonised_from_Xtest_step = colonised_from_Xtest_rate * step_width,
      colonised_from_FN_step = colonised_from_FN_rate * step_width,
      colonised_from_Xtreat_step = colonised_from_Xtreat_rate * step_width,
      colonised_patient_to_patient_step = colonised_patient_to_patient_rate * step_width,
      decolonised_in_hospital_step = decolonised_in_hospital_rate * step_width,
      total_colonised_incidence_step = total_colonised_incidence_rate * step_width,
      
      # Infection
      infected_in_hospital_step = infected_in_hospital_rate * step_width,
      total_infected_incidence_step = total_infected_incidence_rate * step_width,
      
      # Case identification and recovery
      suspected_cases_step = suspected_cases_rate * step_width,
      confirmed_cases_step = confirmed_cases_rate * step_width,
      recovered_cases_step = recovered_cases_rate * step_width,
      
      # Testing
      first_test_step = first_test_rate * step_width,
      retest_false_negatives_step = retest_false_negatives_rate * step_width,
      total_test_step = total_test_rate * step_width,
      false_negatives_step = false_negatives_rate * step_width,
      
      # Missed cases
      missed_untested_cases_step = missed_untested_cases_rate * step_width,
      missed_false_negative_cases_step = missed_false_negative_cases_rate * step_width,
      total_missed_cases_step = total_missed_cases_rate * step_width,
      
      # Discharges
      discharges_susceptible_step = discharges_susceptible_rate * step_width,
      discharges_colonised_step = discharges_colonised_rate * step_width,
      discharges_infected_step = discharges_infected_rate * step_width,
      discharges_FN_step = discharges_FN_rate * step_width,
      discharges_recovered_step = discharges_recovered_rate * step_width,
      discharges_while_infected_step = discharges_while_infected_rate * step_width,
      total_discharges_step = total_discharges_rate * step_width,
      
      # Deaths
      deaths_susceptible_step = deaths_susceptible_rate * step_width,
      deaths_colonised_step = deaths_colonised_rate * step_width,
      deaths_infected_step = deaths_infected_rate * step_width,
      deaths_suspected_cases_step = deaths_suspected_cases_rate * step_width,
      deaths_FN_step = deaths_FN_rate * step_width,
      deaths_confirmed_cases_step = deaths_confirmed_cases_rate * step_width,
      deaths_post_infection_step = deaths_post_infection_rate * step_width,
      total_deaths_step = total_deaths_rate * step_width,
      
      # Exits
      total_exits_step = total_exits_rate * step_width,
      
      # Cumulative flows
      
      # Admissions
      cum_admitted_susceptible = cumsum(admitted_susceptible_step),
      cum_admitted_colonised = cumsum(admitted_colonised_step),
      cum_admitted_infected = cumsum(admitted_infected_step),
      cum_total_admissions = cumsum(total_admissions_step),
      
      # Colonisation and decolonisation
      cum_colonised_in_hospital = cumsum(colonised_in_hospital_step),
      cum_colonised_background = cumsum(colonised_background_step),
      cum_colonised_from_C = cumsum(colonised_from_C_step),
      cum_colonised_from_I = cumsum(colonised_from_I_step),
      cum_colonised_from_Xtest = cumsum(colonised_from_Xtest_step),
      cum_colonised_from_FN = cumsum(colonised_from_FN_step),
      cum_colonised_from_Xtreat = cumsum(colonised_from_Xtreat_step),
      cum_colonised_patient_to_patient = cumsum(colonised_patient_to_patient_step),
      cum_decolonised_in_hospital = cumsum(decolonised_in_hospital_step),
      cum_total_colonised_incidence = cumsum(total_colonised_incidence_step),
      
      # Infections
      cum_infected_in_hospital = cumsum(infected_in_hospital_step),
      cum_total_infected_incidence = cumsum(total_infected_incidence_step),
      
      # Case identification and recovery
      cum_suspected_cases = cumsum(suspected_cases_step),
      cum_confirmed_cases = cumsum(confirmed_cases_step),
      cum_recovered_cases = cumsum(recovered_cases_step),
      
      # Testing
      cum_first_test = cumsum(first_test_step),
      cum_retest_false_negatives = cumsum(retest_false_negatives_step),
      cum_total_test = cumsum(total_test_step),
      cum_false_negatives = cumsum(false_negatives_rate * step_width),
      
      # Missing cases
      cum_missed_untested_cases = cumsum(missed_untested_cases_step),
      cum_missed_false_negative_cases = cumsum(missed_false_negative_cases_step),
      cum_total_missed_cases = cumsum(total_missed_cases_step),
      
      # Discharges   
      cum_discharges_susceptible = cumsum(discharges_susceptible_step),
      cum_discharges_colonised = cumsum(discharges_colonised_step),
      cum_discharges_infected = cumsum(discharges_infected_step),
      cum_discharges_FN = cumsum(discharges_FN_step),
      cum_discharges_recovered = cumsum(discharges_recovered_step),
      cum_discharges_while_infected = cumsum(discharges_while_infected_step),
      cum_total_discharges = cumsum(total_discharges_step),
      
      # Deaths
      cum_deaths_susceptible = cumsum(deaths_susceptible_step),
      cum_deaths_colonised = cumsum(deaths_colonised_step),
      cum_deaths_infected = cumsum(deaths_infected_step),
      cum_deaths_suspected_cases = cumsum(deaths_suspected_cases_step),
      cum_deaths_FN = cumsum(deaths_FN_step),
      cum_deaths_confirmed_cases = cumsum(deaths_confirmed_cases_step),
      cum_deaths_post_infection = cumsum(deaths_post_infection_step),
      cum_total_deaths = cumsum(total_deaths_step),
      
      # Exits
      cum_total_exits = cumsum(total_exits_step)
    ) %>%
    select(-step_width) %>%
    relocate(year, .before = time)
}

## b) Calculate total patient-days per scenario ####

calculate_total_patient_days <- function(df) {
  summarise(df,
            S_patient_days = sum(S_days, na.rm = TRUE),
            C_patient_days = sum(C_days, na.rm = TRUE),
            I_patient_days = sum(I_days, na.rm = TRUE),
            Xtest_patient_days = sum(Xtest_days, na.rm = TRUE),
            FN_patient_days    = sum(FN_days, na.rm = TRUE),
            Xtreat_patient_days = sum(Xtreat_days, na.rm = TRUE),
            general_ward_bed_days = sum(general_ward_bed_days, na.rm = TRUE),
            side_room_bed_days    = sum(side_room_bed_days, na.rm = TRUE),
            side_room_Xtest_bed_days  = sum(side_room_Xtest_bed_days, na.rm = TRUE),
            side_room_Xtreat_bed_days = sum(side_room_Xtreat_bed_days, na.rm = TRUE),
            hospital_bed_days     = sum(hospital_bed_days, na.rm = TRUE)
  )
}

## c) Identify any rows with negative values ####

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
  }
  
  # 3. Explicitly return the filtered dataframe so lapply() captures it
  return(neg_df)
}

## d) Find the day when the model states stabilise ####

# Defined as the first day when all specified states change by less than a given tolerance for a
# specified number of consecutive time steps.
#
# Note: despite its name, consecutive_days counts consecutive rows of the model output (time steps),
# not days. With time_step = 0.5, consecutive_days = 10 means 10 half-day steps = 5 days. Likewise
# tol is the maximum absolute change per time step, not per day.
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
    message(
      "No stabilisation period of ",
      consecutive_days,
      " days found within the model run time frame of ",
      max(df$time),
      " days."
    )
    return(NA_real_)
  }
  
  # If stabilisation is found, save the time, print message, and explicitly RETURN it
  stab_time <- x$time[starts[idx]]
  
  message(
    "Stabilisation achieved in ",
    consecutive_days,
    " consecutive steps starting at day ",
    stab_time,
    " (tolerance = ",
    tol,
    ")."
  )
  
  return(stab_time)
}

## e) Summary table of model states and flows at stabilisation ####

# Along with the stabilisation day and parameters used for the stabilisation assessment.
# consecutive_days counts time steps, not days (see find_stabilisation_day()): the default of 10
# is 5 days at the half-day time step. The consecutive_days column in the output is this step count.

stabilisation_summary <- function(df,
                                  scenario_label,
                                  state_vars = c(
                                    "S",
                                    "C",
                                    "I",
                                    "Xtest",
                                    "FN",
                                    "Xtreat",
                                    "D",
                                    "hospital_bed_occupancy",
                                    "general_ward_bed_occupancy",
                                    "side_room_bed_occupancy"
                                  ),
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
                                    
                                    "first_test_rate",
                                    "retest_false_negatives_rate",
                                    "total_test_rate",
                                    "false_negatives_rate",
                                    
                                    "missed_untested_cases_rate",
                                    "missed_false_negative_cases_rate",
                                    "total_missed_cases_rate",
                                    
                                    "discharges_susceptible_rate",
                                    "discharges_colonised_rate",
                                    "discharges_infected_rate",
                                    "discharges_recovered_rate",
                                    "discharges_while_infected_rate",
                                    "total_discharges_rate",
                                    
                                    "deaths_susceptible_rate",
                                    "deaths_colonised_rate",
                                    "deaths_infected_rate",
                                    "deaths_suspected_cases_rate",
                                    "deaths_confirmed_cases_rate",
                                    "deaths_post_infection_rate",
                                    "total_deaths_rate",
                                    
                                    "total_exits_rate"
                                  ),
                                  extra_vars = c("N", "lambda"),
                                  tol = 1e-4,
                                  consecutive_days = 10) {
  # Stabilisation assessed on living/transmission states (i.e. not "Dead")
  stab_states <- c("S", "C", "I", "Xtest", "FN", "Xtreat")
  
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

## f) Create annual summary table ####

# Years of the analysis period only (year >= 1); the burn-in (year 0) is excluded.
make_annual_summary <- function(df, scenario_label) {
  df %>%
    filter(year >= 1) %>%
    group_by(year) %>%
    summarise(
      # Bed days
      annual_S_patient_days   = sum(S_days, na.rm = TRUE),
      annual_C_patient_days   = sum(C_days, na.rm = TRUE),
      annual_I_patient_days   = sum(I_days, na.rm = TRUE),
      annual_Xtest_patient_days = sum(Xtest_days, na.rm = TRUE),
      annual_FN_patient_days    = sum(FN_days, na.rm = TRUE),
      annual_Xtreat_patient_days = sum(Xtreat_days, na.rm = TRUE),
      annual_hospital_bed_days = sum(hospital_bed_days, na.rm = TRUE),
      annual_general_ward_bed_days = sum(general_ward_bed_days, na.rm = TRUE),
      annual_side_room_bed_days = sum(side_room_bed_days, na.rm = TRUE),
      annual_side_room_Xtest_bed_days  = sum(side_room_Xtest_bed_days, na.rm = TRUE),
      annual_side_room_Xtreat_bed_days = sum(side_room_Xtreat_bed_days, na.rm = TRUE),
      
      # Admissions
      annual_admitted_susceptible = sum(admitted_susceptible_step, na.rm = TRUE),
      annual_admitted_colonised = sum(admitted_colonised_step, na.rm = TRUE),
      annual_admitted_infected = sum(admitted_infected_step, na.rm = TRUE),
      annual_total_admissions = sum(total_admissions_step, na.rm = TRUE),
      
      # Colonisation and decolonisation
      annual_colonised_in_hospital = sum(colonised_in_hospital_step, na.rm = TRUE),
      annual_colonised_background = sum(colonised_background_step, na.rm = TRUE),
      annual_colonised_patient_to_patient = sum(colonised_patient_to_patient_step, na.rm = TRUE),
      annual_decolonised_in_hospital = sum(decolonised_in_hospital_step, na.rm = TRUE),
      annual_total_colonised_incidence = sum(total_colonised_incidence_step, na.rm = TRUE),
      
      # Infection
      annual_infected_in_hospital = sum(infected_in_hospital_step, na.rm = TRUE),
      annual_total_infected_incidence = sum(total_infected_incidence_step, na.rm = TRUE),
      
      # Case identification and recovery
      annual_suspected_cases = sum(suspected_cases_step, na.rm = TRUE),
      annual_confirmed_cases = sum(confirmed_cases_step, na.rm = TRUE),
      annual_recovered_cases = sum(recovered_cases_step, na.rm = TRUE),
      
      # Testing
      annual_first_test = sum(first_test_step, na.rm = TRUE),
      annual_retest_false_negatives = sum(retest_false_negatives_step, na.rm = TRUE),
      annual_total_test = sum(total_test_step, na.rm = TRUE),
      annual_false_negatives = sum(false_negatives_step, na.rm = TRUE),
      
      # Missed cases
      annual_missed_untested_cases = sum(missed_untested_cases_step, na.rm = TRUE),
      annual_missed_false_negative_cases = sum(missed_false_negative_cases_step, na.rm = TRUE),
      annual_total_missed_cases = sum(total_missed_cases_step, na.rm = TRUE),
      
      # Discharges
      annual_discharges_susceptible = sum(discharges_susceptible_step, na.rm = TRUE),
      annual_discharges_colonised = sum(discharges_colonised_step, na.rm = TRUE),
      annual_discharges_infected = sum(discharges_infected_step, na.rm = TRUE),
      annual_discharges_FN = sum(discharges_FN_step, na.rm = TRUE),
      annual_discharges_recovered = sum(discharges_recovered_step, na.rm = TRUE),
      annual_discharges_while_infected = sum(discharges_while_infected_step, na.rm = TRUE),
      annual_total_discharges = sum(total_discharges_step, na.rm = TRUE),
      
      # Deaths
      annual_deaths_susceptible = sum(deaths_susceptible_step, na.rm = TRUE),
      annual_deaths_colonised = sum(deaths_colonised_step, na.rm = TRUE),
      annual_deaths_infected = sum(deaths_infected_step, na.rm = TRUE),
      annual_deaths_suspected_cases = sum(deaths_suspected_cases_step, na.rm = TRUE),
      annual_deaths_FN = sum(deaths_FN_step, na.rm = TRUE),
      annual_deaths_confirmed_cases = sum(deaths_confirmed_cases_step, na.rm = TRUE),
      annual_deaths_post_infection = sum(deaths_post_infection_step, na.rm = TRUE),
      annual_total_deaths = sum(total_deaths_step, na.rm = TRUE),
      
      # Exits
      annual_total_exits = sum(total_exits_step, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(scenario = scenario_label) %>%
    select(scenario, year, everything())
}

# ================================================================================================ #

# 6) Output summary functions ####

## a) Final cumulative outcomes table ####

# Extracts final state values and all cumulative outcome columns from a named
# list of post-processed scenario data frames (output of add_flow_quantities()).
# Returns a one-row-per-scenario tibble.

make_final_outcomes <- function(scenario_post) {
  bind_rows(lapply(names(scenario_post), function(nm) {
    last <- scenario_post[[nm]] %>% filter(time == max(time)) %>% slice(1)
    
    # Post-hoc split of cumulative confirmed cases between first test and retest
    # Based on the proportion of the inflow to Xtest from I and FN
    prop_first_test <- if (last$cum_total_test > 0) {
      last$cum_first_test / last$cum_total_test
    } else {
      0
    }
  
    # Cumulative confirmed cases via first test vs retest
    # Based on the proportion of the inflow to Xtest
    cum_confirmed_via_first_test <- last$cum_confirmed_cases * prop_first_test
    cum_confirmed_via_retest <- last$cum_confirmed_cases * (1 - prop_first_test)
   
    # Hospital-attributable infections
    # Proportion of C pool that is hospital-acquired, using gross inflow ratio.
    # Decolonisation (mu) applies equally to both hospital-acquired and imported C,
    # so it cancels in the ratio and does not need to be netted out.
    prop_C_hospital_acquired <- if (last$cum_total_colonised_incidence > 0) {
      last$cum_colonised_in_hospital / last$cum_total_colonised_incidence
    } else {
      0
    }
    cum_hospital_attributable_I  <- last$cum_infected_in_hospital * prop_C_hospital_acquired
    prop_I_hospital_attributable <- if (last$cum_total_infected_incidence > 0) {
      cum_hospital_attributable_I / last$cum_total_infected_incidence
    } else {
      0
    }
    
    # Split hospital-acquired colonisations and CDI into background acquisition (beta0) and
    # patient-to-patient transmission. Hospital-attributable CDI is apportioned by the share of
    # in-hospital colonisations from each route (same gross-inflow logic as above).
    share_patient_to_patient <- if (last$cum_colonised_in_hospital > 0) {
      last$cum_colonised_patient_to_patient / last$cum_colonised_in_hospital
    } else {
      0
    }
    safe_prop <- function(num, den) if (den > 0) num / den else 0
    cum_hospital_attributable_I_patient_to_patient <- cum_hospital_attributable_I * share_patient_to_patient
    cum_hospital_attributable_I_background         <- cum_hospital_attributable_I * (1 - share_patient_to_patient)
   
    # Add these to the final outcomes tibble 
    tibble(
      scenario     = nm,
      
      # Final state occupancy
      final_S      = last$S,
      final_C      = last$C,
      final_I      = last$I,
      final_Xtest  = last$Xtest,
      final_FN     = last$FN,
      final_Xtreat = last$Xtreat,
      final_D      = last$D,
      final_lambda = last$lambda,
      
      # Final bed occupancy
      final_hospital_bed_occupancy = last$hospital_bed_occupancy,
      final_general_ward_bed_occupancy = last$general_ward_bed_occupancy,
      final_side_room_bed_occupancy = last$side_room_bed_occupancy,
      
      # Admissions
      cum_admitted_susceptible = last$cum_admitted_susceptible,
      cum_admitted_colonised   = last$cum_admitted_colonised,
      cum_admitted_infected    = last$cum_admitted_infected,
      cum_total_admissions     = last$cum_total_admissions,
      
      # Colonisation
      cum_colonised_in_hospital     = last$cum_colonised_in_hospital,
      cum_decolonised_in_hospital   = last$cum_decolonised_in_hospital,
      cum_total_colonised_incidence = last$cum_total_colonised_incidence,
      
      # Colonisation in hospital by source
      cum_colonised_background = last$cum_colonised_background,
      cum_colonised_from_C = last$cum_colonised_from_C,
      cum_colonised_from_I = last$cum_colonised_from_I,
      cum_colonised_from_Xtest = last$cum_colonised_from_Xtest,
      cum_colonised_from_FN = last$cum_colonised_from_FN,
      cum_colonised_from_Xtreat = last$cum_colonised_from_Xtreat,
      cum_colonised_patient_to_patient = last$cum_colonised_patient_to_patient,
      prop_C_background         = safe_prop(last$cum_colonised_background, last$cum_total_colonised_incidence),
      prop_C_patient_to_patient = safe_prop(last$cum_colonised_patient_to_patient, last$cum_total_colonised_incidence),
      prop_p2p_from_C      = safe_prop(last$cum_colonised_from_C,      last$cum_colonised_patient_to_patient),
      prop_p2p_from_I      = safe_prop(last$cum_colonised_from_I,      last$cum_colonised_patient_to_patient),
      prop_p2p_from_Xtest  = safe_prop(last$cum_colonised_from_Xtest,  last$cum_colonised_patient_to_patient),
      prop_p2p_from_FN     = safe_prop(last$cum_colonised_from_FN,     last$cum_colonised_patient_to_patient),
      prop_p2p_from_Xtreat = safe_prop(last$cum_colonised_from_Xtreat, last$cum_colonised_patient_to_patient),
      
      # Infection
      cum_infected_in_hospital     = last$cum_infected_in_hospital,
      cum_total_infected_incidence = last$cum_total_infected_incidence,
      
      # Case identification and recovery
      cum_suspected_cases = last$cum_suspected_cases,
      cum_confirmed_cases = last$cum_confirmed_cases,
      cum_recovered_cases = last$cum_recovered_cases,
      
      # Testing
      cum_first_test             = last$cum_first_test,
      cum_retest_false_negatives = last$cum_retest_false_negatives,
      cum_total_test             = last$cum_total_test,
      cum_false_negatives        = last$cum_false_negatives,
      
      # Missed cases
      cum_missed_untested_cases       = last$cum_missed_untested_cases,
      cum_missed_false_negative_cases = last$cum_missed_false_negative_cases,
      cum_total_missed_cases          = last$cum_total_missed_cases,
      
      # Confirmation via first test vs retest
      cum_confirmed_via_first_test = cum_confirmed_via_first_test,
      cum_confirmed_via_retest = cum_confirmed_via_retest,
     
      # Hospital-attributable infections
      prop_C_hospital_acquired     = prop_C_hospital_acquired,
      prop_C_hospital_attributable = prop_C_hospital_acquired,
      cum_hospital_attributable_I  = cum_hospital_attributable_I,
      prop_I_hospital_attributable = prop_I_hospital_attributable,
      cum_hospital_attributable_I_patient_to_patient = cum_hospital_attributable_I_patient_to_patient,
      cum_hospital_attributable_I_background         = cum_hospital_attributable_I_background,
      prop_I_hospital_attributable_patient_to_patient =
        safe_prop(cum_hospital_attributable_I_patient_to_patient, last$cum_total_infected_incidence),
      prop_I_hospital_attributable_background =
        safe_prop(cum_hospital_attributable_I_background, last$cum_total_infected_incidence),
     
      # Discharges
      cum_discharges_susceptible    = last$cum_discharges_susceptible,
      cum_discharges_colonised      = last$cum_discharges_colonised,
      cum_discharges_infected       = last$cum_discharges_infected,
      cum_discharges_FN             = last$cum_discharges_FN,
      cum_discharges_recovered      = last$cum_discharges_recovered,
      cum_discharges_while_infected = last$cum_discharges_while_infected,
      cum_total_discharges          = last$cum_total_discharges,
      
      # Deaths
      cum_deaths_susceptible     = last$cum_deaths_susceptible,
      cum_deaths_colonised       = last$cum_deaths_colonised,
      cum_deaths_infected        = last$cum_deaths_infected,
      cum_deaths_suspected_cases = last$cum_deaths_suspected_cases,
      cum_deaths_FN              = last$cum_deaths_FN,
      cum_deaths_confirmed_cases = last$cum_deaths_confirmed_cases,
      cum_deaths_post_infection  = last$cum_deaths_post_infection,
      cum_total_deaths           = last$cum_total_deaths,
      
      # Exits
      cum_total_exits             = last$cum_total_exits
    )
  }))
}

## b) Research question outcomes table ####

# Takes final_outcomes (from make_final_outcomes()) and patient_days_summary
# (from calculate_total_patient_days()) and returns a single wide table of all
# research question outcomes, one row per outcome and one column per scenario.
#
# Footnotes in the Outcome column:
#   [1] prop_confirmed + prop_missed do not sum to 1; the gap reflects deaths
#       in the Xtest state (tested but unconfirmed) and unresolved cases still
#       in hospital at the end of the model run.
#   [2] Discharge minus admission prevalence, in percentage points (pp).
#       A positive value indicates the hospital is a net source.
#   [3] Missed cases by origin are apportioned by each origin's share of CDI cases: all cases enter
#       the same I compartment and follow the same rates, so the model cannot track the origin of
#       each missed case, and the proportion missed is the same for every origin.
#
# per_year = TRUE labels counts "(n per year)": use it when the count and patient-day columns have
# been divided by the length of the analysis period. Rates per 10,000 patient bed-days and
# proportions are the same either way.

make_research_outcomes_table <- function(final_outcomes, patient_days_summary, per_year = FALSE) {
  
  safe_div <- function(num, den) ifelse(den > 0, num / den, NA_real_)
  
  fo <- final_outcomes |>
    left_join(
      patient_days_summary |> select(scenario, hospital_bed_days, side_room_bed_days,
                                    side_room_Xtest_bed_days, side_room_Xtreat_bed_days),
      by = "scenario"
    ) |>
    mutate(
      # Case ascertainment proportions
      prop_confirmed            = safe_div(cum_confirmed_cases,        cum_total_infected_incidence),
      prop_missed               = safe_div(cum_total_missed_cases,     cum_total_infected_incidence),
      prop_confirmed_via_retest = safe_div(cum_confirmed_via_retest,   cum_total_infected_incidence),

      # Origin of CDI cases: admitted with CDI or hospital-onset (all in-hospital C -> I)
      prop_I_hospital_onset     = safe_div(cum_infected_in_hospital,   cum_total_infected_incidence),

      # Rates per 10,000 patient bed-days
      rate_cdi                  = 1e4 * safe_div(cum_total_infected_incidence, hospital_bed_days),
      rate_ho                   = 1e4 * safe_div(cum_infected_in_hospital,     hospital_bed_days),
      rate_ha                   = 1e4 * safe_div(cum_hospital_attributable_I,  hospital_bed_days),
      rate_tests                = 1e4 * safe_div(cum_total_test,               hospital_bed_days),

      # Missed cases by origin, apportioned by each origin's share of CDI cases [3]
      missed_admitted           = cum_total_missed_cases * safe_div(cum_admitted_infected,       cum_total_infected_incidence),
      missed_ho                 = cum_total_missed_cases * safe_div(cum_infected_in_hospital,    cum_total_infected_incidence),
      missed_ha                 = cum_total_missed_cases * safe_div(cum_hospital_attributable_I, cum_total_infected_incidence),

      # Discharge and admission prevalence
      prop_admitted_C           = safe_div(cum_admitted_colonised,     cum_total_admissions),
      prop_admitted_I           = safe_div(cum_admitted_infected,      cum_total_admissions),
      prop_discharged_C         = safe_div(cum_discharges_colonised,   cum_total_discharges),
      prop_discharged_CDI       = safe_div(cum_discharges_while_infected, cum_total_discharges),
      diff_prev_C               = prop_discharged_C   - prop_admitted_C,
      diff_prev_CDI             = prop_discharged_CDI - prop_admitted_I,
      
      # Resource use
      avg_LOS                   = safe_div(hospital_bed_days,          cum_total_admissions),
      prop_side_room_bed_days   = safe_div(side_room_bed_days,         hospital_bed_days),
      prop_side_room_Xtest_bed_days  = safe_div(side_room_Xtest_bed_days,  hospital_bed_days),
      prop_side_room_Xtreat_bed_days = safe_div(side_room_Xtreat_bed_days, hospital_bed_days),
      
      # Mortality
      prop_deaths_post_infection = safe_div(cum_deaths_post_infection, cum_total_deaths),
      cum_deaths_missed_cases    = cum_deaths_infected + cum_deaths_FN,
      prop_deaths_missed_cases   = safe_div(cum_deaths_infected + cum_deaths_FN, cum_total_deaths)
    )
  
  outcomes_spec <- tribble(
    ~group,                   ~label,                                                           ~variable,                     ~fmt,
    
    "Case ascertainment",     "(i) Total incident CDI cases (n)",                               "cum_total_infected_incidence", "n",
    "Case ascertainment",     "(i) - admitted with CDI (n)",                                    "cum_admitted_infected",        "n",
    "Case ascertainment",     "(i) - hospital-onset (n)",                                       "cum_infected_in_hospital",     "n",
    "Case ascertainment",     "(i) - hospital-onset (% of all CDI)",                            "prop_I_hospital_onset",        "pct",
    "Case ascertainment",     "(i) CDI cases per 10,000 patient bed-days",                      "rate_cdi",                     "dec",
    "Case ascertainment",     "(i) Confirmed (n)",                                              "cum_confirmed_cases",          "n",
    "Case ascertainment",     "(i) Confirmed (%) [1]",                                          "prop_confirmed",               "pct",
    "Case ascertainment",     "(ii) Missed (n)",                                                "cum_total_missed_cases",       "n",
    "Case ascertainment",     "(ii) Missed (%) [1]",                                            "prop_missed",                  "pct",
    "Case ascertainment",     "(ii) - missed, admitted with CDI (n) [3]",                      "missed_admitted",              "n",
    "Case ascertainment",     "(ii) - missed, hospital-onset (n) [3]",                         "missed_ho",                    "n",
    "Case ascertainment",     "(ii) - missed, hospital-acquired (n) [3]",                      "missed_ha",                    "n",
    "Case ascertainment",     "(iii) Confirmed via retest (n)",                                 "cum_confirmed_via_retest",     "n",
    "Case ascertainment",     "(iii) Confirmed via retest (%) [1]",                             "prop_confirmed_via_retest",    "pct",
    
    "Transmission",           "(i) Hospital-acquired CDI (n)",                                  "cum_hospital_attributable_I",  "n",
    "Transmission",           "(i) Hospital-acquired CDI (% of all CDI)",                       "prop_I_hospital_attributable", "pct",
    "Transmission",           "(i) Hospital-onset CDI per 10,000 patient bed-days",             "rate_ho",                      "dec",
    "Transmission",           "(i) Hospital-acquired CDI per 10,000 patient bed-days",          "rate_ha",                      "dec",
    "Transmission",           "(i) - via patient-to-patient transmission (n)",                 "cum_hospital_attributable_I_patient_to_patient", "n",
    "Transmission",           "(i) - via patient-to-patient transmission (% of all CDI)",      "prop_I_hospital_attributable_patient_to_patient", "pct",
    "Transmission",           "(i) - via background acquisition (n)",                          "cum_hospital_attributable_I_background", "n",
    "Transmission",           "(i) - via background acquisition (% of all CDI)",               "prop_I_hospital_attributable_background", "pct",
    "Transmission",           "(ii) Hospital-acquired colonisations (n)",                       "cum_colonised_in_hospital",    "n",
    "Transmission",           "(ii) Hospital-acquired colonisations (% of all colonisations)",  "prop_C_hospital_attributable", "pct",
    "Transmission",           "(ii) - via patient-to-patient transmission (n)",                "cum_colonised_patient_to_patient", "n",
    "Transmission",           "(ii) - via patient-to-patient transmission (% of all colonisations)", "prop_C_patient_to_patient", "pct",
    "Transmission",           "(ii) - via background acquisition (n)",                         "cum_colonised_background",     "n",
    "Transmission",           "(ii) - via background acquisition (% of all colonisations)",    "prop_C_background",            "pct",
    "Transmission",           "(ii) Patient-to-patient colonisations from C (n)",              "cum_colonised_from_C",         "n",
    "Transmission",           "(ii) Patient-to-patient colonisations from C (%)",              "prop_p2p_from_C",              "pct",
    "Transmission",           "(ii) Patient-to-patient colonisations from I (n)",              "cum_colonised_from_I",         "n",
    "Transmission",           "(ii) Patient-to-patient colonisations from I (%)",              "prop_p2p_from_I",              "pct",
    "Transmission",           "(ii) Patient-to-patient colonisations from Xtest (n)",          "cum_colonised_from_Xtest",     "n",
    "Transmission",           "(ii) Patient-to-patient colonisations from Xtest (%)",          "prop_p2p_from_Xtest",          "pct",
    "Transmission",           "(ii) Patient-to-patient colonisations from FN (n)",             "cum_colonised_from_FN",        "n",
    "Transmission",           "(ii) Patient-to-patient colonisations from FN (%)",             "prop_p2p_from_FN",             "pct",
    "Transmission",           "(ii) Patient-to-patient colonisations from Xtreat (n)",         "cum_colonised_from_Xtreat",    "n",
    "Transmission",           "(ii) Patient-to-patient colonisations from Xtreat (%)",         "prop_p2p_from_Xtreat",         "pct",
    "Transmission",           "(ii) Colonisation prevalence at admission (%)",                 "prop_admitted_C",              "pct",
    "Transmission",           "(iii) Colonised patients discharged (n)",                        "cum_discharges_colonised",     "n",
    "Transmission",           "(iii) Colonisation prevalence at discharge (%)",                 "prop_discharged_C",            "pct",
    "Transmission",           "(iii) Difference in colonisation prevalence (pp) [2]",           "diff_prev_C",                  "pp",
    "Transmission",           "(iii) CDI prevalence at admission (%)",                          "prop_admitted_I",              "pct",
    "Transmission",           "(iii) CDI patients discharged untreated (n)",                    "cum_discharges_while_infected", "n",
    "Transmission",           "(iii) CDI prevalence at discharge (%)",                          "prop_discharged_CDI",          "pct",
    "Transmission",           "(iii) Difference in CDI prevalence (pp) [2]",                   "diff_prev_CDI",                "pp",
    
    "Hospital resource use",  "(i) Total hospital bed-days (n)",                                "hospital_bed_days",            "n",
    "Hospital resource use",  "(ii) Average length of stay (days)",                              "avg_LOS",                      "dec",
    "Hospital resource use",  "(ii) Side-room bed-days (n)",                                    "side_room_bed_days",           "n",
    "Hospital resource use",  "(ii) Side-room bed-days (% of total)",                           "prop_side_room_bed_days",      "pct",
    "Hospital resource use",  "(ii) Side-room bed-days: suspected cases (n)",                   "side_room_Xtest_bed_days",     "n",
    "Hospital resource use",  "(ii) Side-room bed-days: suspected cases (% of total)",          "prop_side_room_Xtest_bed_days", "pct",
    "Hospital resource use",  "(ii) Side-room bed-days: confirmed cases (n)",                   "side_room_Xtreat_bed_days",    "n",
    "Hospital resource use",  "(ii) Side-room bed-days: confirmed cases (% of total)",          "prop_side_room_Xtreat_bed_days", "pct",
    "Hospital resource use",  "(iii) Faecal specimens tested (n)",                              "cum_total_test",               "n",
    "Hospital resource use",  "(iii) Faecal specimens tested per 10,000 patient bed-days",      "rate_tests",                   "dec",
    
    "Mortality",              "(i) Total hospital deaths (n)",                                  "cum_total_deaths",             "n",
    "Mortality",              "(ii) Deaths following CDI (n)",                                  "cum_deaths_post_infection",    "n",
    "Mortality",              "(ii) Deaths following CDI (% of all deaths)",                    "prop_deaths_post_infection",   "pct",
    "Mortality",              "(iii) Deaths among missed CDI cases (n)",                        "cum_deaths_missed_cases",      "n",
    "Mortality",              "(iii) Deaths among missed CDI cases (% of all deaths)",          "prop_deaths_missed_cases",     "pct"
  )
  if (per_year) {
    outcomes_spec <- outcomes_spec |> mutate(label = sub("(n)", "(n per year)", label, fixed = TRUE))
  }

  fmt_val <- function(x, fmt) {
    dplyr::case_when(
      fmt == "pct" ~ scales::percent(x, accuracy = 0.1),
      fmt == "pp"  ~ paste0(sprintf("%+.1f", x * 100), " pp"),
      fmt == "n"   ~ scales::comma(round(x)),
      fmt == "dec" ~ sprintf("%.1f", x),
      TRUE         ~ as.character(round(x, 3))
    )
  }
  
  scenario_cols <- fo |>
    select(scenario, all_of(outcomes_spec$variable)) |>
    pivot_longer(-scenario, names_to = "variable", values_to = "value")
  
  outcomes_spec |>
    left_join(scenario_cols, by = "variable") |>
    mutate(value_fmt = fmt_val(value, fmt)) |>
    select(group, label, scenario, value_fmt) |>
    pivot_wider(names_from = scenario, values_from = value_fmt) |>
    rename("Research question" = group, "Outcome" = label)
}

## c) Final outcomes vs base scenario ####

# Takes the final_outcomes tibble and returns absolute differences (scenario minus base) vs a
# chosen base scenario for all numeric columns. base_scenario is a scenario label and must be
# supplied, e.g. "P1_none_cov050_sens050".
#
# Column names are prefixed with "d_" and suffixed with "_vs_{base_scenario}".
# The "final_" and "cum_" prefixes are stripped for brevity, e.g. with
# base_scenario = "P1_none_cov050_sens050":
#   final_S             -> d_S_vs_P1_none_cov050_sens050
#   cum_confirmed_cases -> d_confirmed_cases_vs_P1_none_cov050_sens050
# Section 11 of the QMD then renames the suffix to the comparison, e.g. "_P2_vs_P1".

make_final_outcomes_vs_base <- function(final_outcomes, base_scenario) {
  base_row <- final_outcomes %>% filter(scenario == base_scenario) %>% slice(1)
  if (nrow(base_row) == 0) {
    stop("base_scenario '",
         base_scenario,
         "' not found in final_outcomes.")
  }
  suffix   <- paste0("_vs_", base_scenario)
  num_cols <- setdiff(names(final_outcomes), "scenario")
  
  final_outcomes %>%
    mutate(across(all_of(num_cols), ~ .x - base_row[[cur_column()]])) %>%
    rename_with(~ paste0("d_", sub("^(final_|cum_)", "", .x), suffix), all_of(num_cols)) %>%
    select(scenario, starts_with("d_"))
}

# ================================================================================================ #

# 7) Plot data preparation functions ####

## a) Long-format state and bed occupancy data for state trajectory line plots ####

# Reshapes scenario_post into long format for state trajectory plots.
# Filtered to times strictly before longest_stabilisation.
make_plot_long_states <- function(scenario_post, longest_stabilisation) {
  bind_rows(lapply(names(scenario_post), function(nm) {
    scenario_post[[nm]] %>%
      select(time, S, C, I, Xtest, FN, Xtreat) %>%
      filter(time < longest_stabilisation) %>%
      pivot_longer(
        cols = c(S, C, I, Xtest, FN, Xtreat),
        names_to = "state",
        values_to = "count"
      ) %>%
      mutate(scenario = nm, 
             state = factor(state, levels = c("S", "C", "I", "Xtest", "FN", "Xtreat")))
  })) %>%
    add_scenario_facets()
}

## b) Long-format bed occupancy data for bed occupancy line plots ####

make_plot_long_bed_occupancy <- function(scenario_post, longest_stabilisation) {
  bind_rows(lapply(names(scenario_post), function(nm) {
    scenario_post[[nm]] %>%
      select(time, hospital_bed_occupancy, general_ward_bed_occupancy, side_room_bed_occupancy) %>%
      filter(time < longest_stabilisation) %>%
      pivot_longer(
        cols = c(hospital_bed_occupancy, general_ward_bed_occupancy, side_room_bed_occupancy),
        names_to = "bed_type",
        values_to = "occupancy"
      ) %>%
      mutate(scenario = nm)
  })) %>%
    add_scenario_facets()
}

## c) Long-format flow data for (instantaneous) flows and cumulative plots ####

# Collects lambda, all instantaneous rate outputs and all cumulative outputs from scenario_post into 
# a single combined data frame.
# Uses tidy-select helpers so new ODE outputs are picked up automatically.

make_plot_long_flows <- function(scenario_post) {
  bind_rows(lapply(names(scenario_post), function(nm) {
    scenario_post[[nm]] %>%
      select(time,
             lambda,
             ends_with("_rate"),
             starts_with("cum_")) %>%
      mutate(scenario = nm)
  })) %>%
    add_scenario_facets()
}

# ================================================================================================ #

# 8) Colour palette and pretty labels ####

# Named colour vector used consistently across all model plots.
# Based on RColorBrewer Set2 and Set3 palettes.
# Plot functions accept a colour_palette argument defaulting to this object,
# so it can be overridden without editing the functions.

model_colour_palette <- c(
  
  # States 
  "S"      = "#66C2A5FF",
  "C"      = "#FC8D62FF",
  "I"      = "#8DA0CBFF",
  "Xtest"  = "#E78AC3FF",
  "FN"     = "#A6D854FF",
  "Xtreat" = "#FFD92FFF",
  "D"      = "#B3B3B3FF",
  
  # Bed occupancy
  "hospital_bed_occupancy"      = "#66C2A5FF",
  "general_ward_bed_occupancy" = "#8DA0CBFF",
  "side_room_bed_occupancy"  = "#FFD92FFF",
  
  # Flows into C and I
  "colonised_in_hospital_rate" = "#FC8D62FF",
  "admitted_colonised_rate"    = "#FC8D62FF",
  "infected_in_hospital_rate" = "#8DA0CBFF",
  "admitted_infected_rate"    = "#8DA0CBFF",
  
  # Case ascertainment flows
  "total_infected_incidence_rate"   = "#8DA0CBFF",
  "suspected_cases_rate"            = "#E78AC3FF",
  "false_negatives_rate"            = "#A6D854FF",
  "confirmed_cases_rate"            = "#FFD92FFF",
  
  # Testing
  "first_test_rate"             = "grey80",
  "retest_false_negatives_rate" = "grey60",
  "total_test_rate"             = "grey40",
  
  # Missed cases flows
  "missed_untested_cases_rate"       = "#8DA0CBFF",
  "missed_false_negative_cases_rate" = "#A6D854FF",
  "total_missed_cases_rate"          = "#B3B3B3FF",
  
  # Deaths by state
  "deaths_susceptible_rate"     = "#66C2A5FF",
  "deaths_colonised_rate"       = "#FC8D62FF",
  "deaths_infected_rate"        = "#8DA0CBFF",
  "deaths_suspected_cases_rate" = "#E78AC3FF",
  "deaths_FN_rate"              = "#A6D854FF",
  "deaths_confirmed_cases_rate" = "#FFD92FFF",
  
  # Discharges by state
  "discharges_susceptible_rate" = "#66C2A5FF",
  "discharges_colonised_rate"   = "#FC8D62FF",
  "discharges_infected_rate"    = "#8DA0CBFF",
  "discharges_FN_rate"          = "#A6D854FF",
  "discharges_recovered_rate"   = "#FFD92FFF",
  
  # Case outcome flows
  "discharges_while_infected_rate"  = "#8DA0CBFF",
  "recovered_cases_rate"            = "#FFD92FFF",
  "deaths_post_infection_rate"      = "#B3B3B3FF"
)

model_pretty_labels <- c(
  
  # States (line and area plots)
  "S"      = "Susceptible",
  "C"      = "Colonised",
  "I"      = "Infected",
  "Xtest"  = "Suspected case",
  "FN"     = "False negative",
  "Xtreat" = "Confirmed case",
  "D"      = "Dead",
  
  # Bed occupancy (line plot)
  "hospital_bed_occupancy"      = "All beds",
  "general_ward_bed_occupancy" = "General ward beds",
  "side_room_bed_occupancy"  = "Side-room beds",
  
  # Flows into C and I (line plot)
  "colonised_in_hospital_rate" = "Colonisations in hospital",
  "admitted_colonised_rate"    = "Colonised admissions",
  "infected_in_hospital_rate"  = "Progressions to infection in hospital",
  "admitted_infected_rate"     = "Infected admissions",
  
  # Case ascertainment flows (line plot)
  "total_infected_incidence_rate" = "Infected cases",
  "suspected_cases_rate"          = "Suspected cases",
  "false_negatives_rate"          = "False negatives",
  "retest_false_negatives_rate"   = "Retested false negatives",
  "confirmed_cases_rate"          = "Confirmed cases",
  
  # Testing (line plot)
  "first_test_rate"             = "First tests",
  "retest_false_negatives_rate" = "Retests of false negatives",
  "total_test_rate"             = "Total tests",
  
  # Missed cases flows (line plot)
  "missed_untested_cases_rate"       = "Missed cases (untested)",
  "missed_false_negative_cases_rate" = "Missed cases (false negatives)",
  "total_missed_cases_rate"          = "Total missed cases",
  
  # Deaths by state (area plot)
  "deaths_susceptible_rate"     = "Susceptible",
  "deaths_colonised_rate"       = "Colonised",
  "deaths_infected_rate"        = "Infected",
  "deaths_suspected_cases_rate" = "Suspected cases",
  "deaths_FN_rate"              = "False negatives",
  "deaths_confirmed_cases_rate" = "Confirmed cases",
  
  # Discharges by state (area plot)
  "discharges_susceptible_rate" = "Susceptible",
  "discharges_colonised_rate"   = "Colonised",
  "discharges_infected_rate"    = "Infected",
  "discharges_FN_rate"          = "False negatives",
  "discharges_recovered_rate"   = "Recovered",
  
  # Case outcome flows (area plot)
  "discharges_while_infected_rate"  = "Discharges from infected",
  "recovered_cases_rate"            = "Recovered cases",
  "deaths_post_infection_rate"      = "Deaths post-infection",
 
  # Cumulative
  
  # Cumulative bed occupancy (bar plots)
  "hospital_bed_days"        = "Total patient bed-days",
  "general_ward_bed_days" = "General ward patient bed-days",
  "side_room_bed_days"    = "Side-room patient bed-days",
  "side_room_Xtest_bed_days"  = "Side-room bed-days: suspected cases",
  "side_room_Xtreat_bed_days" = "Side-room bed-days: confirmed cases",
  
  # Cumulative flows into C and I (bar plots)
  "cum_colonised_in_hospital"     = "Colonisations (in hospital)",
  "admitted_colonised"            = "Colonised admissions",
  "cum_infected_in_hospital"      = "Progressions to infection (in hospital)",
  "admitted_infected"             = "Infected admissions",
  
  # Cumulative case ascertainment (bar plots)
  "cum_total_infected_incidence"   = "Infected cases",
  "cum_suspected_cases"            = "Suspected cases",
  "cum_false_negatives"            = "False negatives",
  "cum_retest_false_negatives"     = "Retested false negatives",
  "cum_confirmed_cases"            = "Confirmed cases",
  
  # Cumulative testing (bar plots)
  "cum_first_test"             = "First tests",
  "cum_retest_false_negatives" = "Retests of false negatives",
  "cum_total_test"             = "Total tests",
  
  # Cumulative missed cases (bar plots)
  "cum_missed_untested_cases"       = "Untested cases",
  "cum_missed_false_negative_cases" = "False negatives",
  "cum_total_missed_cases"          = "Total missed cases",
  
  # Cumulative confirmed cases via first test vs retest (bar plots)
  "cum_confirmed_via_first_test"    = "Cases confirmed via first test",
  "cum_confirmed_via_retest"        = "Cases confirmed via retest",
  
  # Cumulative deaths by state (bar plots)
  "cum_deaths_susceptible"     = "Susceptible",
  "cum_deaths_colonised"       = "Colonised",
  "cum_deaths_infected"        = "Infected",
  "cum_deaths_suspected_cases" = "Suspected cases",
  "cum_deaths_FN"              = "False negatives",
  "cum_deaths_confirmed_cases" = "Confirmed cases",

  # Cumulative discharges by state (bar plots)
  "cum_discharges_susceptible" = "Susceptible",
  "cum_discharges_colonised"   = "Colonised",
  "cum_discharges_infected"    = "Infected",
  "cum_discharges_FN"          = "False negatives",
  "cum_discharges_recovered"   = "Recovered",
  
  # Cumulative case outcomes (bar plots)
  "cum_discharges_while_infected" = "Discharges while infected",
  "cum_recovered_cases"           = "Recovered cases",
  "cum_deaths_post_infection"     = "Deaths post infection"
)

# ================================================================================================ #

# 9) Plot functions ####

# Each function accepts pre-built long-format data (from make_plot_long_states
# or make_plot_long_flows) and returns a ggplot object.
# Saving to file is handled separately in the calling script via ggsave().

font_base_size <- 4

my_theme <- theme(
  plot.title      = element_text(size = font_base_size),
  strip.text.x    = element_text(size = font_base_size),
  strip.text.y    = element_text(size = font_base_size),
  axis.title      = element_text(size = font_base_size),
  axis.text       = element_text(size = font_base_size),
  legend.text     = element_text(size = font_base_size),
  legend.title    = element_text(size = font_base_size),
  legend.key.size = unit(0.3, "lines")
)

## a) State trajectories - line plot ####

plot_states_line <- function(plot_long_states,
                             colour_palette = model_colour_palette,
                             pretty_labels  = model_pretty_labels) {
  ggplot(plot_long_states, aes(time, count, color = state)) +
    geom_line(linewidth = 0.5) +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_color_manual(values = colour_palette, labels = pretty_labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Hospital patients, trajectories by state",
      x     = "Days (up to model stabilisation)",
      y     = "Count",
      color = "State"
    )
}

## b) State trajectories - line plot without the S state ####

plot_states_line_noS <- function(plot_long_states,
                             colour_palette = model_colour_palette,
                             pretty_labels  = model_pretty_labels) {
  plot_long_states %>%
    filter(state != "S") %>%
    ggplot(aes(time, count, color = state)) +
    geom_line(linewidth = 0.5) +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_color_manual(values = colour_palette, labels = pretty_labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Hospital patients (excluding susceptible), trajectories by state",
      x     = "Days (up to model stabilisation)",
      y     = "Count",
      color = "State"
    )
}

## c) State trajectories - area plot ####

plot_states_area <- function(plot_long_states,
                             colour_palette = model_colour_palette,
                             pretty_labels  = model_pretty_labels) {
  ggplot(plot_long_states, aes(time, count, fill = state)) +
    geom_area() +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_fill_manual(values = colour_palette, labels = pretty_labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Hospital patients, composition by state",
      x     = "Days (up to model stabilisation)",
      y     = "Count",
      fill  = "State"
    )
}

## d) State trajectories - area plot without the S state ####

plot_states_area_noS <- function(plot_long_states,
                             colour_palette = model_colour_palette,
                             pretty_labels  = model_pretty_labels) {
  plot_long_states %>%
    filter(state != "S") %>%
    ggplot(aes(time, count, fill = state)) +
    geom_area() +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_fill_manual(values = colour_palette, labels = pretty_labels) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Hospital patients (excluding susceptible), composition by state",
      x     = "Days (up to model stabilisation)",
      y     = "Count",
      fill  = "State"
    )
}

## e) Bed occupancy - area plot ####

# General ward and side-room beds stacked, so the top of each area is total bed occupancy.
plot_bed_occupancy <- function(plot_long_bed_occupancy,
                               longest_stabilisation,
                               colour_palette = model_colour_palette,
                               pretty_labels  = model_pretty_labels) {
    plot_long_bed_occupancy %>%
    filter(time < longest_stabilisation,
           bed_type %in% c("general_ward_bed_occupancy", "side_room_bed_occupancy")) %>%
    mutate(resource = factor(
      bed_type,
      levels = c("side_room_bed_occupancy",
                 "general_ward_bed_occupancy"))) %>%
    ggplot(aes(x = time, y = occupancy, fill = resource)) +
    geom_area() +
    example_facets +
    scale_fill_manual(values = colour_palette, labels = pretty_labels) +
    scale_y_continuous(labels = scales::comma) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    labs(
      title = "Bed occupancy (general ward and side-room beds; total = top of the stack)",
      x     = "Days (up to model stabilisation)",
      y     = "Occupied beds",
      fill  = NULL
    )
}
## g) New colonisations and infections - line plot ####

model_linetype_palette_C_and_I <- c(
  "colonised_in_hospital_rate" = "solid",
  "admitted_colonised_rate"    = "dotted",
  "infected_in_hospital_rate"  = "solid",
  "admitted_infected_rate"     = "dotted"
)

plot_flows_to_C_and_I <- function(plot_long_flows,
                                  longest_stabilisation,
                                  final_outcomes    = NULL,
                                  colour_palette   = model_colour_palette,
                                  linetype_palette = model_linetype_palette_C_and_I,
                                  pretty_labels    = model_pretty_labels) {
  # Compute y position from actual data range to avoid panel clipping
  flow_data <- plot_long_flows %>%
    filter(time < longest_stabilisation) %>%
    select(colonised_in_hospital_rate, admitted_colonised_rate,
           infected_in_hospital_rate, admitted_infected_rate)
  y_min <- min(unlist(flow_data), na.rm = TRUE)
  
  label_df <- if (!is.null(final_outcomes)) {
    final_outcomes %>%
      select(scenario, prop_I_hospital_attributable) %>%
      add_scenario_facets() %>%
      mutate(
        label = paste0("Hosp.-attributable CDI: ",
                       scales::percent(prop_I_hospital_attributable, accuracy = 0.1)),
        x_pos = 0,
        y_pos = y_min + 2  
## !! May need to change y_pos ####
      )
  } else {
    NULL
  }
  
  p <- plot_long_flows %>%
    filter(time < longest_stabilisation) %>%
    select(
      time, scenario, policy, setting,
      colonised_in_hospital_rate, admitted_colonised_rate,
      infected_in_hospital_rate,  admitted_infected_rate
    ) %>%
    pivot_longer(cols = -c(time, scenario, policy, setting), names_to = "flow", values_to = "value") %>%
    mutate(flow = factor(flow, levels = c(
      "colonised_in_hospital_rate", "admitted_colonised_rate",
      "infected_in_hospital_rate",  "admitted_infected_rate"))) %>%
    ggplot(aes(time, value, color = flow, linetype = flow)) +
    geom_line(linewidth = 0.5) +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_color_manual(values = colour_palette, labels = pretty_labels) +
    scale_linetype_manual(values = linetype_palette, labels = pretty_labels) +
    labs(
      title    = "New colonisations and infections",
      x        = "Days (up to model stabilisation)",
      y        = "Rate (per day)",
      color    = "Flow type",
      linetype = "Flow type"
    )
  
  if (!is.null(label_df)) {
    p <- p + geom_text(data = label_df,
                       aes(x = x_pos, y = y_pos, label = label),
                       inherit.aes = FALSE,
                       hjust = 0, vjust = 1,
                       size = font_base_size / .pt)
  }
  p
}

## h) Case incidence and ascertainment - line plot ####

plot_case_ascertainment_line <- function(plot_long_flows,
                             longest_stabilisation,
                             colour_palette = model_colour_palette,
                             pretty_labels  = model_pretty_labels) {
  plot_long_flows %>%
    filter(time < longest_stabilisation) %>%
    select(
      time,
      scenario, policy, setting,
      total_infected_incidence_rate,
      suspected_cases_rate,
      confirmed_cases_rate,
      false_negatives_rate,
      retest_false_negatives_rate
      ) %>%
    pivot_longer(
      cols = -c(time, scenario, policy, setting),
      names_to = "flow",
      values_to = "value"
    ) %>%
    mutate(flow = factor(
      flow,
      levels = c(
        "total_infected_incidence_rate",
        "suspected_cases_rate",
        "confirmed_cases_rate",
        "false_negatives_rate",
        "retest_false_negatives_rate"
      )
    )) %>%
    ggplot(aes(time, value, color = flow)) +
    geom_line(linewidth = 0.5) +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_color_manual(values = colour_palette, labels = pretty_labels) +
    labs(
      title = "Case incidence and ascertainment rates",
      x     = "Days (up to model stabilisation)",
      y     = "Rate (count per day)",
      color = "Rate"
    )
}

## j) Faecal specimens tested per year: first tests and retests - stacked bar plot ####

# Stacked bar plots (sections j, m, n, s, u) show the example scenarios grouped by coverage and
# sensitivity, with P1, P2 and P3 side by side (see plot_example_stacked() in section 10).
plot_testing_bar <- function(final_outcomes) {
  plot_example_stacked(final_outcomes,
                       components = c("First tests"                = "cum_first_test",
                                      "Retests of false negatives" = "cum_retest_false_negatives"),
                       title   = "Faecal specimens tested per year: first tests and retests",
                       y_label = "Faecal specimens tested per year",
                       colours = c("grey70", "grey40"))
}

## m) Missed cases per year: untested and false negatives - stacked bar plot ####

plot_missed_cases_bar <- function(final_outcomes, colour_palette = model_colour_palette) {
  plot_example_stacked(final_outcomes,
                       components = c("Untested (discharged or died in I)"         = "cum_missed_untested_cases",
                                      "False negatives (discharged or died in FN)" = "cum_missed_false_negative_cases"),
                       title   = "Missed CDI cases per year: untested and false negatives",
                       y_label = "Missed CDI cases per year",
                       colours = unname(colour_palette[c("I", "FN")]))
}

## n) Confirmed cases per year: via first test and retest - stacked bar plot ####

plot_case_confirmation_bar <- function(final_outcomes) {
  plot_example_stacked(final_outcomes,
                       components = c("Confirmed via first test" = "cum_confirmed_via_first_test",
                                      "Confirmed via retest"     = "cum_confirmed_via_retest"),
                       title   = "Confirmed CDI cases per year: via first test and retest",
                       y_label = "Confirmed CDI cases per year",
                       colours = c("#FFD92FFF", "#E5A800FF"),
                       caption = "Retests are apportioned post hoc by the share of all tests that were retests.")
}

## p) Deaths by state - area plot (omitting susceptible) ####

plot_deaths_by_state_area_noS <- function(plot_long_flows,
                             longest_stabilisation,
                             colour_palette = model_colour_palette,
                             pretty_labels  = model_pretty_labels) {
  plot_long_flows %>%
    filter(time < longest_stabilisation) %>%
    select(
      time,
      scenario, policy, setting,
      deaths_colonised_rate,
      deaths_infected_rate,
      deaths_suspected_cases_rate,
      deaths_FN_rate,
      deaths_confirmed_cases_rate
    ) %>%
    pivot_longer(
      cols = -c(time, scenario, policy, setting),
      names_to = "flow",
      values_to = "value"
    ) %>%
    mutate(flow = factor(
      flow,
      levels = c(
        "deaths_colonised_rate",
        "deaths_infected_rate",
        "deaths_suspected_cases_rate",
        "deaths_FN_rate",
        "deaths_confirmed_cases_rate"
      )
    )) %>%
    ggplot(aes(time, value, fill = flow)) +
    geom_area() +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_fill_manual(values = colour_palette, labels = pretty_labels) +
    labs(
      title = "Death rates by state (excluding susceptible)",
      x     = "Days (up to model stabilisation)",
      y     = "Death rate (per day)",
      fill  = "State"
    )
}

## r) Discharges by state (except susceptible) - area plot ####

plot_discharges_by_state_area_noS <- function(plot_long_flows,
                             longest_stabilisation,
                             colour_palette = model_colour_palette,
                             pretty_labels  = model_pretty_labels) {
  plot_long_flows %>%
    filter(time < longest_stabilisation) %>%
    select(
      time,
      scenario, policy, setting,
      discharges_colonised_rate,
      discharges_infected_rate,
      discharges_FN_rate,
      discharges_recovered_rate
    ) %>%
    pivot_longer(
      cols = -c(time, scenario, policy, setting),
      names_to = "flow",
      values_to = "value"
    ) %>%
    mutate(flow = factor(
      flow,
      levels = c(
        "discharges_colonised_rate",
        "discharges_infected_rate",
        "discharges_FN_rate",
        "discharges_recovered_rate"
      )
    )) %>%
    ggplot(aes(time, value, fill = flow)) +
    geom_area() +
    example_facets +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_fill_manual(values = colour_palette, labels = pretty_labels) +
    labs(
      title = "Discharge rates by state (excluding susceptible)",
      x     = "Days (up to model stabilisation)",
      y     = "Discharge rate (per day)",
      fill  = "State"
    )
}

## s) Outcomes of CDI cases per year - stacked bar plot ####

# Exits of CDI patients: discharged untreated (from I and FN), recovered after treatment (returned to
# a general ward or discharged from Xtreat) and deaths following CDI (from I, Xtest, FN and Xtreat).
plot_case_outcomes_bar <- function(final_outcomes) {
  plot_example_stacked(final_outcomes,
                       components = c("Discharged untreated (from I and FN)" = "cum_discharges_while_infected",
                                      "Recovered after treatment"             = "cum_recovered_cases",
                                      "Died following CDI"                    = "cum_deaths_post_infection"),
                       title   = "Outcomes of CDI cases per year",
                       y_label = "CDI cases per year",
                       colours = c("#8DA0CBFF", "#FFD92FFF", "#B3B3B3FF"))
}

## u) Hospital-acquired colonisations by source per year - stacked bar plot ####

# Stacks hospital-acquired colonisations (S -> C) from background acquisition (beta0) and
# patient-to-patient transmission (secondary colonisations) from each state. The stacks sum to
# cum_colonised_in_hospital. These are colonisations, not CDI cases. For modelling simplicity the
# background rate does not depend on the number of patients in each state.

plot_colonisation_sources_bar <- function(final_outcomes, colour_palette = model_colour_palette) {
  plot_example_stacked(final_outcomes,
                       components = c("Background"              = "cum_colonised_background",
                                      "Colonised (C)"           = "cum_colonised_from_C",
                                      "Infected (I)"            = "cum_colonised_from_I",
                                      "Suspected case (Xtest)"  = "cum_colonised_from_Xtest",
                                      "False negative (FN)"     = "cum_colonised_from_FN",
                                      "Confirmed case (Xtreat)" = "cum_colonised_from_Xtreat"),
                       title   = "Hospital-acquired colonisations per year by source",
                       y_label = "Hospital-acquired colonisations per year",
                       colours = c("#B3B3B3FF", unname(colour_palette[c("C", "I", "Xtest", "FN", "Xtreat")])),
                       caption = "Background acquisition does not depend on the number of patients in each state (modelling simplification).")
}

## t) Control parameter values - scatter plot ####

plot_scenario_params <- function(control_parms) {
  control_parms |>
    mutate(isolation_policy = factor(isolation_policy, levels = isolation_policy_levels)) |>
    ggplot(aes(test_sens, prop_I_suspected)) +
    geom_point() +
    geom_text(aes(label = scenario), vjust = -0.7, size = font_base_size / .pt) +
    facet_wrap(~ isolation_policy, labeller = as_labeller(isolation_policy_labels)) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_x_continuous(limits = c(0, 1.1), labels = scales::percent) +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
    labs(
      title = "Values of control parameters",
      x     = "Test sensitivity",
      y     = coverage_axis_label_wrapped
    )
}

# ================================================================================================ #

# 10) Scenario grid ####

# The scenario grid crosses every isolation policy with every testing coverage (prop_I_suspected)
# and test sensitivity value. Outcomes are totals over the analysis period after the burn-in, shown
# in the QMD as annual totals.

## a) Isolation policy levels, labels and comparisons ####

# Policies are numbered P1-P3 in this order, from least to most isolation; plots, tables and
# comparisons follow it. P1 (no isolation) is the reference policy.
isolation_policy_levels <- c("none", "confirmed", "suspected_and_confirmed")

isolation_policy_labels <- c(
  none                    = "P1: No isolation",
  confirmed               = "P2: Confirmed cases isolated",
  suspected_and_confirmed = "P3: Suspected and confirmed cases isolated"
)

# Short policy codes used in scenario labels
isolation_policy_codes <- c(
  none                    = "P1_none",
  confirmed               = "P2_confirmed",
  suspected_and_confirmed = "P3_susp_conf"
)

# Reference policy (P1, no isolation)
reference_policy <- "none"

# Policy comparisons: difference = policy minus comparator, at the same testing coverage and test
# sensitivity (negative = less of the outcome than under the comparator). P2 and P3 are compared
# with P1 (no isolation), and P3 with P2 (the incremental effect of also isolating suspected cases).
policy_comparisons <- tribble(
  ~policy,                   ~comparator,
  "confirmed",               "none",
  "suspected_and_confirmed", "none",
  "suspected_and_confirmed", "confirmed"
) |>
  mutate(label = paste(sub(":.*", "", isolation_policy_labels[policy]), "vs",
                       sub(":.*", "", isolation_policy_labels[comparator])))

## b) Build the scenario grid ####

# Returns one row per combination, with a scenario label such as "P2_confirmed_cov050_sens100".
make_scenario_grid <- function(isolation_policy, prop_I_suspected, test_sens) {
  bad <- setdiff(isolation_policy, isolation_policy_levels)
  if (length(bad) > 0) {
    stop("Unknown isolation_policy value(s): ", paste(bad, collapse = ", "), ".")
  }
  expand_grid(
    isolation_policy = isolation_policy,
    prop_I_suspected = prop_I_suspected,
    test_sens        = test_sens
  ) |>
    mutate(
      scenario = sprintf("%s_cov%03d_sens%03d",
                         isolation_policy_codes[isolation_policy],
                         round(100 * prop_I_suspected),
                         round(100 * test_sens)),
      .before = 1
    )
}

## c) Run every grid scenario, keeping summaries only ####

# A full post-processed time series is thousands of rows x ~170 columns per scenario, so keeping one for
# every grid scenario would use a lot of memory. Each run is summarised as soon as it finishes;
# re-run selected scenarios with run_scenario() when time series are needed.
#
# n_cores > 1 runs scenarios in parallel on a local cluster of R worker processes (base R
# 'parallel'; works on Windows). Each worker gets this session's library paths, loads tidyverse and
# deSolve, and receives a copy of every object in the environment where the model functions were
# sourced (the functions, fixed, mod.init and mod.t), so the workers run the same code and inputs.
# Results are identical to a serial run.
#
# burn_in (days) is passed to add_quantities(), so cumulative outcomes cover the period after it.
#
# Returns a list of combined tibbles: parms, final, patient_days, stability, negatives.
run_scenario_grid <- function(grid, time_step, burn_in = 0, n_cores = 1) {
  run_one <- function(scenario, isolation_policy, prop_I_suspected, test_sens) {
    parms <- suppressMessages(make_mod_parms(
      isolation_policy = isolation_policy,
      prop_I_suspected = prop_I_suspected,
      test_sens        = test_sens,
      scenario_label   = scenario
    ))
    post <- run_scenario(parms) |> add_quantities(time_step = time_step, burn_in = burn_in)
    list(
      parms        = as_tibble(parms),
      final        = make_final_outcomes(set_names(list(post), scenario)),
      patient_days = calculate_total_patient_days(post) |> mutate(scenario = scenario, .before = 1),
      stability    = suppressMessages(stabilisation_summary(post, scenario)),
      negatives    = suppressMessages(negative_value_rows(post)) |>
        mutate(scenario = rep(scenario, n()), .before = 1)
    )
  }

  if (n_cores > 1) {
    model_env <- environment(Model1A_ODE)
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(lib) {
      .libPaths(lib)
      suppressPackageStartupMessages({
        library(tidyverse)
        library(deSolve)
      })
      NULL
    }, .libPaths())
    parallel::clusterExport(cl, ls(envir = model_env), envir = model_env)
    # Load-balanced, as run times vary between scenarios
    runs <- parallel::parLapplyLB(cl, seq_len(nrow(grid)), \(i) do.call(run_one, as.list(grid[i, ])))
  } else {
    runs <- pmap(grid, run_one, .progress = "Running scenario grid")
  }

  parts <- c("parms", "final", "patient_days", "stability", "negatives")
  set_names(lapply(parts, \(x) bind_rows(map(runs, x))), parts)
}

## d) Combine grid inputs and cumulative outcomes ####

# One row per scenario: control values, then final cumulative outcomes and patient-days.
make_grid_results <- function(grid, grid_runs) {
  grid |>
    left_join(grid_runs$final,        by = "scenario") |>
    left_join(grid_runs$patient_days, by = "scenario")
}

## e) Outcomes compared across the grid ####

# Shown as annual totals (mean per year of the analysis period after the burn-in). Definitions:
#   - CDI cases: all incident cases, admitted with CDI plus hospital-onset.
#   - Hospital-onset CDI: all in-hospital progressions from C to I (patients admitted colonised plus
#     patients colonised in hospital).
#   - Hospital-acquired CDI: the post-hoc attribution in make_final_outcomes() (in-hospital
#     progressions x share of colonisations acquired in hospital), a subset of hospital-onset CDI.
#   - Missed CDI cases are exits from I and FN by discharge or death (never diagnosed).
# fmt sets how values are labelled: "n" (counts), "rate" (one decimal place) or "pct" (percent).
# plot = TRUE: plotted in the QMD's tabs of outcomes by testing coverage and test sensitivity.
# Outcomes with plot = FALSE follow the same pattern as a plotted one (e.g. hospital-onset and
# hospital-acquired CDI follow all CDI cases), so they are only tabulated.
grid_outcome_spec <- tribble(
  ~variable,                       ~label,                               ~fmt, ~plot,
  "cum_total_infected_incidence",  "CDI cases",                          "n",  TRUE,
  "cum_infected_in_hospital",      "Hospital-onset CDI cases",           "n",  FALSE,
  "cum_hospital_attributable_I",   "Hospital-acquired CDI cases",        "n",  FALSE,
  "cum_total_missed_cases",        "Missed CDI cases",                   "n",  TRUE,
  "cum_deaths_post_infection",     "Deaths following CDI",               "n",  TRUE,
  "cum_colonised_in_hospital",     "Hospital-acquired colonisations",    "n",  FALSE,
  "cum_discharges_colonised",      "Discharged colonisations",           "n",  FALSE,
  "hospital_bed_days",             "Patient bed-days",                   "n",  TRUE,
  "side_room_bed_days",            "Side-room patient bed-days",         "n",  TRUE,
  "cum_total_test",                "Faecal specimens tested",            "n",  TRUE
)

## f) Long-format grid outcomes and policy differences ####

# One row per scenario and outcome. value is divided by `years` (the length of the analysis period)
# to give annual totals; use years = 1 for outcomes that are already rates or proportions.
make_grid_outcomes_long <- function(grid_results, spec = grid_outcome_spec, years = 1) {
  grid_results |>
    select(scenario, isolation_policy, prop_I_suspected, test_sens, all_of(spec$variable)) |>
    pivot_longer(all_of(spec$variable), names_to = "outcome", values_to = "value") |>
    mutate(
      value            = value / years,
      outcome          = factor(outcome, levels = spec$variable),
      isolation_policy = factor(isolation_policy, levels = isolation_policy_levels)
    )
}

# Differences between policies at the same coverage and sensitivity, for each comparison in
# policy_comparisons (diff = policy minus comparator; comparisons whose policies are not both in
# the grid are skipped).
make_grid_diffs <- function(grid_long, comparisons = policy_comparisons) {
  wide <- grid_long |>
    select(outcome, prop_I_suspected, test_sens, isolation_policy, value) |>
    mutate(isolation_policy = as.character(isolation_policy)) |>
    pivot_wider(names_from = isolation_policy, values_from = value)
  comparisons |>
    filter(policy %in% names(wide), comparator %in% names(wide)) |>
    pmap(\(policy, comparator, label) wide |>
           transmute(outcome, prop_I_suspected, test_sens, comparison = label,
                     diff = .data[[policy]] - .data[[comparator]])) |>
    bind_rows() |>
    mutate(comparison = factor(comparison, levels = comparisons$label))
}

## g) Grid plots ####

# Larger text than the time-series plots above, as these are viewed one outcome at a time.
grid_font_size <- 9

grid_theme <- theme_classic(base_size = grid_font_size) +
  theme(
    legend.position       = "bottom",
    legend.title.position = "top",
    strip.background      = element_blank()
  )

# Axis/legend label for prop_I_suspected ("testing coverage"). The model uses this proportion to
# calculate the case identification rate gamma (see make_mod_parms()).
coverage_axis_label         <- "Proportion of CDI cases tested for CDI"
coverage_axis_label_wrapped <- "Proportion of CDI cases\ntested for CDI"

# Testing coverage values drawn as lines on the line plots (a subset of the grid, for legibility)
grid_line_coverage <- seq(0, 100, by = 20) / 100

# Categorical colour scale for coverage on the line plots: 0% yellow to 100% purple
grid_coverage_colour <- scale_colour_viridis_d(
  coverage_axis_label,
  end       = 0.95,
  direction = -1,
  guide     = guide_legend(nrow = 1)
)

# Keep only the coverage values in grid_line_coverage (rounded to avoid floating-point mismatches)
filter_line_coverage <- function(df, coverage = grid_line_coverage) {
  df |> filter(round(prop_I_suspected, 6) %in% round(coverage, 6))
}

# Coverage as an ordered factor of percentages, for the categorical legend
coverage_factor <- function(x) {
  factor(scales::percent(x, accuracy = 1), levels = scales::percent(sort(unique(x)), accuracy = 1))
}

# Axis and legend label functions by value format; diff = TRUE for differences between policies
# (percentages are then differences in percentage points)
fmt_labels <- function(fmt, diff = FALSE) {
  switch(fmt,
    pct  = if (diff) scales::label_number(scale = 100, accuracy = 0.1, suffix = " pp", style_positive = "plus")
           else scales::label_percent(accuracy = 1),
    rate = scales::label_number(accuracy = 0.1, big.mark = ",", style_positive = if (diff) "plus" else "none"),
    scales::label_comma(style_positive = if (diff) "plus" else "none"))
}

spec_row <- function(spec, outcome_var) spec[spec$variable == outcome_var, ]

### i) Line plot: outcome vs test sensitivity, one line per testing coverage ####

plot_grid_outcome_lines <- function(grid_long,
                                    outcome_var,
                                    coverage = grid_line_coverage,
                                    spec     = grid_outcome_spec,
                                    y_label  = "Per year") {
  s <- spec_row(spec, outcome_var)
  grid_long |>
    filter(outcome == outcome_var) |>
    filter_line_coverage(coverage) |>
    mutate(coverage = coverage_factor(prop_I_suspected)) |>
    ggplot(aes(test_sens, value, colour = coverage, group = coverage)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.7) +
    facet_wrap(~ isolation_policy, labeller = as_labeller(isolation_policy_labels)) +
    grid_coverage_colour +
    scale_x_continuous("Test sensitivity", labels = scales::percent) +
    scale_y_continuous(y_label, labels = fmt_labels(s$fmt)) +
    labs(title = s$label) +
    grid_theme
}

### ii) Line plot of differences between policies ####

# One panel per comparison with the reference policy (P2 vs P1, P3 vs P1), placed under P2 and P3 in
# the line plot of totals: the first panel (where P1 vs P1 would be) is left blank. Negative = less
# of the outcome than under P1.
plot_grid_diff_lines <- function(grid_diffs,
                                 outcome_var,
                                 coverage    = grid_line_coverage,
                                 spec        = grid_outcome_spec,
                                 comparisons = policy_comparisons$label[policy_comparisons$comparator == reference_policy]) {
  s <- spec_row(spec, outcome_var)
  grid_diffs |>
    filter(outcome == outcome_var, comparison %in% comparisons) |>
    filter_line_coverage(coverage) |>
    mutate(coverage   = coverage_factor(prop_I_suspected),
           comparison = factor(as.character(comparison), levels = c(" ", comparisons))) |>
    ggplot(aes(test_sens, diff, colour = coverage, group = coverage)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.7) +
    facet_wrap(~ comparison, drop = FALSE) +
    grid_coverage_colour +
    scale_x_continuous("Test sensitivity", labels = scales::percent) +
    scale_y_continuous("Difference (negative = fewer)", labels = fmt_labels(s$fmt, diff = TRUE)) +
    labs(title = paste0(s$label, ": differences from P1 (no isolation)")) +
    grid_theme
}

### iii) Heatmap: test sensitivity (x) by testing coverage (y) ####

# One fill scale shared across the policy panels, so colours are comparable between policies.
plot_grid_outcome_heatmap <- function(grid_long, outcome_var, spec = grid_outcome_spec,
                                      fill_label = NULL) {
  s <- spec_row(spec, outcome_var)
  grid_long |>
    filter(outcome == outcome_var) |>
    ggplot(aes(test_sens, prop_I_suspected, fill = value)) +
    geom_tile() +
    facet_wrap(~ isolation_policy, labeller = as_labeller(isolation_policy_labels)) +
    scale_fill_viridis_c(
      fill_label %||% s$label,
      labels    = fmt_labels(s$fmt),
      direction = -1,   # low values yellow, high values purple
      guide     = guide_colourbar(barwidth = unit(6, "cm"), barheight = unit(0.3, "cm"))
    ) +
    scale_x_continuous("Test sensitivity", labels = scales::percent, expand = c(0, 0)) +
    scale_y_continuous(coverage_axis_label_wrapped, labels = scales::percent, expand = c(0, 0)) +
    coord_fixed() +
    labs(title = s$label) +
    grid_theme
}

### iv) Heatmap of differences between policies ####

# Yellow = largest negative difference (most reduction), purple = largest positive difference.
plot_grid_diff_heatmap <- function(grid_diffs, outcome_var, spec = grid_outcome_spec) {
  s <- spec_row(spec, outcome_var)
  grid_diffs |>
    filter(outcome == outcome_var) |>
    ggplot(aes(test_sens, prop_I_suspected, fill = diff)) +
    geom_tile() +
    facet_wrap(~ comparison) +
    scale_fill_viridis_c(
      "Difference (negative = fewer)",
      labels    = fmt_labels(s$fmt, diff = TRUE),
      direction = -1,
      guide     = guide_colourbar(barwidth = unit(6, "cm"), barheight = unit(0.3, "cm"))
    ) +
    scale_x_continuous("Test sensitivity", labels = scales::percent, expand = c(0, 0)) +
    scale_y_continuous(coverage_axis_label_wrapped, labels = scales::percent, expand = c(0, 0)) +
    coord_fixed() +
    labs(title = paste0(s$label, ": differences between isolation policies")) +
    grid_theme
}

## h) Outcomes avoided by improving test sensitivity ####

# How the change is reported for each outcome. direction = "avoided": value at sens_from minus value
# at sens_to (positive = fewer with the better test). direction = "additional": value at sens_to
# minus value at sens_from (positive = more with the better test), used for side-room bed-days,
# which rise as more cases are confirmed and isolated.
sensitivity_gain_spec <- grid_outcome_spec |>
  mutate(
    direction = if_else(variable == "side_room_bed_days", "additional", "avoided"),
    title     = case_when(
      variable == "side_room_bed_days" ~ "Additional side-room patient bed-days",
      variable == "cum_total_test"     ~ "Faecal specimen tests avoided",
      TRUE                             ~ paste(label, "avoided")
    )
  )

# For each outcome, isolation policy and testing coverage: values at sens_from and sens_to (annual,
# from the grid_long values), and the change reported as set in sensitivity_gain_spec (column gain).
make_sensitivity_gain <- function(grid_long, sens_from, sens_to, spec = sensitivity_gain_spec) {
  available <- unique(round(grid_long$test_sens, 6))
  missing   <- setdiff(round(c(sens_from, sens_to), 6), available)
  if (length(missing) > 0) {
    stop("Test sensitivity value(s) not in the scenario grid: ", paste(missing, collapse = ", "), ".")
  }
  grid_long |>
    mutate(test_sens = round(test_sens, 6)) |>
    filter(test_sens %in% round(c(sens_from, sens_to), 6)) |>
    mutate(sens = if_else(test_sens == round(sens_from, 6), "value_from", "value_to")) |>
    select(outcome, isolation_policy, prop_I_suspected, sens, value) |>
    pivot_wider(names_from = sens, values_from = value) |>
    left_join(spec |> select(variable, direction), by = c("outcome" = "variable")) |>
    mutate(gain = if_else(direction == "additional", value_to - value_from, value_from - value_to)) |>
    arrange(outcome, isolation_policy, prop_I_suspected)
}

# Title for one outcome, e.g. "CDI cases avoided per year by increasing test sensitivity from 50% to
# 80%"
sensitivity_gain_title <- function(outcome_var, sens_from, sens_to, spec = sensitivity_gain_spec) {
  paste0(spec$title[spec$variable == outcome_var], " per year by increasing test sensitivity from ",
         scales::percent(sens_from), " to ", scales::percent(sens_to))
}

# Table for one outcome: testing coverage (rows) by isolation policy (columns)
make_sensitivity_gain_table <- function(gain, outcome_var, coverage = grid_line_coverage) {
  gain |>
    filter(outcome == outcome_var) |>
    filter_line_coverage(coverage) |>
    mutate(coverage = scales::percent(prop_I_suspected),
           policy   = isolation_policy_labels[as.character(isolation_policy)],
           gain     = scales::comma(round(gain))) |>
    select(coverage, policy, gain) |>
    pivot_wider(names_from = policy, values_from = gain) |>
    rename(!!coverage_axis_label := coverage)
}

# Chart for one outcome: change against testing coverage, one line per isolation policy
plot_sensitivity_gain <- function(gain, outcome_var, sens_from, sens_to, spec = sensitivity_gain_spec) {
  y_lab <- if (spec$direction[spec$variable == outcome_var] == "additional") {
    "Additional per year (negative = decrease)"
  } else {
    "Avoided per year (negative = increase)"
  }
  gain |>
    filter(outcome == outcome_var) |>
    ggplot(aes(prop_I_suspected, gain, colour = isolation_policy)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.2) +
    scale_colour_viridis_d("Isolation policy", labels = \(x) isolation_policy_labels[x], end = 0.8) +
    scale_x_continuous(coverage_axis_label, labels = scales::percent, breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(y_lab, labels = scales::comma) +
    labs(title = sensitivity_gain_title(outcome_var, sens_from, sens_to, spec)) +
    grid_theme +
    guides(colour = guide_legend(ncol = 1))
}

## i) Rates per 10,000 patient bed-days, side-room use and missed cases (section 12 tabs) ####

# Rates per 10,000 patient bed-days use cases and bed-days over the same analysis period, so they
# are the same whether calculated per year or over the whole period.
#
# Missed cases by origin: all CDI cases enter the same I compartment and then follow the same
# rates, whatever their origin, so the model cannot track which missed cases were admitted with CDI,
# hospital-onset or hospital-acquired. Missed cases are apportioned by each origin's share of CDI
# cases; the proportion missed is therefore the same for every origin.
rate_spec <- tribble(
  ~variable,          ~label,                                                         ~fmt,
  "rate_cdi",         "CDI cases per 10,000 patient bed-days",                        "rate",
  "rate_ho",          "Hospital-onset CDI cases per 10,000 patient bed-days",         "rate",
  "rate_ha",          "Hospital-acquired CDI cases per 10,000 patient bed-days",      "rate",
  "rate_tests",       "Faecal specimens tested per 10,000 patient bed-days",          "rate",
  "prop_side_room",   "Patient bed-days in a side room (%)",                          "pct",
  "prop_missed",      "CDI cases missed (%)",                                         "pct",
  "missed_admitted",  "Missed CDI cases per year: admitted with CDI",                 "n",
  "missed_ho",        "Missed CDI cases per year: hospital-onset",                    "n",
  "missed_ha",        "Missed CDI cases per year: hospital-acquired (part of hospital-onset)", "n"
)

# Adds the rate_spec columns to grid_results. years = length of the analysis period, used for the
# missed-case counts per year.
make_grid_rates <- function(grid_results, years) {
  per_10k <- function(n, bed_days) if_else(bed_days > 0, 1e4 * n / bed_days, NA_real_)
  share   <- function(n, d) if_else(d > 0, n / d, NA_real_)
  grid_results |>
    mutate(
      rate_cdi        = per_10k(cum_total_infected_incidence, hospital_bed_days),
      rate_ho         = per_10k(cum_infected_in_hospital,     hospital_bed_days),
      rate_ha         = per_10k(cum_hospital_attributable_I,  hospital_bed_days),
      rate_tests      = per_10k(cum_total_test,               hospital_bed_days),
      prop_side_room  = share(side_room_bed_days, hospital_bed_days),
      prop_missed     = share(cum_total_missed_cases, cum_total_infected_incidence),
      missed_admitted = cum_total_missed_cases * share(cum_admitted_infected,       cum_total_infected_incidence) / years,
      missed_ho       = cum_total_missed_cases * share(cum_infected_in_hospital,    cum_total_infected_incidence) / years,
      missed_ha       = cum_total_missed_cases * share(cum_hospital_attributable_I, cum_total_infected_incidence) / years
    )
}

## j) Example scenario facets (section 11) ####

# Adds `policy` (P1-P3 labels) and `setting` (coverage and sensitivity) columns to a data frame with
# a scenario column, parsed from scenario labels such as "P2_confirmed_cov050_sens100". `settings` is
# a data frame of prop_I_suspected and test_sens giving the order of the settings (by default the
# example_settings defined in the QMD; otherwise sorted by coverage then sensitivity).
example_setting_label <- function(cov, sens) {
  paste0("Coverage ", round(100 * cov), "%\nsensitivity ", round(100 * sens), "%")
}

add_scenario_facets <- function(df,
                                settings = get0("example_settings", envir = globalenv(),
                                                ifnotfound = NULL)) {
  parts <- regmatches(df$scenario, regexec("^(P\\d)_.+_cov(\\d{3})_sens(\\d{3})$", df$scenario))
  code  <- vapply(parts, \(x) if (length(x)) x[2] else NA_character_, "")
  cov   <- vapply(parts, \(x) if (length(x)) as.numeric(x[3]) / 100 else NA_real_, 0)
  sens  <- vapply(parts, \(x) if (length(x)) as.numeric(x[4]) / 100 else NA_real_, 0)
  if (is.null(settings)) {
    settings <- distinct(tibble(prop_I_suspected = cov, test_sens = sens)) |>
      arrange(prop_I_suspected, test_sens)
  }
  setting_levels <- example_setting_label(settings$prop_I_suspected, settings$test_sens)
  policy <- names(isolation_policy_codes)[match(code, sub("_.*", "", isolation_policy_codes))]
  df |>
    mutate(
      policy      = factor(isolation_policy_labels[policy], levels = isolation_policy_labels[isolation_policy_levels]),
      policy_code = factor(code, levels = sub("_.*", "", isolation_policy_codes[isolation_policy_levels])),
      setting     = factor(example_setting_label(cov, sens), levels = setting_levels)
    )
}

# Facet layout for the example time-series plots: settings as rows, policies as columns
example_facets <- facet_grid(setting ~ policy)

# Stacked bar plots of example outcomes: one stacked bar per policy (P1, P2, P3 side by side) within
# each setting, with the setting label beneath each group of three bars. `components` is a named
# vector of columns to stack (names are the legend labels, in stacking order from the bottom);
# `colours` an optional vector of fill colours in the same order.
plot_example_stacked <- function(df, components, title, y_label, colours = NULL, caption = NULL) {
  p <- df |>
    add_scenario_facets() |>
    select(policy_code, setting, all_of(unname(components))) |>
    pivot_longer(-c(policy_code, setting), names_to = "component", values_to = "value") |>
    mutate(component = factor(names(components)[match(component, components)],
                              levels = rev(names(components)))) |>
    ggplot(aes(policy_code, value, fill = component)) +
    geom_col() +
    facet_grid(~ setting, switch = "x") +
    scale_y_continuous(labels = scales::comma) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(strip.placement = "outside", strip.background.x = element_blank()) +
    labs(title = title, x = NULL, y = y_label, fill = NULL,
         caption = paste(c("P1: no isolation; P2: confirmed cases isolated; P3: suspected and confirmed cases isolated.",
                           caption), collapse = " "))
  if (is.null(colours)) {
    p + scale_fill_viridis_d(end = 0.85, breaks = names(components))
  } else {
    p + scale_fill_manual(values = set_names(colours, names(components)), breaks = names(components))
  }
}

# ================================================================================================ #
# End of script ####
