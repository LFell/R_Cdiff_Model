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
    
    # Number of patients in Xtest and/or Xtreat at each time step
    # Include Xtest and Xtreat if isolate_before_test = TRUE
    # Include only Xtreat when isolate_before_test = FALSE
    
    if (isolate_before_test) {
      side_room_bed_occupancy <- Xtest + Xtreat
    } else {
      side_room_bed_occupancy <- Xtreat
    }
    
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
        general_ward_bed_occupancy = general_ward_bed_occupancy,
        
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

add_quantities <- function(df, time_step) {
  df %>%
    mutate(
      # If time is 0, assign year 1; otherwise calculate year based on time/365
      # using ceiling() to bin days into years: (0-365 -> Year 1, 365.5-730 -> Year 2)
      year = if_else(time == 0, 1, ceiling(time / 365)),
      
      # Patient-days per state
      S_days     = S * time_step,
      C_days     = C * time_step,
      I_days     = I * time_step,
      Xtest_days = Xtest * time_step,
      FN_days    = FN * time_step,
      Xtreat_days= Xtreat * time_step,
      general_ward_bed_days = general_ward_bed_occupancy * time_step,
      side_room_bed_days    = side_room_bed_occupancy * time_step,
      hospital_bed_days     = hospital_bed_occupancy * time_step,

      # Per-step flows
      
      # Admissions
      admitted_susceptible_step = admitted_susceptible_rate * time_step,
      admitted_colonised_step = admitted_colonised_rate * time_step,
      admitted_infected_step = admitted_infected_rate * time_step,
      total_admissions_step = total_admissions_rate * time_step,
      
      # Colonisation and decolonisation
      colonised_in_hospital_step = colonised_in_hospital_rate * time_step,
      decolonised_in_hospital_step = decolonised_in_hospital_rate * time_step,
      total_colonised_incidence_step = total_colonised_incidence_rate * time_step,
      
      # Infection
      infected_in_hospital_step = infected_in_hospital_rate * time_step,
      total_infected_incidence_step = total_infected_incidence_rate * time_step,
      
      # Case identification and recovery
      suspected_cases_step = suspected_cases_rate * time_step,
      confirmed_cases_step = confirmed_cases_rate * time_step,
      recovered_cases_step = recovered_cases_rate * time_step,
      
      # Testing
      first_test_step = first_test_rate * time_step,
      retest_false_negatives_step = retest_false_negatives_rate * time_step,
      total_test_step = total_test_rate * time_step,
      false_negatives_step = false_negatives_rate * time_step,
      
      # Missed cases
      missed_untested_cases_step = missed_untested_cases_rate * time_step,
      missed_false_negative_cases_step = missed_false_negative_cases_rate * time_step,
      total_missed_cases_step = total_missed_cases_rate * time_step,
      
      # Discharges
      discharges_susceptible_step = discharges_susceptible_rate * time_step,
      discharges_colonised_step = discharges_colonised_rate * time_step,
      discharges_infected_step = discharges_infected_rate * time_step,
      discharges_FN_step = discharges_FN_rate * time_step,
      discharges_recovered_step = discharges_recovered_rate * time_step,
      discharges_while_infected_step = discharges_while_infected_rate * time_step,
      total_discharges_step = total_discharges_rate * time_step,
      
      # Deaths
      deaths_susceptible_step = deaths_susceptible_rate * time_step,
      deaths_colonised_step = deaths_colonised_rate * time_step,
      deaths_infected_step = deaths_infected_rate * time_step,
      deaths_suspected_cases_step = deaths_suspected_cases_rate * time_step,
      deaths_FN_step = deaths_FN_rate * time_step,
      deaths_confirmed_cases_step = deaths_confirmed_cases_rate * time_step,
      deaths_post_infection_step = deaths_post_infection_rate * time_step,
      total_deaths_step = total_deaths_rate * time_step,
      
      # Exits
      total_exits_step = total_exits_rate * time_step,
      
      # Cumulative flows
      
      # Admissions
      cum_admitted_susceptible = cumsum(admitted_susceptible_step),
      cum_admitted_colonised = cumsum(admitted_colonised_step),
      cum_admitted_infected = cumsum(admitted_infected_step),
      cum_total_admissions = cumsum(total_admissions_step),
      
      # Colonisation and decolonisation
      cum_colonised_in_hospital = cumsum(colonised_in_hospital_step),
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
      cum_false_negatives = cumsum(false_negatives_rate * time_step),
      
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

make_annual_summary <- function(df, scenario_label) {
  df %>%
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
      
      # Admissions
      annual_admitted_susceptible = sum(admitted_susceptible_step, na.rm = TRUE),
      annual_admitted_colonised = sum(admitted_colonised_step, na.rm = TRUE),
      annual_admitted_infected = sum(admitted_infected_step, na.rm = TRUE),
      annual_total_admissions = sum(total_admissions_step, na.rm = TRUE),
      
      # Colonisation and decolonisation
      annual_colonised_in_hospital = sum(colonised_in_hospital_step, na.rm = TRUE),
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
      cum_hospital_attributable_I  = cum_hospital_attributable_I,
      prop_I_hospital_attributable = prop_I_hospital_attributable,
     
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

## b) Final outcomes vs base scenario ####

# Takes the final_outcomes tibble and returns absolute differences vs a chosen
# base scenario for all numeric columns.
#
# Column names are prefixed with "d_" and suffixed with "_vs_{base_scenario}".
# The "final_" and "cum_" prefixes are stripped for brevity:
#   final_S                  -> d_S_vs_S1
#   cum_confirmed_cases      -> d_confirmed_cases_vs_S1
#   cum_side_room_bed_occupancy -> d_side_room_bed_occupancy_vs_S1

make_final_outcomes_vs_base <- function(final_outcomes, base_scenario = "S1") {
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
  }))
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
  }))
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
  }))
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
  "Xtreat_patient_days"     = "Confirmed cases (side-room pbds)",
  
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
    facet_wrap( ~ scenario) +
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
    facet_wrap( ~ scenario) +
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
    facet_wrap( ~ scenario) +
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
    facet_wrap( ~ scenario) +
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

## e) Bed occupancy - line plot ####

plot_bed_occupancy <- function(plot_long_bed_occupancy,
                               longest_stabilisation,
                               colour_palette = model_colour_palette,
                               pretty_labels  = model_pretty_labels) {
    plot_long_bed_occupancy %>%
    filter(time < longest_stabilisation) %>%
    mutate(resource = factor(
      bed_type,
      levels = c("hospital_bed_occupancy",
                 "general_ward_bed_occupancy",
                 "side_room_bed_occupancy"))) %>%
    ggplot(aes(x = time, y = occupancy, color = resource)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~scenario) +
    scale_color_manual(values = colour_palette, labels = pretty_labels) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    labs(
      title = "Bed occupancy rates",
      x     = "Days (model run period)",
      y     = "Bed occupancy per day",
      color = NULL
    )
}
## f) Cumulative bed occupancy - bar plot ####

plot_cumulative_bed_occupancy <- function(patient_days_summary, pretty_labels = model_pretty_labels) {
  patient_days_summary %>%
    select(scenario, 
           hospital_bed_days, 
           general_ward_bed_days,
           side_room_bed_days,
           Xtreat_patient_days) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c(
        "hospital_bed_days",
        "general_ward_bed_days",
        "side_room_bed_days",
        "Xtreat_patient_days"
      )
    )) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      scales = "free_y",
      ncol = 2,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Cumulative bed occupancy by scenario",
      x     = "",
      y     = "Cumulative patient bed-days"
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
      time, scenario,
      colonised_in_hospital_rate, admitted_colonised_rate,
      infected_in_hospital_rate,  admitted_infected_rate
    ) %>%
    pivot_longer(cols = -c(time, scenario), names_to = "flow", values_to = "value") %>%
    mutate(flow = factor(flow, levels = c(
      "colonised_in_hospital_rate", "admitted_colonised_rate",
      "infected_in_hospital_rate",  "admitted_infected_rate"))) %>%
    ggplot(aes(time, value, color = flow, linetype = flow)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ scenario) +
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
      scenario,
      total_infected_incidence_rate,
      suspected_cases_rate,
      confirmed_cases_rate,
      false_negatives_rate,
      retest_false_negatives_rate
      ) %>%
    pivot_longer(
      cols = -c(time, scenario),
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
    facet_wrap( ~ scenario) +
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

## i) Cumulative case incidence and ascertainment - bar plot ####

plot_case_ascertainment_bar <- function(final_outcomes, pretty_labels = model_pretty_labels) {
  final_outcomes %>%
    select(
      scenario,
      cum_total_infected_incidence,
      cum_suspected_cases,
      cum_confirmed_cases,
      cum_false_negatives,
      cum_retest_false_negatives
    ) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c(
        "cum_total_infected_incidence",
        "cum_suspected_cases",
        "cum_confirmed_cases",
        "cum_false_negatives",
        "cum_retest_false_negatives"
      )
    )) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      #scales = "free_y",
      ncol = 3,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Cumulative case incidence and ascertainment by scenario",
      x     = "",
      y     = "Cumulative count"
    )
}

## j) Cumulative testing - bar plot ####

plot_testing_bar <- function(final_outcomes, pretty_labels = model_pretty_labels) {
  final_outcomes %>%
    select(
      scenario,
      cum_first_test,
      cum_retest_false_negatives,
      cum_total_test
    ) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c(
        "cum_first_test",
        "cum_retest_false_negatives",
        "cum_total_test"
      )
    )) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      #scales = "free_y",
      ncol = 3,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Cumulative testing by scenario",
      x     = "",
      y     = "Cumulative count"
    )
}

## m) Cumulative missed cases - bar plot ####

plot_missed_cases_bar <- function(final_outcomes,
                                  pretty_labels = model_pretty_labels) {
  final_outcomes %>%
    select(
      scenario,
      cum_missed_untested_cases,
      cum_missed_false_negative_cases,
      cum_total_missed_cases
    ) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c(
        "cum_missed_untested_cases",
        "cum_missed_false_negative_cases",
        "cum_total_missed_cases"
      )
    )) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      #scales = "free_y",
      ncol = 3,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Breakdown of missed cases scenario",
      x     = "",
      y     = "Cumulative count"
    )
}

## n) Confirmation via first test vs retest - cumulative bar plot ####

plot_case_confirmation_bar <- function(final_outcomes,
                                  pretty_labels = model_pretty_labels) {
  final_outcomes %>%
    select(
      scenario,
      cum_confirmed_via_first_test,
      cum_confirmed_via_retest,
      cum_confirmed_cases
    ) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c(
        "cum_confirmed_via_first_test",
        "cum_confirmed_via_retest",
        "cum_confirmed_cases"
      )
    )) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      #scales = "free_y",
      ncol = 3,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Cumulative confirmed cases via first test vs retest",
      x     = "",
      y     = "Cumulative count"
    )
}

## o) Cumulative deaths by state - bar plot ####

plot_deaths_by_state_bar <- function(final_outcomes, pretty_labels = model_pretty_labels) {
  final_outcomes %>%
    select(
      scenario,
      cum_deaths_susceptible,
      cum_deaths_colonised,
      cum_deaths_infected,
      cum_deaths_suspected_cases,
      cum_deaths_FN,
      cum_deaths_confirmed_cases
    ) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c(
        "cum_deaths_susceptible",
        "cum_deaths_colonised",
        "cum_deaths_infected",
        "cum_deaths_suspected_cases",
        "cum_deaths_FN",
        "cum_deaths_confirmed_cases"
      )
    )) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      scales = "free_y",
      ncol = 3,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Cumulative deaths by state",
      x     = "",
      y     = "Cumulative count"
    )
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
      scenario,
      deaths_colonised_rate,
      deaths_infected_rate,
      deaths_suspected_cases_rate,
      deaths_FN_rate,
      deaths_confirmed_cases_rate
    ) %>%
    pivot_longer(
      cols = -c(time, scenario),
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
    facet_wrap( ~ scenario) +
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

## q) Cumulative discharges by state - bar plot ####

plot_discharges_by_state_bar <- function(final_outcomes, pretty_labels = model_pretty_labels) {
  final_outcomes %>%
    select(
      scenario,
      cum_discharges_susceptible,
      cum_discharges_colonised,
      cum_discharges_infected,
      cum_discharges_FN,
      cum_discharges_recovered
    ) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c(
        "cum_discharges_susceptible",
        "cum_discharges_colonised",
        "cum_discharges_infected",
        "cum_discharges_FN",
        "cum_discharges_recovered"
      )
    )) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      scales = "free_y",
      ncol = 3,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Cumulative discharges by state",
      x     = "",
      y     = "Cumulative count"
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
      scenario,
      discharges_colonised_rate,
      discharges_infected_rate,
      discharges_FN_rate,
      discharges_recovered_rate
    ) %>%
    pivot_longer(
      cols = -c(time, scenario),
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
    facet_wrap( ~ scenario) +
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

## s) Cumulative case outcomes - bar plot ####

plot_case_outcomes_bar <- function(final_outcomes, pretty_labels = model_pretty_labels) {
  final_outcomes %>%
    select(
      scenario,
      cum_discharges_while_infected,
      cum_recovered_cases,
      cum_deaths_post_infection
    ) %>%
    pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
    # Set the order of the metric factor levels to control the order of bars and facets
    mutate(metric = factor(
      metric,
      levels = c(
       "cum_discharges_while_infected",
       "cum_recovered_cases",
       "cum_deaths_post_infection"
      ))) %>%
    ggplot(aes(scenario, value, fill = scenario)) +
    geom_col() +
    facet_wrap(
      ~ metric,
      scales = "free_y",
      ncol = 3,
      strip.position = "top",
      labeller = labeller(metric = as_labeller(pretty_labels))
    ) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Cumulative outcomes by scenario", x     = "", y     = "Count")
}

## t) Control parameter values - scatter plot ####

plot_scenario_params <- function(control_parms) {
  isolate_value <- unique(control_parms$isolate_before_test)
  # If >1 value, this will display all (comma separated), otherwise just TRUE/FALSE.
  iso_text <- paste("Isolate before testing:", paste(isolate_value, collapse = ", "))
  
  ggplot(control_parms, aes(test_sens, prop_I_suspected)) +
    geom_point() +
    geom_text(aes(label = scenario), vjust = -0.7, size = font_base_size / .pt) +
    theme_classic(base_size = font_base_size) +
    my_theme +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Values of control parameters",
      x     = "Test sensitivity",
      y     = "Proportion of infections identified\nas suspected cases"
    ) +
    annotate("text", x = 0.00, y = 0.05, label = iso_text, size = font_base_size / .pt, hjust = 0)
}

# ================================================================================================ #
# End of script ####
