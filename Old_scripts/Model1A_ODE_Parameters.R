# Model 1A Parameters ####

# 1) Introduction ####

# This script defines the parameters for the ODE model of hospital transmission of a pathogen
# with a focus on the impact of identification of suspected cases, test turnaround time and 
# test/algorithm sensitivity on outcomes.

## a) Model structure ####

# There is no parameter for test specificity as it is assumed only patients with the pathogen of
# interest are tested, so there are no false positives.

# The model structure is specified in the accompanying ODE function script (Model1A_ODE_Functions.R).
# The parameters defined here are used to run the model and generate outputs in the main quarto
# document (Model1A_ODE_Model_Run.qmd).

# The model compartments are:
# S: Susceptible (not colonised or infected)
# C: Colonised (carrying the pathogen but not infected)
# I: Infected (showing symptoms caused by the pathogen)
# Xtest: Suspected cases: Isolated or on a general ward, awaiting test result
# Xtreat: Confirmed cases: Isolated and receiving treatment, awaiting recovery & discharge or return 
# to a general ward
# D: Dead (cumulative count of deaths, absorbing state)

## b) Time units ####

# The model time unit is 1 day, which means that all rates are defined per day.
# However, the ODE solver is set to evaluate the derivatives every 0.5 days (day, set in
# Model1A_ODE_Model_Run.qmd) to allow for more granular outputs and to capture transitions that may
# occur within a day (e.g., identification and isolation of infected patients).

# IMPORTANT:
# - All transition parameters in this file (e.g., beta*, alpha, mu, gamma, theta, pi, dis*, mort*)
#   are rates per DAY (i.e., per model time unit), not per solver step.
# - The ODE solver evaluates derivatives at each time point in 'mod.t'; with time_step = 0.5,
#   outputs are produced every half-day.
# - To convert an output rate to a quantity accumulated over one solver step, multiply by time_step
#   in the run script (e.g., flow_per_step = flow_rate * time_step).

## c) Formulae for rates ####

# Need to check evidence sources for parameter values, in particular to confirm the specification of
# rates. E.g., this coding assumes progression is reported as 'x% of colonised progress to infected,
# and the mean time to progression is y days', so we can calculate the rate as
# (proportion progressing) / (mean time to progression).

# If progression is reported as 'x% of colonised progress to infected per day', then this is already
# a rate and we can use it directly.

# If progression is reported as 'x% of colonised become infected within y days', then this is a
# cumulative risk and we need  to convert to a rate using the prob2rate function.

# =========================================================== #

# 2) Initial setup ####

## a) Scientific notation ####
options(scipen = 1000)

## b) Rate function ####

# This function converts a cumulative probability of an event occurring by
# time t (p) into a constant hazard rate (rate) that would yield that cumulative
# probability over the time period t.

# p = 1 - exp(-rate * t)
# Therefore, rate = -log(1 - p) / t

# This is useful for parameters like mortality or progression where we have a cumulative risk
# over a known time period, but the model requires a constant rate.

# The function includes input validation to ensure p is in [0,1] and t is > 0,
# and will throw informative errors if the inputs are invalid.

prob2rate <- function(p, t) {
  if (!is.numeric(p) ||
      length(p) != 1 || is.na(p) || p < 0 || p >= 1) {
    stop("prob2rate(): 'p' must be a single numeric value in [0, 1).")
  }
  if (!is.numeric(t) || length(t) != 1 || is.na(t) || t <= 0) {
    stop("prob2rate(): 't' must be a single numeric value > 0.")
  }
  - (1 / t) * log(1 - p)
}

## c) Function to validate fixed parameters ####

# The validate_fixed_parameters function checks the fixed parameters for validity and plausibility.
# It should be called once after defining the fixed parameters to catch any issues early.

# It checks that admission prevalences are non-negative and sum to <= 1, that proportions are in
# [0,1], that LOS and time parameters are > 0, and it issues warnings for potentially unrealistic
# values (e.g., very high admission prevalence or very short LOS).

# If any critical issues are found, it will stop execution with an informative error message.

# If only warnings are found, it will allow execution to continue but alert the user to check the
# assumptions.

validate_fixed_parameters <- function(prop_adm_C,
                                      prop_adm_I,
                                      beta0,
                                      betaC,
                                      betaI,
                                      betaXtreat,
                                      prop_progress_C_to_I,
                                      abx_prop,
                                      abx_mult,
                                      t_progress_C_to_I,
                                      prop_decolonise_C_to_S,
                                      t_decolonise_C_to_S,
                                      losS,
                                      losC,
                                      losI,
                                      losXtreat,
                                      prop_Xtreat_to_S) {
  # Admission prevalence checks
  if (prop_adm_C < 0 ||
      prop_adm_I < 0 || (prop_adm_C + prop_adm_I) > 1) {
    stop(
      "Invalid admission prevalence: require prop_adm_C >= 0, prop_adm_I >= 0, and prop_adm_C + prop_adm_I <= 1."
    )
  }
  
  # Transmission parameters must be non-negative
  if (beta0 < 0 || betaC < 0 || betaI < 0 || betaXtreat < 0) {
    stop("Transmission parameters (beta0, betaC, betaI, betaXtreat) must be >= 0.")
  }
  
  # Proportions in [0,1]
  if (prop_progress_C_to_I < 0 ||
      prop_progress_C_to_I > 1)
    stop("prop_progress_C_to_I must be in [0,1].")
  if (abx_prop < 0 ||
      abx_prop > 1)
    stop("abx_prop must be in [0,1].")
  if (prop_decolonise_C_to_S < 0 ||
      prop_decolonise_C_to_S > 1)
    stop("prop_decolonise_C_to_S must be in [0,1].")
  if (prop_Xtreat_to_S < 0 ||
      prop_Xtreat_to_S > 1)
    stop("prop_Xtreat_to_S must be in [0,1].")
  
  # Positive multipliers/times/LOS
  if (abx_mult <= 0)
    stop("abx_mult must be > 0.")
  if (t_progress_C_to_I <= 0 ||
      t_decolonise_C_to_S <= 0)
    stop("Progression/decolonisation times must be > 0.")
  if (losS <= 0 ||
      losC <= 0 ||
      losI <= 0 ||
      losXtreat <= 0)
    stop("All LOS parameters must be > 0.")
  
  # Non-fatal plausibility warnings
  if ((prop_adm_C + prop_adm_I) > 0.5) {
    warning("High imported prevalence (prop_adm_C + prop_adm_I > 0.5). Check assumptions.")
  }
  if (losS < 1 || losC < 1 || losI < 1 || losXtreat < 1) {
    warning("Some LOS values are < 1 day. Check units.")
  }
  
  message("All fixed parameters passed validation.")
}

## d) Function to validate scenario-specific inputs ####

# The validate_scenario_inputs() function is called within make_mod_parms() to validate the scenario-specific 
# parameters (isolate_before_test, prop_I_suspected, test_sens, test_turnaround) for each scenario 
# as they are created.

# It checks that the inputs are numeric, of length 1, not NA, and within the expected ranges (e.g., 
# proportions in [0,1], test_turnaround > 0).

# It will stop execution with an informative error message if any of the inputs are invalid, and
# it will issue warnings for potentially unrealistic values (e.g., very low test sensitivity or very 
# short turnaround time).


validate_scenario_inputs <- function(isolate_before_test,
                                     prop_I_suspected,
                                     test_sens,
                                     test_turnaround,
                                     scenario_label = NA_character_,
                                     verbose = TRUE) {
  nm <- if (!is.na(scenario_label))
    paste0("[", scenario_label, "] ")
  else
    ""
  
  if (!is.logical(isolate_before_test) ||
      length(isolate_before_test) != 1 ||
      is.na(isolate_before_test)) {
    stop(nm, "isolate_before_test must be TRUE or FALSE.")
  }
  if (!is.numeric(prop_I_suspected) ||
      length(prop_I_suspected) != 1 || is.na(prop_I_suspected) ||
      prop_I_suspected < 0 || prop_I_suspected > 1) {
    stop(nm, "prop_I_suspected must be a single numeric value in [0,1].")
  }
  if (!is.numeric(test_sens) ||
      length(test_sens) != 1 || is.na(test_sens) ||
      test_sens < 0 || test_sens > 1) {
    stop(nm, "test_sens must be a single numeric value in [0,1].")
  }
  if (!is.numeric(test_turnaround) ||
      length(test_turnaround) != 1 || is.na(test_turnaround) ||
      test_turnaround <= 0) {
    stop(nm, "test_turnaround must be a single numeric value > 0.")
  }
  
  # Non-fatal warnings
  if (test_turnaround < 1)
    warning(nm, "test_turnaround < 1 ; check units/realism.")
  
  # Message to confirm not issues found with scenario-specific parameters
  message(
    "All scenario-specific parameters have been defined, are positive and proportions are between 0 and 1."
  )
}

## e) Function to validate rate parameters  ####

# The validate_rate_parameters() function is called within make_mod_parms() to validate the calculated 
# transition rate parameters for each scenario, ensuring that they are finite and non-negative, and 
# to check for any potential issues with the combination of parameters (e.g., theta + pi > 1).

# It takes all the calculated rate parameters as inputs, along with the fixed rates and the scenario 
# label for informative error messages.

# If any of the rate parameters are non-finite or negative, it will stop execution with an 
# informative error message listing the problematic rate parameters.

# It also checks if theta + pi > 1, which could indicate an issue with the test sensitivity and 
# turnaround time leading to unrealistically fast transitions from Xtest. If any issues are found, 
# it will alert the user to check the scenario parameters and assumptions.

validate_rate_parameters <- function(alpha,
                           mu,
                           gamma,
                           gamma0,
                           pi,
                           theta,
                           delta,
                           disS,
                           disC,
                           disI,
                           disXtreat,
                           mortS,
                           mortC,
                           mortI,
                           mortXtest,
                           mortXtreat,
                           scenario_label = NA_character_,
                           verbose = TRUE) {
  nm <- if (!is.na(scenario_label))
    paste0("[", scenario_label, "] ")
  else
    ""
  
  vals <- c(
    alpha,
    mu,
    gamma,
    gamma0,
    pi,
    theta,
    delta,
    disS,
    disC,
    disI,
    disXtreat,
    mortS,
    mortC,
    mortI,
    mortXtest,
    mortXtreat
  )
  
  nms  <- c(
    "alpha",
    "mu",
    "gamma",
    "gamma0",
    "pi",
    "theta",
    "delta",
    "disS",
    "disC",
    "disI",
    "disXtreat",
    "mortS",
    "mortC",
    "mortI",
    "mortXtest",
    "mortXtreat"
  )
  
  if (any(!is.finite(vals))) {
    bad <- paste(nms[!is.finite(vals)], collapse = ", ")
    stop(nm, "Non-finite rate parameter(s): ", bad, ".")
  }
  if (any(vals < 0)) {
    bad <- paste(nms[vals < 0], collapse = ", ")
    stop(nm, "Negative rate parameter(s): ", bad, ".")
  }
  
  # Non-fatal warnings
  if (theta + pi > 1)
    warning(nm, "theta + pi > 1; transitions from Xtest may be very fast.")
  
  # Message to confirm not issues found with rates
  message("All rate parameters have been defined, are finite and positive.")
}

# =========================================================== #

# 3) Fixed assumptions ####

# Note: these are fixed across all scenarios, so we can define them once here.

## a) Prevalence of C & I on admission ####

# Proportion of admissions that are colonised
prop_adm_C <- 0.11 

# Proportion of admissions that are infected
prop_adm_I <- 0.002 

# The functions script assumes the rest of admissions are to S, so admS <- exits_total - admC - admI.

## b) Transmission parameters ####

# Rates of effective contact with individuals in S, per occupancy of each state per day.

# The transmission parameters are used in the force of colonisation calculation in the ODE function, 
# which determines the flow from S to C.

# Background transmission (rate of effective exposure from the hospital environment)
beta0 <- 0.3 

# Colonised transmission (rate of effective contact between individuals in S and C)
betaC <- 0.50 

# Infected transmission (rate of effective contact between individuals in S and I)
betaI <- 0.80  

# Treated transmission (rate of effective contact between individuals in S and Xtreat)
betaXtreat <- 0.50 

# betaXtest is defined in the make_mod_parms() function as it depends on whether patients are 
# isolated before testing or not.

## c) Progression C -> I ####

# Proportion of colonised that progress from colonised to infected
prop_progress_C_to_I <- 0.50 

# Proportion of colonised patients receiving antibiotics
abx_prop <- 0.50 

# Multiplier for progression due to antibiotics (relative risk)
abx_mult <- 2.0  

# Mean time from colonisation to infection (days)
t_progress_C_to_I <- 5.0 

# Parameter for rate of progression from colonised to infected (per occupancy of C per day)
alpha <- (prop_progress_C_to_I * ((1 - abx_prop) + (abx_prop * abx_mult)) / t_progress_C_to_I) 

## d) Decolonisation without treatment (C -> S) ####

# Proportion that decolonise without treatment
prop_decolonise_C_to_S <- 0.05 

# Mean time required for decolonisation (days)
t_decolonise_C_to_S <- 5.0 

# Parameter for rate of recovery from colonised to susceptible without treatment (per occupancy of 
# C day)
mu <- prop_decolonise_C_to_S / t_decolonise_C_to_S 

## e) Length of stay (days) ####

# Mean length of stay for susceptible patients
losS <- 7 
# Mean length of stay for colonised patients
losC <- 7 
# Mean length of stay for infected patients
losI <- 14 
# Mean length of stay for confirmed cases in isolation for treatment and recovery
losXtreat <- 14

# Length of stay in Xtest is determined by the test turnaround time, so is set in the 
# make_mod_parms() function.

## f) Recovery Xtreat -> S ####

# Proportion of treated cases that recover and are returned to a general ward
prop_Xtreat_to_S <- 0.80  

# Parameter for rate of recovery for treated cases and return to a general ward per occupancy of 
# Xtreat per day
delta <- prop_Xtreat_to_S / losXtreat 

## g) Discharge rate parameters ####

# Per occupancy of each state per day

# Discharge of susceptible patients
disS <- 1 / losS 

# Discharge of colonised patients
disC <- 1 / losC 

# Discharge of infected patients
disI <- 1 / losI 

# Discharge of treated, confirmed cases
disXtreat <- (1 - prop_Xtreat_to_S) / losXtreat 

# Assume no discharge from Xtest as patients are awaiting test results and thus remain in hospital
# until they are either treated or returned to a general ward.

## h) Mortality rate parameters ####

# Per occupancy of each state per day

# Mortality among susceptible patients
mortS <- prob2rate(0.001, losS) 

# Mortality rate for colonised patients
# Assume same as susceptible, as colonisation is asymptomatic and not associated with increased 
# mortality risk.
mortC <- mortS 

# Mortality among infected patients
mortI <- prob2rate(0.010, losI) 

# Mortality among suspected cases in Xtest
# Assume same as infected as they are still infected and awaiting test results, but not yet being 
# treated.
mortXtest <- mortI 

# Mortality among confirmed cases in Xtreat
mortXtreat <- prob2rate(0.005, losXtreat) 

## i) Validate fixed parameters ####

# Use the validate_fixed_parameters here to catch any issues before we start building scenarios.
validate_fixed_parameters(
  prop_adm_C = prop_adm_C,
  prop_adm_I = prop_adm_I,
  beta0 = beta0,
  betaC = betaC,
  betaI = betaI,
  betaXtreat = betaXtreat,
  prop_progress_C_to_I = prop_progress_C_to_I,
  abx_prop = abx_prop,
  abx_mult = abx_mult,
  t_progress_C_to_I = t_progress_C_to_I,
  prop_decolonise_C_to_S = prop_decolonise_C_to_S,
  t_decolonise_C_to_S = t_decolonise_C_to_S,
  losS = losS,
  losC = losC,
  losI = losI,
  losXtreat = losXtreat,
  prop_Xtreat_to_S = prop_Xtreat_to_S
)

# =========================================================== #

# 4) Scenario builder function ####

make_mod_parms <- function(isolate_before_test,
                           prop_I_suspected,
                           test_sens,
                           test_turnaround,
                           scenario_label) {
  
# a) Validate scenario-specific parameters ####
  
  validate_scenario_inputs(
    isolate_before_test = isolate_before_test,
    prop_I_suspected = prop_I_suspected,
    test_sens = test_sens,
    test_turnaround = test_turnaround,
    scenario_label = scenario_label
  )
  
# b) Set betaXtest ####
  
  # If patients are isolated before testing, BetaXtest is the same as BetaXtreat.
  
  # If not, BetaXtest = BetaI, as they are still infected and transmitting at the same rate and not 
  # in isolation.
  
  betaXtest <- if (isolate_before_test)
    betaXtreat
  else
    betaI
  
# c) Calculate gamma (I -> Xtest) ####
  
  # Mean time from symptom onset to suspected case identification and isolation (days)
  t_isolate_I_to_Xtest <- 4
  
  # Parameter for rate of suspected case identification and isolation (per day)
  gamma <- prop_I_suspected * (1 / t_isolate_I_to_Xtest) 
  
  # Parameter for rate of missed untested cases (per day)
  gamma0 <- (1 - prop_I_suspected) * (1 / t_isolate_I_to_Xtest)
  
# d) Calculate theta (Xtest -> Xtreat)####
  
  # Parameter for rate of confirmation of suspected cases (true positives) 
  theta <- test_sens / test_turnaround
  
# e) Calculate pi (Xtest -> I) ####
  
  # Parameter for rate of return to I for false negatives. Assumes that false negatives are returned to a general 
  # ward after the test turnaround time, so the rate is determined by the test sensitivity and 
  # turnaround time.
  pi <- (1 - test_sens) / test_turnaround     
  
# f) Validate rate parameters ####
  
  validate_rate_parameters(
    alpha = alpha,
    mu = mu,
    gamma = gamma,
    gamma0 = gamma0,
    pi = pi,
    theta = theta,
    delta = delta,
    disS = disS,
    disC = disC,
    disI = disI,
    disXtreat = disXtreat,
    mortS = mortS,
    mortC = mortC,
    mortI = mortI,
    mortXtest = mortXtest,
    mortXtreat = mortXtreat,
    scenario_label = scenario_label
  )

# g) Return list of parameters for the ODE model ####
  
  list(
    # Scenario controls 
    scenario_label = scenario_label,
    isolate_before_test = isolate_before_test,
    prop_I_suspected = prop_I_suspected,
    test_sens = test_sens,
    test_turnaround = test_turnaround,
    
    # Allocation of admissions to C, I and S
    prop_adm_C = prop_adm_C,
    prop_adm_I = prop_adm_I,
    
    # Transmission 
    beta0 = beta0,
    betaC = betaC,
    betaI = betaI,
    betaXtest = betaXtest,
    betaXtreat = betaXtreat,
    
    # Transition between states 
    # C -> I (Progression to infection)
    alpha = alpha,
    # C -> S (Decolonisation without treatment)
    mu = mu,
    # I -> Xtest (Identification of suspected cases)
    gamma = gamma,
    # I -> I (Missed cases that are not tested)
    gamma0 = gamma0,
    # Xtest -> I (Return of false negatives)
    pi = pi,
    # Xtest -> Xtreat (Confirmation of true positives)
    theta = theta,
    # Xtreat -> S (Recovery of confirmed cases and return to general ward)
    delta = delta, 
    
    # Discharges 
    disS = disS,
    disC = disC,
    disI = disI,
    # Xtreat -> discharge (Recovery of confirmed cases and discharge from hospital)
    disXtreat = disXtreat,
    
    # Mortality
    mortS = mortS,
    mortC = mortC,
    mortI = mortI,
    mortXtest = mortXtest,
    mortXtreat = mortXtreat
  )
}

# =========================================================== #

# 5) Build lists of parameter values for scenarios ####

## a) S1 ####
mod.parms.S1 <- make_mod_parms(
  isolate_before_test = TRUE,
  prop_I_suspected = 1.00,
  test_sens = 1.00,
  test_turnaround = 2,
  scenario_label = "S1"
)

## b) S2 ####
mod.parms.S2 <- make_mod_parms(
  isolate_before_test = TRUE,
  prop_I_suspected = 1.00,
  test_sens = 0.50,
  test_turnaround = 2,
  scenario_label = "S2"
)

## c) S3 ####
mod.parms.S3 <- make_mod_parms(
  isolate_before_test = TRUE,
  prop_I_suspected = 0.50,
  test_sens = 1.00,
  test_turnaround = 2,
  scenario_label = "S3"
)

## d) S4 ####
mod.parms.S4 <- make_mod_parms(
  isolate_before_test = TRUE,
  prop_I_suspected = 0.50,
  test_sens = 0.50,
  test_turnaround = 2,
  scenario_label = "S4"
)

# =========================================================== #
# END OF SCRIPT ####
