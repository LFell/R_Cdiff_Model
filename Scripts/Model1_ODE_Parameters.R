# Model 1 Parameters ####

# 1) Introduction ####
# This script defines the parameters for the ODE model of hospital transmission of a pathogen
# with a focus on the impact of isolation policy, identification of suspected cases (testing
# coverage) and test/algorithm sensitivity on outcomes. Test turnaround time is a fixed parameter.
#
# The script is called by the main model script (Model1_ODE_Model_Run.R) to create a list of 
# parameters for each scenario.
#
# IMPORTANT:
# Need to check evidence sources for parameter values, in particular to confirm the specification of
# rates. E.g., this coding assumes progression is reported as 'x% of colonised progress to infected,
# and the mean time to progression is y days', so we can calculate the rate as
# (proportion progressing) / (mean time to progression).
#
# If progression is reported as 'x% of colonised progress to infected per day', then this is already
# a rate and we can use it directly.
#
# If progression is reported as 'x% of colonised become infected within y days', then this is a
# cumulative risk and we need to convert to a rate using the prob2rate function.


# =========================================================== #

# 2) Initial setup ####

## a) Scientific notation ####
options(scipen = 1000)

## b) Rate function ####
prob2rate <- function(p, t) {
  if (!is.numeric(p) || length(p) != 1 || is.na(p) || p < 0 || p >= 1) {
    stop("prob2rate(): 'p' must be a single numeric value in [0, 1).")
  }
  if (!is.numeric(t) || length(t) != 1 || is.na(t) || t <= 0) {
    stop("prob2rate(): 't' must be a single numeric value > 0.")
  }
  - (1 / t) * log(1 - p)
}

## c) Function to validate fixed parameters ####
validate_fixed_parameters <- function(adm_rate,
                                      prop_adm_C,
                                      prop_adm_I,
                                      beta0,
                                      betaC,
                                      betaI,
                                      betaXtreat,
                                      beta_mult_isolated,
                                      prop_progress_C_to_I,
                                      abx_prop,
                                      abx_mult,
                                      t_progress_C_to_I,
                                      prop_decolonise_C_to_S,
                                      t_decolonise_C_to_S,
                                      losS,
                                      losC,
                                      losI,
                                      losFN,
                                      losXtreat,
                                      t_patient_transfer,
                                      t_identify_suspected,
                                      t_test_turnaround,
                                      t_wait_retest,
                                      prop_Xtreat_to_S,
                                      prob_mort_S,
                                      mortI_mult,
                                      mortXtreat_mult) {
  # ---- Must be > 0 ----
  if (adm_rate <= 0)
    stop("Admission rate (adm_rate) must be > 0.")
  if (abx_mult <= 0)
    stop("abx_mult must be > 0.")
  if (t_progress_C_to_I <= 0 || t_decolonise_C_to_S <= 0)
    stop("Progression/decolonisation times must be > 0.")
  if (losS <= 0 || losC <= 0 || losI <= 0 || losFN <= 0 || losXtreat <= 0)
    stop("All LOS parameters must be > 0.")
  if (t_patient_transfer <= 0 || t_identify_suspected <= 0 || t_test_turnaround <= 0 ||
      t_wait_retest <= 0)
    stop("Time parameters must be > 0.")
  
  # ---- Must be >= 0 ----
  if (beta0 < 0 || betaC < 0 || betaI < 0 || betaXtreat < 0)
    stop("Transmission parameters (beta0, betaC, betaI, betaXtreat) must be >= 0.")
  if (mortI_mult < 0 || mortXtreat_mult < 0)
    stop("mort*_mult values must be >= 0.")
  
  # ---- Must be in [0,1] ----
  if (beta_mult_isolated < 0 || beta_mult_isolated > 1)
    stop("beta_mult_isolated must be in [0,1] (isolation should not increase transmission).")
  if (prop_adm_C < 0 || prop_adm_I < 0 || (prop_adm_C + prop_adm_I) > 1)
    stop("Invalid admission prevalence: require prop_adm_C >= 0, prop_adm_I >= 0, and prop_adm_C + prop_adm_I <= 1.")
  if (prop_progress_C_to_I < 0 || prop_progress_C_to_I > 1)
    stop("prop_progress_C_to_I must be in [0,1].")
  if (abx_prop < 0 || abx_prop > 1)
    stop("abx_prop must be in [0,1].")
  if (prop_decolonise_C_to_S < 0 || prop_decolonise_C_to_S > 1)
    stop("prop_decolonise_C_to_S must be in [0,1].")
  if (prop_Xtreat_to_S < 0 || prop_Xtreat_to_S > 1)
    stop("prop_Xtreat_to_S must be in [0,1].")
  if (prob_mort_S < 0 || prob_mort_S >= 1)
    stop("prob_mort_S must be in [0,1).")
  
  # ---- Warnings (plausibility checks) ----
  if ((prop_adm_C + prop_adm_I) > 0.5)
    warning("High imported prevalence (prop_adm_C + prop_adm_I > 0.5). Check assumptions.")
  if (losS < 1 || losC < 1 || losI < 1 || losFN < 1 || losXtreat < 1)
    warning("Some LOS values are < 1 day. Check units.")
  if (t_test_turnaround < 1)
    warning("t_test_turnaround < 1 ; check units/realism.")

  message("All fixed parameters passed validation.")
}

## d) Function to validate scenario control parameters ####
validate_scenario_inputs <- function(isolation_policy,
                                     prop_I_suspected,
                                     test_sens,
                                     scenario_label = NA_character_,
                                     verbose = TRUE) {
  nm <- if (!is.na(scenario_label))
    paste0("[", scenario_label, "] ")
  else
    ""
  
  # ---- Must be one of the allowed isolation policies ----
  isolation_levels <- c("none", "confirmed", "suspected_and_confirmed")
  if (!is.character(isolation_policy) ||
      length(isolation_policy) != 1 ||
      !isolation_policy %in% isolation_levels) {
    stop(nm, "isolation_policy must be one of: ",
         paste(isolation_levels, collapse = ", "), ".")
  }
  
  # ---- Must be in [0,1] ----
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

  message(
    "All scenario-specific parameters have been defined, are positive and proportions are between 0 and 1."
  )
}
## e) Function to validate rate parameters  ####
validate_rate_parameters <- function(alpha,
                                     mu,
                                     gamma,
                                     pi,
                                     theta,
                                     sigma,
                                     delta,
                                     disS,
                                     disC,
                                     disI,
                                     disFN,
                                     disXtreat,
                                     mortS,
                                     mortC,
                                     mortI,
                                     mortXtest,
                                     mortFN,
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
    pi,
    theta,
    sigma,
    delta,
    disS,
    disC,
    disI,
    disFN,
    disXtreat,
    mortS,
    mortC,
    mortI,
    mortXtest,
    mortFN,
    mortXtreat
  )
  
  nms  <- c(
    "alpha",
    "mu",
    "gamma",
    "pi",
    "theta",
    "sigma",
    "delta",
    "disS",
    "disC",
    "disI",
    "disFN",
    "disXtreat",
    "mortS",
    "mortC",
    "mortI",
    "mortXtest",
    "mortFN",
    "mortXtreat"
  )
  
  # ---- Must be finite ----
  if (any(!is.finite(vals))) {
    bad <- paste(nms[!is.finite(vals)], collapse = ", ")
    stop(nm, "Non-finite rate parameter(s): ", bad, ".")
  }
  
  # ---- Must be >= 0 ----
  if (any(vals < 0)) {
    bad <- paste(nms[vals < 0], collapse = ", ")
    stop(nm, "Negative rate parameter(s): ", bad, ".")
  }
  
  # ---- Warnings ----
  if (theta + pi > 1)
    warning(nm, "theta + pi > 1; transitions from Xtest may be very fast.")
  
  message("All rate parameters have been defined, are finite and positive.")
}

# =========================================================== #

# 3) Fixed assumptions ####
# Single source of truth for fixed parameters.
fixed <- list(
  # Admissions
  adm_rate = 133,
  prop_adm_C = 0.086,
  prop_adm_I = 0.002,
  
  # Transmission parameters
  beta0 = 0.002,
  betaC = 0.008,
  betaI = 0.040,
  # Transmission from confirmed cases on treatment (Xtreat) on a general ward. Treatment is assumed
  # to reduce shedding from betaI (0.04) at the start of treatment to 0 at the end; the model uses a
  # constant rate equal to the average across the treatment period (0.02).
  betaXtreat = 0.020,
  # Multiplier on betaI (Xtest) and betaXtreat (Xtreat) for patients isolated in a side room
  beta_mult_isolated = 0.050,
  
  # Progression C -> I
  prop_progress_C_to_I = 0.50,
  abx_prop = 0.373,
  abx_mult = 2.0,
  t_progress_C_to_I = 5.0,
  
  # Decolonisation without treatment (C -> S)
  prop_decolonise_C_to_S = 0.05,
  t_decolonise_C_to_S = 5.0,
  
  # Length of stay (days)
  losS = 7,
  losC = 7,
  # mean LOS is assumed to be higher for infected than susceptible and colonised,
  # as infected patients are symptomatic and less likely to be discharged until symptoms resolve
  losI = 14, 
  # mean LOS for false negatives is assumed to be the same as for infected, as they are also 
  # symptomatic and less likely to be discharged until symptoms resolve.
  # Note: losFN governs the discharge rate (disFN = 1 / losFN), not the
  # actual mean time spent in the FN state. The true mean sojourn in FN is shorter when patients
  # leave via retesting (sigma) before being discharged, since sigma and disFN are competing exits.
  losFN = 14,
  # mean time to recovery and return to general ward for treated cases
  losXtreat = 14, 
  
  # Patient transfer time (days)
  t_patient_transfer = 0.5,
  
  # Time from symptom onset to case identification
  t_identify_suspected = 2.0,

  # Time from test administration to result (days). Fixed across scenarios; scenarios vary
  # isolation_policy, prop_I_suspected and test_sens.
  t_test_turnaround = 2.0,

  # Time from FN result to retest
  # Patients in Xtest are assumed to remain in hospital awaiting their result, so there is no
  # discharge from Xtest; all patients exit via a positive result (theta) or false negative (pi).
  t_wait_retest = 7.0,
  
  # Proportion recover Xtreat -> S
  prop_Xtreat_to_S = 0.80,
  
  # Mortality: cumulative probability of death over LOS for S only
  prob_mort_S = 0.001,
  
  # Daily mortality multipliers (relative to mortS)
  mortI_mult = 2.0,
  mortXtreat_mult = 1.2
)

# Validate fixed parameters once
validate_fixed_parameters(
  adm_rate = fixed$adm_rate,
  prop_adm_C = fixed$prop_adm_C,
  prop_adm_I = fixed$prop_adm_I,
  beta0 = fixed$beta0,
  betaC = fixed$betaC,
  betaI = fixed$betaI,
  betaXtreat = fixed$betaXtreat,
  beta_mult_isolated = fixed$beta_mult_isolated,
  prop_progress_C_to_I = fixed$prop_progress_C_to_I,
  abx_prop = fixed$abx_prop,
  abx_mult = fixed$abx_mult,
  t_progress_C_to_I = fixed$t_progress_C_to_I,
  prop_decolonise_C_to_S = fixed$prop_decolonise_C_to_S,
  t_decolonise_C_to_S = fixed$t_decolonise_C_to_S,
  losS = fixed$losS,
  losC = fixed$losC,
  losI = fixed$losI,
  losFN = fixed$losFN,
  losXtreat = fixed$losXtreat,
  t_patient_transfer = fixed$t_patient_transfer,
  t_identify_suspected = fixed$t_identify_suspected,
  t_test_turnaround = fixed$t_test_turnaround,
  t_wait_retest = fixed$t_wait_retest,
  prop_Xtreat_to_S = fixed$prop_Xtreat_to_S,
  prob_mort_S = fixed$prob_mort_S,
  mortI_mult = fixed$mortI_mult,
  mortXtreat_mult = fixed$mortXtreat_mult
)

# ================================================================================================ #

# 4) Compose a list of parameters for scenarios ####

make_mod_parms <- function(
    # ---- Scenario-specific parameters (must always be supplied) ----
    # isolation_policy: "none", "confirmed" or "suspected_and_confirmed"
    isolation_policy,
    prop_I_suspected,
    test_sens,
    scenario_label,

    # ---- Fixed parameters (NULL = use fixed list; supply a value to override) ----
    adm_rate               = NULL,
    prop_adm_C             = NULL,
    prop_adm_I             = NULL,
    beta0                  = NULL,
    betaC                  = NULL,
    betaI                  = NULL,
    betaXtreat             = NULL,
    beta_mult_isolated    = NULL,
    prop_progress_C_to_I   = NULL,
    abx_prop               = NULL,
    abx_mult               = NULL,
    t_progress_C_to_I      = NULL,
    prop_decolonise_C_to_S = NULL,
    t_decolonise_C_to_S    = NULL,
    losS                   = NULL,
    losC                   = NULL,
    losI                   = NULL,
    losFN                  = NULL,
    losXtreat              = NULL,
    t_patient_transfer     = NULL,
    t_identify_suspected   = NULL,
    t_test_turnaround      = NULL,
    t_wait_retest          = NULL,
    prop_Xtreat_to_S       = NULL,
    prob_mort_S            = NULL,
    mortI_mult             = NULL,
    mortXtreat_mult        = NULL,
    fixed                  = NULL
) {
  # Resolve fixed parameter list: use supplied value, or fall back to calling environment
  if (is.null(fixed)) {
    fixed <- tryCatch(
      get("fixed", envir = parent.frame(), inherits = TRUE),
      error = function(e) stop("'fixed' must be supplied or available in the calling environment.")
    )
  }
  
  # a) Resolve NULL fixed parameters from fixed list ####
  if (is.null(adm_rate))               adm_rate               <- fixed$adm_rate
  if (is.null(prop_adm_C))             prop_adm_C             <- fixed$prop_adm_C
  if (is.null(prop_adm_I))             prop_adm_I             <- fixed$prop_adm_I
  if (is.null(beta0))                  beta0                  <- fixed$beta0
  if (is.null(betaC))                  betaC                  <- fixed$betaC
  if (is.null(betaI))                  betaI                  <- fixed$betaI
  if (is.null(betaXtreat))             betaXtreat             <- fixed$betaXtreat
  if (is.null(beta_mult_isolated))    beta_mult_isolated    <- fixed$beta_mult_isolated
  if (is.null(prop_progress_C_to_I))   prop_progress_C_to_I   <- fixed$prop_progress_C_to_I
  if (is.null(abx_prop))               abx_prop               <- fixed$abx_prop
  if (is.null(abx_mult))               abx_mult               <- fixed$abx_mult
  if (is.null(t_progress_C_to_I))      t_progress_C_to_I      <- fixed$t_progress_C_to_I
  if (is.null(prop_decolonise_C_to_S)) prop_decolonise_C_to_S <- fixed$prop_decolonise_C_to_S
  if (is.null(t_decolonise_C_to_S))    t_decolonise_C_to_S    <- fixed$t_decolonise_C_to_S
  if (is.null(losS))                   losS                   <- fixed$losS
  if (is.null(losC))                   losC                   <- fixed$losC
  if (is.null(losI))                   losI                   <- fixed$losI
  if (is.null(losFN))                  losFN                  <- fixed$losFN
  if (is.null(losXtreat))              losXtreat              <- fixed$losXtreat
  if (is.null(t_patient_transfer))     t_patient_transfer     <- fixed$t_patient_transfer
  if (is.null(t_identify_suspected))   t_identify_suspected   <- fixed$t_identify_suspected
  if (is.null(t_test_turnaround))      t_test_turnaround      <- fixed$t_test_turnaround
  if (is.null(t_wait_retest))          t_wait_retest          <- fixed$t_wait_retest
  if (is.null(prop_Xtreat_to_S))       prop_Xtreat_to_S       <- fixed$prop_Xtreat_to_S
  if (is.null(prob_mort_S))            prob_mort_S            <- fixed$prob_mort_S
  if (is.null(mortI_mult))             mortI_mult             <- fixed$mortI_mult
  if (is.null(mortXtreat_mult))        mortXtreat_mult        <- fixed$mortXtreat_mult
  
  # b) Validate fixed parameters (including overrides) ####
  validate_fixed_parameters(
    adm_rate = adm_rate,
    prop_adm_C = prop_adm_C,
    prop_adm_I = prop_adm_I,
    beta0 = beta0,
    betaC = betaC,
    betaI = betaI,
    betaXtreat = betaXtreat,
    beta_mult_isolated = beta_mult_isolated,
    prop_progress_C_to_I = prop_progress_C_to_I,
    abx_prop = abx_prop,
    abx_mult = abx_mult,
    t_progress_C_to_I = t_progress_C_to_I,
    prop_decolonise_C_to_S = prop_decolonise_C_to_S,
    t_decolonise_C_to_S = t_decolonise_C_to_S,
    losS = losS,
    losC = losC,
    losI = losI,
    losFN = losFN,
    losXtreat = losXtreat,
    t_patient_transfer = t_patient_transfer,
    t_identify_suspected = t_identify_suspected,
    t_test_turnaround = t_test_turnaround,
    t_wait_retest = t_wait_retest,
    prop_Xtreat_to_S = prop_Xtreat_to_S,
    prob_mort_S = prob_mort_S,
    mortI_mult = mortI_mult,
    mortXtreat_mult = mortXtreat_mult
  )
  
  # c) Validate scenario-specific parameters ####
  validate_scenario_inputs(
    isolation_policy    = isolation_policy,
    prop_I_suspected    = prop_I_suspected,
    test_sens           = test_sens,
    scenario_label      = scenario_label
  )
  
  # d) Calculate derived rate parameters ####
  
  # isolation_policy determines which patients are in a side room:
  #   "none"                    - neither suspected (Xtest) nor confirmed (Xtreat) cases
  #   "confirmed"               - confirmed cases (Xtreat) only
  #   "suspected_and_confirmed" - both suspected (Xtest) and confirmed (Xtreat) cases
  # Patient transfer time is added to a rate whenever the transition involves a change of location
  # (general ward <-> side room).
  
  # Transmission parameters for suspected (Xtest) and confirmed (Xtreat) cases
  
  # Suspected cases (Xtest) transmit at betaI on a general ward. Confirmed cases (Xtreat) transmit at
  # the fixed betaXtreat (on treatment) on a general ward. Patients in a side room transmit at the
  # general ward value * beta_mult_isolated. The returned betaXtreat is this policy-adjusted value.
  betaXtest <- switch(isolation_policy,
    none = , confirmed      = betaI,
    suspected_and_confirmed = betaI * beta_mult_isolated)
  
  betaXtreat <- switch(isolation_policy,
    none                                  = betaXtreat,
    confirmed = , suspected_and_confirmed = betaXtreat * beta_mult_isolated)
  
  # Transmission parameter for false negatives 
  
  # Assume the same as for infected (i.e. symptomatic but not isolated)
  betaFN <- betaI
  
  # Progression from colonised to infected (C -> I) 
  
  # Progression is modified by antibiotic use, which enhances progression by a factor of abx_mult in 
  # the proportion of colonised receiving antibiotics (abx_prop).
  alpha <- (prop_progress_C_to_I * ((1 - abx_prop) + (abx_prop * abx_mult)) /
              t_progress_C_to_I)
  
  # Decolonisation from C -> S 
  
  # Decolonisation occurs without treatment is assumed to occur at a constant rate, so we can 
  # calculate the rate as the proportion decolonising divided by the mean time to decolonisation.
  mu <- prop_decolonise_C_to_S / t_decolonise_C_to_S
  
  # Recovery from treatment and return to a general ward (Xtreat -> S) 
  
  # Depends on the proportion recovering and returning to S divided by the mean time to recovery.
  # If confirmed cases are in a side room, the time to patient transfer back to a general ward is
  # added to the length of stay in treatment.
  delta <- switch(isolation_policy,
    none                                  = prop_Xtreat_to_S / losXtreat,
    confirmed = , suspected_and_confirmed = prop_Xtreat_to_S / (losXtreat + t_patient_transfer))
  
  # Discharges
  
  # Discharges occur at a constant rate, the inverse of the length of stay,
  # which represents the average time until discharge.
  disS      <- 1 / losS
  disC      <- 1 / losC
  disI      <- 1 / losI
  disFN     <- 1 / losFN
  disXtreat <- (1 - prop_Xtreat_to_S) / (losXtreat)
  
  # Mortality 
  
  # Mortality rates are calculated from the baseline mortality rate for S, which is derived from the
  # cumulative probability of death over the length of stay for S using the prob2rate function.
  # A multiplier is applied to the baseline mortality rate to calculate the mortality rates for C, I 
  # Xtest, FN and Xtreat, which reflects the increased risk of death associated with infection.
  mortS      <- prob2rate(prob_mort_S, losS)
  mortC      <- mortS
  mortI      <- mortI_mult * mortS
  mortXtest  <- mortI
  mortFN     <- mortI
  mortXtreat <- mortXtreat_mult * mortS
  
  # Rate of identification of suspected cases (I -> Xtest) 
  
  # Calculated as the proportion of I that are suspected divided by the mean time to identification.
  # Note: gamma competes with discharge and death from I, so the share of infected patients who
  # actually reach Xtest is gamma / (gamma + disI + mortI), lower than prop_I_suspected (e.g. about
  # 88% at prop_I_suspected = 1 and 78% at 0.5 under the confirmed and none policies).
  # If suspected cases are isolated, the time to identification includes both the time to identify
  # as suspected and the time to transfer to a side room. Otherwise it is just the time to identify
  # as suspected.
  gamma <- switch(isolation_policy,
    none = , confirmed      = prop_I_suspected / t_identify_suspected,
    suspected_and_confirmed = prop_I_suspected / (t_identify_suspected + t_patient_transfer))
  
  # Rate of confirmation of true positives (Xtest -> Xtreat)
  
  # Transfer time is added only when the patient moves from a general ward (Xtest) to a side room
  # (Xtreat), i.e. when only confirmed cases are isolated.
  theta <- switch(isolation_policy,
    none = , suspected_and_confirmed = test_sens / t_test_turnaround,
    confirmed                        = test_sens / (t_test_turnaround + t_patient_transfer))
  
  # Rate of false negatives (Xtest -> FN)
  
  # False negatives are on a general ward, so transfer time is added if suspected cases are isolated.
  # Note: as transfer time is added to theta (confirmed policy) or pi (suspected_and_confirmed
  # policy) but not both, the share of results that are positive, theta / (theta + pi), differs
  # from test_sens under those policies (e.g. at test_sens = 0.5: 44% confirmed, 56%
  # suspected_and_confirmed, 50% none).
  pi <- switch(isolation_policy,
    none = , confirmed      = (1 - test_sens) / t_test_turnaround,
    suspected_and_confirmed = (1 - test_sens) / (t_test_turnaround + t_patient_transfer))
  
  # Rate of retesting after false negative (FN -> Xtest)
  
  # Inverse of the mean time to retest after a false negative, plus patient transfer time if
  # suspected cases are isolated.
  sigma <- switch(isolation_policy,
    none = , confirmed      = 1 / t_wait_retest,
    suspected_and_confirmed = 1 / (t_wait_retest + t_patient_transfer))
  
  # e) Validate rate parameters ####
  validate_rate_parameters(
    alpha      = alpha,
    mu         = mu,
    gamma      = gamma,
    theta      = theta,
    pi         = pi,
    sigma      = sigma,
    delta      = delta,
    disS       = disS,
    disC       = disC,
    disI       = disI,
    disFN      = disFN,
    disXtreat  = disXtreat,
    mortS      = mortS,
    mortC      = mortC,
    mortI      = mortI,
    mortXtest  = mortXtest,
    mortFN     = mortFN,
    mortXtreat = mortXtreat,
    scenario_label = scenario_label
  )
  
  # f) Return list of parameters for the ODE model ####
  list(
    # Scenario controls
    scenario_label       = scenario_label,
    isolation_policy     = isolation_policy,
    prop_I_suspected     = prop_I_suspected,
    test_sens            = test_sens,

    # Fixed test turnaround time (reported for reference; enters the model via theta and pi)
    t_test_turnaround    = t_test_turnaround,

    # Allocation of admissions to C, I and S
    adm_rate = adm_rate,
    prop_adm_C = prop_adm_C,
    prop_adm_I = prop_adm_I,
    
    # Transmission
    beta0      = beta0,
    betaC      = betaC,
    betaI      = betaI,
    beta_mult_isolated = beta_mult_isolated,
    betaFN     = betaFN,
    betaXtest  = betaXtest,
    betaXtreat = betaXtreat,
    
    # Transition rates between states
    alpha  = alpha,   # C -> I  (progression to infection)
    mu     = mu,      # C -> S  (decolonisation without treatment)
    gamma  = gamma,   # I -> Xtest  (identification of suspected cases)
    theta  = theta,   # Xtest -> Xtreat  (confirmation of true positives)
    pi     = pi,      # Xtest -> FN  (false negatives)
    sigma  = sigma,   # FN -> Xtest  (retesting after false negative)
    delta  = delta,   # Xtreat -> S  (recovery and return to general ward)
    
    # Discharge rates
    disS      = disS,
    disC      = disC,
    disI      = disI,
    disFN     = disFN,
    disXtreat = disXtreat,
    
    # Mortality rates
    mortS      = mortS,
    mortC      = mortC,
    mortI      = mortI,
    mortXtest  = mortXtest,
    mortFN     = mortFN,
    mortXtreat = mortXtreat
  )
}

# =========================================================== #
# END OF SCRIPT ####
