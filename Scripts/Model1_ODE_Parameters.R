# Model 1 Parameters ####

# 1) Introduction ####
# This script defines the parameters for the ODE model of hospital transmission of a pathogen
# with a focus on the impact of identification of suspected cases, test turnaround time and 
# test/algorithm sensitivity on outcomes.
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
  if (t_patient_transfer <= 0 || t_identify_suspected <= 0 || t_wait_retest <= 0)
    stop("Time parameters must be > 0.")
  
  # ---- Must be >= 0 ----
  if (beta0 < 0 || betaC < 0 || betaI < 0 || betaXtreat < 0)
    stop("Transmission parameters (beta0, betaC, betaI, betaXtreat) must be >= 0.")
  if (mortI_mult < 0 || mortXtreat_mult < 0)
    stop("mort*_mult values must be >= 0.")
  
  # ---- Must be in [0,1] ----
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
  
  message("All fixed parameters passed validation.")
}

## d) Function to validate scenario control parameters ####
validate_scenario_inputs <- function(isolate_before_test,
                                     prop_I_suspected,
                                     test_sens,
                                     t_test_turnaround,
                                     scenario_label = NA_character_,
                                     verbose = TRUE) {
  nm <- if (!is.na(scenario_label))
    paste0("[", scenario_label, "] ")
  else
    ""
  
  # ---- Must be TRUE/FALSE ----
  if (!is.logical(isolate_before_test) ||
      length(isolate_before_test) != 1 ||
      is.na(isolate_before_test)) {
    stop(nm, "isolate_before_test must be TRUE or FALSE.")
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
  
  # ---- Must be > 0 ----
  if (!is.numeric(t_test_turnaround) ||
      length(t_test_turnaround) != 1 || is.na(t_test_turnaround) ||
      t_test_turnaround <= 0) {
    stop(nm, "t_test_turnaround must be a single numeric value > 0.")
  }
  
  # ---- Warnings ----
  if (t_test_turnaround < 1)
    warning(nm, "t_test_turnaround < 1 ; check units/realism.")
  
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
  beta0 = 0.005,
  betaC = 0.008,
  betaI = 0.040,
  betaXtreat = 0.004,
  
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
  # Note: losFN governs the discharge rate (disFN = 1 / (losFN + t_patient_transfer)), not the
  # actual mean time spent in the FN state. The true mean sojourn in FN is shorter when patients
  # leave via retesting (sigma) before being discharged, since sigma and disFN are competing exits.
  losFN = 14,
  # mean time to recovery and return to general ward for treated cases
  losXtreat = 14, 
  
  # Patient transfer time (days)
  t_patient_transfer = 0.5,
  
  # Time from symptom onset to case identification
  t_identify_suspected = 2.0,
  
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
    isolate_before_test,
    prop_I_suspected,
    test_sens,
    t_test_turnaround,
    scenario_label,
    
    # ---- Fixed parameters (NULL = use fixed list; supply a value to override) ----
    adm_rate               = NULL,
    prop_adm_C             = NULL,
    prop_adm_I             = NULL,
    beta0                  = NULL,
    betaC                  = NULL,
    betaI                  = NULL,
    betaXtreat             = NULL,
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
    t_wait_retest = t_wait_retest,
    prop_Xtreat_to_S = prop_Xtreat_to_S,
    prob_mort_S = prob_mort_S,
    mortI_mult = mortI_mult,
    mortXtreat_mult = mortXtreat_mult
  )
  
  # c) Validate scenario-specific parameters ####
  validate_scenario_inputs(
    isolate_before_test = isolate_before_test,
    prop_I_suspected    = prop_I_suspected,
    test_sens           = test_sens,
    t_test_turnaround   = t_test_turnaround,
    scenario_label      = scenario_label
  )
  
  # d) Calculate derived rate parameters ####
  
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
  
  # Depends on the proportion recovering and returning to S divided by the mean time to recovery and 
  # return to a general ward
  delta <- prop_Xtreat_to_S / (losXtreat + t_patient_transfer)
  
  # Discharges
  
  # Discharges occur at a constant rate, the inverse of the length of stay plus patient transfer time,
  # which represents the average time until discharge.
  disS      <- 1 / (losS + t_patient_transfer)
  disC      <- 1 / (losC + t_patient_transfer)
  disI      <- 1 / (losI + t_patient_transfer)
  disFN     <- 1 / (losFN + t_patient_transfer)
  disXtreat <- (1 - prop_Xtreat_to_S) / (losXtreat + t_patient_transfer)
  
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
  
  # Transmission parameter for suspected cases (Xtest) 
  
  # Depends on whether they are isolated before tested. 
  # If they are isolated, we assume the same transmission parameter as treated cases (betaXtreat). 
  # If not, we assume the same transmission parameter as infected cases (betaI).
  betaXtest <- if (isolate_before_test) betaXtreat else betaI
  
  # Rate of identification of suspected cases (I -> Xtest) 
  
  # Depends on whether they are isolated before tested.
  # Calculated as the proportion of I that are suspected divided by the mean time to identification.
  # If isolation occurs before testing, we assume that the time to identification includes both the 
  # time to identify as suspected and the time to transfer to isolation, since they would be 
  # isolated as soon as they are identified as suspected. 
  # If isolation does not occur before testing, then the time to identification is just the time
  # to identify as suspected.
  if (isolate_before_test) {
    gamma <- prop_I_suspected * (1 / (t_identify_suspected + t_patient_transfer))
  } else {
    gamma <- prop_I_suspected * (1 / t_identify_suspected)
  }
  
  # Rate of confirmation of true positives (Xtest -> Xtreat)
  if (isolate_before_test) {
    theta <- test_sens / t_test_turnaround
  } else {
    theta <- test_sens / (t_test_turnaround + t_patient_transfer)
  }
  
  # Rate of false negatives (Xtest -> FN)
  
  if (isolate_before_test) {
    pi <- (1 - test_sens) / (t_test_turnaround + t_patient_transfer)
  } else {
    pi <- (1 - test_sens) / t_test_turnaround
  }
  
  # Rate of retesting after false negative (FN -> Xtest)
  
  # Assumed to be the inverse of the mean time to retest after a false negative, plus patient transfer 
  # time if isolation occurs before testing.
  if (isolate_before_test) {
    sigma <- 1 / (t_wait_retest + t_patient_transfer)
  } else {
    sigma <- 1 / t_wait_retest
  }
 
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
    isolate_before_test  = isolate_before_test,
    prop_I_suspected     = prop_I_suspected,
    test_sens            = test_sens,
    t_test_turnaround    = t_test_turnaround,
    
    # Allocation of admissions to C, I and S
    adm_rate = adm_rate,
    prop_adm_C = prop_adm_C,
    prop_adm_I = prop_adm_I,
    
    # Transmission
    beta0      = beta0,
    betaC      = betaC,
    betaI      = betaI,
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
    