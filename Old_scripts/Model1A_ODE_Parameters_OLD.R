# Model 1A Parameters ####

# 1) Introduction ####

# This script defines the parameters for the ODE model of hospital transmission of a pathogen (e.g., C difficile) with a focus on 
# the impact of identification of suspected cases, test turnaround time and test sensitivity on outcomes.

# There is no parameter for test specificity as it is assumed only patients with the pathogen of interest are tested, so there are no false positives.

# The model structure is specified in the accompanying ODE function script (Model1A_ODE_Functions.R).
# The parameters defined here are used to run the model and generate outputs in the main script (Model1A_ODE_Model_Run.qmd).

# The model compartments are:
# S: Susceptible (not colonised or infected)
# C: Colonised (carrying the pathogen but not infected)
# I: Infected (showing symptoms and can transmit the pathogen)
# Xtest: Isolated awaiting test result (identified as suspected cases and isolated for testing)
# Xtreat: Isolated awaiting treatment (identified as confirmed cases and receiving treatment in isolation)
# D: Dead

# The time steps are in half-days to allow for transition steps that occur within a day (e.g., identification and isolation of infected patients).

# !! Take note ---- 

# Check evidence sources for parameter values and to confirm the specification of rates.
# E.g., this coding assumes progression is reported as 'x% of colonised progress to infected, and the mean time to progression is y half-days', 
# so we can calculate the rate as (proportion progressing) / (mean time to progression).
# If progression is reported as 'x% of colonised progress to infected per half-day', then this is already a rate and we can use it directly.
# If progression is reported as 'x% of colonised become infected within y half-days', then this is a cumulative risk and we need to convert to a rate using the prob2rate function. 

# ========================================================== #

# 2) Initial setup ####

# Use options to avoid scientific notation for small numbers (e.g., transmission rates)
options(scipen = 1000)

# Helper function to convert probability of event over a time interval (cumulative probability or interval risk) into 
# a continuous-time hazard (e.g. daily rate), assuming an exponential distribution of time to event.
prob2rate <- function(p, t) {
  -(1 / t) * log(1 - p)
}
# This means the rate is the negative of the natural log of (1 - probability) divided by the time interval.

# ========================================================== #

# 3) Fixed assumptions ####

# The values of these parameters are fixed across all scenarios unless explicitly changed in the make_mod_parms function.

# Prevalence of C and I on admission to hospital
prop_adm_C <- 0.03 # Proportion of admissions that are colonised
prop_adm_I <- 0.002 # Proportion of admissions that are infected

# Transmission
beta0 <- 0.00020 # Background transmission (rate of effective exposure from the environment per half-day)
betaC <- 0.00035 # Colonised transmission (rate of effective contact between individuals in S and C per half-day)
betaI <- 0.00060 # Infected transmission (rate of effective contact between individuals in S and I per half-day)
betaXtest <- 0.00010 # Isolated transmission (rate of effective contact between individuals in S and Xtest per half-day)
betaXtreat <- 0.00005 # Treated transmission (rate of effective contact between individuals in S and Xtreat per half-day)

# Progression from colonisation to infection (C to I)
prop_progress_C_to_I <- 0.010  # Proportion of colonised that progress from colonised to infected 
abx_mult <- 1.30  # Multiplier for progression due to antibiotics (relative risk)
abx_prop <- 0.20 # Proportion of colonised patients receiving antibiotics 
t_progress_C_to_I <- 10.0 # Mean time from colonisation to infection (half-days)
alpha <- (prop_progress_C_to_I * ((1 - abx_prop) + (abx_prop * abx_mult)) / t_progress_C_to_I) # Rate of progression from colonised to infected (per half-day)

# Decolonisation without treatment (C to S)
prop_decolonise_C_to_S <- 0.05 ## Proportion that decolonise without treatment
t_decolonise_C_to_S <- 5.0 # Mean time required for decolonisation (half-days)
mu <- prop_decolonise_C_to_S / t_decolonise_C_to_S # Rate of recovery from colonised to susceptible without treatment (per half-day)

# Length of stay (half-days)
losS <- 7 # Mean length of stay for susceptible patients
losC <- 10 # Mean length of stay for colonised patients
losI <- 12 # Mean length of stay for infected patients
losXtreat <- 14 # Mean length of stay for confirmed cases in isolation for treatment and recovery before return to S or discharge 
# (Treatment and recovery duration + 72 hours symptom-free observation period + time to move to a general ward or discharge from hospital)

# Recovery of confirmed cases (Xtreat to S and Xtreat to discharge)
prop_Xtreat_to_S <- 0.80  # Proportion of treated confirmed cases that recover to susceptible state and return to a general ward (Xtreat to S)
delta <- prop_Xtreat_to_S / losXtreat # Rate of recovery for confirmed cases and return to a general ward (Xtreat to S, per half-day)

# Discharge rates (per half-day)
disS <- 1 / losS # Discharge rate for susceptible patients (per half-day)
disC <- 1 / losC # Discharge rate for colonised patients (per half-day)
disI <- 1 / losI # Discharge rate for infected patients (per half-day)
disXtreat <- (1 - prop_Xtreat_to_S) / losXtreat # Discharge rate for treated patients in isolation (per half-day)
# Assume no discharge from Xtest as patients are awaiting test results and thus remain in hospital 
# until they are either treated or returned to a general ward.

# Mortality rates (per half-day)
mortS <- prob2rate(0.02, losS) # Mortality rate for susceptible patients (per half-day) 
mortC <- prob2rate(0.02, losC) # Mortality rate for colonised patients (per half-day) - same as susceptible
mortI <- prob2rate(0.20, losI) # Mortality rate for infected patients (per half-day)
mortXtreat <- prob2rate(0.10, losXtreat) # Mortality rate for confirmed cases in isolation (per half-day)
# Length of stay in Xtest is determined by the test turnaround time, so is not fixed across scenarios.

# Safety checks for fixed parameters

stopifnot(
  prop_adm_C >= 0, prop_adm_I >= 0, prop_adm_C + prop_adm_I <= 1,
  abx_prop >= 0, abx_prop <= 1,
  prop_progress_C_to_I >= 0, prop_progress_C_to_I <= 1,
  prop_decolonise_C_to_S >= 0, prop_decolonise_C_to_S <= 1,
  prop_Xtreat_to_S >= 0, prop_Xtreat_to_S <= 1,
  losS > 0, losC > 0, losI > 0, losXtreat > 0,
  t_progress_C_to_I > 0, t_decolonise_C_to_S > 0
)

# ========================================================== #

# 4) Create function to build parameter list for each scenario ####

make_mod_parms <- function(
    prop_I_suspected,
    test_sens,
    test_turnaround,
    scenario_label
) {
  # Validate inputs
  stopifnot(
    length(prop_I_suspected) == 1, prop_I_suspected >= 0, prop_I_suspected <= 1,
    length(test_sens) == 1, test_sens >= 0, test_sens <= 1,
    length(test_turnaround) == 1, test_turnaround > 0
  )
  
  # I -> Xtest (suspected case identification and isolation)
  t_isolate_I_to_Xtest <- 4
  gamma0 <- 1 / t_isolate_I_to_Xtest # Rate of isolation if all infected patients are identified as suspected cases (per half-day)
  gamma <- prop_I_suspected * gamma0 # Rate of isolation if a proportion of infected patients are identified as suspected cases (per half-day)
  
  # Xtest -> Xtreat (true positives) and Xtest -> I (false negatives)
  theta <- test_sens / test_turnaround        # Xtest -> Xtreat (true positive)
  pi <- (1 - test_sens) / test_turnaround     # Xtest -> I (false negative)
  
  # Mortality in Xtest (depends on duration there)
  mortXtest <- prob2rate(0.05, test_turnaround)
  
  # Safety checks after computing rates
  stopifnot(gamma0 >= 0, gamma >= 0, theta >= 0, pi >= 0, delta >= 0)
  stopifnot(disS >= 0, disC >= 0, disI >= 0, disXtreat >= 0)
  stopifnot(mortS >= 0, mortC >= 0, mortI >= 0, mortXtest >= 0, mortXtreat >= 0)
  
  # Create and return parameter list for the model
  list(
    scenario = scenario_label,
    
    # Admissions/transmission/natural history
    prop_adm_C = prop_adm_C, prop_adm_I = prop_adm_I,
    beta0 = beta0, betaC = betaC, betaI = betaI,
    betaXtest = betaXtest, betaXtreat = betaXtreat,
    alpha = alpha, abx_prop = abx_prop, abx_mult = abx_mult, 
    mu = mu,
    
    # Discharge rates
    disS = disS, disC = disC, disI = disI, disXtreat = disXtreat,
    
    # Mortality rates
    mortS = mortS, mortC = mortC, mortI = mortI, mortXtest = mortXtest, mortXtreat = mortXtreat,
    
    # Scenario controls
    prop_I_suspected = prop_I_suspected,
    test_sens = test_sens,
    test_turnaround = test_turnaround,
    
    # Transition rates
    gamma0 = gamma0,     # I -> Xtest if all infected patients are identified as suspected cases
    gamma = gamma,       # I -> Xtest if a proportion of infected patients are identified as suspected cases
    theta = theta,       # Xtest -> Xtreat (true positives)
    pi = pi,             # Xtest -> I (false negatives)
    delta = delta        # Xtreat -> S (recovery and return to general ward)
  )
}

# =========================================================== #

# 5) Build lists of parameter values for scenarios ####

# S1 is base case 
mod.parms.S1 <- make_mod_parms(
  prop_I_suspected = 0.50, test_sens = 0.80, test_turnaround = 5, scenario_label = "S1_midSus_highSe_slow"
)

# S2 is same as S1 but with lower test sensitivity and thus more false negatives, so more patients remain in the I compartment and are not isolated for treatment.
mod.parms.S2 <- make_mod_parms(
  prop_I_suspected = 0.50, test_sens = 0.50, test_turnaround = 5, scenario_label = "S2_midSus_lowSe_slow"
  )

# S3 is same as S1 but with faster test turnaround time, so patients spend less time in the Xtest compartment awaiting results and thus have less opportunity to transmit while awaiting results.
mod.parms.S3 <- make_mod_parms(
  prop_I_suspected = 0.50, test_sens = 0.80, test_turnaround = 3, scenario_label = "S3_midSus_highSe_fast"
)

# S4 is same as S1 but with both lower test sensitivity and faster turnaround time, so patients spend less time in the Xtest compartment awaiting results but more patients remain in the I compartment and are not isolated for treatment.
mod.parms.S4 <- make_mod_parms(
   prop_I_suspected = 0.50, test_sens = 0.50, test_turnaround = 3, scenario_label = "S4_midSus_low_Se_fast"
  
)

# ========================================================== #

# END OF SCRIPT ####