# Helper function to convert probability of event over time t to a daily rate
prob2rate <- function(p, t) {
  -(1 / t) * log(1 - p)
}

# ---- Fixed assumptions ----
# Transmission (per day)
beta0 <- 0.00020
betaC <- 0.00035
betaI <- 0.00060
betaX <- 0.00010

# Natural history
# Colonisation -> infection progression

alpha0 <- 0.010      # baseline progression RATE per day (if this is a probability, see note below)
abx_mult <- 1.30
abx_prop <- 0.20

# If alpha0 is a daily RATE, use:
alpha <- alpha0 * ((1 - abx_prop) + (abx_prop * abx_mult))

# If instead alpha0 is a PROBABILITY over t_progress_C_to_I days, use this instead:
# t_progress_C_to_I <- 10
# alpha_base <- prob2rate(alpha0, t_progress_C_to_I)
# alpha <- alpha_base * ((1 - abx_prop) + (abx_prop * abx_mult))

# Natural decolonisation C -> S
mu0 <- 0.004  # baseline decolonisation RATE per day (if probability over period, use prob2rate)
mu <- mu0

# If mu0 is a probability over 30 days, do:
# t_recover_C_to_S <- 30
# mu <- prob2rate(mu0, t_recover_C_to_S)

# Length of stay (days)
losS <- 7
losC <- 10
losI <- 12
losX <- 6

# Discharge rates (per day)
dS <- 1 / losS
dC <- 1 / losC
dI <- 1 / losI
dX <- 1 / losX

# Mortality rates (per day)
mS <- prob2rate(0.02, losS)
mC <- prob2rate(0.02, losC)
mI <- prob2rate(0.20, losI)  # infected
mX <- prob2rate(0.20, losX)  # isolated

# Admission prevalence
prop_adm_C <- 0.03
prop_adm_I <- 0.002

# Suspected-case identification
p_I_suspected <- 0.85
t_identify_isolate_I_to_X <- 1
gamma <- p_I_suspected / t_identify_isolate_I_to_X
p_not_suspected <- 1 - p_I_suspected