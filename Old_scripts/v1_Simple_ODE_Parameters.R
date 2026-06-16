# Simple ODE model parameters

# Fixed values
births <- 0        # set >0 if needed
mu <- 0            # background mortality
mu_i <- 0.001      # infection mortality
gamma <- 1/7       # recovery rate (7-day infectious duration)

# Costs
c_treat_day <- 25
c_test <- 8

make_mod_parms <- function(
  R0 = 1.8,
  infectious_days = 7,
  tests_per_day = 0
) {
  gamma_local <- 1 / infectious_days
  beta_local <- R0 * gamma_local

  list(
    beta = beta_local,
    gamma = gamma_local,
    mu = mu,
    mu_i = mu_i,
    births = births,
    c_treat_day = c_treat_day,
    c_test = c_test,
    tests_per_day = tests_per_day
  )
}

# Example scenarios
mod.parms.base <- make_mod_parms(R0 = 1.8, infectious_days = 7, tests_per_day = 0)
mod.parms.int1 <- make_mod_parms(R0 = 1.4, infectious_days = 7, tests_per_day = 50)  # intervention lowers transmission
mod.parms.int2 <- make_mod_parms(R0 = 1.8, infectious_days = 5, tests_per_day = 80)  # intervention speeds recovery