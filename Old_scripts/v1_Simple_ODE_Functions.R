# Simple infectious disease ODE functions

# Compartments: S, I, R
# Optional outcomes: incidence, prevalence, deaths, costs

simple_ode_model <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    N <- S + I + R

    # Force of infection
    lambda <- beta * I / N

    # ODE system
    dS <- births - lambda * S - mu * S
    dI <- lambda * S - gamma * I - mu * I - mu_i * I
    dR <- gamma * I - mu * R

    # Epidemiologic outputs
    new_infections <- lambda * S
    new_recoveries <- gamma * I
    new_deaths_infection <- mu_i * I
    prevalence <- ifelse(N > 0, I / N, 0)

    # Cost outputs (optional)
    daily_treatment_cost <- c_treat_day * I
    daily_testing_cost <- c_test * tests_per_day
    daily_total_cost <- daily_treatment_cost + daily_testing_cost

    list(
      c(dS, dI, dR),
      c(
        N = N,
        lambda = lambda,
        new_infections = new_infections,
        new_recoveries = new_recoveries,
        new_deaths_infection = new_deaths_infection,
        prevalence = prevalence,
        daily_treatment_cost = daily_treatment_cost,
        daily_testing_cost = daily_testing_cost,
        daily_total_cost = daily_total_cost
      )
    )
  })
}

# Helper: probability to rate
prob2rate <- function(p, t) {
  -(1 / t) * log(1 - p)
}