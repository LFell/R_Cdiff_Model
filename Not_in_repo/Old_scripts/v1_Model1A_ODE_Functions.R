Model1A_ODE <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    N <- S + C + I + X
    lambda <- beta0 + (betaC * C + betaI * I + betaX * X) / N
    alpha_eff <- alpha * abx_mult

    dS <- aS + mu * C + delta * X - lambda * S - (dS_rate + mS) * S
    dC <- aC + lambda * S - (alpha_eff + mu + dC_rate + mC) * C
    dI <- aI + alpha_eff * C + pi * X - (gamma + dI_rate + mI) * I
    dX <- gamma * I - (delta + pi + dX_rate + mX) * X

    # Outputs for validation
    total_cases_flow <- alpha_eff * C + aI
    identified_cases_flow <- delta * X + dX_rate * X
    false_neg_flow <- pi * X
    discharged_colonised_flow <- dC_rate * C
    discharged_infected_flow <- dI_rate * I

    # Optional side-room patient-days proxy
    side_room_occupancy <- X

    list(
      c(dS, dC, dI, dX),
      c(
        N = N,
        lambda = lambda,
        alpha_eff = alpha_eff,
        total_cases_flow = total_cases_flow,
        identified_cases_flow = identified_cases_flow,
        missed_cases_flow = total_cases_flow - identified_cases_flow,
        false_neg_flow = false_neg_flow,
        discharged_colonised_flow = discharged_colonised_flow,
        discharged_infected_flow = discharged_infected_flow,
        side_room_occupancy = side_room_occupancy
      )
    )
  })
}