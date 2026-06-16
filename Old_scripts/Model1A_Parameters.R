# Model1A_Parameters.R
# Baseline parameters and scenario controls for Model 1A.

# Build parameters for Model 1A

build_model1a_parameters <- function(
  scenario_id = "001",
  scenario_label = "Baseline",
  max_time = 365,
  dt = 0.5,
  initial_state = NULL
) {
  if (is.null(initial_state)) {
    initial_state <- c(
      S = 950,
      C = 30,
      I = 15,
      X = 5,
      E = 0,
      CumCases = 0,
      CumIdentified = 0,
      CumFalseNeg = 0,
      CumDischC = 0,
      CumDischI = 0,
      CumXDays = 0,
      CumDeaths = 0,
      CumCI = 0,
      CumAdmI = 0,
      CumXS = 0,
      CumXI = 0,
      CumXDisch = 0
    )
  }

  list(
    scenario_id = scenario_id,
    scenario_label = scenario_label,
    max_time = max_time,
    dt = dt,
    initial_state = initial_state,

    # Admissions
    A = 20,
    pS_adm = 0.90,
    pC_adm = 0.07,
    pI_adm = 0.03,

    # Time-to-rate inputs (days)
    t_incubation = 5,
    t_decolonisation = 20,
    t_identify = 1,
    t_transfer = 0.0833333,    # 2 hours
    t_test = 1,
    t_treat = 1,
    t_recover = 3,             # 72 hours
    t_transfer_back = 0.5,
    t_discharge = 1,

    # Length of stay inputs (days)
    LOS_S = 5,
    LOS_C = 6,
    LOS_I = 8,
    LOS_X = 12,

    # Transmission / environment
    betaC = 0.002,
    betaI = 0.005,
    betaX = 0.001,
    beta0 = 0.0,
    etaC = 0.001,
    etaI = 0.002,
    etaX = 0.0005,
    rhoE = 0.2,

    # Antibiotics effect on progression C -> I
    abx_mult = 1.5,

    # Mortality rates
    mS = 0.001,
    mC = 0.001,
    mI = 0.004,
    mX = 0.003,

    # Numerical settings
    stabilisation_tol = 1e-6,
    stabilisation_window = 10
  )
}
