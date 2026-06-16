# Model1A_Functions.R
# ODE model, derived outputs, stabilisation detection, and plotting helpers.

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
})

make_model1a_times <- function(max_time, dt = 0.5) {
  seq(0, max_time, by = dt)
}

model1a_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    N <- pmax(S + C + I + X, 1e-12)

    # Convert times to rates.
    alpha <- 1 / t_incubation
    mu <- 1 / t_decolonisation
    gamma <- 1 / (t_identify + t_transfer)
    delta <- 1 / (t_test + t_treat + t_recover)
    pi <- 1 / (t_test + t_transfer_back)
    dXA <- 1 / (t_test + t_treat + t_recover + t_discharge)
    rateS <- 1 / LOS_S
    rateC <- 1 / LOS_C
    rateI <- 1 / LOS_I
    rateX <- 1 / LOS_X

    alpha_eff <- alpha * abx_mult

    lambda <- betaC * (C / N) + betaI * (I / N) + betaX * (X / N) + beta0 * E

    aS <- A * pS_adm
    aC <- A * pC_adm
    aI <- A * pI_adm

    # Core compartments
    dS_dt <- aS - lambda * S + mu * C + delta * X - rateS * S - mS * S
    dC_dt <- aC + lambda * S - (alpha_eff + mu + rateC + mC) * C
    dI_dt <- aI + alpha_eff * C - (gamma + rateI + mI) * I + pi * X
    dX_dt <- gamma * I - (delta + pi + dXA + rateX + mX) * X
    dE_dt <- etaC * C + etaI * I + etaX * X - rhoE * E

    # Cumulative outcomes
    dCumCases_dt <- alpha_eff * C + aI
    dCumIdentified_dt <- delta * X + dXA * X
    dCumFalseNeg_dt <- pi * X
    dCumDischC_dt <- rateC * C
    dCumDischI_dt <- rateI * I
    dCumXDays_dt <- X
    dCumDeaths_dt <- mS * S + mC * C + mI * I + mX * X
    dCumCI_dt <- alpha_eff * C
    dCumAdmI_dt <- aI
    dCumXS_dt <- delta * X
    dCumXI_dt <- pi * X
    dCumXDisch_dt <- dXA * X

    list(
      c(
        dS_dt, dC_dt, dI_dt, dX_dt, dE_dt,
        dCumCases_dt, dCumIdentified_dt, dCumFalseNeg_dt,
        dCumDischC_dt, dCumDischI_dt, dCumXDays_dt, dCumDeaths_dt,
        dCumCI_dt, dCumAdmI_dt, dCumXS_dt, dCumXI_dt, dCumXDisch_dt
      ),
      c(
        N = N,
        lambda = lambda,
        alpha_eff = alpha_eff,
        gamma = gamma,
        delta = delta,
        pi = pi,
        dXA = dXA,
        rateS = rateS,
        rateC = rateC,
        rateI = rateI,
        rateX = rateX
      )
    )
  })
}

calculate_model1a_derived_outputs <- function(df) {
  df$TotalCases <- df$CumCI + df$CumAdmI
  df$IdentifiedCases <- df$CumXS + df$CumXDisch
  df$MissedCases <- df$TotalCases - df$IdentifiedCases
  df$FalseNegatives <- df$CumFalseNeg
  df$InvasiveNotSuspected <- df$MissedCases - df$FalseNegatives
  dt <- if (nrow(df) > 1) df$time[2] - df$time[1] else 0
  df$SideRoomUse_PatientDays <- cumsum(c(0, head(df$X, -1))) * dt
  df
}

find_model1a_stabilisation <- function(df, tol = 1e-6, window = 10) {
  core <- c("S", "C", "I", "X")
  if (nrow(df) < (window + 2)) return(nrow(df))

  deltas <- as.matrix(abs(df[-1, core] - df[-nrow(df), core]))
  max_step_change <- apply(deltas, 1, max)
  stable_flag <- max_step_change < tol

  if (sum(stable_flag) < window) return(nrow(df))

  run_sum <- stats::filter(as.numeric(stable_flag), rep(1, window), sides = 1)
  idx <- which(!is.na(run_sum) & run_sum >= window)[1]
  if (is.na(idx)) nrow(df) else idx + 1
}

summarise_model1a_stabilisation <- function(df, stab_row) {
  stable <- df[stab_row, ]
  data.frame(
    time = stable$time,
    S = stable$S,
    C = stable$C,
    I = stable$I,
    X = stable$X,
    CumCases = stable$CumCases,
    CumIdentified = stable$CumIdentified,
    CumFalseNeg = stable$CumFalseNeg,
    CumDischC = stable$CumDischC,
    CumDischI = stable$CumDischI,
    CumXDays = stable$CumXDays,
    CumDeaths = stable$CumDeaths,
    CumCI = stable$CumCI,
    CumAdmI = stable$CumAdmI,
    CumXS = stable$CumXS,
    CumXI = stable$CumXI,
    CumXDisch = stable$CumXDisch,
    TotalCases = stable$TotalCases,
    IdentifiedCases = stable$IdentifiedCases,
    MissedCases = stable$MissedCases,
    FalseNegatives = stable$FalseNegatives,
    InvasiveNotSuspected = stable$InvasiveNotSuspected,
    SideRoomUse_PatientDays = stable$SideRoomUse_PatientDays,
    row.names = NULL,
    check.names = FALSE
  )
}

save_model1a_plot <- function(df, out_path, title, y_cols, y_label) {
  plot_df <- df[, c("time", y_cols)]
  long_df <- reshape(
    plot_df,
    varying = y_cols,
    v.names = "value",
    timevar = "series",
    times = y_cols,
    direction = "long"
  )
  long_df$series <- factor(long_df$series, levels = y_cols)

  p <- ggplot(long_df, aes(x = time, y = value, colour = series)) +
    geom_line(linewidth = 0.8) +
    labs(title = title, x = "Time (days)", y = y_label, colour = NULL) +
    theme_minimal(base_size = 12)

  ggsave(filename = out_path, plot = p, width = 10, height = 6, dpi = 300)
  invisible(p)
}

save_model1a_outputs <- function(df, stab_row, output_dir, scenario_tag) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  state_file <- file.path(output_dir, paste0(scenario_tag, "_state_trajectories.csv"))
  cum_file <- file.path(output_dir, paste0(scenario_tag, "_cumulative_outcomes.csv"))
  stab_file <- file.path(output_dir, paste0(scenario_tag, "_stabilisation_summary.csv"))
  plot_states_file <- file.path(output_dir, paste0(scenario_tag, "_states.png"))
  plot_cum_file <- file.path(output_dir, paste0(scenario_tag, "_cumulative_outcomes.png"))
  plot_outcomes_file <- file.path(output_dir, paste0(scenario_tag, "_key_outcomes.png"))

  write.csv(df[, c("time", "S", "C", "I", "X", "E")], state_file, row.names = FALSE)
  write.csv(
    df[, c(
      "time", "CumCases", "CumIdentified", "CumFalseNeg", "CumDischC", "CumDischI",
      "CumXDays", "CumDeaths", "CumCI", "CumAdmI", "CumXS", "CumXI", "CumXDisch",
      "TotalCases", "IdentifiedCases", "MissedCases", "FalseNegatives", "InvasiveNotSuspected",
      "SideRoomUse_PatientDays"
    )],
    cum_file,
    row.names = FALSE
  )
  write.csv(summarise_model1a_stabilisation(df, stab_row), stab_file, row.names = FALSE)

  save_model1a_plot(df, plot_states_file, "Model 1A: State trajectories", c("S", "C", "I", "X"), "People")
  save_model1a_plot(df, plot_cum_file, "Model 1A: Cumulative outcomes", c("CumCases", "CumIdentified", "CumFalseNeg", "CumDischC", "CumDischI"), "Cumulative count")
  save_model1a_plot(df, plot_outcomes_file, "Model 1A: Outcomes of interest", c("MissedCases", "FalseNegatives", "InvasiveNotSuspected", "SideRoomUse_PatientDays"), "Count / patient-days")

  list(
    state_file = state_file,
    cum_file = cum_file,
    stab_file = stab_file,
    plot_states_file = plot_states_file,
    plot_cum_file = plot_cum_file,
    plot_outcomes_file = plot_outcomes_file
  )
}
