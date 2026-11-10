# ================================================================================================ #
# Model 4 flow diagrams
#
# Builds a Graphviz flow diagram of the model for one isolation policy from a table of the model's
# flows (model_flows), drawn with DiagrammeR, with the states at fixed positions in three rows
# (diagram_pos) and discharge and death drawn from each state. Arrows are labelled with the parameters
# only (no values), so the diagram is general; a key box gives the formula for each rate shown (mu,
# alpha, lambda, gamma, sigma, tau, and the dis and mort rates).
#
# Flows that leave a state at one rate and are then split into outcomes by proportions are drawn as a
# forked arrow: an arrow from the state to a junction dot, labelled with the time (TIME), and arrows
# from the dot to each outcome, labelled with the proportion (PROPORTION). In Model 4 these are:
#   admissions:          adm_rate, split into S, C and I by prop_adm_C and prop_adm_I
#   test result:         1 / t_test_turnaround, split into positive (test_sens) and negative results
#   end of treatment:    1 / t_treatment, split into recovered (prop_Xtreat_to_S) and discharged
# Every other flow is a single arrow labelled with its rate.
#
# check_model_flows() checks that the table gives the same derivatives as model_ode() for a scenario,
# so a change to model_ode() that is not copied into the table is caught.
#
# Requires DiagrammeR (and DiagrammeRsvg and rsvg to save PNG files), and Model4_Parameters.R and
# Model4_Functions.R to be sourced first.
# ================================================================================================ #

# 1) Flow table ####

# One row per flow. Columns:
#   from, to: states ("Admissions" is the source of admissions; "Exit" is discharge or death, drawn as
#             a small grey node beside each state; "D" is death, counted in the model's D state)
#   route:    expression giving 1 if the flow is used under the scenario's policy, 0 if not (from the
#             via_T* flags in make_mod_parms())
#   rate:     per-capita rate for a single flow (multiplied by the size of `from`); for admissions, the
#             number per day
#   fork:     name of the fork the flow belongs to (NA for a single flow); its rate is stem * branch
#   stem:     the fork's exit rate (the same for every branch)
#   branch:   the proportion taking this branch
#   label:    edge label (single flows), or branch label (forks)
#   type:     "move" (between states), "exit_dis" (discharge) or "exit_mort" (death)
model_flows <- tibble::tribble(
  ~from,        ~to,          ~route,                   ~rate,         ~fork,       ~stem,                    ~branch,                         ~label,                                       ~type,
  "Admissions", "S",          "1",                      NA,            "adm",       "adm_rate",               "1 - prop_adm_C - prop_adm_I",   "1 - prop_adm_C - prop_adm_I (to S)",         "move",
  "Admissions", "C",          "1",                      NA,            "adm",       "adm_rate",               "prop_adm_C",                    "prop_adm_C (to C)",                          "move",
  "Admissions", "I",          "1",                      NA,            "adm",       "adm_rate",               "prop_adm_I",                    "prop_adm_I (to I)",                          "move",
  "S",          "C",          "1",                      "lambda",      NA,          NA,                       NA,                              "lambda (Force of colonisation)",             "move",
  "C",          "I",          "1",                      "alpha",       NA,          NA,                       NA,                              "alpha (Progression to CDI)",                 "move",
  "C",          "S",          "1",                      "mu",          NA,          NA,                       NA,                              "mu (Decolonisation)",                        "move",
  "I",          "Xtest",      "1 - via_T3",             "gamma",       NA,          NA,                       NA,                              "gamma (Sent for testing)",                   "move",
  "I",          "T3",         "via_T3",                 "gamma",       NA,          NA,                       NA,                              "gamma (Sent for testing)",                   "move",
  "T3",         "Xtest",      "1",                      "tau",         NA,          NA,                       NA,                              "tau (Patient transfer)",                     "move",
  "Xtest",      "Xtreat",     "1 - via_T1",             NA,            "result",    "1 / t_test_turnaround",  "test_sens",                     "test_sens (True positive)",                  "move",
  "Xtest",      "T1",         "via_T1",                 NA,            "result",    "1 / t_test_turnaround",  "test_sens",                     "test_sens (True positive)",                  "move",
  "Xtest",      "FN",         "1 - via_T4",             NA,            "result",    "1 / t_test_turnaround",  "1 - test_sens",                 "1 - test_sens (False negative)",             "move",
  "Xtest",      "T4",         "via_T4",                 NA,            "result",    "1 / t_test_turnaround",  "1 - test_sens",                 "1 - test_sens (False negative)",             "move",
  "T4",         "FN",         "1",                      "tau",         NA,          NA,                       NA,                              "tau (Patient transfer)",                     "move",
  "FN",         "Xtest",      "1 - via_T3",             "sigma",       NA,          NA,                       NA,                              "sigma (Sent for retesting)",                 "move",
  "FN",         "T3",         "via_T3",                 "sigma",       NA,          NA,                       NA,                              "sigma (Sent for retesting)",                 "move",
  "T1",         "Xtreat",     "1",                      "tau",         NA,          NA,                       NA,                              "tau (Patient transfer)",                     "move",
  "Xtreat",     "S",          "1 - via_T2 - via_T5",    NA,            "treat",     "1 / t_treatment",        "prop_Xtreat_to_S",              "prop_Xtreat_to_S (Recovered)",               "move",
  "Xtreat",     "T2",         "via_T2",                 NA,            "treat",     "1 / t_treatment",        "prop_Xtreat_to_S",              "prop_Xtreat_to_S (Recovered)",               "move",
  "Xtreat",     "T5",         "via_T5",                 NA,            "treat",     "1 / t_treatment",        "prop_Xtreat_to_S",              "prop_Xtreat_to_S (Recovered)",               "move",
  "Xtreat",     "Exit",       "1",                      NA,            "treat",     "1 / t_treatment",        "1 - prop_Xtreat_to_S",          "1 - prop_Xtreat_to_S (Discharged)",          "exit_dis",
  "T2",         "S",          "1",                      "tau",         NA,          NA,                       NA,                              "tau (Patient transfer)",                     "move",
  "T5",         "S",          "1",                      "tau",         NA,          NA,                       NA,                              "tau (Patient transfer)",                     "move",
  # Discharges (not a model state) and deaths (to D)
  "S",          "Exit",       "1",                      "disS",        NA,          NA,                       NA,                              "disS",                                       "exit_dis",
  "C",          "Exit",       "1",                      "disC",        NA,          NA,                       NA,                              "disC",                                       "exit_dis",
  "I",          "Exit",       "1",                      "disI",        NA,          NA,                       NA,                              "disI",                                       "exit_dis",
  "FN",         "Exit",       "1",                      "disFN",       NA,          NA,                       NA,                              "disFN",                                      "exit_dis",
  "T4",         "Exit",       "1",                      "disT4",       NA,          NA,                       NA,                              "disT4",                                      "exit_dis",
  "S",          "D",          "1",                      "mortS",       NA,          NA,                       NA,                              "mortS",                                      "exit_mort",
  "C",          "D",          "1",                      "mortC",       NA,          NA,                       NA,                              "mortC",                                      "exit_mort",
  "I",          "D",          "1",                      "mortI",       NA,          NA,                       NA,                              "mortI",                                      "exit_mort",
  "Xtest",      "D",          "1",                      "mortXtest",   NA,          NA,                       NA,                              "mortXtest",                                  "exit_mort",
  "FN",         "D",          "1",                      "mortFN",      NA,          NA,                       NA,                              "mortFN",                                     "exit_mort",
  "Xtreat",     "D",          "1",                      "mortXtreat",  NA,          NA,                       NA,                              "mortXtreat",                                 "exit_mort",
  "T1",         "D",          "1",                      "mortT1",      NA,          NA,                       NA,                              "mortT1",                                     "exit_mort",
  "T2",         "D",          "1",                      "mortT2",      NA,          NA,                       NA,                              "mortT2",                                     "exit_mort",
  "T3",         "D",          "1",                      "mortT3",      NA,          NA,                       NA,                              "mortT3",                                     "exit_mort",
  "T4",         "D",          "1",                      "mortT4",      NA,          NA,                       NA,                              "mortT4",                                     "exit_mort",
  "T5",         "D",          "1",                      "mortT5",      NA,          NA,                       NA,                              "mortT5",                                     "exit_mort"
)

# States: label, whether infectious (has a transmission term in lambda) and where the patient is.
# location: "ward" (general ward), "side" (side room) or "policy" (Xtest and Xtreat: depends on the
# isolation policy).
model_states <- tibble::tribble(
  ~state,   ~label,                                         ~infectious, ~location,
  "S",      "S\nSusceptible",                               FALSE,       "ward",
  "C",      "C\nColonised",                                 TRUE,        "ward",
  "I",      "I\nInfected\nUntested CDI cases",              TRUE,        "ward",
  "Xtest",  "Xtest\nSuspected CDI cases\nbeing tested",     TRUE,        "policy",
  "FN",     "FN\nFalse negatives",                          TRUE,        "ward",
  "Xtreat", "Xtreat\nConfirmed CDI cases\nbeing treated",   TRUE,        "policy",
  "T1",     "T1\nConfirmed CDI cases\nwaiting side room",   TRUE,        "ward",
  "T2",     "T2\nRecovered\nAwaiting ward bed",             FALSE,       "side",
  "T3",     "T3\nSuspected CDI cases\nawaiting side room",  TRUE,        "ward",
  "T4",     "T4\nFalse negatives\nawaiting ward bed",       TRUE,        "side",
  "T5",     "T5\nRecovered CDI cases\nawaiting ward bed",   FALSE,       "side"
)

# 2) Evaluate the flow table for one scenario ####

# Environment holding the scenario's parameters (make_mod_parms() output plus the fixed parameters
# it does not return, e.g. t_treatment), the state sizes and the force of colonisation (lambda).
flow_env <- function(parms, state = NULL, fixed_parms = fixed) {
  vals <- c(parms, fixed_parms[setdiff(names(fixed_parms), names(parms))])
  if (!is.null(state)) {
    vals <- c(vals, as.list(state))
    vals$N <- with(vals, S + C + I + Xtest + FN + Xtreat + T1 + T2 + T3 + T4 + T5)
    vals$lambda <- with(vals, beta0 + (betaC * C + betaI * I + betaXtest * Xtest + betaFN * FN +
                                         betaXtreat * Xtreat + betaT1 * T1 + betaT3 * T3 + betaT4 * T4) / N)
  }
  list2env(vals)
}

ev <- function(expr, env) if (is.na(expr)) NA_real_ else eval(parse(text = expr), envir = env)

# The flows used under the scenario's policy (route = 1), with the numeric route, rate, stem and branch.
# For forks, rate = stem * branch.
scenario_flows <- function(parms, state = NULL, flows = model_flows) {
  env <- flow_env(parms, state)
  flows |>
    dplyr::mutate(
      route_val  = vapply(route, ev, numeric(1), env = env),
      stem_val   = vapply(stem, ev, numeric(1), env = env),
      branch_val = vapply(branch, ev, numeric(1), env = env),
      rate_val   = dplyr::if_else(is.na(fork), vapply(dplyr::coalesce(rate, NA_character_), \(r)
                     if (is.na(r) || r == "lambda" && is.null(state)) NA_real_ else ev(r, env), numeric(1)),
                     stem_val * branch_val)
    ) |>
    dplyr::filter(route_val > 0)
}

# 3) Check the table against model_ode() ####

# For a scenario's parameters and a random (positive) state, compares the derivative of every state
# from the flow table (inflows minus outflows) with model_ode(). Stops if any differ.
check_model_flows <- function(parms, seed = 1, tol = 1e-9) {
  set.seed(seed)
  states <- c("S", "C", "I", "Xtest", "FN", "Xtreat", "T1", "T2", "T3", "T4", "T5", "D")
  state  <- setNames(runif(length(states), 1, 100), states)
  f <- scenario_flows(parms, state)
  size <- c(state, Admissions = 1)
  f$flow <- f$rate_val * size[f$from]
  d_table <- vapply(states, \(s) sum(f$flow[f$to == s]) - sum(f$flow[f$from == s]), numeric(1))
  d_model <- setNames(model_ode(0, state, parms)[[1]], states)
  diff <- abs(d_table - d_model)
  if (any(diff > tol * pmax(1, abs(d_model)))) {
    stop("model_flows does not match model_ode() for ", parms$isolation_policy, ": ",
         paste(sprintf("%s (table %.6g, model %.6g)", states, d_table, d_model)[diff > tol * pmax(1, abs(d_model))],
               collapse = "; "))
  }
  invisible(TRUE)
}

# 4) Build the diagram ####

# Splits a label "name (description)" into its two parts
label_name <- function(x) sub(" \\(.*$", "", x)
label_desc <- function(x) ifelse(grepl("\\(", x), sub("^.*\\((.*)\\)$", "\\1", x), "")

# Fixed positions (inches, x to the right, y up) for the three-row layout:
#   top:    admission, colonisation and testing, left to right
#   middle: false negatives, right to left, with the retest arrow back up to T3 / Xtest
#   bottom: treatment and recovery, right to left, with the recovery arrow back up to S
diagram_pos <- tibble::tribble(
  ~node,          ~x,    ~y,
  "Admissions",   0.0,   7.6,
  "fork_adm",     2.6,   7.6,
  "S",            1.5,   5.6,
  "C",            4.5,   5.6,
  "I",            7.5,   5.6,
  "T3",          11.0,   5.6,
  "Xtest",       14.5,   5.6,
  "fork_result", 14.5,   3.0,
  "T4",          11.0,   3.0,
  "FN",           7.5,   3.0,
  "T1",          14.5,   0.4,
  "Xtreat",      11.0,   0.4,
  "fork_treat",   8.3,   0.4,
  "T2",           3.6,   0.4,
  "T5",           3.6,   0.4,
  "S_exit_dis",   0.5,   4.5,  "S_exit_mort",   1.4,   4.4,
  "C_exit_dis",   3.9,   4.5,  "C_exit_mort",   5.1,   4.5,
  "I_exit_dis",   6.9,   4.5,  "I_exit_mort",   8.1,   4.5,
  "T3_exit_mort", 11.0,  4.3,
  "Xtest_exit_mort", 15.9, 4.6,
  "T4_exit_dis",  10.4,  1.8,  "T4_exit_mort",  11.6,  1.8,
  "FN_exit_dis",  6.9,   1.8,  "FN_exit_mort",  8.1,   1.8,
  "T1_exit_mort", 14.5, -0.9,
  "Xtreat_exit_dis", 7.5, -0.9, "Xtreat_exit_mort", 11.0, -0.9,
  "T2_exit_mort", 3.6,  -0.9,
  "T5_exit_mort", 3.6,  -0.9
)

# Position of the formula key box in the three-row layout (the empty space above T3 and Xtest)
diagram_pos <- dplyr::bind_rows(diagram_pos, tibble::tibble(node = "formula_key", x = 12.3, y = 7.7))

# Edge labels drawn at fixed positions in the three-row layout (as text placed below each arrow),
# instead of being placed by Graphviz. Node id "lbl_<from>_<to>".
#   C -> S (decolonisation, upper arrow): just below the arrow, between the two S-C arrows
#   S -> C (force of colonisation, lower arrow): below the arrow
fixed_edge_labels <- tibble::tribble(
  ~from, ~to,  ~x,   ~y,
  "C",   "S",  3.0,  5.62,
  "S",   "C",  3.0,  4.95
)
diagram_pos <- dplyr::bind_rows(
  diagram_pos,
  fixed_edge_labels |> dplyr::transmute(node = paste0("lbl_", from, "_", to), x, y))

# pos attributes for the nodes that are drawn (`present`)
pos_lines <- function(present, pos = diagram_pos) {
  p <- pos[pos$node %in% present, ]
  sprintf("  %s [pos=\"%.2f,%.2f!\"];", p$node, p$x, p$y)   # neato reads pos in inches
}

# Readable text for each fork stem (shown on the arrow into the junction dot)
fork_stems <- c(
  adm    = "adm_rate",
  result = "TIME: Test turnaround",
  treat  = "TIME: Treatment\nrecovery & clearance"
)

# Formula key: one line per rate drawn on the diagram (`rates`: the rate names used by the drawn flows;
# `states`: the states drawn, for the terms of lambda), left-justified ("\l" ends each line in DOT)
formula_key <- function(rates, states, parms) {
  pol <- parms$isolation_policy
  # lambda: environmental exposure plus one term per infectious state drawn
  beta_of <- c(C = "betaC", I = "betaI", Xtest = "betaXtest", FN = "betaFN", Xtreat = "betaXtreat",
               T1 = "betaT1", T3 = "betaT3", T4 = "betaT4")
  inf <- intersect(names(beta_of), states)
  terms <- paste0(beta_of[inf], " × ", inf)
  lambda_lines <- c(
    paste0("lambda = beta0 + (", paste(head(terms, 4), collapse = " + "), if (length(terms) > 4) " +" else ") / N"),
    if (length(terms) > 4) paste0("              ", paste(terms[-(1:4)], collapse = " + "), ") / N"),
    "    N = number of patients in hospital (all states except D)",
    paste0("    betaXtest = ", if (pol == "suspected_and_confirmed") "betaI × beta_mult_isolated (side room)" else "betaI",
           "; betaFN = betaI"),
    if (pol != "none") "    betaXtreat = betaXtreat × beta_mult_isolated (side room)",
    if (any(c("T1", "T3") %in% states))
      paste0("    ", paste(c("betaT1", "betaT3")[c("T1", "T3") %in% states], collapse = " = "), " = betaI"),
    if ("T4" %in% states) "    betaT4 = betaI × beta_mult_isolated"
  )
  lines <- c(
    "KEY TO RATES (per day)",
    if ("mu" %in% rates)     "mu = prop_decolonise_C_to_S / t_decolonise_C_to_S",
    if ("alpha" %in% rates)  "alpha = prop_progress_C_to_I × ((1 - abx_prop) + abx_prop × abx_mult) / t_progress_C_to_I",
    if ("lambda" %in% rates) lambda_lines,
    if ("gamma" %in% rates)  "gamma = -log(1 - prob_test_daily)",
    if ("sigma" %in% rates)  "sigma = 1 / t_wait_retest",
    if ("tau" %in% rates)    "tau = 1 / t_patient_transfer",
    if (any(grepl("^dis", rates)))  "disS = 1/losS, disC = 1/losC, disI = 1/losI, disFN = 1/losFN; disT4 = disFN",
    if (any(grepl("^mort", rates))) c(
      "mortS = -log(1 - prob_mort_S)/losS; mortC, mortT2, mortT5 = mortS",
      "mortI, mortXtest, mortFN, mortT1, mortT3, mortT4 = mortI_mult × mortS; mortXtreat = mortXtreat_mult × mortS")
  )
  # Graphviz drops leading spaces, so indents use non-breaking spaces
  lines <- sub("^ +", strrep(" ", 6), lines)
  paste0(paste(lines, collapse = "\\l"), "\\l")
}

# DOT text for one isolation policy. parms: make_mod_parms() output for that policy (only
# isolation_policy affects the diagram; it sets which transfer states are used).
build_flow_dot <- function(parms) {
  f <- scenario_flows(parms)
  # Keep only states patients can reach under this policy (e.g. no T3-T5 under P2)
  active <- "Admissions"
  repeat {
    new <- unique(c(active, f$to[f$type == "move" & f$from %in% active]))
    if (length(new) == length(active)) break
    active <- new
  }
  f <- f |> dplyr::filter(from %in% active)
  q <- \(x) paste0("\"", gsub("\"", "\\\\\"", x), "\"")

  # Node styles: infectious states black (thick) outline, non-infectious grey; side rooms blue fill
  side_room <- function(state) {
    loc <- model_states$location[match(state, model_states$state)]
    pol <- parms$isolation_policy
    ifelse(loc == "side", TRUE,
      ifelse(loc == "policy",
             (state == "Xtest" & pol == "suspected_and_confirmed") |
             (state == "Xtreat" & pol %in% c("confirmed", "suspected_and_confirmed")),
             FALSE))
  }
  used <- unique(c(f$from, f$to))
  st <- model_states |> dplyr::filter(state %in% used)
  node_lines <- sprintf(
    "  %s [label=%s, color=%s, penwidth=%s, fillcolor=%s, style=\"filled%s\"];",
    st$state, q(st$label),
    ifelse(st$infectious, "black", "grey55"),
    ifelse(st$infectious, "2", "1"),
    ifelse(side_room(st$state), "\"#CFE3F4\"", "white"),
    ifelse(grepl("^T", st$state), ",dashed", ""))

  # Labels: description on the first line (long ones wrapped onto two lines, so labels fit between
  # nearby states), then the parameter
  lab <- function(text) {
    desc <- paste(strwrap(label_desc(text), 20), collapse = "\n")
    paste0(if (desc != "") paste0(desc, "\n") else "", label_name(text))
  }

  edge_lines <- character(0)
  extra_nodes <- character(0)
  exit_node <- function(s, type) {
    id <- paste0(s, "_", type)
    extra_nodes <<- unique(c(extra_nodes, sprintf(
      "  %s [label=%s, shape=plaintext, fontcolor=grey40, fontsize=9];", id,
      q(if (type == "exit_mort") "Died" else "Discharged"))))
    id
  }
  target <- function(row) {
    if (row$type == "move") row$to else exit_node(row$from, row$type)
  }

  # Single flows
  singles <- f |> dplyr::filter(is.na(fork))
  for (i in seq_len(nrow(singles))) {
    r <- singles[i, ]
    is_exit <- r$type != "move"
    text <- lab(r$label)
    # Two states with arrows both ways (S and C): draw one above and one below, so they don't overlap
    two_way <- r$type == "move" && any(singles$from == r$to & singles$to == r$from)
    ports <- if (two_way) {
      if (r$from < r$to) ", tailport=n, headport=n" else ", tailport=s, headport=s"
    } else if (r$to == "S") ", headport=s" else ""   # recovered patients return up into S
    # Labels with a fixed position (fixed_edge_labels): a separate text node, and the arrow unlabelled
    if (any(fixed_edge_labels$from == r$from & fixed_edge_labels$to == r$to)) {
      id <- paste0("lbl_", r$from, "_", r$to)
      extra_nodes <- c(extra_nodes, sprintf("  %s [shape=plaintext, fontsize=10, margin=0, label=%s];", id, q(text)))
      edge_lines <- c(edge_lines, sprintf("  %s -> %s [%s];", r$from, target(r), sub("^, ", "", ports)))
      next
    }
    edge_lines <- c(edge_lines, sprintf(
      "  %s -> %s [label=%s%s%s];", r$from, target(r), q(text), ports,
      if (is_exit) ", color=grey60, fontcolor=grey40, fontsize=8, arrowsize=0.6" else
        ""))
  }

  # Forks: state -> junction dot (TIME), dot -> each outcome (PROPORTION)
  forks <- f |> dplyr::filter(!is.na(fork))
  for (k in unique(forks$fork)) {
    fk <- forks |> dplyr::filter(fork == k)
    from <- fk$from[1]
    dot  <- paste0("fork_", k)
    stem_text <- fork_stems[[k]]
    extra_nodes <- c(extra_nodes, sprintf("  %s [shape=point, width=0.09, label=\"\"];", dot))
    if (from == "Admissions") {
      extra_nodes <- c(extra_nodes, "  Admissions [shape=plaintext, fontsize=11];")
    }
    edge_lines <- c(edge_lines, sprintf("  %s -> %s [label=%s, fontcolor=\"#B00000\", arrowhead=none];",
                                        from, dot, q(stem_text)))
    for (i in seq_len(nrow(fk))) {
      r <- fk[i, ]
      text <- paste0("PROPORTION: ", lab(r$label))
      edge_lines <- c(edge_lines, sprintf(
        "  %s -> %s [label=%s, fontcolor=\"#B00000\"%s];", dot, target(r), q(text),
        if (r$type != "move") ", color=grey60, arrowsize=0.6" else
          if (r$to == "S" && k == "treat") ", headport=s" else ""))
    }
  }

  policy_name <- c(none = "P1: No isolation", confirmed = "P2: Confirmed cases isolated",
                   suspected_and_confirmed = "P3: Suspected and confirmed cases isolated")[parms$isolation_policy]
  legend <- paste("Thick black outline: infectious. Grey outline: not infectious. Blue fill: in a side room.",
                  "Dashed: awaiting transfer. Dot: one exit time split into outcomes by proportion.")
  extra_nodes <- c(extra_nodes, sprintf(
    "  formula_key [shape=box, style=rounded, color=grey60, fontsize=10, margin=\"0.15,0.1\", label=%s];",
    q(formula_key(unique(singles$rate), st$state, parms))))

  paste(c(
    "digraph model {",
    sprintf("  graph [layout=neato, splines=true, forcelabels=true, fontname=Arial, fontsize=14, labelloc=t, label=%s];",
            q(paste0(policy_name, "\n", legend, "\n\n"))),
    "  node [shape=ellipse, fontname=Arial, fontsize=12, margin=\"0.04,0.02\"];",
    "  edge [fontname=Arial, fontsize=10, arrowsize=0.7];",
    node_lines, extra_nodes,
    # Edge labels are external labels (xlabel), which Graphviz places beside the arrow, clear of other
    # labels and nodes, rather than centred on it
    sub(" \\[label=", " [xlabel=", edge_lines),
    pos_lines(c(st$state, sub(" .*", "", trimws(extra_nodes)))),
    "}"), collapse = "\n")
}

# Draw (DiagrammeR htmlwidget, for the Quarto report or the RStudio viewer)
plot_flow_diagram <- function(parms) DiagrammeR::grViz(build_flow_dot(parms))

# Save as PNG (width in pixels)
save_flow_diagram <- function(parms, path, width = 2600) {
  svg <- DiagrammeRsvg::export_svg(plot_flow_diagram(parms))
  rsvg::rsvg_png(charToRaw(svg), path, width = width)
  invisible(path)
}
