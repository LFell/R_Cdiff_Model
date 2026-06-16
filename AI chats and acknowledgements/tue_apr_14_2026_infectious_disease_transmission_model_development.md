# Copilot Chat Conversation Export: Infectious disease transmission model development

**User:** @LFell
**Thread URL:** https://github.com/copilot/c/21289557-e173-4b11-9d81-7eb25d0b6165

## @LFell

I want to create a, ODE infectious  disease transmission model, using a similar script structure as on the attached. My model is simpler, but has a similar aim. 

## Copilot

Great idea — and your attached structure is a very strong template to reuse.

For a **simpler infectious disease ODE model** with a similar workflow, I recommend you keep the same 3-file pattern:

1. `Model_Functions.R` (ODE + helper functions)  
2. `Model_Parameters.R` (fixed values + scenario builder)  
3. `Model_Run.qmd` (run scenarios, summaries, plots, CEA-like outputs)

Here’s a clean starter you can adapt.

---

### 1) `Simple_ODE_Functions.R`

```r name=Simple_ODE_Functions.R
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
```

---

### 2) `Simple_ODE_Parameters.R`

```r name=Simple_ODE_Parameters.R
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
```

---

### 3) `Simple_ODE_Model_Run.qmd`

```r name=Simple_ODE_Model_Run.qmd
---
title: "Simple ODE Model Run"
format: html
execute:
  echo: true
  warning: false
  message: false
---

```{r}
rm(list = ls())
library(deSolve)
library(tidyverse)

dir.create("Outputs", showWarnings = FALSE, recursive = TRUE)

source("Scripts/Simple_ODE_Functions.R")
source("Scripts/Simple_ODE_Parameters.R")
```

```{r}
# Time settings
model.time <- 365
mod.t <- seq(0, model.time, by = 1)

# Initial states
init <- c(
  S = 9990,
  I = 10,
  R = 0
)

run_scenario <- function(parms) {
  as.data.frame(ode(y = init, times = mod.t, func = simple_ode_model, parms = parms))
}

add_derived <- function(df) {
  df %>%
    mutate(
      cum_infections = cumsum(new_infections),
      cum_deaths_inf = cumsum(new_deaths_infection),
      cum_cost = cumsum(daily_total_cost)
    )
}
```

```{r}
# Run scenarios
scenarios <- list(
  base = mod.parms.base,
  int1 = mod.parms.int1,
  int2 = mod.parms.int2
)

results <- lapply(scenarios, run_scenario)
results <- lapply(results, add_derived)

# Save outputs
walk2(results, names(results), ~ write_csv(.x, paste0("Outputs/solution_", .y, ".csv")))
```

```{r}
# Summary table
summarise_outcomes <- function(df) {
  tibble(
    peak_I = max(df$I),
    total_infections = max(df$cum_infections),
    total_inf_deaths = max(df$cum_deaths_inf),
    total_cost = max(df$cum_cost)
  )
}

summary_tbl <- bind_rows(
  base = summarise_outcomes(results$base),
  int1 = summarise_outcomes(results$int1),
  int2 = summarise_outcomes(results$int2),
  .id = "scenario"
)

write_csv(summary_tbl, "Outputs/summary_table.csv")
knitr::kable(summary_tbl)
```

```{r}
# Plot I over time
plot_df <- bind_rows(
  lapply(names(results), function(nm) results[[nm]] %>% mutate(scenario = nm))
)

ggplot(plot_df, aes(time, I, color = scenario)) +
  geom_line(size = 1) +
  theme_classic() +
  labs(title = "Infectious population over time", y = "I", x = "Day")

ggsave("Outputs/I_over_time.png", width = 8, height = 5, dpi = 300)
```
```

---

## Why this maps well to your current structure

- Your current model separates **equations**, **parameters**, and **analysis** — keep that.
- You can scale from SIR to SEIR or add testing/isolation later without changing architecture.
- You keep reproducible outputs (`Outputs/*.csv`, plots, summary tables), just like your attached scripts.

---

If you want, next I can give you a **drop-in SEIR version** (still simple) with:
- `S, E, I, R`
- imported infections
- intervention toggles (reduce beta, increase detection)
- cost-effectiveness summary vs baseline exactly like your i0/i1/i2 layout.

## @LFell

Can I give you some text that describes the model

## Copilot

Yes — please do.  

Send whatever you have (even rough notes), and I’ll turn it into:

1. **Compartment structure**  
2. **ODE equations**  
3. **Parameter list + assumptions**  
4. **R scripts in the same structure** (`Functions`, `Parameters`, `Run.qmd`)

If helpful, use this quick template:

- Population/setting:
- Compartments (e.g., S, E, I, R):
- Transitions (who moves where, and why):
- Key rates/probabilities:
- Interventions/scenarios to compare:
- Outputs you want (cases, deaths, costs, ICER, etc.):
- Time horizon + time step:
- Any fixed numbers you already have:

Even partial info is fine.

## @LFell

Model 1A states: 
Susceptible (S) – all patients not currently colonised, infected or in isolation. Assumed to have no immunity to C. difficile.
Colonised (C) – patients that are colonised with toxigenic C. difficile, that would test positive for organism but negative for toxins by gold standard refence methods. It includes patients that become colonised within hospital, or that have been colonised for some time (e.g., colonised on admission to hospital). 
Infected (I) – patients with an invasive infection, that would test positive for organism and toxins by gold standard reference methods. For modelling simplicity, it is assumed that there is no recovery from an invasive infection without specific treatment for CDI, as toxin positivity is associated with severe symptoms. In reality, mild infections may be resolved by suspending treatment with antibiotics and proton pump inhibitors, and may not test positive for toxins.
Enhanced IPC (X) – patients transferred to a side room and transmission-based precautions applied. In reality, transmission-based precautions may be applied to patients on a general ward, or cohorted with other patients with the same infection. 
There is no ‘Recovered’ or ‘Immune’ state. Once treated,  patients in X revert to ‘Susceptible’ on recovery and are either transferred to a general ward or discharged into the community. Recovered cases can be counted by adding up the X to S transfers and discharges from X into the community. 
Transitions between states:
The force of colonisation (λ) is a function of the numbers of patients in each of the states S, C, I and X, and the transmission parameters for C, I and X. The model could differentiate transmission parameters, e.g. βI > βC > βX, and include an environmental transmission parameter β0, which could vary with the numbers in C, I and X in a previous time period (so environmental transmission can increase due to increased contamination with spores).
The rate of progression from colonisation to invasive infection (α) is an average, covering patients that go through a colonisation stage before infection, and those that immediately become infected.
Colonised patients can become de-colonised without intervention, reverting to S at a rate μ.
The default transitions between I and X are described in Model 1A, where suspected cases are isolated before testing. Model 1B represents circumstances where only confirmed cases are isolated.
Antibiotic prescribing:
Antibiotic prescribing can be included in this model in two ways, which should produce the same outcomes:
enhancing the rate of progression from C to I, on the basis that antibiotic prescribing increases the risk of invasive infection among those already colonised; or
enhancing the force of colonisation, on the basis that antibiotic prescribing increases susceptibility to colonisation. 
This model will use option i).
Model 1A is consistent with the recommended ‘SIGHT’ Protocol, which requires transfer of suspected cases to a side room within 2 hours of an episode of potentially infectious diarrhoea. Patients are released from isolation and either discharged from hospital or transferred back to a general ward once they have been symptom-free for 72 hours.
Testing occurs during isolation in a side room, and results either in treatment following a positive test result, or returning the patient to a general ward following a negative test result.
The suspected case isolation rate (γ) is a function of the proportion of I that are identified as a suspected case, and the time between symptom onset and transfer to a side room.
Patients with an invasive infection that are not identified as suspected cases are retained in I, and may be discharged into the community at a rate dI.
The recovery rate (δ) is a function of test sensitivity and proportion of treated cases that are transferred to a general ward after treatment, as well as the time required to test, treat, ensure absence of symptoms for 72 hours and transfer into a general ward. 
False negatives are returned to a general ward, reverting to I at a rate π which is a function of test sensitivity and time required to test and transfer a patient back to a general ward. They may be discharged at a rate dI.
Treated cases can also be discharged from hospital, so the discharge rate from X (dXA) is a function of test sensitivity and proportion of treated cases that are discharged, as well as the time required to test, treat, ensure absence of symptoms for 72 hours and discharge from hospital. 
Control variables:
Proportion of invasive infections that are identified as suspected cases
Sensitivity of testing methods 

Parameters:
Prevalence of colonised and infected in patients admitted to hospital.
Length of stay for S, C, I and X.
Mortality rates for S, C, I and X.
Force of colonisation: transmission parameters for C, I, X and environmental transmission and duration of latent period.
Rate of progression from colonisation to an invasive infection: proportion that progress and duration of incubation period.
Antibiotic prescribing rates and multiplier for the rate of progression from colonisation to an invasive infection.
Time required to:
identify a suspected case, e.g. observe 3 episodes over 24 hours;
transfer to a side room (I to X) (Model A: SIGHT recommends transfer within 2 hours of diarrhoeal episode);
obtain a test result;
treat;
recover (to include 72 hours symptom-free);
discharge to the community from X; and
return to a general ward (X to I or X to S).
Total cases = Σ flows from C to I (α) + admissions to I (aI)
Identified cases = Σ flows from X to S (δ or ε) + Σ discharges from X (dXA or dXB)
Missed cases = Total cases – Identified cases
Discharged colonised = Σ discharges from C (dC)
Discharged infected = Σ discharges from I (dI)

Model 1A: 
False negatives = Σ flows from X to I (π)
Invasive infections not identified as suspected cases = Missed cases – False negatives
Side room use (patient days) = 
      Σ flows from X to S (δ) * days spent in X (including testing, treatment, recovery and transfer) +
      Σ flows from X to I (π) * days spent in X (including testing and transfer)
      Σ discharges from X (dXA) * days spent in X (including testing, treatment, recovery and transfer)


## Copilot

Yes — this is enough to define a clean **Model 1A ODE**.  
Below is the direct model formulation + an R-ready script structure matching your existing style.

---

## 1) Proposed Model 1A compartments and flows

You have 4 states:

- **S**: susceptible
- **C**: colonised
- **I**: invasive infection (not isolated)
- **X**: isolated suspected/confirmed cases in side room

No R state (recovered go back to S).

### Force of colonisation
Use:

\[
\lambda = \beta_0 + \frac{\beta_C C + \beta_I I + \beta_X X}{N}
\]

where \(N=S+C+I+X\).  
(You can later add lagged environmental contamination if needed.)

### Core transitions
- \(S \to C\): \(\lambda S\)
- \(C \to I\): \(\alpha_{eff} C\), where \(\alpha_{eff} = \alpha \times abx\_mult\)
- \(C \to S\): \(\mu C\) (decolonisation)
- \(I \to X\): \(\gamma I\) (suspected-case isolation)
- \(X \to S\): \(\delta X\) (true positives treated/recovered OR true negatives clinically resolved)
- \(X \to I\): \(\pi X\) (false negatives returned to general ward)
- discharges: \(d_S S, d_C C, d_I I, d_X X\)
- mortality: \(m_S S, m_C C, m_I I, m_X X\)
- admissions into S/C/I: \(a_S, a_C, a_I\) (optional dynamic balancing can be added later)

---

## 2) ODE system (Model 1A)

\[
\frac{dS}{dt}=a_S+\mu C+\delta X-\lambda S-(d_S+m_S)S
\]

\[
\frac{dC}{dt}=a_C+\lambda S-(\alpha_{eff}+\mu+d_C+m_C)C
\]

\[
\frac{dI}{dt}=a_I+\alpha_{eff}C+\pi X-(\gamma+d_I+m_I)I
\]

\[
\frac{dX}{dt}=\gamma I-(\delta+\pi+d_X+m_X)X
\]

---

## 3) R function file (drop-in skeleton)

```r name=Model1A_ODE_Functions.R
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
```

---

## 4) Key point on \(\gamma, \delta, \pi\)

You described these as functions of detection proportion, sensitivity, and timing.  
Best practice: compute them in your **parameter script** from primitive inputs, e.g.:

- `gamma = p_suspected / t_identify_isolate`
- `pi = (1 - sens) / t_test_return`
- `delta = sens * p_treated_recover / t_test_treat_recover`

This keeps assumptions transparent.

---

If you want, next message I can produce the full three files (`Functions.R`, `Parameters.R`, `Run.qmd`) in your exact style with:
- base + 2 intervention scenarios,
- annual summaries,
- missed vs identified cases,
- side-room patient-days,
- and CSV outputs.

## @LFell

can you create base + 3 scenarios; identify the model states and flows at stabilisation; break missed cases down into untested and false negatives; include all plots and tables for results

## Copilot

Absolutely — yes.  

I can build that full package for you in the same structure as your current scripts, with:

- **Base + 3 scenarios** (4 total runs)
- **Stabilisation detection** and **state/flow values at stabilisation**
- **Missed cases split** into:
  - **untested/unsuspected** (I not moved to X)
  - **false negatives** (X → I)
- **Full results tables + plots + CSV outputs**

## What I’ll include

### States
- `S, C, I, X`

### Core flow outputs (daily + cumulative)
- `S_to_C`, `C_to_I`, `C_to_S`, `I_to_X`, `X_to_S`, `X_to_I`
- admissions/discharges/deaths by state
- `total_cases`, `identified_cases`, `missed_cases`
- `missed_untested`, `missed_false_negative`

### Stabilisation outputs
- day of stabilisation per scenario
- values of `S,C,I,X` at stabilisation
- key flows at stabilisation (`C_to_I`, `I_to_X`, `X_to_I`, `X_to_S`, etc.)

### Plots
- state trajectories (line and stacked area)
- key daily flows over time
- cumulative cases (total/identified/missed split)
- scenario comparison bar charts (totals, missed split, side-room use)
- optional faceted plots across scenarios

### Tables
- annual summaries (post burn-in)
- final cumulative outcomes
- stabilisation summary
- model validity checks (negatives, mass balance indicators)
- scenario comparison including deltas vs base

---

## One quick thing I need from you before I generate the scripts

Please confirm these 3 choices (I can default if you say “use defaults”):

1. **Time horizon / burn-in**: keep `6 years` with `1-year burn-in`?
2. **Hospital population handling**:
   - A) fixed total census (admissions auto-balance exits), or  
   - B) explicit fixed admission rates (`aS, aC, aI`)?
3. **Scenario definitions** (example):
   - Base
   - S1: improved suspected-case detection (`p_suspected ↑`)
   - S2: higher test sensitivity (`sens ↑`)
   - S3: combined detection + sensitivity + faster turnaround (`t_test ↓`)

If you want, I can proceed immediately with sensible defaults and deliver all 3 files in one go.

## @LFell

1 OK; 2A OK; 3 S1 High sensitivity & High specificity S2 Low sensitivity & high specificity S3 High sensitivity and low specificity S4 Low sensitivity low specificity. I might add other scenarios later to change the turnaround time

## Copilot

Perfect — great spec.  
One correction: you asked for “base + 3 scenarios” earlier, but you’ve now defined **4 scenarios (S1–S4)**. I’ll provide **Base + S1 + S2 + S3 + S4** (5 total), with easy expansion later (e.g., turnaround time).

Use the three files below.

```r name=Scripts/Model1A_ODE_Functions.R
# Model 1A ODE Functions ----
# States: S, C, I, X
# S = susceptible
# C = colonised
# I = invasive infection (not isolated)
# X = enhanced IPC/isolation side-room

Model1A_ODE <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    # Total in-hospital census
    N <- S + C + I + X

    # Force of colonisation
    lambda <- beta0 + (betaC * C + betaI * I + betaX * X) / N

    # Antibiotic-adjusted progression C -> I
    alpha_eff <- alpha * abx_mult

    # Dynamic admissions to keep census approximately constant (Option 2A)
    # Replace all exits immediately.
    exits_total <- (dS_rate + mS) * S + (dC_rate + mC) * C + (dI_rate + mI) * I + (dX_rate + mX) * X

    # Allocate admissions by preset proportions and import prevalences
    # p_adm_C = admitted colonised prevalence
    # p_adm_I = admitted infected prevalence
    aI <- exits_total * p_adm_I
    aC <- exits_total * p_adm_C
    aS <- exits_total - aC - aI
    if (aS < 0) aS <- 0

    # ODEs
    dS <- aS + mu * C + delta * X - lambda * S - (dS_rate + mS) * S
    dC <- aC + lambda * S - (alpha_eff + mu + dC_rate + mC) * C
    dI <- aI + alpha_eff * C + pi * X - (gamma + dI_rate + mI) * I
    dX <- gamma * I - (delta + pi + dX_rate + mX) * X

    # Daily flows (for outputs)
    flow_S_to_C <- lambda * S
    flow_C_to_I <- alpha_eff * C
    flow_C_to_S <- mu * C
    flow_I_to_X <- gamma * I
    flow_X_to_S <- delta * X
    flow_X_to_I <- pi * X

    flow_dis_S <- dS_rate * S
    flow_dis_C <- dC_rate * C
    flow_dis_I <- dI_rate * I
    flow_dis_X <- dX_rate * X

    flow_death_S <- mS * S
    flow_death_C <- mC * C
    flow_death_I <- mI * I
    flow_death_X <- mX * X

    # Case definitions
    total_cases <- flow_C_to_I + aI
    identified_cases <- flow_X_to_S + flow_dis_X
    missed_false_negative <- flow_X_to_I
    missed_untested <- p_not_suspected * I
    missed_cases <- missed_untested + missed_false_negative

    # Side-room occupancy (patient-days per day)
    side_room_patient_days <- X

    list(
      c(dS, dC, dI, dX),
      c(
        N = N, lambda = lambda, alpha_eff = alpha_eff,
        aS = aS, aC = aC, aI = aI,

        flow_S_to_C = flow_S_to_C,
        flow_C_to_I = flow_C_to_I,
        flow_C_to_S = flow_C_to_S,
        flow_I_to_X = flow_I_to_X,
        flow_X_to_S = flow_X_to_S,
        flow_X_to_I = flow_X_to_I,

        flow_dis_S = flow_dis_S,
        flow_dis_C = flow_dis_C,
        flow_dis_I = flow_dis_I,
        flow_dis_X = flow_dis_X,

        flow_death_S = flow_death_S,
        flow_death_C = flow_death_C,
        flow_death_I = flow_death_I,
        flow_death_X = flow_death_X,

        total_cases = total_cases,
        identified_cases = identified_cases,
        missed_cases = missed_cases,
        missed_untested = missed_untested,
        missed_false_negative = missed_false_negative,

        side_room_patient_days = side_room_patient_days
      )
    )
  })
}
```

```r name=Scripts/Model1A_ODE_Parameters.R
# Model 1A Parameters ----

options(scipen = 1000)

prob2rate <- function(p, t) {
  -(1 / t) * log(1 - p)
}

# ---- Fixed baseline assumptions ----

# Transmission
beta0 <- 0.00020
betaC <- 0.00035
betaI <- 0.00060
betaX <- 0.00010

# Natural history
alpha <- 0.010   # baseline C -> I
abx_mult <- 1.30 # antibiotic multiplier for progression C -> I
mu <- 0.004      # decolonisation C -> S

# LOS assumptions (days)
losS <- 7
losC <- 10
losI <- 12
losX <- 6

# Discharge rates
dS_rate <- 1 / losS
dC_rate <- 1 / losC
dI_rate <- 1 / losI
dX_rate <- 1 / losX

# Mortality (per-day rates)
mS

## @LFell

S1 is the base

## Copilot

Great — understood: **S1 is the base case**.

Also, I noticed the parameter file got cut off at `mS`.  
If you want, I’ll now send a **complete, runnable set of all 3 files** with:

- **S1 = Base**
- **S2 = Low sensitivity, high specificity**
- **S3 = High sensitivity, low specificity**
- **S4 = Low sensitivity, low specificity**
- fixed census admissions (your option 2A)
- stabilisation states + flows
- missed cases split: **untested** and **false negatives**
- all tables + plots + CSV exports

I’ll keep turnaround time fixed for now, but structure parameters so you can add turnaround scenarios later in one line.

If you reply **“go ahead”**, I’ll post the full finished scripts end-to-end.

## @LFell

go ahead

## Copilot

```r name=Scripts/Model1A_ODE_Functions.R
# Model 1A ODE Functions ----
# States:
# S = susceptible
# C = colonised
# I = infected (not isolated)
# X = enhanced IPC / side-room isolation

Model1A_ODE <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    N <- S + C + I + X

    # FOI
    lambda <- beta0 + (betaC * C + betaI * I + betaX * X) / N

    # Antibiotic-adjusted progression
    alpha_eff <- alpha * abx_mult

    # Fixed census admissions (Option 2A): replace all exits
    exits_total <- (dS_rate + mS) * S +
      (dC_rate + mC) * C +
      (dI_rate + mI) * I +
      (dX_rate + mX) * X

    aI <- exits_total * p_adm_I
    aC <- exits_total * p_adm_C
    aS <- exits_total - aC - aI
    if (aS < 0) aS <- 0

    # ODEs
    dS <- aS + mu * C + delta * X - lambda * S - (dS_rate + mS) * S
    dC <- aC + lambda * S - (alpha_eff + mu + dC_rate + mC) * C
    dI <- aI + alpha_eff * C + pi * X - (gamma + dI_rate + mI) * I
    dX <- gamma * I - (delta + pi + dX_rate + mX) * X

    # Daily flows
    flow_S_to_C <- lambda * S
    flow_C_to_I <- alpha_eff * C
    flow_C_to_S <- mu * C
    flow_I_to_X <- gamma * I
    flow_X_to_S <- delta * X
    flow_X_to_I <- pi * X

    flow_dis_S <- dS_rate * S
    flow_dis_C <- dC_rate * C
    flow_dis_I <- dI_rate * I
    flow_dis_X <- dX_rate * X

    flow_death_S <- mS * S
    flow_death_C <- mC * C
    flow_death_I <- mI * I
    flow_death_X <- mX * X

    # Cases
    total_cases <- flow_C_to_I + aI
    identified_cases <- flow_X_to_S + flow_dis_X
    missed_false_negative <- flow_X_to_I
    missed_untested <- p_not_suspected * I
    missed_cases <- missed_untested + missed_false_negative

    # Side-room use
    side_room_patient_days <- X

    # Return
    list(
      c(dS, dC, dI, dX),
      c(
        N = N, lambda = lambda, alpha_eff = alpha_eff,
        aS = aS, aC = aC, aI = aI,

        flow_S_to_C = flow_S_to_C,
        flow_C_to_I = flow_C_to_I,
        flow_C_to_S = flow_C_to_S,
        flow_I_to_X = flow_I_to_X,
        flow_X_to_S = flow_X_to_S,
        flow_X_to_I = flow_X_to_I,

        flow_dis_S = flow_dis_S,
        flow_dis_C = flow_dis_C,
        flow_dis_I = flow_dis_I,
        flow_dis_X = flow_dis_X,

        flow_death_S = flow_death_S,
        flow_death_C = flow_death_C,
        flow_death_I = flow_death_I,
        flow_death_X = flow_death_X,

        total_cases = total_cases,
        identified_cases = identified_cases,
        missed_cases = missed_cases,
        missed_untested = missed_untested,
        missed_false_negative = missed_false_negative,

        side_room_patient_days = side_room_patient_days
      )
    )
  })
}
```

```r name=Scripts/Model1A_ODE_Parameters.R
# Model 1A Parameters ----
options(scipen = 1000)

# Helper
prob2rate <- function(p, t) {
  -(1 / t) * log(1 - p)
}

# ---- Fixed assumptions ----
# Transmission
beta0 <- 0.00020
betaC <- 0.00035
betaI <- 0.00060
betaX <- 0.00010

# Natural history
alpha <- 0.010
abx_mult <- 1.30
mu <- 0.004

# Length of stay (days)
losS <- 7
losC <- 10
losI <- 12
losX <- 6

# Discharge rates
dS_rate <- 1 / losS
dC_rate <- 1 / losC
dI_rate <- 1 / losI
dX_rate <- 1 / losX

# Mortality rates (per day)
mS <- prob2rate(0.02, losS)
mC <- mS
mI <- prob2rate(0.20, losI)
mX <- mI

# Admission prevalence
p_adm_C <- 0.03
p_adm_I <- 0.002

# Suspected-case identification process (fixed across scenarios unless changed)
p_suspected <- 0.85
t_identify_isolate <- 1.0  # days
gamma_base <- p_suspected / t_identify_isolate
p_not_suspected <- 1 - p_suspected

# Recovery process components
p_treated_recover_to_S <- 0.80
t_test_treat_recover <- 5.0  # days

# False negative return timing
t_test_return_false_neg <- 1.0  # days

# Scenario builder from Se/Sp
# Se affects delta and pi
# Sp included for outputs/extension (e.g., false positives from other structures)
make_mod_parms <- function(
    sens,
    spec,
    scenario_label,
    turnaround_mult = 1
) {
  stopifnot(length(sens) == 1, sens >= 0, sens <= 1)
  stopifnot(length(spec) == 1, spec >= 0, spec <= 1)
  stopifnot(turnaround_mult > 0)

  gamma <- gamma_base
  delta <- (sens * p_treated_recover_to_S) / (t_test_treat_recover * turnaround_mult)
  pi <- (1 - sens) / (t_test_return_false_neg * turnaround_mult)

  list(
    scenario = scenario_label,

    beta0 = beta0, betaC = betaC, betaI = betaI, betaX = betaX,
    alpha = alpha, abx_mult = abx_mult, mu = mu,

    dS_rate = dS_rate, dC_rate = dC_rate, dI_rate = dI_rate, dX_rate = dX_rate,
    mS = mS, mC = mC, mI = mI, mX = mX,

    p_adm_C = p_adm_C, p_adm_I = p_adm_I,

    p_suspected = p_suspected,
    p_not_suspected = p_not_suspected,

    sens = sens, spec = spec,
    gamma = gamma, delta = delta, pi = pi
  )
}

# ---- Scenarios ----
# S1 is base case (as requested)
mod.parms.S1 <- make_mod_parms(
  sens = 0.95, spec = 0.95, scenario_label = "S1_base_highSe_highSp"
)

mod.parms.S2 <- make_mod_parms(
  sens = 0.70, spec = 0.95, scenario_label = "S2_lowSe_highSp"
)

mod.parms.S3 <- make_mod_parms(
  sens = 0.95, spec = 0.70, scenario_label = "S3_highSe_lowSp"
)

mod.parms.S4 <- make_mod_parms(
  sens = 0.70, spec = 0.70, scenario_label = "S4_lowSe_lowSp"
)
```

````markdown name=Model1A_ODE_Model_Run.qmd
---
title: "Model 1A ODE Run"
format:
  html:
    toc: true
    code-fold: false
execute:
  echo: true
  warning: false
  message: false
---

## 1) Setup

```{r}
rm(list = ls())

library(tidyverse)
library(deSolve)
library(scales)

options(scipen = 1000)

dir.create("Outputs", showWarnings = FALSE, recursive = TRUE)

source("Scripts/Model1A_ODE_Functions.R")
source("Scripts/Model1A_ODE_Parameters.R")
```

## 2) Model settings

```{r}
model.time <- 365 * 6
time.step <- 1
burn_in_days <- 365
mod.t <- seq(0, model.time, by = time.step)
```

## 3) Initial conditions

```{r}
S0 <- 620
C0 <- 15
I0 <- 5
X0 <- 0

mod.init <- c(S = S0, C = C0, I = I0, X = X0)
```

## 4) Utility functions

```{r}
run_scenario <- function(parms) {
  as.data.frame(ode(y = mod.init, times = mod.t, func = Model1A_ODE, parms = parms))
}

apply_burn_in <- function(df, burn_days = 365) {
  df %>% filter(time > burn_days) %>% mutate(time = time - burn_days)
}

add_derived_outputs <- function(df) {
  df %>%
    mutate(year = if_else(time == 0, 1, ((time - 1) %/% 365) + 1)) %>%
    relocate(year, .before = time) %>%
    mutate(
      flow_dis_total = flow_dis_S + flow_dis_C + flow_dis_I + flow_dis_X,
      flow_death_total = flow_death_S + flow_death_C + flow_death_I + flow_death_X
    ) %>%
    group_by(year) %>%
    mutate(
      cum_total_cases = cumsum(total_cases),
      cum_identified_cases = cumsum(identified_cases),
      cum_missed_cases = cumsum(missed_cases),
      cum_missed_untested = cumsum(missed_untested),
      cum_missed_false_negative = cumsum(missed_false_negative),
      cum_side_room_patient_days = cumsum(side_room_patient_days),
      cum_S_to_C = cumsum(flow_S_to_C),
      cum_C_to_I = cumsum(flow_C_to_I),
      cum_I_to_X = cumsum(flow_I_to_X),
      cum_X_to_S = cumsum(flow_X_to_S),
      cum_X_to_I = cumsum(flow_X_to_I),
      cum_dis_total = cumsum(flow_dis_total),
      cum_death_total = cumsum(flow_death_total)
    ) %>%
    ungroup()
}

negative_value_rows <- function(df) {
  df %>% filter(if_any(where(is.numeric), ~ .x < 0))
}

negative_value_summary <- function(df) {
  df %>%
    summarise(across(where(is.numeric), ~ sum(.x < 0, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "num_negative_values") %>%
    filter(num_negative_values > 0)
}

find_stabilisation_day <- function(df, states, tol = 1e-4, consecutive_days = 30) {
  stopifnot(all(states %in% names(df)))

  x <- df %>% arrange(time) %>% select(time, all_of(states))
  d <- x %>% mutate(across(all_of(states), ~ abs(.x - lag(.x))))
  stable_flag <- d %>%
    mutate(stable = if_all(all_of(states), ~ !is.na(.x) & .x <= tol)) %>%
    pull(stable)

  r <- rle(stable_flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  idx <- which(r$values & r$lengths >= consecutive_days)[1]

  if (is.na(idx)) return(NA_real_)
  x$time[starts[idx]]
}

stabilisation_summary <- function(df, scenario_label,
                                  state_vars = c("S","C","I","X"),
                                  flow_vars = c("flow_S_to_C","flow_C_to_I","flow_I_to_X","flow_X_to_S","flow_X_to_I",
                                                "total_cases","identified_cases","missed_cases","missed_untested","missed_false_negative"),
                                  tol = 1e-4, consecutive_days = 30) {
  t_star <- find_stabilisation_day(df, state_vars, tol = tol, consecutive_days = consecutive_days)

  if (is.na(t_star)) {
    return(tibble(
      scenario = scenario_label,
      stabilisation_day = NA_real_,
      tolerance = tol,
      consecutive_days = consecutive_days,
      stabilised = FALSE
    ))
  }

  row_star <- df %>%
    filter(time == t_star) %>%
    slice(1) %>%
    select(all_of(c(state_vars, flow_vars)))

  tibble(
    scenario = scenario_label,
    stabilisation_day = t_star,
    tolerance = tol,
    consecutive_days = consecutive_days,
    stabilised = TRUE
  ) %>% bind_cols(row_star)
}
```

## 5) Run scenarios

```{r}
scenario_parms <- list(
  S1_base = mod.parms.S1,
  S2_lowSe_highSp = mod.parms.S2,
  S3_highSe_lowSp = mod.parms.S3,
  S4_lowSe_lowSp = mod.parms.S4
)

scenario_raw <- lapply(scenario_parms, run_scenario)
scenario_abso <- lapply(scenario_raw, function(x) x %>% apply_burn_in(burn_in_days) %>% add_derived_outputs())

walk2(scenario_abso, names(scenario_abso), ~ write_csv(.x, paste0("Outputs/solution_", .y, ".csv")))
```

## 6) Annual summary tables

```{r}
make_annual_summary <- function(df, scenario_label) {
  df %>%
    filter(time %% 365 == 0 & time != 0) %>%
    mutate(year_after_burnin = time / 365) %>%
    transmute(
      scenario = scenario_label,
      year_after_burnin,
      S, C, I, X,
      cum_total_cases,
      cum_identified_cases,
      cum_missed_cases,
      cum_missed_untested,
      cum_missed_false_negative,
      missed_share_untested = if_else(cum_missed_cases > 0, cum_missed_untested / cum_missed_cases, NA_real_),
      missed_share_falseneg = if_else(cum_missed_cases > 0, cum_missed_false_negative / cum_missed_cases, NA_real_),
      cum_side_room_patient_days
    )
}

annual_summary_all <- bind_rows(
  make_annual_summary(scenario_abso$S1_base, "S1_base"),
  make_annual_summary(scenario_abso$S2_lowSe_highSp, "S2_lowSe_highSp"),
  make_annual_summary(scenario_abso$S3_highSe_lowSp, "S3_highSe_lowSp"),
  make_annual_summary(scenario_abso$S4_lowSe_lowSp, "S4_lowSe_lowSp")
)

write_csv(annual_summary_all, "Outputs/annual_summary_all.csv")
knitr::kable(annual_summary_all, caption = "Annual summary (post burn-in)")
```

## 7) Stabilisation summary (states + flows)

```{r}
stability_all <- bind_rows(
  stabilisation_summary(scenario_abso$S1_base, "S1_base"),
  stabilisation_summary(scenario_abso$S2_lowSe_highSp, "S2_lowSe_highSp"),
  stabilisation_summary(scenario_abso$S3_highSe_lowSp, "S3_highSe_lowSp"),
  stabilisation_summary(scenario_abso$S4_lowSe_lowSp, "S4_lowSe_lowSp")
)

write_csv(stability_all, "Outputs/stabilisation_summary_states_flows.csv")
knitr::kable(stability_all, caption = "Stabilisation: model states and flows")
```

## 8) Negative values checks

```{r}
neg_tables <- lapply(scenario_abso, negative_value_rows)
neg_summaries <- lapply(scenario_abso, negative_value_summary)

walk2(neg_tables, names(neg_tables), ~ write_csv(.x, paste0("Outputs/negative_rows_", .y, ".csv")))
walk2(neg_summaries, names(neg_summaries), ~ write_csv(.x, paste0("Outputs/negative_summary_", .y, ".csv")))
```

## 9) Final cumulative outcomes table

```{r}
final_outcomes <- bind_rows(lapply(names(scenario_abso), function(nm) {
  d <- scenario_abso[[nm]]
  last <- d %>% filter(time == max(time)) %>% slice(1)

  tibble(
    scenario = nm,
    final_S = last$S,
    final_C = last$C,
    final_I = last$I,
    final_X = last$X,
    cum_total_cases = last$cum_total_cases,
    cum_identified_cases = last$cum_identified_cases,
    cum_missed_cases = last$cum_missed_cases,
    cum_missed_untested = last$cum_missed_untested,
    cum_missed_false_negative = last$cum_missed_false_negative,
    cum_side_room_patient_days = last$cum_side_room_patient_days
  )
}))

# Incrementals vs base S1
base_row <- final_outcomes %>% filter(scenario == "S1_base")
final_outcomes_vs_base <- final_outcomes %>%
  mutate(
    d_total_cases_vs_S1 = cum_total_cases - base_row$cum_total_cases,
    d_identified_vs_S1 = cum_identified_cases - base_row$cum_identified_cases,
    d_missed_vs_S1 = cum_missed_cases - base_row$cum_missed_cases,
    d_side_room_days_vs_S1 = cum_side_room_patient_days - base_row$cum_side_room_patient_days
  )

write_csv(final_outcomes_vs_base, "Outputs/final_outcomes_vs_base.csv")
knitr::kable(final_outcomes_vs_base, caption = "Final cumulative outcomes (with deltas vs S1 base)")
```

## 10) Plots

```{r}
plot_long_states <- bind_rows(lapply(names(scenario_abso), function(nm) {
  scenario_abso[[nm]] %>%
    select(time, S, C, I, X) %>%
    pivot_longer(cols = c(S,C,I,X), names_to = "state", values_to = "count") %>%
    mutate(scenario = nm)
}))

# 10.1 State lines by scenario (faceted)
p_states_line <- ggplot(plot_long_states, aes(time, count, color = state)) +
  geom_line() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "State trajectories by scenario", x = "Days post burn-in", y = "Count")
ggsave("Outputs/plot_states_lines_facet.png", p_states_line, width = 12, height = 8, dpi = 300)

# 10.2 State stacked area by scenario
p_states_area <- ggplot(plot_long_states, aes(time, count, fill = state)) +
  geom_area() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "State composition over time", x = "Days post burn-in", y = "Count")
ggsave("Outputs/plot_states_area_facet.png", p_states_area, width = 12, height = 8, dpi = 300)

# Flows for plotting
plot_long_flows <- bind_rows(lapply(names(scenario_abso), function(nm) {
  scenario_abso[[nm]] %>%
    select(time, flow_S_to_C, flow_C_to_I, flow_I_to_X, flow_X_to_S, flow_X_to_I,
           total_cases, identified_cases, missed_cases, missed_untested, missed_false_negative,
           side_room_patient_days,
           cum_total_cases, cum_identified_cases, cum_missed_cases,
           cum_missed_untested, cum_missed_false_negative, cum_side_room_patient_days) %>%
    mutate(scenario = nm)
}))

# 10.3 Key daily flows
p_daily_flows <- plot_long_flows %>%
  select(time, scenario, flow_S_to_C, flow_C_to_I, flow_I_to_X, flow_X_to_S, flow_X_to_I) %>%
  pivot_longer(cols = -c(time, scenario), names_to = "flow", values_to = "value") %>%
  ggplot(aes(time, value, color = flow)) +
  geom_line() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "Key daily flows", x = "Days post burn-in", y = "Flow/day")
ggsave("Outputs/plot_daily_flows_facet.png", p_daily_flows, width = 12, height = 8, dpi = 300)

# 10.4 Cumulative total/identified/missed
p_cum_cases <- plot_long_flows %>%
  select(time, scenario, cum_total_cases, cum_identified_cases, cum_missed_cases) %>%
  pivot_longer(cols = -c(time, scenario), names_to = "series", values_to = "value") %>%
  ggplot(aes(time, value, color = series)) +
  geom_line(size = 1) +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "Cumulative cases", x = "Days post burn-in", y = "Cumulative count")
ggsave("Outputs/plot_cumulative_cases_facet.png", p_cum_cases, width = 12, height = 8, dpi = 300)

# 10.5 Missed cases breakdown
p_missed_breakdown <- plot_long_flows %>%
  select(time, scenario, cum_missed_untested, cum_missed_false_negative) %>%
  pivot_longer(cols = -c(time, scenario), names_to = "missed_type", values_to = "value") %>%
  ggplot(aes(time, value, color = missed_type)) +
  geom_line(size = 1) +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "Cumulative missed cases breakdown", x = "Days post burn-in", y = "Cumulative count")
ggsave("Outputs/plot_cumulative_missed_breakdown_facet.png", p_missed_breakdown, width = 12, height = 8, dpi = 300)

# 10.6 Side-room patient-days over time
p_side_room <- plot_long_flows %>%
  ggplot(aes(time, side_room_patient_days, color = scenario)) +
  geom_line(size = 1) +
  theme_classic() +
  labs(title = "Side-room occupancy over time", x = "Days post burn-in", y = "Patient-days/day")
ggsave("Outputs/plot_side_room_occupancy.png", p_side_room, width = 10, height = 6, dpi = 300)

# 10.7 Final outcomes comparison bars
final_long <- final_outcomes_vs_base %>%
  select(scenario, cum_total_cases, cum_identified_cases, cum_missed_cases, cum_side_room_patient_days) %>%
  pivot_longer(-scenario, names_to = "metric", values_to = "value")

p_final_bars <- ggplot(final_long, aes(scenario, value, fill = scenario)) +
  geom_col() +
  facet_wrap(~metric, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "Final cumulative outcomes by scenario", x = "", y = "Value")
ggsave("Outputs/plot_final_outcomes_bars.png", p_final_bars, width = 12, height = 8, dpi = 300)

p_states_line
p_states_area
p_daily_flows
p_cum_cases
p_missed_breakdown
p_side_room
p_final_bars
```

## 11) Save parameter table

```{r}
all_parms <- bind_rows(
  S1_base = as_tibble(mod.parms.S1),
  S2_lowSe_highSp = as_tibble(mod.parms.S2),
  S3_highSe_lowSp = as_tibble(mod.parms.S3),
  S4_lowSe_lowSp = as_tibble(mod.parms.S4),
  .id = "scenario"
)
write_csv(all_parms, "Outputs/model_parameters_all_scenarios.csv")
knitr::kable(all_parms, caption = "Model parameters for all scenarios")
```
````

This is ready to run as-is.  
If you want, next I can give you a **“turnaround-time extension block”** so you can add `S5/S6` by only editing `turnaround_mult`.

## @LFell

Can you check this? # Natural history
alpha0 <- 0.010  # Proportion of colonised that progress from colonised to infected (per day)
abx_mult <- 1.30  # Multiplier for progression due to antibiotics (relative risk)
abx_prop <- 0.20 # Proportion of colonised patients receiving antibiotics (per day)
alpha <- alpha0 * (1 + (abx_mult - 1) * abx_prop) # Adjusted progression rate from colonised to infected (per day)
alpha <- alpha0 * ((1 - abx_prop) + (abx_prop * abx_mult)) # Adjusted progression rate from colonised to infected (per day)
alpha <- 0.010 * ((1 - 0.20) + (0.2*1.3)) # Check calculation of alpha

## Copilot

Yes — this is correct, and your two formulas are algebraically identical.

- `alpha <- alpha0 * (1 + (abx_mult - 1) * abx_prop)`
- `alpha <- alpha0 * ((1 - abx_prop) + (abx_prop * abx_mult))`

Both give:

\[
\alpha = 0.010 \times ((1-0.20) + 0.20 \times 1.30)
= 0.010 \times (0.80 + 0.26)
= 0.010 \times 1.06
= 0.0106
\]

So adjusted progression is **`alpha = 0.0106` per day**.

### One cleanup suggestion
In code, keep only one assignment to avoid accidental overwrite confusion:

```r name=alpha_adjustment.R
alpha0 <- 0.010
abx_mult <- 1.30
abx_prop <- 0.20

alpha <- alpha0 * (1 + (abx_mult - 1) * abx_prop)  # 0.0106
```

And in the ODE, use this `alpha` directly (or call it `alpha_eff` once in parameters, not both places).

## @LFell

Does this look correct so far? 
# Helper function to convert probability of event over time t to a daily rate
prob2rate <- function(p, t) {
  -(1 / t) * log(1 - p)
}

# ---- Fixed assumptions ----
# Transmission
beta0 <- 0.00020 # Background transmission (per day)
betaC <- 0.00035 # Colonised transmission (per day)
betaI <- 0.00060 # Infected transmission (per day)
betaX <- 0.00010 # Isolated transmission (per day)

# Natural history

# Colonisation to infection progression

alpha0 <- 0.010  # Proportion of colonised that progress from colonised to infected 
abx_mult <- 1.30  # Multiplier for progression due to antibiotics (relative risk)
abx_prop <- 0.20 # Proportion of colonised patients receiving antibiotics 
t_progress_C_to_I <- 10.0 # Time from colonisation to infection (days)
alpha <- (alpha0 * ((1 - abx_prop) + (abx_prop * abx_mult)) / t_progress_to_I) # Rate of progression from colonised to infected (per day)

# Recovery from colonisation to susceptible without treatment

mu0 <- 0.004 ## Proportion that recover from colonisation to susceptible without treatment
t_recover_C_to_S <- 30.0 # Time from colonisation to natural recovery (days)
mu <- mu0 / t_recover_C_to_S # Rate of recovery from colonised to susceptible without treatment (per day)

# Length of stay (days)
losS <- 7 # Average length of stay for susceptible patients
losC <- 10 # Average length of stay for colonised patients
losI <- 12 # Average length of stay for infected patients
losX <- 6 # Average length of stay for isolated patients

# Discharge rates
dS <- 1 / losS # Discharge rate for susceptible patients (per day)
dC <- 1 / losC # Discharge rate for colonised patients (per day)
dI <- 1 / losI # Discharge rate for infected patients (per day)
dX <- 1 / losX # Discharge rate for isolated patients (per day)

# Mortality rates (per day)
mS <- prob2rate(0.02, losS) # Mortality rate for susceptible patients (per day) 
mC <- prob2rate(0.02, losC) # Mortality rate for colonised patients (per day) - same as susceptible
mI <- prob2rate(0.20, losI) # Mortality rate for infected patients (per day # 
mX <- prob2rate(0.20, losX) # Mortality rate for isolated patients (per day) - same as infected

# Prevalence of C and I on admission to hospital
prop_adm_C <- 0.03 # Proportion of admissions that are colonised
prop_adm_I <- 0.002 # Proportion of admissions that are infected

# Suspected-case identification process (fixed across scenarios unless changed)
p_I_suspected <- 0.85 # Proportion of infected patients that are suspected and thus isolated and tested 
t_identify_isolate_I_to_X <- 1  # Time from infection onset to identification and isolation (days)
gamma <- p_suspected / t_identify_isolate # Isolation rate for infected patients (per day)
p_not_suspected <- 1 - p_suspected # Proportion of infected patients that are not suspected and thus remain in general wards (per day)


## Copilot

Almost there — good structure, but there are a few **important variable-name and rate/proportion issues** to fix.

## Key fixes

1. **Name mismatch / typo**
   - You define `t_progress_C_to_I` but use `t_progress_to_I`.
2. **`alpha0` meaning is inconsistent**
   - If `alpha0` is already a **rate per day**, don’t divide by time again.
   - If `alpha0` is a **proportion over a period**, convert with `prob2rate(alpha0, t_progress_C_to_I)`.
3. **Same issue for `mu0`**
   - If `mu0` is proportion over 30 days, use `prob2rate(mu0, 30)`, not `mu0 / 30`.
4. **Another name mismatch**
   - You define `p_I_suspected` and `t_identify_isolate_I_to_X` but use `p_suspected` and `t_identify_isolate`.
5. **Comment typo**
   - `mI` line has an unclosed comment parenthesis.

---

Here is a corrected version, keeping your intent and naming consistent:

```r name=parameters_snippet_corrected.R
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
```

If you want, I can now do a **full consistency pass** on your entire parameter file so every symbol matches the ODE function exactly.

## @LFell

can i savethis chat?

## Copilot

Yes — definitely.

Quick ways:

1. **Copy/paste** into a document (Word/Google Docs/Notes).
2. **Download/export chat** (if your interface has an export/share option).
3. **Save code blocks directly** into files:
   - `Scripts/Model1A_ODE_Functions.R`
   - `Scripts/Model1A_ODE_Parameters.R`
   - `Model1A_ODE_Model_Run.qmd`

If you want, I can also give you a **single consolidated “final version”** of all scripts in one message so it’s easier to save.
