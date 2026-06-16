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

## @LFell

I am not sure when it is best to code a rate as proportion/duration and when to use the prob2rate function

## Copilot

Great question. Short rule:

- Use **`1 / duration`** when you’re modeling a **continuous transition with constant hazard** and you only know the **mean time** in a state.
- Use **`prob2rate(p, t)`** when you know a **cumulative probability** over a period \(t\) (e.g., “20% die within 12 days”).

### Why
`prob2rate(p,t)` gives hazard \(r\) such that:
\[
p = 1 - e^{-rt}
\]
So \(r = -\ln(1-p)/t\).

`p/t` is only a rough approximation (okay only when \(p\) is very small).

### In your model
- **LOS/discharge**: usually `d = 1/LOS` ✅
- **If paper says “X% mortality over LOS days”**: use `prob2rate(X, LOS)` ✅
- **If you have “probability of progression over T days”**: use `prob2rate(p, T)` ✅
- **If parameter is already a daily rate** from literature: use it directly ✅

### Your specific snippet
If `alpha0 <- 0.010` means **daily rate**, do **not** divide by duration again.  
If it means “1% progress over 10 days,” then:
```r
alpha_base <- prob2rate(0.01, 10)
```
then apply your antibiotic multiplier.

If you want, I can label each parameter in your file as “rate” vs “probability over time” so it’s unambiguous.

## @LFell

so prob rate is good when the event could happen at any time during the x days and 1/duration is good when its ok to assume a smooth distribution?

## Copilot

Exactly — that’s a good intuition.  

Slightly sharper wording:

- **`prob2rate(p, t)`**: use when you know **risk over an interval** (“by day \(t\), probability \(p\)”). It converts that interval risk into a continuous-time hazard.
- **`1/duration`**: use when you know the **mean waiting time** and assume an **exponential (memoryless) process** with a constant hazard.

In fact, `1/duration` is also a constant-hazard assumption — just parameterized by mean time instead of interval probability.

So they’re consistent frameworks, just different inputs:
- input = probability over time ⟶ `prob2rate`
- input = mean time ⟶ `1/mean_time`

Your “smooth distribution” phrasing is close; technically it’s “constant hazard / exponential waiting time.”

## @LFell

in R, Can a qmd call and run a qmd to load functions and parameter values, or does it need to be a script

## Copilot

Yes — a `.qmd` **can** run another `.qmd`, but for loading functions/parameters it’s usually better to use `.R` scripts.

## Best practice
Use:

- `Scripts/functions.R`
- `Scripts/parameters.R`

and in your main `.qmd`:

```r
source("Scripts/functions.R")
source("Scripts/parameters.R")
```

## Can a qmd call another qmd?
Yes, via rendering:

```r
quarto::quarto_render("other_file.qmd")
```

But that **renders a document**, not just “loads objects into current session” in a clean/simple way. It’s heavier and less predictable for modular model code.

## Recommendation for your model
Keep reusable code in `.R` files; keep `.qmd` as orchestration + results/tables/plots.

## @LFell

Once in isolation, there are four competing risks: a false negative test result (X to I), a positive test result + treatment + recovery + return to a general ward within the hospital (X to S), a positive test result + treatment + recovery + discharge from hospital (X to d) and death. The length of stay in X is dependent on the test turnaround time so needs to be included in the scenario builder if test turnaround is to be varied. The time to test, treat and recover is the same for X to S and X to d. So I am wondering if it would be better to have two compartments Xtest and Xtreat, so I can have a mean length of stay for X treat and a constant mortality assumption, and divide X to S and X to d using a fixed proportion. Would that be the simplest way to do it or are there other solutions that modellers use?

## Copilot

Yes — **splitting into `Xtest` and `Xtreat` is the cleanest and most standard approach** for what you want.

## Recommendation
Use:

- `Xtest`: in side room awaiting test result
- `Xtreat`: confirmed positive, on treatment/recovery pathway in side room

Then model:

- `I -> Xtest` (suspected cases isolated)
- `Xtest -> I` (false negatives)
- `Xtest -> Xtreat` (true positives)
- `Xtreat -> S` (recover and return ward)
- `Xtreat -> discharge` (recover and leave hospital)
- mortality from both `Xtest` and `Xtreat` (or one shared rate if you want simplicity)

This solves your exact issue:
- **turnaround time** affects only `Xtest` dwell time
- **treat/recover duration** affects only `Xtreat`
- you can split `Xtreat` exits to `S` vs discharge with a fixed proportion

## Why this is better than one X compartment
With one `X`, you’re forcing several processes (testing, treatment, discharge decision) into one lumped rate set, which gets messy when turnaround changes by scenario.

## Common alternatives modellers use
1. **Single X with scenario-specific effective rates** (simpler but less transparent)
2. **Delay differential equations** for fixed delays (more complex)
3. **Stage/Erlang compartments** (`Xtest1->Xtest2...`) if you want less variable waiting times

For your project, **`Xtest + Xtreat` is the best simplicity/fidelity balance**.

If you want, I can refactor your current model equations and scenario builder to this 2-isolation-compartment structure next.

## @LFell

I have modified the parameter script, but need some help with the scenario builder   #### Model 1A Parameters ----

# Introduction ----

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

# Initial setup ----

# Use options to avoid scientific notation for small numbers (e.g., transmission rates)
options(scipen = 1000)

# Helper function to convert probability of event over a time interval (cumulative probability or interval risk) into 
# a continuous-time hazard (e.g. daily rate), assuming an exponential distribution of time to event.
prob2rate <- function(p, t) {
  -(1 / t) * log(1 - p)
}
# This means the rate is the negative of the natural log of (1 - probability) divided by the time interval.

# ========================================================== #

# ---- Fixed assumptions ----

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

# Recovery of confirmed cases (X to S and X to discharge)
prop_Xtreat_to_S <- 0.80  # Proportion of treated confirmed cases that recover to susceptible state and return to a general ward (X to S)
delta <- prop_Xtreat_to_S / losXtreat # Rate of recovery for confirmed cases and return to a general ward (X to S, per half-day)

# Discharge rates (per half-day)
dS <- 1 / losS # Discharge rate for susceptible patients (per half-day)
dC <- 1 / losC # Discharge rate for colonised patients (per half-day)
dI <- 1 / losI # Discharge rate for infected patients (per half-day)
dXtreat <- -prop_Xtreat_to_S / losXtreat # Discharge rate for treated patients in isolation (per half-day)
# Assume no discharge from Xtest as patients are awaiting test results and thus remain in hospital 
# until they are either treated or returned to a general ward.

# Mortality rates (per half-day)
mS <- prob2rate(0.02, losS) # Mortality rate for susceptible patients (per half-day) 
mC <- prob2rate(0.02, losC) # Mortality rate for colonised patients (per half-day) - same as susceptible
mI <- prob2rate(0.20, losI) # Mortality rate for infected patients (per half-day)
mXtreat <- prob2rate(0.10, losXtreat) # Mortality rate for confirmed cases in isolation (per half-day)
# Length of stay in Xtest is determined by the test turnaround time, so is not fixed across scenarios.

# ========================================================== #

# Function to build list of parameters for ODE model for each scenario ----

make_mod_parms <- function(
    prop_I_suspected,
    test_sens,
    test_turnaround,
    scenario_label
) {
  stopifnot(length(sens) == 1, sens >= 0, sens <= 1)
  stopifnot(length(spec) == 1, spec >= 0, spec <= 1)
  stopifnot(test_turnaround > 0)
  
  # Suspected-case identification (I to Xtest)
  prop_I_suspected <- 0.50 # Proportion of infected patients that are suspected and thus isolated and tested 
  t_isolate_I_to_Xtest <- 4  # Time from infection onset to identification and isolation (half-days)
  gamma <- prop_I_suspected / t_isolate_I_to_Xtest # Isolation rate for suspected cases (per half-day)
  
  # Testing (Xtest to I and Xtest to Xtreat)
  test_sens <- 0.95 # Test sensitivity (proportion of true positives correctly identified)
  test_turnaround <- 4 # Time to test result and return to a general ward or remain in isolation for treatment (half-days)
  theta <- test_sens / test_turnaround # Rate of progression to Xtreat for true positives (per half-day)
  pi <- (1 - test_sens) / test_turnaround # Rate of false negative return from isolation to general ward (per half-day)
  mXtest <- prob2rate(0.05, test_turnaround) # Mortality rate suspected cases in isolation awaiting test results (per half-day)
  
  list(
    scenario = scenario_label,
    
    beta0 = beta0, betaC = betaC, betaI = betaI, betaXtest = betaXtest, betaXtreat = betaXtreat,
    alpha = alpha, abx_prop = abx_prop, abx_mult = abx_mult, mu = mu,
    
    dS_rate = dS_rate, dC_rate = dC_rate, dI_rate = dI_rate, dXtreat_rate = dXtreat_rate,
    mS = mS, mC = mC, mI = mI, mXtest = mXtest, mXtreat = mXtreat,
    
    p_adm_C = p_adm_C, p_adm_I = p_adm_I,
    
    prop_I_suspected = prop_I_suspected,
    
    sens = sens, spec = spec,
    gamma = gamma, delta = delta, pi = pi, theta = theta
  )
}

# =========================================================== #

# Scenarios ----
# S1 is base case 
mod.parms.S1 <- make_mod_parms(
  prop_I_suspected = 0.50, test_sens = 0.80, test_turnaround = 5, scenario_label = "S1_midSus_highSe_slow"
)

mod.parms.S2 <- make_mod_parms(
  prop_I_suspected = 0.50, test_sens = 0.50, test_turnaround = 5, scenario_label = "S2_midSus_lowSe_slow"
  )

mod.parms.S3 <- make_mod_parms(
  prop_I_suspected = 0.50, test_sens = 0.80, test_turnaround = 3, scenario_label = "S3_midSus_highSe_fast"
)

mod.parms.S4 <- make_mod_parms(
   prop_I_suspected = 0.50, test_sens = 0.50, test_turnaround = 3, scenario_label = "S4_midSus_low_Se_fast"
  
)

## Copilot

You’re very close. Main issue: your scenario builder currently **overwrites inputs** and has several **name mismatches** (`sens/spec`, `dS_rate`, etc.). Also `dXtreat` is incorrectly negative.

Here’s a corrected builder section you can paste in.

```r name=make_mod_parms_fixed.R
# Function to build parameter list for each scenario ----
make_mod_parms <- function(
  prop_I_suspected,
  test_sens,
  test_turnaround,
  scenario_label
) {
  # Validate inputs
  stopifnot(length(prop_I_suspected) == 1, prop_I_suspected >= 0, prop_I_suspected <= 1)
  stopifnot(length(test_sens) == 1, test_sens >= 0, test_sens <= 1)
  stopifnot(length(test_turnaround) == 1, test_turnaround > 0)

  # I -> Xtest
  t_isolate_I_to_Xtest <- 4
  gamma <- prop_I_suspected / t_isolate_I_to_Xtest

  # Xtest split by sensitivity over turnaround time
  theta <- test_sens / test_turnaround        # Xtest -> Xtreat (true positive)
  pi <- (1 - test_sens) / test_turnaround     # Xtest -> I (false negative)

  # Mortality in Xtest (depends on duration there)
  mXtest <- prob2rate(0.05, test_turnaround)

  list(
    scenario = scenario_label,

    # Admissions/transmission/natural history
    prop_adm_C = prop_adm_C, prop_adm_I = prop_adm_I,
    beta0 = beta0, betaC = betaC, betaI = betaI,
    betaXtest = betaXtest, betaXtreat = betaXtreat,
    alpha = alpha, abx_prop = abx_prop, abx_mult = abx_mult, mu = mu,

    # Discharge rates
    dS = dS, dC = dC, dI = dI, dXtreat = dXtreat,

    # Mortality rates
    mS = mS, mC = mC, mI = mI, mXtest = mXtest, mXtreat = mXtreat,

    # Scenario controls
    prop_I_suspected = prop_I_suspected,
    test_sens = test_sens,
    test_turnaround = test_turnaround,

    # Transition rates
    gamma = gamma,       # I -> Xtest
    theta = theta,       # Xtest -> Xtreat
    pi = pi,             # Xtest -> I
    delta = delta        # Xtreat -> S
  )
}
```

Also fix these earlier lines in your parameter file:

```r name=key_fixes.R
# discharge from Xtreat should be positive
dXtreat <- (1 - prop_Xtreat_to_S) / losXtreat

# keep names consistent with builder/list:
# use prop_adm_C / prop_adm_I everywhere (not p_adm_C / p_adm_I)
```

And remove these from builder (they should NOT be hardcoded inside):
- `prop_I_suspected <- 0.50`
- `test_sens <- 0.95`
- `test_turnaround <- 4`

Those must come from scenario inputs (which you already pass).

If you want, I can next give you the matching **ODE function equations for S/C/I/Xtest/Xtreat/D** with these exact parameter names so it runs first try.

## @LFell

Does this look better?  #### Model 1A Parameters ####

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

# Recovery of confirmed cases (X to S and X to discharge)
prop_Xtreat_to_S <- 0.80  # Proportion of treated confirmed cases that recover to susceptible state and return to a general ward (X to S)
delta <- prop_Xtreat_to_S / losXtreat # Rate of recovery for confirmed cases and return to a general ward (X to S, per half-day)

# Discharge rates (per half-day)
dS <- 1 / losS # Discharge rate for susceptible patients (per half-day)
dC <- 1 / losC # Discharge rate for colonised patients (per half-day)
dI <- 1 / losI # Discharge rate for infected patients (per half-day)
dXtreat <- (1 - prop_Xtreat_to_S) / losXtreat # Discharge rate for treated patients in isolation (per half-day)
# Assume no discharge from Xtest as patients are awaiting test results and thus remain in hospital 
# until they are either treated or returned to a general ward.

# Mortality rates (per half-day)
mS <- prob2rate(0.02, losS) # Mortality rate for susceptible patients (per half-day) 
mC <- prob2rate(0.02, losC) # Mortality rate for colonised patients (per half-day) - same as susceptible
mI <- prob2rate(0.20, losI) # Mortality rate for infected patients (per half-day)
mXtreat <- prob2rate(0.10, losXtreat) # Mortality rate for confirmed cases in isolation (per half-day)
# Length of stay in Xtest is determined by the test turnaround time, so is not fixed across scenarios.

# ========================================================== #

# 4) Create function to build parameter list for each scenario ####

make_mod_parms <- function(
    prop_I_suspected,
    test_sens,
    test_turnaround,
    scenario_label
) {
  # Validate inputs
  stopifnot(length(prop_I_suspected) == 1, prop_I_suspected >= 0, prop_I_suspected <= 1)
  stopifnot(length(test_sens) == 1, test_sens >= 0, test_sens <= 1)
  stopifnot(length(test_turnaround) == 1, test_turnaround > 0)
  
  # I -> Xtest (suspected case identification and isolation)
  t_isolate_I_to_Xtest <- 4
  gamma <- prop_I_suspected / t_isolate_I_to_Xtest
  
  # Xtest -> Xtreat (true positives) and Xtest -> I (false negatives)
  theta <- test_sens / test_turnaround        # Xtest -> Xtreat (true positive)
  pi <- (1 - test_sens) / test_turnaround     # Xtest -> I (false negative)
  
  # Mortality in Xtest (depends on duration there)
  mXtest <- prob2rate(0.05, test_turnaround)
  
  list(
    scenario = scenario_label,
    
    # Admissions/transmission/natural history
    prop_adm_C = prop_adm_C, prop_adm_I = prop_adm_I,
    beta0 = beta0, betaC = betaC, betaI = betaI,
    betaXtest = betaXtest, betaXtreat = betaXtreat,
    alpha = alpha, abx_prop = abx_prop, abx_mult = abx_mult, mu = mu,
    
    # Discharge rates
    dS = dS, dC = dC, dI = dI, dXtreat = dXtreat,
    
    # Mortality rates
    mS = mS, mC = mC, mI = mI, mXtest = mXtest, mXtreat = mXtreat,
    
    # Scenario controls
    prop_I_suspected = prop_I_suspected,
    test_sens = test_sens,
    test_turnaround = test_turnaround,
    
    # Transition rates
    gamma = gamma,       # I -> Xtest
    theta = theta,       # Xtest -> Xtreat
    pi = pi,             # Xtest -> I
    delta = delta        # Xtreat -> S
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

## Copilot

Yes — **this is much better** and mostly consistent now.  
You’re very close to production-ready.

Two quick suggestions to make it robust:

1. Ensure your ODE function uses exactly these names (`dS`, `dC`, `dI`, `dXtreat`, `gamma`, `theta`, `pi`, `delta`, etc.).
2. Your comments say “half-days,” but LOS values like `losS <- 7` look like days. That’s fine, just be consistent with `times` in the solver.

Here’s your scenario builder snippet (as requested) in markdown:

```r
# Function to build parameter list for each scenario ----
make_mod_parms <- function(
  prop_I_suspected,
  test_sens,
  test_turnaround,
  scenario_label
) {
  # Validate inputs
  stopifnot(length(prop_I_suspected) == 1, prop_I_suspected >= 0, prop_I_suspected <= 1)
  stopifnot(length(test_sens) == 1, test_sens >= 0, test_sens <= 1)
  stopifnot(length(test_turnaround) == 1, test_turnaround > 0)

  # I -> Xtest
  t_isolate_I_to_Xtest <- 4
  gamma <- prop_I_suspected / t_isolate_I_to_Xtest

  # Xtest split by sensitivity over turnaround time
  theta <- test_sens / test_turnaround        # Xtest -> Xtreat (true positive)
  pi <- (1 - test_sens) / test_turnaround     # Xtest -> I (false negative)

  # Mortality in Xtest (depends on duration there)
  mXtest <- prob2rate(0.05, test_turnaround)

  list(
    scenario = scenario_label,

    # Admissions/transmission/natural history
    prop_adm_C = prop_adm_C, prop_adm_I = prop_adm_I,
    beta0 = beta0, betaC = betaC, betaI = betaI,
    betaXtest = betaXtest, betaXtreat = betaXtreat,
    alpha = alpha, abx_prop = abx_prop, abx_mult = abx_mult, mu = mu,

    # Discharge rates
    dS = dS, dC = dC, dI = dI, dXtreat = dXtreat,

    # Mortality rates
    mS = mS, mC = mC, mI = mI, mXtest = mXtest, mXtreat = mXtreat,

    # Scenario controls
    prop_I_suspected = prop_I_suspected,
    test_sens = test_sens,
    test_turnaround = test_turnaround,

    # Transition rates
    gamma = gamma,       # I -> Xtest
    theta = theta,       # Xtest -> Xtreat
    pi = pi,             # Xtest -> I
    delta = delta        # Xtreat -> S
  )
}
```

If you want, next I can do a **line-by-line ODE compatibility check** against this parameter list so you can run without name errors on first attempt.

## @LFell

Could you redraft the R script Model1A_ODE_Functions.R to take these changes into account? 

## Copilot

Absolutely. Your current function still uses the old single isolation compartment (`X`) and old parameter names.  
Here is a redraft aligned to your updated parameter script with **`Xtest` + `Xtreat`** and matching names.

```r name=Model1A_ODE_Functions.R
# Model 1A ODE Functions ----
# States:
# S      = susceptible
# C      = colonised
# I      = infected (not isolated)
# Xtest  = isolated awaiting test result
# Xtreat = isolated confirmed positive, on treatment/recovery
# D      = cumulative deaths (absorbing count state)

Model1A_ODE <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    # Total in-hospital living census
    N <- S + C + I + Xtest + Xtreat

    # Force of colonisation
    lambda <- beta0 + (betaC * C + betaI * I + betaXtest * Xtest + betaXtreat * Xtreat) / N

    # Fixed census admissions (Option 2A): replace all exits from living states
    exits_total <- (dS + mS) * S +
      (dC + mC) * C +
      (dI + mI) * I +
      (mXtest) * Xtest +              # no discharge from Xtest
      (dXtreat + mXtreat) * Xtreat

    # Admission split by import prevalence
    aI <- exits_total * prop_adm_I
    aC <- exits_total * prop_adm_C
    aS <- exits_total - aC - aI
    if (aS < 0) aS <- 0

    # ODEs
    dS <- aS + mu * C + delta * Xtreat - lambda * S - (dS + mS) * S
    dC <- aC + lambda * S - (alpha + mu + dC + mC) * C
    dI <- aI + alpha * C + pi * Xtest - (gamma + dI + mI) * I
    dXtest <- gamma * I - (theta + pi + mXtest) * Xtest
    dXtreat <- theta * Xtest - (delta + dXtreat + mXtreat) * Xtreat
    dD <- mS * S + mC * C + mI * I + mXtest * Xtest + mXtreat * Xtreat

    # Flows (for outputs)
    flow_S_to_C <- lambda * S
    flow_C_to_I <- alpha * C
    flow_C_to_S <- mu * C
    flow_I_to_Xtest <- gamma * I
    flow_Xtest_to_I <- pi * Xtest
    flow_Xtest_to_Xtreat <- theta * Xtest
    flow_Xtreat_to_S <- delta * Xtreat

    flow_dis_S <- dS * S
    flow_dis_C <- dC * C
    flow_dis_I <- dI * I
    flow_dis_Xtreat <- dXtreat * Xtreat
    flow_dis_total <- flow_dis_S + flow_dis_C + flow_dis_I + flow_dis_Xtreat

    flow_death_S <- mS * S
    flow_death_C <- mC * C
    flow_death_I <- mI * I
    flow_death_Xtest <- mXtest * Xtest
    flow_death_Xtreat <- mXtreat * Xtreat
    flow_death_total <- flow_death_S + flow_death_C + flow_death_I + flow_death_Xtest + flow_death_Xtreat

    # Cases
    total_cases <- flow_C_to_I + aI
    identified_cases <- flow_Xtreat_to_S + flow_dis_Xtreat
    missed_false_negative <- flow_Xtest_to_I
    missed_untested <- (1 - prop_I_suspected) * I
    missed_cases <- missed_untested + missed_false_negative

    # Side-room occupancy (patient-days per model time step)
    side_room_patient_days <- Xtest + Xtreat

    # Return
    list(
      c(dS, dC, dI, dXtest, dXtreat, dD),
      c(
        N = N, lambda = lambda,
        aS = aS, aC = aC, aI = aI,

        flow_S_to_C = flow_S_to_C,
        flow_C_to_I = flow_C_to_I,
        flow_C_to_S = flow_C_to_S,
        flow_I_to_Xtest = flow_I_to_Xtest,
        flow_Xtest_to_I = flow_Xtest_to_I,
        flow_Xtest_to_Xtreat = flow_Xtest_to_Xtreat,
        flow_Xtreat_to_S = flow_Xtreat_to_S,

        flow_dis_S = flow_dis_S,
        flow_dis_C = flow_dis_C,
        flow_dis_I = flow_dis_I,
        flow_dis_Xtreat = flow_dis_Xtreat,
        flow_dis_total = flow_dis_total,

        flow_death_S = flow_death_S,
        flow_death_C = flow_death_C,
        flow_death_I = flow_death_I,
        flow_death_Xtest = flow_death_Xtest,
        flow_death_Xtreat = flow_death_Xtreat,
        flow_death_total = flow_death_total,

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

### Important follow-up
In your run script, update initial states to include:
- `Xtest`, `Xtreat`, and `D`

Example:
```r
mod.init <- c(S = 620, C = 15, I = 5, Xtest = 0, Xtreat = 0, D = 0)
```

If you want, I can now redraft your `Run.qmd` chunks so all summaries/plots match these new variable names automatically.

## @LFell

I think in an earlier script I also had the cumulative totals for total_cases, identified_cases, missed_cases, missed_untested, missed_false_negative, with the option for a burn in period to exclude from the results

## Copilot

Yes — that’s a good setup, and you should definitely keep it.

Here’s a clean chunk you can drop into your run `.qmd` to do exactly that (with burn-in exclusion + cumulative totals).

```r name=burnin_and_cumulative_outputs.R
# Apply burn-in and reset time
apply_burn_in <- function(df, burn_days = 365) {
  df %>%
    dplyr::filter(time > burn_days) %>%
    dplyr::mutate(time = time - burn_days)
}

# Add cumulative outputs after burn-in
add_case_cumulatives <- function(df) {
  df %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(
      cum_total_cases = cumsum(total_cases),
      cum_identified_cases = cumsum(identified_cases),
      cum_missed_cases = cumsum(missed_cases),
      cum_missed_untested = cumsum(missed_untested),
      cum_missed_false_negative = cumsum(missed_false_negative)
    )
}

# Optional: annual cumulatives (reset each year)
add_annual_case_cumulatives <- function(df) {
  df %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(year = if_else(time == 0, 1L, ((time - 1) %/% 365) + 1L)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(
      annual_cum_total_cases = cumsum(total_cases),
      annual_cum_identified_cases = cumsum(identified_cases),
      annual_cum_missed_cases = cumsum(missed_cases),
      annual_cum_missed_untested = cumsum(missed_untested),
      annual_cum_missed_false_negative = cumsum(missed_false_negative)
    ) %>%
    dplyr::ungroup()
}
```

Example usage:

```r name=usage_example.R
scenario_raw <- lapply(scenario_parms, run_scenario)

scenario_postburn <- lapply(scenario_raw, function(x) {
  x %>%
    apply_burn_in(burn_days = burn_in_days) %>%
    add_case_cumulatives() %>%
    add_annual_case_cumulatives()
})
```

If you want, I can also give you a matching `summarise_outcomes()` that returns final (post-burn-in) cumulative totals per scenario in one table.

## @LFell

Just checking that calculating cumulative cases doesn't require any changes to the ODE function script

## Copilot

Correct — **no changes needed in the ODE function**.

As long as your ODE function outputs the per-time-step flows (`total_cases`, `identified_cases`, `missed_cases`, etc.), cumulative values should be calculated **after solving** (in your run/analysis script) using `cumsum()`.

So:
- ODE script: keep **instantaneous rates/flows**
- Run script: compute **cumulative totals**, burn-in filtering, annual summaries

## @LFell

Is it a problem that I have discharge rates dS dC dI dXtreat in the paramenters script and changes in S, C labellled as dS, dC etc in the functions script? shall I rename the dischages

## Copilot

Yes — that **is a problem** (or at least a major source of bugs/confusion).  
You should rename one set.

In ODE code, `dS`, `dC`, etc. conventionally mean **state derivatives**.  
If you also use `dS` as a **discharge rate parameter**, names collide conceptually (and can collide in `with(as.list(...))` logic).

## Recommended fix
Rename discharge-rate parameters to something explicit, e.g.:

- `disS`, `disC`, `disI`, `disXtreat`
  or
- `r_dis_S`, `r_dis_C`, `r_dis_I`, `r_dis_Xtreat`

Then keep derivatives as:
- `dS`, `dC`, `dI`, `dXtest`, `dXtreat`, `dD`

### Example
```r
# parameters
disS <- 1 / losS
disC <- 1 / losC
disI <- 1 / losI
disXtreat <- (1 - prop_Xtreat_to_S) / losXtreat
```

```r
# ODE
dS <- aS + mu * C + delta * Xtreat - lambda * S - (disS + mS) * S
```

So yes — **rename discharge rates now**; it will save debugging time later.

## @LFell

I have corrected the parameter and functions scripts to make the names clearer. I am not sure about the method for calculating missed untested cases. Either there could be a gamma0 = 1 /  t_isolate_I_to_Xtest and missed_untested <- gamma0*I or, I could run two scenarios with prop_I_suspected = 1 and prop_I_suspected = 0.5, and calculate the difference in total missed cases

## Copilot

Use **neither** of those as your primary definition.

Best definition (instantaneous flow) is:

\[
\text{missed\_untested\_flow} = (1-\text{prop\_I\_suspected}) \times \gamma_0 \times I
\]
where \(\gamma_0 = 1/t_{\text{isolate}}\).

Why:
- `gamma0 * I` alone is the **total potential suspicion flow** if everyone were suspectable.
- Multiplying by `(1 - prop_I_suspected)` gives the portion that is **not identified** at that same decision point.
- Using `(1 - prop_I_suspected) * I` (without gamma0) is not a flow (units mismatch).

Scenario-difference (prop=1 vs 0.5) can be a **sensitivity check**, but not a clean decomposition because system dynamics feedback changes everything else too.

So in code:

```r name=missed_untested_flow.R
gamma0 <- 1 / t_isolate_I_to_Xtest
gamma  <- prop_I_suspected * gamma0

flow_I_to_Xtest <- gamma * I
missed_untested <- (1 - prop_I_suspected) * gamma0 * I
```

Then:
```r
missed_cases <- missed_untested + missed_false_negative
```

This gives a consistent flow-based accounting.

## @LFell

Could you also have # I -> Xtest (suspected case identification and isolation)
  t_isolate_I_to_Xtest <- 4
  gamma0 <- 1 / t_isolate_I_to_Xtest # Rate of isolation if all infected patients are identified as suspected cases (per half-day)
  gamma <- prop_I_suspected * gamma0 # Rate of isolation if a proportion of infected patients are identified as suspected cases (per half-day)   ## flow_I_to_Xtest <- gamma * I ## missed_untested <- gamma0 * I - flow_I_to_Xtest
  

## Copilot

Yes — that formulation is correct, and it’s a good way to keep units consistent.

```r
# I -> Xtest (suspected case identification and isolation)
t_isolate_I_to_Xtest <- 4
gamma0 <- 1 / t_isolate_I_to_Xtest              # if 100% suspected
gamma  <- prop_I_suspected * gamma0             # realised isolation rate

flow_I_to_Xtest <- gamma * I
missed_untested <- gamma0 * I - flow_I_to_Xtest
# equivalent: missed_untested <- (1 - prop_I_suspected) * gamma0 * I
```

That gives:
- `flow_I_to_Xtest`: identified/suspected flow
- `missed_untested`: not-identified-at-screening flow
- and `missed_cases = missed_untested + missed_false_negative` as desired.

One tiny safety guard you can add:
```r
missed_untested <- pmax(0, gamma0 * I - flow_I_to_Xtest)
```

## @LFell

Here are the current versions. can you check for naming consistency and suggest any other safeguards to avoid negative values

## Copilot

Yes — this is very close. You have a few naming mismatches that will break execution.

## Critical naming fixes

1) In `Functions.R`:
- You compute `admS/admC/admI` but check `aS`:
```r
if (aS < 0) aS <- 0
```
Should be:
```r
if (admS < 0) admS <- 0
```

2) In `Functions.R`:
- `total_cases <- flow_C_to_I + aI` should use `admI`:
```r
total_cases <- flow_C_to_I + admI
```

3) In `Parameters.R` list:
```r
disS = disS, disC = dC, disI = disI, ...
```
`dC` doesn’t exist now. Should be:
```r
disS = disS, disC = disC, disI = disI, ...
```

4) ODE state comments still mention old `X` in header earlier versions; your latest uses `Xtest/Xtreat` correctly—good.

---

## Safeguards for negative values (recommended)

Add these in ODE function after `N` and before calculations:

```r
# Guard against zero/negative census division
N <- max(N, 1e-12)
```

Clamp potentially tiny negative states due to solver tolerance:

```r
S <- max(S, 0); C <- max(C, 0); I <- max(I, 0); Xtest <- max(Xtest, 0); Xtreat <- max(Xtreat, 0)
```

Keep your admission guard:

```r
admS <- exits_total - admC - admI
admS <- max(admS, 0)
```

For missed untested (you already did this well):

```r
missed_untested <- pmax(0, gamma0 * I - flow_I_to_Xtest)
```

Also enforce parameter sanity in `make_mod_parms()`:

```r
stopifnot(prop_adm_C + prop_adm_I <= 1)
stopifnot(delta >= 0, disXtreat >= 0, theta >= 0, pi >= 0, gamma >= 0)
```

---

## One modeling consistency note
You currently apply antibiotic effect in parameters (`alpha` already adjusted), and in ODE you use `alpha` directly — that is consistent ✅ (don’t multiply by `abx_mult` again in the ODE).

If you want, I can provide a **fully corrected final pair** of files with just these fixes applied verbatim.

## @LFell

Is this ok for parameters

## Copilot

Yes — this parameters file now looks **good and internally consistent** ✅

Only small safeguards I’d add:

```r name=Model1A_ODE_Parameters.R
# after fixed admission prevalence
stopifnot(prop_adm_C >= 0, prop_adm_I >= 0, (prop_adm_C + prop_adm_I) <= 1)

# inside make_mod_parms(), after existing stopifnot()
stopifnot(prop_I_suspected <= 1, prop_I_suspected >= 0)
stopifnot(test_sens <= 1, test_sens >= 0)
stopifnot(test_turnaround > 0)

# optional numeric safety checks after computing rates
stopifnot(gamma0 >= 0, gamma >= 0, theta >= 0, pi >= 0, delta >= 0)
stopifnot(disS >= 0, disC >= 0, disI >= 0, disXtreat >= 0)
stopifnot(mortS >= 0, mortC >= 0, mortI >= 0, mortXtest >= 0, mortXtreat >= 0)
```

### One conceptual note (not an error)
You currently calculate rates like:
- `alpha <- proportion / time`
- `mu <- proportion / time`

That is fine **if** those are interpreted as “fraction transitioning per mean time.”  
If instead those are cumulative risks over `t`, use `prob2rate()`.

But structurally, your file is now in good shape.

## @LFell

Have i added the safety checks in the right places

## Copilot

Yes — mostly right, and nicely done. ✅  
Just two tweaks:

1. You duplicated some input checks in `make_mod_parms` (harmless, but redundant).  
2. Add a few **fixed-parameter sanity checks** once in section 3 (for LOS, proportions, etc.), not only inside `make_mod_parms`.

Here’s a compact cleaned pattern:

```r
# Fixed-parameter checks (Section 3, after fixed values defined)
stopifnot(
  prop_adm_C >= 0, prop_adm_I >= 0, prop_adm_C + prop_adm_I <= 1,
  abx_prop >= 0, abx_prop <= 1,
  prop_progress_C_to_I >= 0, prop_progress_C_to_I <= 1,
  prop_decolonise_C_to_S >= 0, prop_decolonise_C_to_S <= 1,
  prop_Xtreat_to_S >= 0, prop_Xtreat_to_S <= 1,
  losS > 0, losC > 0, losI > 0, losXtreat > 0,
  t_progress_C_to_I > 0, t_decolonise_C_to_S > 0
)
```

```r
# In make_mod_parms: keep one set of input checks only
stopifnot(
  length(prop_I_suspected) == 1, prop_I_suspected >= 0, prop_I_suspected <= 1,
  length(test_sens) == 1, test_sens >= 0, test_sens <= 1,
  length(test_turnaround) == 1, test_turnaround > 0
)
```

Your post-rate checks (`gamma >= 0`, `theta >= 0`, etc.) are in a good place.  
So: **yes, correct placement overall** — just simplify duplicates and add the fixed-value block.

## @LFell

If i am getting it to stop should I also get it to create warnings so i know why its stopped

## Copilot

Yes — great idea.  
Best practice is:

- use `stop()` for **fatal invalid inputs**
- use `warning()` for **non-fatal but suspicious** values

`stopifnot()` is quick, but `stop("...")` gives clearer messages.

```r name=validation_helpers.R
validate_inputs <- function(prop_I_suspected, test_sens, test_turnaround) {
  if (length(prop_I_suspected) != 1 || is.na(prop_I_suspected) ||
      prop_I_suspected < 0 || prop_I_suspected > 1) {
    stop("prop_I_suspected must be a single number in [0,1].")
  }

  if (length(test_sens) != 1 || is.na(test_sens) ||
      test_sens < 0 || test_sens > 1) {
    stop("test_sens must be a single number in [0,1].")
  }

  if (length(test_turnaround) != 1 || is.na(test_turnaround) ||
      test_turnaround <= 0) {
    stop("test_turnaround must be a single number > 0.")
  }

  # non-fatal warnings
  if (test_turnaround < 1) {
    warning("test_turnaround < 1 half-day: check units/realism.")
  }
  if (test_sens < 0.3) {
    warning("Very low test_sens (<0.3): check evidence source.")
  }
}
```

Then call at top of `make_mod_parms()`:

```r
validate_inputs(prop_I_suspected, test_sens, test_turnaround)
```

If you want, I can convert your current checks into this style for the full parameter script.

## @LFell

OK can you use this version 

## Copilot

Perfect — yes, use this version. It looks solid. ✅

I’d proceed with it as your active `Model1A_ODE_Parameters.R`.

Only optional polish:
- if you want clearer error messages later, replace some `stopifnot(...)` with explicit `if (...) stop("...")`.
- keep an eye on units (“half-days”) in the run script time grid so everything stays aligned.

If you want, next I can do a final **paired consistency check** between this parameters file and your current `Model1A_ODE_Functions.R` and give you a “ready-to-run” final functions file with the `admS/admI` naming fixes.

## @LFell

OK this is what I have at the moment for the functions script.

## Copilot

Looks very good — this is essentially ready. ✅

Just do these **2 small cleanups**:

1) Remove the duplicate assignment of `missed_untested` (keep only the `pmax` line).  
2) Optional robustness: clamp states to non-negative before using them.

Here’s the exact edit for the cases block:

```r
# Cases
total_cases <- flow_C_to_I + admI
identified_cases <- flow_Xtreat_to_S + flow_dis_Xtreat
missed_false_negative <- flow_Xtest_to_I
missed_untested <- pmax(0, gamma0 * I - flow_I_to_Xtest)
missed_cases <- missed_untested + missed_false_negative
```

Optional state clamp near top (after `with(...)`):

```r
S <- max(S, 0); C <- max(C, 0); I <- max(I, 0); Xtest <- max(Xtest, 0); Xtreat <- max(Xtreat, 0)
```

Everything else is naming-consistent with your latest parameters file (`dis*`, `mort*`, `adm*`, `gamma0`, etc.).  
So yes — you’re in a strong place to run.

## @LFell

Can you provide clean, corrected versions of both scripts, with your suggestions for using stop() and warning() instead of stopifnot(); validate_inputs, remove the duplicate assignment of missed_untested and clamp states to non-negative before using them. 

## Copilot

Absolutely — here are clean, corrected versions of both scripts with all requested changes.

### `Model1A_ODE_Parameters.R`

```r name=Model1A_ODE_Parameters.R
# Model 1A Parameters ####

# 1) Introduction ####
# This script defines the parameters for the ODE model of hospital transmission of a pathogen
# (e.g., C difficile) with a focus on identification of suspected cases, test turnaround time,
# and test sensitivity.

# 2) Initial setup ####
options(scipen = 1000)

# Convert cumulative probability over interval t to constant hazard rate
prob2rate <- function(p, t) {
  if (!is.numeric(p) || length(p) != 1 || is.na(p) || p < 0 || p >= 1) {
    stop("prob2rate(): 'p' must be a single numeric value in [0, 1).")
  }
  if (!is.numeric(t) || length(t) != 1 || is.na(t) || t <= 0) {
    stop("prob2rate(): 't' must be a single numeric value > 0.")
  }
  -(1 / t) * log(1 - p)
}

# ---- Validation helpers ----

validate_fixed_parameters <- function() {
  if (prop_adm_C < 0 || prop_adm_I < 0 || (prop_adm_C + prop_adm_I) > 1) {
    stop("Invalid admission prevalence: require prop_adm_C >= 0, prop_adm_I >= 0, and prop_adm_C + prop_adm_I <= 1.")
  }

  if (abx_prop < 0 || abx_prop > 1) stop("abx_prop must be in [0,1].")
  if (prop_progress_C_to_I < 0 || prop_progress_C_to_I > 1) stop("prop_progress_C_to_I must be in [0,1].")
  if (prop_decolonise_C_to_S < 0 || prop_decolonise_C_to_S > 1) stop("prop_decolonise_C_to_S must be in [0,1].")
  if (prop_Xtreat_to_S < 0 || prop_Xtreat_to_S > 1) stop("prop_Xtreat_to_S must be in [0,1].")

  if (losS <= 0 || losC <= 0 || losI <= 0 || losXtreat <= 0) stop("All LOS parameters must be > 0.")
  if (t_progress_C_to_I <= 0 || t_decolonise_C_to_S <= 0) stop("Progression/decolonisation times must be > 0.")

  # Non-fatal plausibility warnings
  if (prop_adm_C + prop_adm_I > 0.5) warning("High imported prevalence (prop_adm_C + prop_adm_I > 0.5). Check assumptions.")
  if (losS < 1 || losC < 1 || losI < 1 || losXtreat < 1) warning("Some LOS values are < 1 half-day. Check units.")
}

validate_inputs <- function(prop_I_suspected, test_sens, test_turnaround, scenario_label = NA_character_) {
  nm <- if (!is.na(scenario_label)) paste0("[", scenario_label, "] ") else ""

  if (!is.numeric(prop_I_suspected) || length(prop_I_suspected) != 1 || is.na(prop_I_suspected) ||
      prop_I_suspected < 0 || prop_I_suspected > 1) {
    stop(nm, "prop_I_suspected must be a single numeric value in [0,1].")
  }

  if (!is.numeric(test_sens) || length(test_sens) != 1 || is.na(test_sens) ||
      test_sens < 0 || test_sens > 1) {
    stop(nm, "test_sens must be a single numeric value in [0,1].")
  }

  if (!is.numeric(test_turnaround) || length(test_turnaround) != 1 || is.na(test_turnaround) ||
      test_turnaround <= 0) {
    stop(nm, "test_turnaround must be a single numeric value > 0.")
  }

  # Non-fatal warnings
  if (test_turnaround < 1) warning(nm, "test_turnaround < 1 half-day; check units/realism.")
  if (test_sens < 0.3) warning(nm, "Very low test_sens (<0.3); check evidence source.")
  if (prop_I_suspected < 0.2) warning(nm, "Low prop_I_suspected (<0.2); many infected may remain unisolated.")
}

validate_rate_outputs <- function(gamma0, gamma, theta, pi, delta, disS, disC, disI, disXtreat,
                                  mortS, mortC, mortI, mortXtest, mortXtreat, scenario_label = NA_character_) {
  nm <- if (!is.na(scenario_label)) paste0("[", scenario_label, "] ") else ""

  vals <- c(gamma0, gamma, theta, pi, delta, disS, disC, disI, disXtreat, mortS, mortC, mortI, mortXtest, mortXtreat)
  nms  <- c("gamma0","gamma","theta","pi","delta","disS","disC","disI","disXtreat","mortS","mortC","mortI","mortXtest","mortXtreat")

  if (any(!is.finite(vals))) {
    bad <- paste(nms[!is.finite(vals)], collapse = ", ")
    stop(nm, "Non-finite rate(s): ", bad, ".")
  }
  if (any(vals < 0)) {
    bad <- paste(nms[vals < 0], collapse = ", ")
    stop(nm, "Negative rate(s): ", bad, ".")
  }

  if (theta + pi > 1) warning(nm, "theta + pi > 1 per half-day; transitions from Xtest may be very fast for chosen timestep.")
}

# 3) Fixed assumptions ####

# Admission prevalence
prop_adm_C <- 0.03
prop_adm_I <- 0.002

# Transmission
beta0 <- 0.00020
betaC <- 0.00035
betaI <- 0.00060
betaXtest <- 0.00010
betaXtreat <- 0.00005

# Progression C -> I
prop_progress_C_to_I <- 0.010
abx_mult <- 1.30
abx_prop <- 0.20
t_progress_C_to_I <- 10.0
alpha <- (prop_progress_C_to_I * ((1 - abx_prop) + (abx_prop * abx_mult)) / t_progress_C_to_I)

# Decolonisation C -> S
prop_decolonise_C_to_S <- 0.05
t_decolonise_C_to_S <- 5.0
mu <- prop_decolonise_C_to_S / t_decolonise_C_to_S

# Length of stay
losS <- 7
losC <- 10
losI <- 12
losXtreat <- 14

# Xtreat outcomes
prop_Xtreat_to_S <- 0.80
delta <- prop_Xtreat_to_S / losXtreat

# Discharge rates
disS <- 1 / losS
disC <- 1 / losC
disI <- 1 / losI
disXtreat <- (1 - prop_Xtreat_to_S) / losXtreat

# Mortality rates
mortS <- prob2rate(0.02, losS)
mortC <- prob2rate(0.02, losC)
mortI <- prob2rate(0.20, losI)
mortXtreat <- prob2rate(0.10, losXtreat)

# Validate fixed assumptions once
validate_fixed_parameters()

# 4) Scenario builder ####

make_mod_parms <- function(
  prop_I_suspected,
  test_sens,
  test_turnaround,
  scenario_label
) {
  validate_inputs(prop_I_suspected, test_sens, test_turnaround, scenario_label)

  # I -> Xtest
  t_isolate_I_to_Xtest <- 4
  gamma0 <- 1 / t_isolate_I_to_Xtest
  gamma  <- prop_I_suspected * gamma0

  # Xtest -> Xtreat / I
  theta <- test_sens / test_turnaround
  pi    <- (1 - test_sens) / test_turnaround

  # Mortality in Xtest
  mortXtest <- prob2rate(0.05, test_turnaround)

  validate_rate_outputs(
    gamma0, gamma, theta, pi, delta,
    disS, disC, disI, disXtreat,
    mortS, mortC, mortI, mortXtest, mortXtreat,
    scenario_label
  )

  list(
    scenario = scenario_label,

    # Admissions/transmission/natural history
    prop_adm_C = prop_adm_C, prop_adm_I = prop_adm_I,
    beta0 = beta0, betaC = betaC, betaI = betaI,
    betaXtest = betaXtest, betaXtreat = betaXtreat,
    alpha = alpha, abx_prop = abx_prop, abx_mult = abx_mult, mu = mu,

    # Discharge/mortality
    disS = disS, disC = disC, disI = disI, disXtreat = disXtreat,
    mortS = mortS, mortC = mortC, mortI = mortI, mortXtest = mortXtest, mortXtreat = mortXtreat,

    # Scenario controls
    prop_I_suspected = prop_I_suspected,
    test_sens = test_sens,
    test_turnaround = test_turnaround,

    # Transition rates
    gamma0 = gamma0,
    gamma = gamma,
    theta = theta,
    pi = pi,
    delta = delta
  )
}

# 5) Scenarios ####

mod.parms.S1 <- make_mod_parms(0.50, 0.80, 5, "S1_midSus_highSe_slow")
mod.parms.S2 <- make_mod_parms(0.50, 0.50, 5, "S2_midSus_lowSe_slow")
mod.parms.S3 <- make_mod_parms(0.50, 0.80, 3, "S3_midSus_highSe_fast")
mod.parms.S4 <- make_mod_parms(0.50, 0.50, 3, "S4_midSus_low_Se_fast")
```

---

### `Model1A_ODE_Functions.R`

```r name=Model1A_ODE_Functions.R
# Model 1A ODE Functions ----
# States:
# S      = susceptible
# C      = colonised
# I      = infected (not isolated)
# Xtest  = isolated awaiting test result
# Xtreat = isolated confirmed positive, on treatment/recovery
# D      = cumulative deaths (absorbing count state)

Model1A_ODE <- function(time, state, parms) {
  with(as.list(c(state, parms)), {

    # Clamp states to avoid tiny negative values from numerical solver tolerance
    S <- max(S, 0)
    C <- max(C, 0)
    I <- max(I, 0)
    Xtest <- max(Xtest, 0)
    Xtreat <- max(Xtreat, 0)
    D <- max(D, 0)

    # Total in-hospital living census
    N <- S + C + I + Xtest + Xtreat
    N <- max(N, 1e-12)  # guard division by zero

    # Force of colonisation
    lambda <- beta0 + (betaC * C + betaI * I + betaXtest * Xtest + betaXtreat * Xtreat) / N

    # Fixed census admissions (balance discharges + deaths)
    exits_total <- (disS + mortS) * S +
      (disC + mortC) * C +
      (disI + mortI) * I +
      (mortXtest) * Xtest +                 # no discharge from Xtest
      (disXtreat + mortXtreat) * Xtreat

    # Admissions split by import prevalence
    admI <- exits_total * prop_adm_I
    admC <- exits_total * prop_adm_C
    admS <- exits_total - admC - admI
    admS <- max(admS, 0)                    # non-negative guard

    # ODEs
    dS <- admS + mu * C + delta * Xtreat - lambda * S - (disS + mortS) * S
    dC <- admC + lambda * S - (alpha + mu + disC + mortC) * C
    dI <- admI + alpha * C + pi * Xtest - (gamma + disI + mortI) * I
    dXtest <- gamma * I - (theta + pi + mortXtest) * Xtest
    dXtreat <- theta * Xtest - (delta + disXtreat + mortXtreat) * Xtreat
    dD <- mortS * S + mortC * C + mortI * I + mortXtest * Xtest + mortXtreat * Xtreat

    # Flows for outputs
    flow_S_to_C <- lambda * S
    flow_C_to_I <- alpha * C
    flow_C_to_S <- mu * C
    flow_I_to_Xtest <- gamma * I
    flow_Xtest_to_I <- pi * Xtest
    flow_Xtest_to_Xtreat <- theta * Xtest
    flow_Xtreat_to_S <- delta * Xtreat

    flow_dis_S <- disS * S
    flow_dis_C <- disC * C
    flow_dis_I <- disI * I
    flow_dis_Xtreat <- disXtreat * Xtreat
    flow_dis_total <- flow_dis_S + flow_dis_C + flow_dis_I + flow_dis_Xtreat

    flow_death_S <- mortS * S
    flow_death_C <- mortC * C
    flow_death_I <- mortI * I
    flow_death_Xtest <- mortXtest * Xtest
    flow_death_Xtreat <- mortXtreat * Xtreat
    flow_death_total <- flow_death_S + flow_death_C + flow_death_I + flow_death_Xtest + flow_death_Xtreat

    # Case accounting
    total_cases <- flow_C_to_I + admI
    identified_cases <- flow_Xtreat_to_S + flow_dis_Xtreat
    missed_false_negative <- flow_Xtest_to_I
    missed_untested <- pmax(0, gamma0 * I - flow_I_to_Xtest)  # single assignment only
    missed_cases <- missed_untested + missed_false_negative

    # Side-room occupancy
    side_room_patient_days <- Xtest + Xtreat

    list(
      c(dS, dC, dI, dXtest, dXtreat, dD),
      c(
        N = N, lambda = lambda,
        admS = admS, admC = admC, admI = admI,

        flow_S_to_C = flow_S_to_C,
        flow_C_to_I = flow_C_to_I,
        flow_C_to_S = flow_C_to_S,
        flow_I_to_Xtest = flow_I_to_Xtest,
        flow_Xtest_to_I = flow_Xtest_to_I,
        flow_Xtest_to_Xtreat = flow_Xtest_to_Xtreat,
        flow_Xtreat_to_S = flow_Xtreat_to_S,

        flow_dis_S = flow_dis_S,
        flow_dis_C = flow_dis_C,
        flow_dis_I = flow_dis_I,
        flow_dis_Xtreat = flow_dis_Xtreat,
        flow_dis_total = flow_dis_total,

        flow_death_S = flow_death_S,
        flow_death_C = flow_death_C,
        flow_death_I = flow_death_I,
        flow_death_Xtest = flow_death_Xtest,
        flow_death_Xtreat = flow_death_Xtreat,
        flow_death_total = flow_death_total,

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

If you want, I can also give you a **matching run-script template** (`ode()` call + burn-in + cumulative outputs) aligned to these exact names.

## @LFell

as abx_mult and abx_prop are only used to calculate alpha, do they need to be included in the list of model parameters 

## Copilot

Great question: **no, they don’t need to be in `parms`** if the ODE function never uses them directly.

Since you already pass `alpha`, and ODE uses `alpha`, you can remove:

- `abx_prop`
- `abx_mult`

from the returned list safely.

Keep them in the parameters script for transparency/documentation, just not in the runtime parameter list unless needed for outputs/sensitivity reporting.

## @LFell

are there any others that could be removed?

## Copilot

Yes — from your current ODE function, you can remove any parameter not referenced inside `Model1A_ODE()`.

Safe to remove from `parms` list:

- `scenario` (keep separately as metadata, not needed by ODE math)
- `test_sens`
- `test_turnaround`
- `prop_I_suspected`
- `abx_prop`
- `abx_mult`

Potentially remove `gamma0` **only if** you change missed-untested formula to:
\[
(1-\text{prop\_I\_suspected}) \cdot \gamma0 \cdot I
\]
But since your current ODE directly uses `gamma0`, keep it for now.

So your minimal required `parms` for the current function are:
- `prop_adm_C`, `prop_adm_I`
- `beta0`, `betaC`, `betaI`, `betaXtest`, `betaXtreat`
- `alpha`, `mu`
- `disS`, `disC`, `disI`, `disXtreat`
- `mortS`, `mortC`, `mortI`, `mortXtest`, `mortXtreat`
- `gamma0`, `gamma`, `theta`, `pi`, `delta`

If you want, I can give you a “lean” `make_mod_parms()` that returns only these plus a separate scenario metadata table.

## @LFell

do I need to put anything inside the brackets to call validate_fixed_parameters()?

## Copilot

No — call it with empty parentheses:

```r
validate_fixed_parameters()
```

You only pass arguments if the function is defined to accept them.  
Your current version uses objects from the script environment, so empty call is correct.

## @LFell

Could we add something to confirm it has run and not found any issues

## Copilot

Yes—great idea. Add a `message()` at the end of the validator.

```r
validate_fixed_parameters <- function() {
  # ... all your checks with stop()/warning()

  message("✅ validate_fixed_parameters(): checks passed.")
}
```

Then call:

```r
validate_fixed_parameters()
```

If no fatal issues:
- warnings may print (if any),
- then you’ll see the success message.

If you want it quieter in production, wrap with a flag:

```r
validate_fixed_parameters <- function(verbose = TRUE) {
  # checks...
  if (verbose) message("✅ validate_fixed_parameters(): checks passed.")
}
```

## @LFell

Could we add similar messages for the other two validation functions

## Copilot

Yes — definitely. Add an optional `verbose = TRUE` argument to each and print a success `message()` at the end.

```r name=validation_functions_with_messages.R
validate_inputs <- function(prop_I_suspected, test_sens, test_turnaround,
                            scenario_label = NA_character_, verbose = TRUE) {
  nm <- if (!is.na(scenario_label)) paste0("[", scenario_label, "] ") else ""

  if (!is.numeric(prop_I_suspected) || length(prop_I_suspected) != 1 || is.na(prop_I_suspected) ||
      prop_I_suspected < 0 || prop_I_suspected > 1) {
    stop(nm, "prop_I_suspected must be a single numeric value in [0,1].")
  }
  if (!is.numeric(test_sens) || length(test_sens) != 1 || is.na(test_sens) ||
      test_sens < 0 || test_sens > 1) {
    stop(nm, "test_sens must be a single numeric value in [0,1].")
  }
  if (!is.numeric(test_turnaround) || length(test_turnaround) != 1 || is.na(test_turnaround) ||
      test_turnaround <= 0) {
    stop(nm, "test_turnaround must be a single numeric value > 0.")
  }

  if (test_turnaround < 1) warning(nm, "test_turnaround < 1 half-day; check units/realism.")
  if (test_sens < 0.3) warning(nm, "Very low test_sens (<0.3); check evidence source.")
  if (prop_I_suspected < 0.2) warning(nm, "Low prop_I_suspected (<0.2); many infected may remain unisolated.")

  if (verbose) message("✅ ", nm, "validate_inputs(): checks passed.")
}

validate_rate_outputs <- function(gamma0, gamma, theta, pi, delta, disS, disC, disI, disXtreat,
                                  mortS, mortC, mortI, mortXtest, mortXtreat,
                                  scenario_label = NA_character_, verbose = TRUE) {
  nm <- if (!is.na(scenario_label)) paste0("[", scenario_label, "] ") else ""

  vals <- c(gamma0, gamma, theta, pi, delta, disS, disC, disI, disXtreat,
            mortS, mortC, mortI, mortXtest, mortXtreat)
  nms  <- c("gamma0","gamma","theta","pi","delta","disS","disC","disI","disXtreat",
            "mortS","mortC","mortI","mortXtest","mortXtreat")

  if (any(!is.finite(vals))) {
    bad <- paste(nms[!is.finite(vals)], collapse = ", ")
    stop(nm, "Non-finite rate(s): ", bad, ".")
  }
  if (any(vals < 0)) {
    bad <- paste(nms[vals < 0], collapse = ", ")
    stop(nm, "Negative rate(s): ", bad, ".")
  }

  if (theta + pi > 1) warning(nm, "theta + pi > 1 per half-day; transitions from Xtest may be very fast.")

  if (verbose) message("✅ ", nm, "validate_rate_outputs(): checks passed.")
}
```

And if you also want consistency, do the same in fixed validator:

```r
if (verbose) message("✅ validate_fixed_parameters(): checks passed.")
```

Then in `make_mod_parms()` you can keep them on:
```r
validate_inputs(..., verbose = TRUE)
validate_rate_outputs(..., verbose = TRUE)
```

or silence during batch runs:
```r
validate_inputs(..., verbose = FALSE)
validate_rate_outputs(..., verbose = FALSE)
```

## @LFell

OK this is the final parameter script. I have check and it runs OK 

## Copilot

Excellent — this is in very good shape, and great that it runs ✅

One important cleanup before you lock it in: you currently define `validate_rate_outputs()` **twice**. Keep only one definition (the earlier one with message/verbose, or the later one), otherwise the second silently overwrites the first.

Also minor typo in a message: `"betenn"` → `"between"`.

If you want, this is the safe final check list:
- [x] Single definition per helper function
- [x] `validate_fixed_parameters()` called once after fixed params
- [x] Scenario builder calls `validate_inputs()` and `validate_rate_outputs()`
- [x] No duplicate parameter names
- [x] ODE script names match parameter list exactly

If all that is true, you’re ready to proceed to analysis outputs.

## @LFell

OK this is the parameter script with those changes. 

## Copilot

Looks great overall ✅ — just **one fix needed**:

In `validate_rate_outputs()`, you reference `test_turnaround`, but it is **not an argument** to that function, so that warning line will error.

Remove this line:
```r
if (test_turnaround < 1) warning(nm, "test_turnaround < 1 half-day; check units/realism.")
```
(you already check that in `validate_inputs()`, so it’s covered).

Everything else is consistent and good to go.

## @LFell

Given that parameters script, are there any changes required in the functions script attached? 

## Copilot

Yes — **two changes required** and one recommended.

### Required
1) **Clamp states to non-negative** (you asked for this, but it’s not yet in this version).
2) **Remove duplicate `missed_untested` assignment** (keep only the `pmax(...)` line).

### Recommended
3) Use `admS <- max(admS, 0)` for consistency/readability.

Here is the corrected block pattern:

```r name=Model1A_ODE_Functions.R
# inside with(...), near top
S <- max(S, 0); C <- max(C, 0); I <- max(I, 0); Xtest <- max(Xtest, 0); Xtreat <- max(Xtreat, 0); D <- max(D, 0)

N <- S + C + I + Xtest + Xtreat
N <- max(N, 1e-12)

# ...
admS <- exits_total - admC - admI
admS <- max(admS, 0)

# ...
missed_false_negative <- flow_Xtest_to_I
missed_untested <- pmax(0, gamma0 * I - flow_I_to_Xtest)
missed_cases <- missed_untested + missed_false_negative
```

Everything else matches your current parameter names (`dis*`, `mort*`, `gamma0`, `gamma`, `theta`, `pi`, etc.).

## @LFell

Is this correct

## Copilot

Yes — this is correct and consistent with your current parameter script. ✅

```text
S <- max(S, 0); C <- max(C, 0); I <- max(I, 0); Xtest <- max(Xtest, 0); Xtreat <- max(Xtreat, 0); D <- max(D, 0)
```

Quick final check:
- non-negative state clamp: ✅
- `N <- max(N, 1e-12)` guard: ✅
- `admS <- max(admS, 0)` guard: ✅
- single `missed_untested` assignment using `pmax`: ✅
- names match parameters (`dis*`, `mort*`, `gamma0`, `gamma`, `theta`, `pi`, `delta`): ✅

You’re good to proceed.

## @LFell

Can you provide an update to Model1A_ODE_Model_Run.qmd that is consistent with the parameters and functions scripts?

## Copilot

Yes — your run script needs updates for renamed states/flows (`X` → `Xtest`/`Xtreat`, etc.).  
Here is a fully consistent version.

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
Xtest0 <- 0
Xtreat0 <- 0
D0 <- 0

mod.init <- c(S = S0, C = C0, I = I0, Xtest = Xtest0, Xtreat = Xtreat0, D = D0)
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
      cum_I_to_Xtest = cumsum(flow_I_to_Xtest),
      cum_Xtest_to_I = cumsum(flow_Xtest_to_I),
      cum_Xtest_to_Xtreat = cumsum(flow_Xtest_to_Xtreat),
      cum_Xtreat_to_S = cumsum(flow_Xtreat_to_S),

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
                                  state_vars = c("S","C","I","Xtest","Xtreat"),
                                  flow_vars = c("flow_S_to_C","flow_C_to_I","flow_I_to_Xtest",
                                                "flow_Xtest_to_I","flow_Xtest_to_Xtreat","flow_Xtreat_to_S",
                                                "total_cases","identified_cases","missed_cases",
                                                "missed_untested","missed_false_negative"),
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
  S2_lowSe = mod.parms.S2,
  S3_fastTT = mod.parms.S3,
  S4_lowSe_fastTT = mod.parms.S4
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
      S, C, I, Xtest, Xtreat, D,
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
  make_annual_summary(scenario_abso$S2_lowSe, "S2_lowSe"),
  make_annual_summary(scenario_abso$S3_fastTT, "S3_fastTT"),
  make_annual_summary(scenario_abso$S4_lowSe_fastTT, "S4_lowSe_fastTT")
)

write_csv(annual_summary_all, "Outputs/annual_summary_all.csv")
knitr::kable(annual_summary_all, caption = "Annual summary (post burn-in)")
```

## 7) Stabilisation summary (states + flows)

```{r}
stability_all <- bind_rows(
  stabilisation_summary(scenario_abso$S1_base, "S1_base"),
  stabilisation_summary(scenario_abso$S2_lowSe, "S2_lowSe"),
  stabilisation_summary(scenario_abso$S3_fastTT, "S3_fastTT"),
  stabilisation_summary(scenario_abso$S4_lowSe_fastTT, "S4_lowSe_fastTT")
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
    final_Xtest = last$Xtest,
    final_Xtreat = last$Xtreat,
    final_D = last$D,
    cum_total_cases = last$cum_total_cases,
    cum_identified_cases = last$cum_identified_cases,
    cum_missed_cases = last$cum_missed_cases,
    cum_missed_untested = last$cum_missed_untested,
    cum_missed_false_negative = last$cum_missed_false_negative,
    cum_side_room_patient_days = last$cum_side_room_patient_days
  )
}))

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
    select(time, S, C, I, Xtest, Xtreat) %>%
    pivot_longer(cols = c(S, C, I, Xtest, Xtreat), names_to = "state", values_to = "count") %>%
    mutate(scenario = nm)
}))

p_states_line <- ggplot(plot_long_states, aes(time, count, color = state)) +
  geom_line() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "State trajectories by scenario", x = "Days post burn-in", y = "Count")
ggsave("Outputs/plot_states_lines_facet.png", p_states_line, width = 12, height = 8, dpi = 300)

p_states_area <- ggplot(plot_long_states, aes(time, count, fill = state)) +
  geom_area() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "State composition over time", x = "Days post burn-in", y = "Count")
ggsave("Outputs/plot_states_area_facet.png", p_states_area, width = 12, height = 8, dpi = 300)

plot_long_flows <- bind_rows(lapply(names(scenario_abso), function(nm) {
  scenario_abso[[nm]] %>%
    select(
      time,
      flow_S_to_C, flow_C_to_I, flow_I_to_Xtest, flow_Xtest_to_I, flow_Xtest_to_Xtreat, flow_Xtreat_to_S,
      total_cases, identified_cases, missed_cases, missed_untested, missed_false_negative,
      side_room_patient_days,
      cum_total_cases, cum_identified_cases, cum_missed_cases,
      cum_missed_untested, cum_missed_false_negative, cum_side_room_patient_days
    ) %>%
    mutate(scenario = nm)
}))

p_daily_flows <- plot_long_flows %>%
  select(time, scenario, flow_S_to_C, flow_C_to_I, flow_I_to_Xtest, flow_Xtest_to_I, flow_Xtest_to_Xtreat, flow_Xtreat_to_S) %>%
  pivot_longer(cols = -c(time, scenario), names_to = "flow", values_to = "value") %>%
  ggplot(aes(time, value, color = flow)) +
  geom_line() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "Key daily flows", x = "Days post burn-in", y = "Flow/day")
ggsave("Outputs/plot_daily_flows_facet.png", p_daily_flows, width = 12, height = 8, dpi = 300)

p_cum_cases <- plot_long_flows %>%
  select(time, scenario, cum_total_cases, cum_identified_cases, cum_missed_cases) %>%
  pivot_longer(cols = -c(time, scenario), names_to = "series", values_to = "value") %>%
  ggplot(aes(time, value, color = series)) +
  geom_line(linewidth = 1) +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "Cumulative cases", x = "Days post burn-in", y = "Cumulative count")
ggsave("Outputs/plot_cumulative_cases_facet.png", p_cum_cases, width = 12, height = 8, dpi = 300)

p_missed_breakdown <- plot_long_flows %>%
  select(time, scenario, cum_missed_untested, cum_missed_false_negative) %>%
  pivot_longer(cols = -c(time, scenario), names_to = "missed_type", values_to = "value") %>%
  ggplot(aes(time, value, color = missed_type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~scenario, scales = "free_y") +
  theme_classic() +
  labs(title = "Cumulative missed cases breakdown", x = "Days post burn-in", y = "Cumulative count")
ggsave("Outputs/plot_cumulative_missed_breakdown_facet.png", p_missed_breakdown, width = 12, height = 8, dpi = 300)

p_side_room <- plot_long_flows %>%
  ggplot(aes(time, side_room_patient_days, color = scenario)) +
  geom_line(linewidth = 1) +
  theme_classic() +
  labs(title = "Side-room occupancy over time", x = "Days post burn-in", y = "Patient-days/day")
ggsave("Outputs/plot_side_room_occupancy.png", p_side_room, width = 10, height = 6, dpi = 300)

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
  S2_lowSe = as_tibble(mod.parms.S2),
  S3_fastTT = as_tibble(mod.parms.S3),
  S4_lowSe_fastTT = as_tibble(mod.parms.S4),
  .id = "scenario"
)
write_csv(all_parms, "Outputs/model_parameters_all_scenarios.csv")
knitr::kable(all_parms, caption = "Model parameters for all scenarios")
```
````

If you want, I can also give a **minimal diff list** (old name → new name) so you can patch your current file quickly instead of replacing the whole qmd.

## @LFell

Can you provide a README.md consistent with all three?

## Copilot

Absolutely — here is a README consistent with your three scripts (`Model1A_ODE_Parameters.R`, `Model1A_ODE_Functions.R`, `Model1A_ODE_Model_Run.qmd`).

````markdown name=README.md
# Model 1A: ODE model of hospital transmission with testing and isolation

## Overview

This repository contains a deterministic compartmental ODE model for in-hospital transmission of a pathogen (e.g., *C. difficile*), with explicit representation of:

- suspected case identification and isolation,
- test sensitivity and turnaround time,
- false negatives returning to general care,
- treated confirmed cases in isolation,
- discharge and mortality flows,
- cumulative missed cases and side-room occupancy.

The model is implemented in R with `deSolve` and run through a Quarto workflow.

---

## Model structure

### State variables

- `S`: susceptible (not colonised or infected)
- `C`: colonised
- `I`: infected (not isolated)
- `Xtest`: isolated, awaiting test result
- `Xtreat`: isolated, confirmed positive, on treatment/recovery
- `D`: cumulative deaths (absorbing count state)

Total in-hospital living census is:

\[
N = S + C + I + Xtest + Xtreat
\]

---

## Key transitions

- `S -> C`: force of colonisation (`lambda`)
- `C -> I`: progression (`alpha`)
- `C -> S`: decolonisation (`mu`)
- `I -> Xtest`: suspected-case isolation (`gamma`)
- `Xtest -> Xtreat`: true positives (`theta`)
- `Xtest -> I`: false negatives (`pi`)
- `Xtreat -> S`: recovery to general ward (`delta`)
- Discharge and mortality from relevant compartments
- Deaths accumulate in `D`

---

## Case accounting outputs

From the ODE function:

- `total_cases = flow_C_to_I + admI`
- `identified_cases = flow_Xtreat_to_S + flow_dis_Xtreat`
- `missed_false_negative = flow_Xtest_to_I`
- `missed_untested = pmax(0, gamma0 * I - flow_I_to_Xtest)`
- `missed_cases = missed_untested + missed_false_negative`
- `side_room_patient_days = Xtest + Xtreat`

---

## Scripts

### 1) `Scripts/Model1A_ODE_Parameters.R`

Defines:
- fixed parameters,
- helper `prob2rate()` conversion,
- validation helpers:
  - `validate_fixed_parameters()`
  - `validate_inputs()`
  - `validate_rate_outputs()`
- scenario builder `make_mod_parms(...)`,
- scenario parameter lists:
  - `mod.parms.S1`
  - `mod.parms.S2`
  - `mod.parms.S3`
  - `mod.parms.S4`

### 2) `Scripts/Model1A_ODE_Functions.R`

Defines:
- `Model1A_ODE(time, state, parms)` for `deSolve::ode()`,
- non-negative clamping of state variables to reduce numerical artifacts,
- flow outputs and derived case outputs.

### 3) `Scripts/Model1A_ODE_Model_Run.qmd`

Runs end-to-end workflow:
- setup and sourcing scripts,
- model execution for all scenarios,
- burn-in filtering,
- cumulative derived outputs,
- annual summaries,
- stabilisation checks,
- negative-value diagnostics,
- final outcomes table (including deltas vs base scenario),
- plots,
- CSV export to `Outputs/`.

---

## Requirements

- R (>= 4.1 recommended)
- Packages:
  - `tidyverse`
  - `deSolve`
  - `scales`
  - `knitr` (for rendering tables in Quarto)

Install if needed:

```r
install.packages(c("tidyverse", "deSolve", "scales", "knitr"))
```

For report rendering:
- Quarto installed and available on PATH.

---

## Run instructions

From project root:

1. Ensure file layout matches paths in the run script:
   - `Scripts/Model1A_ODE_Functions.R`
   - `Scripts/Model1A_ODE_Parameters.R`
   - `Model1A_ODE_Model_Run.qmd`

2. Render the model run report:

```bash
quarto render Model1A_ODE_Model_Run.qmd
```

Or run interactively in RStudio by rendering the `.qmd`.

---

## Outputs

All outputs are written to `Outputs/` (created automatically), including:

- `solution_*.csv` (time-series outputs by scenario)
- `annual_summary_all.csv`
- `stabilisation_summary_states_flows.csv`
- `negative_rows_*.csv`
- `negative_summary_*.csv`
- `final_outcomes_vs_base.csv`
- `model_parameters_all_scenarios.csv`
- plot PNG files

---

## Scenario definitions

Current scenarios differ by test sensitivity and turnaround (with fixed `prop_I_suspected = 0.50`):

- `S1`: `test_sens = 0.80`, `test_turnaround = 5`
- `S2`: `test_sens = 0.50`, `test_turnaround = 5`
- `S3`: `test_sens = 0.80`, `test_turnaround = 3`
- `S4`: `test_sens = 0.50`, `test_turnaround = 3`

---

## Validation and safeguards

Implemented safeguards include:

- explicit `stop()` checks for invalid inputs/parameters,
- warnings for potentially unrealistic values,
- non-negative state clamping in ODE function,
- `N <- max(N, 1e-12)` to prevent divide-by-zero,
- `admS <- max(admS, 0)` admission guard,
- non-negative guard for `missed_untested` via `pmax`.

---

## Notes on units

- Model time is in **half-days** in parameter interpretation.
- Run script time grid currently uses daily steps (`time.step = 1`) and yearly summaries at 365-day intervals.
- Keep this convention consistent when changing timescale assumptions.

---

## Citation / reproducibility note

This is a scenario-analysis model intended for structured comparison of assumptions.
Parameter values should be reviewed against current evidence before policy interpretation.
````
