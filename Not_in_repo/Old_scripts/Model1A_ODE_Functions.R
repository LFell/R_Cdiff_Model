# Model 1A ODE Functions ####

# 1) Introduction ####

# This script defines the ODE function for Model 1A, which is a compartmental model of hospital 
# transmission of a pathogen with the following states:

# The model compartments are:
# S: Susceptible, not colonised or infected
# C: Colonised, carrying the pathogen but not infected
# I: Infected, with symptoms caused by the pathogen, not yet suspected cases or false negatives
# Xtest: Suspected cases, isolated or on a general ward, awaiting test result
# Xtreat: Confirmed cases, isolated and receiving treatment, awaiting recovery & discharge or return 
# to a general ward
# D: Dead (cumulative count of deaths, absorbing state)

# The parameters are defined in the accompanying script "Model1A_ODE_Parameters.R". Together, these 
# scripts are used to run the model and generate outputs in the main quarto document 
# (Model1A_ODE_Model_Run.qmd).


# 2) ODE function ####

Model1A_ODE <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    
    # a) Ensure non-negative state variables ####
    
    # Avoid negative values for S, C and I as these are used in the force of colonisation 
    # calculation and negative values can lead to negative force of colonisation and thus 
    # negative flow from S to C.
    
    S <- max(S, 0)
    C <- max(C, 0)
    I <- max(I, 0)
    
    Xtest <- max(Xtest, 0)
    Xtreat <- max(Xtreat, 0)
    
    D <- max(D, 0)
    
    # b) Calculate total in-hospital population (N) ####
    
    # Total in-hospital census - the number of patients in the hospital at a given time (N)
    # is the sum of those in the S, C, I, Xtest and Xtreat states.
    
    N <- S + C + I + Xtest + Xtreat
    
    # c) Ensure non-negative, non-zero N ####
    
    # Ensure the number of patients in the hospital (N) is not negative or zero to avoid 
    # issues with division in the force of colonisation calculation.
    
    N <- max(N, 1e-12)
    
    # d) Calculate the force of colonisation (lambda) ####
    
    # The force of colonisation (lambda) is the rate at which susceptible patients become colonised.
    
    # It is calculated as the sum of colonisation due to environmental exposure (beta0) and 
    # contributions from colonised, infected, and isolated patients, weighted by their respective 
    # transmission rates (betaC, betaI, betaXtest, betaXtreat) and normalized by the total hospital 
    # census (N).
    
    lambda <- beta0 + ((betaC * C + betaI * I + betaXtest * Xtest + betaXtreat * Xtreat) / N)
    
    # e) Calculate total exits ####
    
    # Total exits (discharges + deaths) from the hospital are calculated as the sum of the exit 
    # rates (discharge + mortality) for each state multiplied by the number of patients in that 
    # state.
    
    # It is assumed there are no discharges from Xtest as patients are awaiting test results.
    
    exits_total <- (disS + mortS) * S +
      (disC + mortC) * C +
      (disI + mortI) * I +
      (mortXtest) * Xtest +             
      (disXtreat + mortXtreat) * Xtreat
    
    # f) Assign admissions to S, C and I ####
    
    # The total number of admissions to the hospital is determined by the total exits (discharges + 
    # deaths) to maintain a fixed census. 
    
    # The admissions are split between the S, C and I states, based on the a fixed import prevalence 
    # of colonised and infected patients, with the remaining admissions to the S state.
    
    # Assume there are no admissions direct to Xtest or Xtreat as patients are only moved there 
    # after admission and subsequent identification as suspected/confirmed cases.
    
    admI <- exits_total * prop_adm_I
    admC <- exits_total * prop_adm_C
    admS <- exits_total - admC - admI
    
    # g) Ensure non-negative admissions to S ####
    
    # Ensure that admissions to S are not negative, which can occur if the import prevalence of 
    # colonised and infected patients is very high and/or if discharge and death rates are low, 
    # leading to a situation where the calculated admissions to C and I exceed the total exits.
    
    admS <- max(admS, 0)
  
    # h) ODEs ####
    
    # Formulae for the derivatives of each state variable (dS, dC, dI, dXtest, dXtreat, dD),
    # based on the flows into and out of each state.
    
    dS <- admS + mu * C + delta * Xtreat - lambda * S - (disS + mortS) * S
    dC <- admC + lambda * S - (alpha + mu + disC + mortC) * C
    dI <- admI + alpha * C + pi * Xtest - (gamma + disI + mortI) * I
    dXtest <- gamma * I - (theta + pi + mortXtest) * Xtest
    dXtreat <- theta * Xtest - (delta + disXtreat + mortXtreat) * Xtreat
    dD <- mortS * S + mortC * C + mortI * I + mortXtest * Xtest + mortXtreat * Xtreat
    
    # i) Flows (instantaneous rates per day) ####
    
    # These are rates per day, evaluated at each ODE time step (half-day).
    
    colonised_in_hospital_rate <- lambda * S
    admitted_colonised_rate <- admC
    decolonised_in_hospital_rate <- mu * C
    total_colonised_incidence_rate <- colonised_in_hospital_rate + admitted_colonised_rate
    
    infected_in_hospital_rate <- alpha * C
    admitted_infected_rate <- admI
    total_infected_incidence_rate <- infected_in_hospital_rate + admitted_infected_rate
    
    suspected_cases_rate <- gamma * I
    missed_untested_cases_rate <- gamma0 * I  
    confirmed_cases_rate <- theta * Xtest
    missed_false_negative_cases_rate <- pi * Xtest
    total_missed_cases_rate <- missed_untested_cases_rate + missed_false_negative_cases_rate
    recovered_cases_rate <- (delta * Xtreat) + (disXtreat * Xtreat)
    
    discharges_susceptible_rate <- disS * S
    discharges_colonised_rate <- disC * C
    discharges_infected_rate <- disI * I
    # no discharges from Xtest as patients are awaiting test results
    discharges_recovered_rate <- disXtreat * Xtreat
    
    deaths_susceptible_rate <- mortS * S
    deaths_colonised_rate <- mortC * C
    deaths_infected_rate <- mortI * I
    deaths_suspected_cases_rate <- mortXtest * Xtest
    deaths_confirmed_cases_rate <- mortXtreat * Xtreat
   
    total_admissions_rate <- admS + admC + admI
    total_discharges_rate <- (disS * S) + (disC * C) + (disI * I) + (disXtreat * Xtreat)
    total_deaths_rate <- (mortS * S) + (mortC * C) + (mortI * I) + (mortXtest * Xtest) + (mortXtreat * Xtreat)
    total_exits_rate <- total_discharges_rate + total_deaths_rate
    
    # j) Side-room occupancy ####
    
    # Number of patients in Xtest and/or Xtreat at each time step (half-day)
    # Include Xtest and Xtreat if isolate_before_test = TRUE
    # Include only Xtreat when isolate_before_test = FALSE
    
    if (isolate_before_test) {
      side_room_patient_days <- Xtest + Xtreat
      } else {
      side_room_patient_days <- Xtreat
      }
    
    # k) Return derivatives and outputs as a list ####
    list(
      c(dS, dC, dI, dXtest, dXtreat, dD),
      c(
        N = N, 
        lambda = lambda,
        colonised_in_hospital_rate = colonised_in_hospital_rate,
        admitted_colonised_rate = admitted_colonised_rate,
        decolonised_in_hospital_rate = decolonised_in_hospital_rate,
        total_colonised_incidence_rate = total_colonised_incidence_rate,
        
        infected_in_hospital_rate = infected_in_hospital_rate,
        admitted_infected_rate = admitted_infected_rate,
        total_infected_incidence_rate = total_infected_incidence_rate,
        
        suspected_cases_rate = suspected_cases_rate,
        confirmed_cases_rate = confirmed_cases_rate,
        missed_untested_cases_rate = missed_untested_cases_rate,
        missed_false_negative_cases_rate = missed_false_negative_cases_rate,
        total_missed_cases_rate = total_missed_cases_rate,
        recovered_cases_rate = recovered_cases_rate,
        
        discharges_susceptible_rate = discharges_susceptible_rate,
        discharges_colonised_rate = discharges_colonised_rate,
        discharges_infected_rate = discharges_infected_rate,
        discharges_recovered_rate = discharges_recovered_rate,
        
        deaths_susceptible_rate = deaths_susceptible_rate,
        deaths_colonised_rate = deaths_colonised_rate,
        deaths_infected_rate = deaths_infected_rate,
        deaths_suspected_cases_rate = deaths_suspected_cases_rate,
        deaths_confirmed_cases_rate = deaths_confirmed_cases_rate,
        
        total_admissions_rate = total_admissions_rate,
        total_discharges_rate = total_discharges_rate,
        total_deaths_rate = total_deaths_rate,
        total_exits_rate = total_exits_rate,
        
        side_room_patient_days = side_room_patient_days
        ))
  })
}

# End of script ####