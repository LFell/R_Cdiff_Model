#CPE Compartmental Model Functions ####

# 1) Introduction ####

# This script contains the functions for the CPE compartmental model, including 
# - the ODE function that defines the model structure and equations, 
# - a function to generate model outputs for validation and cost calculations,
# - a function to convert probabilities to rates,
# - a function to generate the beta distribution for a proportion parameter with uncertainty.

# These functions are called in the main model script, CPE_ODE_Model.R, which runs 
# the model and generates outputs for the base case and sensitivity analyses.


# ========================================================== #

# 2) Define ODE model functions ####


#  Ordinary differential equation (ODE)/compartmental model function #

CPE.ODE.model <- function(times, state, parms) {
  #compartmental model with 12 hospital compartments: S, C and I on a general ward; St, Ct and It awaiting test results; Sx, Cx, and Ix tested positive; Sy, Cy, and Iy tested negative
  #with transmission in the hospital
  #and movement of patients to and from the community
  #with 8 community compartments; Sda (S recent discharge), Sdb (S non-recent discharge), Cda (unknown colonised recent discharge), Cdb (unknown colonised non-recent discharge), Sdxa (susceptible recent discharge w/+ve flag), Sdxb (susceptible non-recent discharge w/+ve flag), Cdxa (colonised recent discharge w/+ve flag), Cdxb (colonised non-recent discharge w/+ve flag)
  
# 2.1) Define states ####

  #in-hospital states
  Sg <- state[["Sg"]] #susceptible general hospital population
  Cg <- state[["Cg"]] #colonised general hospital population
  Ig <- state[["Ig"]] #infected general hospital population
  St <- state[["St"]] #susceptible tested
  Ct <- state[["Ct"]] #colonised tested
  It <- state[["It"]] #infected tested
  Sx <- state[["Sx"]] #susceptible positive result
  Cx <- state[["Cx"]] #colonised positive result
  Ix <- state[["Ix"]] #infected positive result
  Sy <- state[["Sy"]] #susceptible negative result
  Cy <- state[["Cy"]] #colonised negative result
  Iy <- state[["Iy"]] #infected negative result
  Tot_hos <- state[["Tot_hos"]] #total hospital patient population sense check
  
  #community states
  Sdxa <- state[["Sdxa"]] #susceptible, positive flag recently discharged
  Cdxa <- state[["Cdxa"]] #colonised, positive flag recently discharged
  Sdxb <- state[["Sdxb"]] #susceptible, positive flag non-recently discharged
  Cdxb <- state[["Cdxb"]] #colonised, positive flag non-recently discharged
  Sda <- state[["Sda"]] #susceptible, recently discharged
  Cda <- state[["Cda"]] #colonised, unknown recently discharged
  Sdb <- state[["Sdb"]] #susceptible, non-recently discharged
  Cdb <- state[["Cdb"]] #colonised, unknown non-recently discharged
  Tot_com <- state[["Tot_com"]] #total community population sense check
  
# 2.2) Extract parameters ####

  #transitions within hospital
  beta0 <- parms[["beta0"]] #transmission rate from other sources, e.g. sinks
  beta1 <- parms[["beta1"]] #transmission rate from C and I patients
  Zg <- parms[["Zg"]] #effectiveness of environmental cleaning in general wards
  Zt <- parms[["Zt"]] #effectiveness of environmental cleaning for patients awaiting results
  Zx <- parms[["Zx"]] #effectiveness of environmental cleaning for patients tested CRO positive
  Zy <- parms[["Zy"]] #effectiveness of environmental cleaning for patients tested CRO negative
  chi <- parms[["chi"]] #effectiveness of IPC measures for patients awaiting test results
  epsilon <- parms[["epsilon"]] #effectiveness of IPC measures for patients tested CRO positive
  alpha <- parms[["alpha"]] #infection progression rate
  gammasg <- parms[["gammasg"]] #testing rate of susceptible patients
  gammacg <- parms[["gammacg"]] #testing rate of colonised patients
  gammaig <- parms[["gammaig"]] #testing rate of infected patients
  gammasx <- parms[["gammasx"]] #testing rate of susceptible tested positive patients
  gammacx <- parms[["gammacx"]] #testing rate of colonised tested positive patients
  gammaix <- parms[["gammaix"]] #testing rate of infected tested positive patients
  gammasy <- parms[["gammasy"]] #testing rate of susceptible tested negative patients
  gammacy <- parms[["gammacy"]] #testing rate of colonised tested negative patients
  gammaiy <- parms[["gammaiy"]] #testing rate of infected tested negative patients
  n <- parms[["n"]] #1/test turnaround time (days)
  theta <- parms[["theta"]] #test sensitivity
  phi <- parms[["phi"]] #test specificity
  
  #discharged from hospital
  dsg <- parms[["dsg"]] #discharge rate of susceptible patients
  dcg <- parms[["dcg"]] #discharge rate of colonised patients
  dig <- parms[["dig"]] #discharge rate of infected patients
  dst <- parms[["dst"]] #discharge rate of susceptible testing patients
  dct <- parms[["dct"]] #discharge rate of colonised testing patients
  dit <- parms[["dit"]] #discharge rate of infected testing patients
  dsx <- parms[["dsx"]] #discharge rate of susceptible tested positive patients
  dcx <- parms[["dcx"]] #discharge rate of colonised tested positive patients
  dix <- parms[["dix"]] #discharge rate of infected tested positive patients
  dsy <- parms[["dsy"]] #discharge rate of susceptible tested negative patients
  dcy <- parms[["dcy"]] #discharge rate of colonised tested negative patients
  diy <- parms[["diy"]] #discharge rate of infected tested negative patients
  
  #transitions within community
  psi <- parms[["psi"]] #1/recently discharged period
  delta <- parms[["delta"]] #1/non-recent discharged period - recent discharge period
  omega <- parms[["omega"]] #1/CRO positive flag period on patient records - recent discharge period
  mu <- parms[["mu"]] #decolonisation rate in the community
  
  #hospital readmission from recently discharged
  rhos <- parms[["rhos"]] #readmission rate amongst recently discharged susceptible patients
  rhoc <- parms[["rhoc"]] #readmission rate amongst recently discharged colonised patients
  rhop <- parms[["rhop"]] #readmission rate among non-recently discharged patients, i.e. general population
  nudx <- parms[["nudx"]] #proportion discharged with positive flag screened on readmission
  nud <- parms[["nud"]] #proportion discharged without positive flag screened on readmission
  
  #hospital admission from wider community
  sigmaec <- parms[["sigmaec"]] #prevalence of colonisation upon admission from England wider community
  sigmanec <- parms[["sigmanec"]] #prevalence of colonisation upon admission from non-England wider community
  tau <- parms[["tau"]] #probability of colonised patient presenting with infection
  zeta <- parms[["zeta"]] #proportion of all admissions from non-England geography
  nuec <- parms[["nuec"]] #proportion of England wider community patients screened on admission
  nunec <- parms[["nunec"]] #proportion of non-England wider community patients screened on admission
  
  #mortality
  msg <- parms[["msg"]] #mortality rate susceptible patients
  mcg <- parms[["mcg"]] #mortality rate colonised patients
  mig <- parms[["mig"]] #mortality rate infected patients
  mst <- parms[["mst"]] #mortality rate of susceptible testing patients
  mct <- parms[["mct"]] #mortality rate of colonised testing patients
  mit <- parms[["mit"]] #mortality rate of infected testing patients
  msx <- parms[["msx"]] #mortality rate of susceptible tested positive patients
  mcx <- parms[["mcx"]] #mortality rate of colonised tested positive patients
  mix <- parms[["mix"]] #mortality rate of infected tested positive patients
  msy <- parms[["msy"]] #mortality rate of susceptible tested negative patients
  mcy <- parms[["mcy"]] #mortality rate of colonised tested negative patients
  miy <- parms[["miy"]] #mortality rate of infected tested negative patients
  msda <- parms[["msda"]] #mortality rate recently discharged susceptible
  mcda <- parms[["mcda"]] #mortality rate recently discharged unknown colonised
  msdb <- parms[["msdb"]] #mortality rate non-recently discharged susceptible
  mcdb <- parms[["mcdb"]] #mortality rate non-recently discharged unknown colonised
  msdxa <- parms[["msdxa"]] #mortality rate recently discharged susceptible positive flag
  mcdxa <- parms[["mcdxa"]] #mortality rate recently discharged colonised positive flag
  msdxb <- parms[["msdxb"]] #mortality rate non-recently discharged susceptible positive flag
  mcdxb <- parms[["mcdxb"]] #mortality rate non-recently discharged colonised positive flag
  
# Hospital costs #
  
  #bed-day unit costs
  cGenWardBed <- parms[["cGenWardBed"]] #general ward bed day cost, per patient per day
  cICUBed <- parms[["cICUBed"]] #ICU bed day cost, per patient per day
  
  #screening test unit costs
  cPCRTest <- parms[["cPCRTest"]] #one-off cost of a PCR test, per test
  cCultureTest <- parms[["cCultureTest"]] #one-off cost of a culture test, per test
  cScreenTest <- parms[["cScreenTest"]] #unit cost for the test of interest (either culture or PCR) to be defined when running model
  #NOTE: ensure cPCRTest incorporates unit cost for culture when running the model, as all PCRs are followed by culture test to determine susceptibility
  
  #screening test staff unit opportunity costs
  cScreenStaff <- parms[["cScreenStaff"]] #unit cost of staff time opportunity cost, per test
  
  #treatment (i.e., abx) unit costs
  cInfectionTreat <- parms[["cInfectionTreat"]] #infection treatment cost, per patient per day
  cInfectionTreatMeet <- parms[["cInfectionTreatMeet"]] #infection treatment staff costs (i.e., multidisciplinary meetings re patient treatment), per new infection
  
  #toxicity test unit costs
  cToxInit <- parms[["cToxInit"]] #initial toxicity test cost, per new infection
  cToxOngo <- parms[["cToxOngo"]] #ongoing toxicity test cost, per patient per day
  
  #IPC costs
  cPPEequip <- parms[["cPPEequip"]] #PPE (gloves and aprons) for contact precautions, per patient per day
  cIPCStaff <- parms[["cIPCStaff"]] #staff opportunity cost of IPC measures, per patient per day
  cRoomClean <- parms[["cRoomClean"]] #cost of room cleaning for patients under IPC measures, per patient per day
  cStockDisp <- parms[["cStockDisp"]] #cost of stock disposal after regular discharge, per patient discharge (dead or alive)
  cStockDispICU <- parms[["cStockDispICU"]] #cost of stock disposal after regular discharge, per patient discharge (dead or alive)
  cInfectWaste <- parms[["cInfectWaste"]] #cost of infectious waste, per isolating patient per day
  
  #mortality review panel
  cMortRevPanel <- parms[["cMortRevPanel"]] #cost of mortality review panel, per positive patient death
  
  #outbreak costs
  cOutbreakpI <- parms[["cOutbreakpI"]] #cost of outbreak, per new infection
  
# 2.3) Calculate number of patients moving in and out of hospital ####

  Nh <- Sg + Cg + Ig + St + Ct + It + Sx + Cx + Ix + Sy + Cy + Iy #total hospital population size, needed for the FOI equations below
  D <- dsg * Sg + dcg * Cg + dig * Ig + dst * St + dct * Ct + dit * It + dsx *
    Sx + dcx * Cx + dix * Ix + dsy * Sy + dcy * Cy + diy * Iy #patients discharged
  M <- msg * Sg + mcg * Cg + mig * Ig + mst * St + mct * Ct + mit * It + msx *
    Sx + mcx * Cx + mix * Ix + msy * Sy + mcy * Cy + miy * Iy #patients who died in hospital
  R <- rhos * Sdxa + rhoc * Cdxa + rhop * Sdxb + rhop * Cdxb + rhos * Sda + rhoc *
    Cda + rhop * Sdb + rhop * Cdb #patients readmitted
  H <- (D + M) - R #patients hospitalised from the wider community
  A <- H + R #total hospital admissions; i.e., admissions from wider community + readmissions
  Anec <- zeta * A #hospital admissions from non-England geography
  Aec <- H - Anec #hospital admissions from England wider community
  Aer <- R #hospital admissions from England readmissions; NOTE: all readmissions assumed to be from England
  Asc <- A - (Anec + Aec + Aer) #sense check for all admissions; Asc should = 0
  patient_movement <- c(D, M, R, H, A, Anec, Aec, Aer, Asc) #collate vars into single object
  names(patient_movement) <- c("D", "M", "R", "H", "A", "Anec", "Aec", "Aer", "Asc") #name vars within object
  
# 2.4) Define force of infection ####
  
  lambdag <- (1 - Zg) * beta0 + (beta1 * (Cg + Ig + Cy + Iy + (1 - chi) *
                                            (Ct + It) + (1 - epsilon) * (Cx + Ix)) / Nh)
  lambdat <- (1 - Zt) * beta0 + ((1 - chi) * beta1 * (Cg + Ig + Cy + Iy +
                                                        Ct + It + (1 - epsilon) * (Cx + Ix)) / Nh)
  lambdax <- (1 - Zx) * beta0 + ((1 - epsilon) * beta1 * (Cg + Ig + Cy +
                                                            Iy + (1 - chi) * (Ct + It) + Cx + Ix) / Nh)
  lambday <- (1 - Zy) * beta0 + (beta1 * (Cg + Ig + Cy + Iy + (1 - chi) *
                                            (Ct + It) + (1 - epsilon) * (Cx + Ix)) / Nh)
  
  #----------------------#
# 2.5) Differential equations ####
  #----------------------#
  
  #in-hospital states
  dSg <- (1 - sigmaec) * (1 - nuec) * Aec + (1 - sigmanec) * (1 - nunec) *
    Anec + rhos * (1 - nud) * Sda + rhop * (1 - nud) * Sdb + rhos * (1 - nudx) *
    Sdxa + rhop * (1 - nudx) * Sdxb - Sg * (lambdag + gammasg + dsg + msg) #admitted from wider community and not screened + readmitted from recent discharge no flag and not screened + readmitted from non-recent discharge no flag and not screened + readmitted from recent discharge with flag and not screened + readmitted from non-recent discharge with flag and not screened – (FOI + in-hospital testing + discharge + mortality)
  dCg <- sigmaec * (1 - tau) * (1 - nuec) * Aec + sigmanec * (1 - tau) *
    (1 - nunec) * Anec + rhoc * (1 - tau) * (1 - nud) * Cda + rhop * (1 - tau) *
    (1 - nud) * Cdb + rhoc * (1 - tau) * (1 - nudx) * Cdxa + rhop * (1 - tau) *
    (1 - nudx) * Cdxb + lambdag * Sg - Cg * (alpha + gammacg + dcg + mcg) #admitted from wider community and not screened + readmitted from recent discharge no flag and not screened + readmitted from non-recent discharge no flag and not screened + readmitted from recent discharge with flag and not screened + readmitted from non-recent discharge with flag and not screened + FOI – (infection progression + in-hospital testing + discharge + mortality)
  dIg <- sigmaec * tau * (1 - nuec) * Aec + sigmanec * tau * (1 - nunec) *
    Anec + rhoc * tau * (1 - nud) * Cda + rhop * tau * (1 - nud) * Cdb + rhoc *
    tau * (1 - nudx) * Cdxa + rhop * tau * (1 - nudx) * Cdxb + alpha * Cg - Ig *
    (gammaig + dig + mig) #admitted from wider community and not screened + readmitted from recent discharge no flag and not screened + readmitted from non-recent discharge no flag and not screened + readmitted from recent discharge with flag and not screened + readmitted from non-recent discharge with flag and not screened + infection progression – (in-hospital testing + discharge + mortality)
  dSt <- (1 - sigmaec) * nuec * Aec + (1 - sigmanec) * nunec * Anec + rhos *
    nud * Sda + rhop * nud * Sdb + rhos * nudx * Sdxa + rhop * nudx * Sdxb + gammasg *
    Sg + gammasx * Sx + gammasy * Sy - St * (lambdat + n + dst + mst) #admitted from wider community and screened + readmitted from recent discharge no flag and screened + readmitted from non-recent discharge no flag and screened + readmitted from recent discharge with flag and screened + readmitted from non-recent discharge with flag and screened + testing – (FOI + test results + discharge + mortality)
  dCt <- sigmaec * (1 - tau) * nuec * Aec + sigmanec * (1 - tau) * nunec *
    Anec + rhoc * (1 - tau) * nud * Cda + rhop * (1 - tau) * nud * Cdb + rhoc *
    (1 - tau) * nudx * Cdxa + rhop * (1 - tau) * nudx * Cdxb + gammacg * Cg + gammacx *
    Cx + gammacy * Cy + lambdat * St - Ct * (alpha + n + dct + mct) #admitted from wider community and screened + readmitted from recent discharge no flag and screened + readmitted from non-recent discharge no flag and screened + readmitted from recent discharge with flag and screened + readmitted from non-recent discharge with flag and screened + testing + FOI – (infection progression + test results + discharge + mortality)
  dIt <- sigmaec * tau * nuec * Aec + sigmanec * tau * nunec * Anec + rhoc *
    tau * nud * Cda + rhop * tau * nud * Cdb + rhoc * tau * nudx * Cdxa + rhop *
    tau * nudx * Cdxb + gammaig * Ig + gammaix * Ix + gammaiy * Iy + alpha *
    Ct - It * (n + dit + mit) #admitted from wider community and screened + readmitted from recent discharge no flag and screened + readmitted from non-recent discharge no flag and screened + readmitted from recent discharge with flag and screened + readmitted from non-recent discharge with flag and screened + testing + infection progression – (test results + discharge + mortality)
  dSx <- (1 - phi) * n * St - Sx * (lambdax + gammasx + dsx + msx) #false positive – (FOI + tested + discharged + mortality)
  dCx <- theta * n * Ct + lambdax * Sx - Cx * (alpha + gammacx + dcx + mcx) #true positive + FOI – (infection progression + tested + discharged + mortality)
  dIx <- theta * n * It + alpha * Cx - Ix * (gammaix + dix + mix) #true positive + infection progression – (tested + discharged + mortality)
  dSy <- phi * n * St - Sy * (lambday + gammasy + dsy + msy) #true negative – (FOI + tested + discharged + mortality)
  dCy <- (1 - theta) * n * Ct + lambday * Sy - Cy * (alpha + gammacy + dcy +
                                                       mcy) #false negative + FOI – (infection progression + tested + discharged + mortality)
  dIy <- (1 - theta) * n * It + alpha * Cy - Iy * (gammaiy + diy + miy) #false negative + infection progression – (tested + discharged + mortality)
  dTot_hos <- dSg + dCg + dIg + dSt + dCt + dIt + dSx + dCx + dIx + dSy + dCy + dIy #total hospital population
  
  #community states
  dSdxa <- dsx * Sx + (1 - phi) * dst * St + mu * Cdxa - Sdxa * (rhos +
                                                                   psi + msdxa) #tested positive discharged + false positive discharged + decolonisation – (readmitted + move to non-recent + mortality)
  dCdxa <- dcx * Cx + dix * Ix + theta * dct * Ct + theta * dit * It - Cdxa *
    (rhoc + mu + psi + mcdxa) #tested positive discharged + true positive discharged – (readmitted + decolonised + move to non-recent + mortality)
  dSdxb <- psi * Sdxa + mu * Cdxb - Sdxb * (rhop + omega + msdxb) #move to non-recent + decolonisation – (readmission + move to wider community + mortality)
  dCdxb <- psi * Cdxa - Cdxb * (rhop + mu + omega + mcdxb) #move to non-recent – (readmission + decolonisation + move to wider community + mortality)
  dSda <- dsg * Sg + phi * dst * St + dsy * Sy + mu * Cda - Sda * (rhos +
                                                                     psi + msda) #not tested discharged + true negative discharged + tested negative discharged + decolonisation – (readmission + move to non-recent + mortality)
  dCda <- dcg * Cg + dig * Ig + (1 - theta) * dct * Ct + (1 - theta) * dit *
    It + dcy * Cy + diy * Iy - Cda * (rhoc + mu + psi + mcda)
  dSdb <- psi * Sda + mu * Cdb - Sdb * (rhop + delta + msdb)
  dCdb <- psi * Cda - Cdb * (rhop + mu + delta + mcdb)
  dTot_com <- dSdxa + dCdxa + dSdxb + dCdxb + dSda + dCda + dSdb + dCdb
  
# 2.6) Create an object with all the necessary compartments ####
  
  compartments <- c(
    dSg,
    dCg,
    dIg,
    dSt,
    dCt,
    dIt,
    dSx,
    dCx,
    dIx,
    dSy,
    dCy,
    dIy,
    dTot_hos,
    #hospital compartments
    dSdxa,
    dCdxa,
    dSdxb,
    dCdxb,
    dSda,
    dCda,
    dSdb,
    dCdb,
    dTot_com #community compartments
  )
  
# 2.7) Generate model outputs for validation ####
  
# New cases
  
  #new susceptible case in hospital (excluding S -> S movements within hospital)
  new_case_s <- (1 - sigmaec) * (1 - nuec) * Aec + (1 - sigmanec) * (1 -
                                                                       nunec) * Anec + rhos * (1 - nud) * Sda + rhop * (1 - nud) * Sdb + rhos *
    (1 - nudx) * Sdxa + rhop * (1 - nudx) * Sdxb + #number moving into Sg from community
    (1 - sigmaec) * nuec * Aec + (1 - sigmanec) * nunec * Anec + rhos *
    nud * Sda + rhop * nud * Sdb + rhos * nudx * Sdxa + rhop * nudx * Sdxb #number moving into St from community
  
  #new colonised case in hospital (excluding C -> C movements within hospital)
  new_case_c <- (lambdag * Sg) + (lambdat * St) + (lambdax * Sx) + (lambday *
                                                                      Sy) + #number moving from S to C in hospital
    sigmaec * (1 - tau) * (1 - nuec) * Aec + sigmanec * (1 - tau) * (1 -
                                                                       nunec) * Anec + rhoc * (1 - tau) * (1 - nud) * Cda + rhop * (1 - tau) *
    (1 - nud) * Cdb + rhoc * (1 - tau) * (1 - nudx) * Cdxa + rhop * (1 - tau) *
    (1 - nudx) * Cdxb + #number moving into Cg from community
    sigmaec * (1 - tau) * nuec * Aec + sigmanec * (1 - tau) * nunec *
    Anec + rhoc * (1 - tau) * nud * Cda + rhop * (1 - tau) * nud * Cdb + rhoc *
    (1 - tau) * nudx * Cdxa + rhop * (1 - tau) * nudx * Cdxb #number moving into Ct from community
  
  #new acquired colonisation
  new_case_acquired_c <- (lambdag * Sg) + (lambdat * St) + (lambdax *
                                                              Sx) + (lambday * Sy) #number moving from S to C in hospital
  
  #new imported colonisation
  new_case_import_c <-
    sigmaec * (1 - tau) * (1 - nuec) * Aec + sigmanec * (1 - tau) * (1 -
                                                                       nunec) * Anec + rhoc * (1 - tau) * (1 - nud) * Cda + rhop * (1 - tau) *
    (1 - nud) * Cdb + rhoc * (1 - tau) * (1 - nudx) * Cdxa + rhop * (1 - tau) *
    (1 - nudx) * Cdxb + #number moving into Cg from community
    sigmaec * (1 - tau) * nuec * Aec + sigmanec * (1 - tau) * nunec *
    Anec + rhoc * (1 - tau) * nud * Cda + rhop * (1 - tau) * nud * Cdb + rhoc *
    (1 - tau) * nudx * Cdxa + rhop * (1 - tau) * nudx * Cdxb #number moving into Ct from community
  
  #new infected case in hospital (excluding I -> I movements within hospital)
  new_case_i <- alpha * (Cg + Ct + Cx + Cy) + #number moving from C to I in hospital
    sigmaec * tau * (1 - nuec) * Aec + sigmanec * tau * (1 - nunec) *
    Anec + rhoc * tau * (1 - nud) * Cda + rhop * tau * (1 - nud) * Cdb + rhoc *
    tau * (1 - nudx) * Cdxa + rhop * tau * (1 - nudx) * Cdxb + #number moving into Ig from community
    sigmaec * tau * nuec * Aec + sigmanec * tau * nunec * Anec + rhoc *
    tau * nud * Cda + rhop * tau * nud * Cdb + rhoc * tau * nudx * Cdxa + rhop *
    tau * nudx * Cdxb #number moving into It from community
  
  #new acquired infection
  new_case_acquired_i <- alpha * (Cg + Ct + Cx + Cy) #number moving from C to I in hospital
  
  #new imported infection
  new_case_import_i <-
    sigmaec * tau * (1 - nuec) * Aec + sigmanec * tau * (1 - nunec) *
    Anec + rhoc * tau * (1 - nud) * Cda + rhop * tau * (1 - nud) * Cdb + rhoc *
    tau * (1 - nudx) * Cdxa + rhop * tau * (1 - nudx) * Cdxb + #number moving into Ig from community
    sigmaec * tau * nuec * Aec + sigmanec * tau * nunec * Anec + rhoc *
    tau * nud * Cda + rhop * tau * nud * Cdb + rhoc * tau * nudx * Cdxa + rhop *
    tau * nudx * Cdxb #number moving into It from community
  
# New known cases
  
  #new known colonised cases
  new_case_known_c <- theta * n * Ct + lambdax * Sx #number of colonised testing positive + number susceptible that have tested +ve and then become colonised
  
  #new known infected cases
  new_case_known_i <- theta * n * It + alpha * Cx #number of infected testing positive + number colonised that have tested +ve and then progress to infection
  
# New discharges
  
  #new susceptible hospital discharge
  new_dis_s <- dsg * Sg + dst * St + dsx * Sx + dsy * Sy
  
  #new susceptible hospital discharge
  new_dis_c <- dcg * Cg + dct * Ct + dcx * Cx + dcy * Cy
  
  #new susceptible hospital discharge
  new_dis_i <- dig * Ig + dit * It + dix * Ix + diy * Iy
  
# New hospital deaths
  
  #new hospital deaths from S patients
  new_mort_hosp_s <- msg * Sg + mst * St + msx * Sx + msy * Sy
  
  #new hospital deaths from C patients
  new_mort_hosp_c <- mcg * Cg + mct * Ct + mcx * Cx + mcy * Cy
  
  #new hospital deaths from I patients
  new_mort_hosp_i <- mig * Ig + mit * It + mix * Ix + miy * Iy
  
# New readmissions
  
  #new readmission within 30 days from S w/+ve flag
  new_readmit30_s_flag <- rhos * Sdxa
  
  #new readmission within 30 days from C w/+ve flag
  new_readmit30_c_flag <- rhoc * Cdxa
  
  #new readmission within 30 days from S w/o +ve flag
  new_readmit30_s <- rhos * Sda
  
  #new readmission within 30 days from C w/o +ve flag
  new_readmit30_c <- rhoc * Cda
  
  #new readmission 31-365 days from S w/+ve flag
  new_readmit365_s_flag <- rhop * Sdxb
  
  #new readmission 31-365 days from C w/+ve flag
  new_readmit365_c_flag <- rhop * Cdxb
  
  #new readmission 31-365 days from S w/o +ve flag
  new_readmit365_s <- rhop * Sdb
  
  #new readmission 31-365 days from C w/o +ve flag
  new_readmit365_c <- rhop * Cdb
  
  #new admissions for wider community
  
  #new S admitted from wider community
  new_admit_com_s <- (1 - sigmaec) * (1 - nuec) * Aec + (1 - sigmanec) *
    (1 - nunec) * Anec + #number moving into Sg from wider community
    (1 - sigmaec) * nuec * Aec + (1 - sigmanec) * nunec * Anec #number moving into St from wider community
  
  #new C admitted from wider community
  new_admit_com_c <- sigmaec * (1 - tau) * (1 - nuec) * Aec + sigmanec *
    (1 - tau) * (1 - nunec) * Anec + #number moving into Cg from wider community
    sigmaec * (1 - tau) * nuec * Aec + sigmanec * (1 - tau) * nunec *
    Anec  #number moving into Ct from wider community
  
  #new I admitted from wider community
  new_admit_com_i <- sigmaec * tau * (1 - nuec) * Aec + sigmanec * tau *
    (1 - nunec) * Anec + #number moving into Ig from wider community
    sigmaec * tau * nuec * Aec + sigmanec * tau * nunec * Anec #number moving into It from wider community
  
  # #exploring readmission among those in community with CPE+ flag
  #
  #   #Sdxa: susceptible recently discharged with +ve flag
  #   sdxa_in_dis <- dsx*Sx + (1-phi)*dst*St #discharged positive from hospital
  #   sdxa_in_decol <- mu*Cdxa #decolonised
  #
  #   sdxa_out_read <- Sdxa*rhos #readmitted to hospital
  #   sdxa_out_nonrec <- Sdxa*psi #moved to non-recent discharge
  #   sdxa_out_mort <- Sdxa*msdxa #died
  #
  #   #Cdxa: colonised recently discharged with +ve flag
  #   cdxa_in_dis <- dcx*Cx + dix*Ix + theta*dct*Ct + theta*dit*It #discharged positive from hospital
  #
  #   cdxa_out_read <- Cdxa*rhoc #readmitted to hospital
  #   cdxa_out_decol <- Cdxa*mu #decolonised
  #   cdxa_out_nonrec <- Cdxa*psi #moved to non-recent discharge
  #   cdxa_out_mort <- Cdxa*mcdxa #died
  #
  #   #Sdxb: susceptible non-recently discharged with +ve flag
  #   sdxb_in_nonrec <- psi*Sdxa #moved to non-recent discharge
  #   sdxb_in_decol <- mu*Cdxb #decolonised
  #
  #   sdxb_out_read <- Sdxb*rhop #readmitted to hospital
  #   sdxb_out_com <- Sdxb*omega #moved to non-recent discharge
  #   sdxb_out_mort <- Sdxb*msdxb #died
  #
  #   #Cdxb: colonised non-recently discharged with +ve flag
  #   cdxb_in_nonrec <- psi*Cdxa #moved to non-recent discharge
  #
  #   cdxb_out_read <- Cdxb*rhop #readmitted to hospital
  #   cdxb_out_decol <- Cdxb*mu #decolonised
  #   cdxb_out_com <- Cdxb*omega #moved to wider community
  #   cdxb_out_mort <- Cdxb*mcdxb #died
  
# 2.8) Create an object with all the necessary variables ####
  
  validation <- c(
    new_case_s,
    new_case_c,
    new_case_acquired_c,
    new_case_import_c,
    new_case_i,
    new_case_acquired_i,
    new_case_import_i,
    #new cases
    new_case_known_c,
    new_case_known_i,
    #new known cases
    new_dis_s,
    new_dis_c,
    new_dis_i,
    #new discharges
    new_mort_hosp_s,
    new_mort_hosp_c,
    new_mort_hosp_i,
    #new hospital deaths
    new_readmit30_s_flag,
    new_readmit30_c_flag,
    new_readmit30_s,
    new_readmit30_c,
    #new hospital readmissions within 30 days
    new_readmit365_s_flag,
    new_readmit365_c_flag,
    new_readmit365_s,
    new_readmit365_c,
    #new hospital readmissions 31-365 days
    new_admit_com_s,
    new_admit_com_c,
    new_admit_com_i #new admissions from wider community
  )
  # sdxa_in_dis, sdxa_in_decol, sdxa_out_read, sdxa_out_nonrec, sdxa_out_mort,
  # cdxa_in_dis, cdxa_out_read, cdxa_out_decol, cdxa_out_nonrec, cdxa_out_mort,
  # sdxb_in_nonrec, sdxb_in_decol, sdxb_out_read, sdxb_out_com, sdxb_out_mort,
  # cdxb_in_nonrec, cdxb_out_read, cdxb_out_decol, cdxb_out_com, cdxb_out_mort
  # )
  
# 2.9) Assign names to the model validation output ####
  
  names(validation) <- c(
    "new_case_s",
    "new_case_c",
    "new_case_acquired_c",
    "new_case_import_c",
    "new_case_i",
    "new_case_acquired_i",
    "new_case_import_i",
    "new_case_known_c",
    "new_case_known_i",
    "new_dis_s",
    "new_dis_c",
    "new_dis_i",
    "new_mort_hosp_s",
    "new_mort_hosp_c",
    "new_mort_hosp_i",
    "new_readmit30_s_flag",
    "new_readmit30_c_flag",
    "new_readmit30_s",
    "new_readmit30_c",
    "new_readmit365_s_flag",
    "new_readmit365_c_flag",
    "new_readmit365_s",
    "new_readmit365_c",
    "new_admit_com_s",
    "new_admit_com_c",
    "new_admit_com_i"
  )
  # "sdxa_in_dis", "sdxa_in_decol", "sdxa_out_read", "sdxa_out_nonrec", "sdxa_out_mort",
  # "cdxa_in_dis", "cdxa_out_read", "cdxa_out_decol", "cdxa_out_nonrec", "cdxa_out_mort",
  # "sdxb_in_nonrec", "sdxb_in_decol", "sdxb_out_read", "sdxb_out_com", "sdxb_out_mort",
  # "cdxb_in_nonrec", "cdxb_out_read", "cdxb_out_decol", "cdxb_out_com", "cdxb_out_mort"
  #  )
 
# 2.10) Hospital cost calculations ####

  #bed-day costs
  bed_c <-
    cGenWardBed * (Sg + Cg + St + Ct + Sx + Cx + Sy + Cy) +
    cICUBed * (Ig + It + Ix + Iy) #infected patients assumed to be in ICU
  
  #screening costs
  screen_tests <- (
    #within hospital screening
    (gammasg * Sg + gammacg * Cg) + #from general ward
      (gammasx * Sx + gammacx * Cx) + #from tested positive group
      (gammasy * Sy + gammacy * Cy) + #from tested negative group
      
      #admission screening
      (nud * (rhos * Sda)) + #from recently discharged S, without positive flag
      (nud * (rhoc * Cda)) + #from recently discharged C, without positive flag
      (nud * (rhop * Sdb)) + #from non-recently discharged S, without positive flag
      (nud * (rhop * Cdb)) + #from non-recently discharged C, without positive flag
      
      (nudx * (rhos * Sdxa)) + #from recently discharged S, with positive flag
      (nudx * (rhoc * Cdxa)) + #from recently discharged C, with positive flag
      (nudx * (rhop * Sdxb)) + #from non-recently discharged S, with positive flag
      (nudx * (rhop * Cdxb)) + #from non-recently discharged C, with positive flag
      
      (nuec * Aec) + #from England wider community
      (nunec * Anec) #from non-England wider community
  )
  
  clinical_tests <- (
    #within hospital clinical testing
    (gammaig * Ig) + #from general ward
      (gammaix * Ix) + #from tested positive group
      (gammaiy * Iy)  #from tested negative group
  )
  
  #screening test costs
  screen_test_c <- cScreenTest * sum(screen_tests, clinical_tests)
  
  #screening test staff opportunity cost
  screen_staff_c <- cScreenStaff * sum(screen_tests, clinical_tests)
  
  #treatment abx costs
  inf_trt_c <- cInfectionTreat * (Ig + It + Ix + Iy) #infected patients assumed to be treated with same abx therapy
  
  #treatment staff costs, i.e., multidisciplinary meetings re patient treatment
  inf_trt_meet_c <- cInfectionTreatMeet * new_case_i #infected patients all assumed to have one hour of multidisciplinary team meeting time
  
  #toxicity test costs
  
  #inital toxicity test cost, on first day of treatment (i.e., when patient become infected)
  tox_initial_test_c <- cToxInit * new_case_i #infected patients all assumed to be tested for toxicity
  
  #ongoing toxicity test cost, on each day patient is infected
  tox_ongoing_test_c <- cToxOngo * (Ig + It + Ix + Iy) #infected patients all assumed to be tested for toxicity
  
  #PPE for contact precaution
  ppe_equip_c <- cPPEequip * (Sx + Cx + Ix) #tested positive assumed to be under IPC measures
  
  #IPC staff opportunity cost (from Manoukian et al. (2022), which appears to combine Otter et al. (2017) IPC staff costs into a daily unit)
  ipc_staff_c <- cIPCStaff * (Sx + Cx + Ix) #tested positive assumed to be under IPC measures
  
  #room cleaning
  room_clean_c <- cRoomClean * (Sx + Cx + Ix) #tested positive assumed to be under IPC measures
  
  #stock disposal
  stock_disp_c <- #stock disposal only assumed to occur among patients that test positive
    
    #stock disposal, discharged alive
    ((dsx * Sx) * cStockDisp) +     #discharged Sx patients
    ((dcx * Cx) * cStockDisp) +     #discharged Cx patients
    ((dix * Ix) * cStockDispICU) +  #discharged Ix patients; infected patients assumed to be in ICU
    
    #stock disposal, after death
    ((msx * Sx) * cStockDisp) +     #dead Sx patients
    ((mcx * Cx) * cStockDisp) +     #dead Cx patients
    ((mix * Ix) * cStockDispICU)    #dead Ix patients; infected patients assumed to be in ICU
  
  #infectious waste stream
  infect_waste_c <- cInfectWaste * (Sx + Cx + Ix) #tested positive assumed to be under IPC measures
  
  #mortality review panel
  mort_rev_c <- #mortality review only assumed to occur among infected patients
    ((mig * Ig) * cMortRevPanel) +     #dead Ig patients
    ((mit * It) * cMortRevPanel) +     #dead It patients
    ((mix * Ix) * cMortRevPanel) +     #dead Ix patients
    ((miy * Iy) * cMortRevPanel)       #dead Iy patients
  
  #outbreak costs
  outbreak_c <- cOutbreakpI * (Ig + It + Ix + Iy) #infected patients each assumed to contribute equally to outbreak costs, based on average outbreak cost per infected patient in Otter et al. (2017)
  
  #collate vars into single object
  hos_costs <- c(
    bed_c,
    screen_tests,
    clinical_tests,
    screen_test_c,
    screen_staff_c,
    inf_trt_c,
    inf_trt_meet_c,
    tox_initial_test_c,
    tox_ongoing_test_c,
    ppe_equip_c,
    ipc_staff_c,
    room_clean_c,
    stock_disp_c,
    infect_waste_c,
    mort_rev_c,
    outbreak_c
  )
  
  names(hos_costs) <- c(
    "bed_c",
    "screen_tests",
    "clinical_tests",
    "screen_test_c",
    "screen_staff_c",
    "inf_trt_c",
    "inf_trt_meet_c",
    "tox_initial_test_c",
    "tox_ongoing_test_c",
    "ppe_equip_c",
    "ipc_staff_c",
    "room_clean_c",
    "stock_disp_c",
    "infect_waste_c",
    "mort_rev_c",
    "outbreak_c"
  )

  
# 2.11) Extract results from the model function#

  #extract results in order of interest
  results <- list(compartments, validation, patient_movement, hos_costs)
  
  #return the results
  return(results)
}

# ========================================================== #

# 3) Method of moment beta function ####

# This function takes the mean and variance of a proportion and returns the alpha and beta parameters of the beta distribution,
# which describes the probability distribution of the proportion. This is useful for sampling the value of a proportion in
# probabilistic sensitivity analysis (PSA).

# The alpha and beta parameters specify the shape of the beta distribution, which is a continuous
# probability distribution defined on the interval [0, 1].

# The function also includes a check to ensure that the variance is less than the maximum
# possible variance for a beta distribution with the given mean, which is mean * (1 - mean).
# If the variance is greater than or equal to this value, an error message is returned.

beta_mom <- function(mean, var) {
  term <- mean * ((1 - mean) / var) - 1
  alpha <- mean * term
  beta <- (1 - mean) * term
  if (var >= mean * (1 - mean))
    stop("var must be less than mean * (1 - mean)")
  return(list(alpha = alpha, beta = beta))
} #see https://devinincerti.com/2018/02/10/psa.html#beta-distribution

# ========================================================== #

# END OF SCRIPT ####

