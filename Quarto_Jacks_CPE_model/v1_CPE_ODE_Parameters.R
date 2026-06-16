# CPE Compartmental Model - parameters for base case and intervention scenarios ####

# 1) Introduction ####

# This script defines the parameters for the CPE compartmental model, which is defined in 
# the 'CPE_ODE_Functions.R' script.

# The parameters are called by the CPE_ODE_Model_Run.qmd, which runs the model and 
# generates outputs for validation and cost-effectiveness analysis.

# ========================================================== #

# 2) Initial setup ####

# Load libraries #

library(dplyr)

# Avoid R using scientific notation for values #

options(scipen = 1000)
# scipen short for scientific penalty, where a high value essentially tells R not
# to use scientific notation unless it's absolutely necessary

# ========================================================== #
# 3) Create function to convert probabilities to rates ####

#prob2rate: function that converts probability to a rate; assume t = days

prob2rate <- function(p, t) {
  r <- -(1 / t) * (log(1 - p))
  return(r)
}

# ========================================================== #

# 4) Set the fixed model parameter values ####

# These do not change between the base case and intervention scenarios, but are used in both.

# Transition within hospital #

beta0 <- 0.00045954 #transmission rate from other sources, i.e., background
beta1 <- 0.00037131 #transmission rate from C and I patients
Zg <- 0 #effectiveness of environmental cleaning in general wards
Zt <- 0 #effectiveness of environmental cleaning for patients awaiting results
Zx <- 0.75 #effectiveness of environmental cleaning for patients tested positive
Zy <- 0 #effectiveness of environmental cleaning for patients tested negative
chi <- 0 #effectiveness of IPC measures for patients awaiting test results
epsilon <- 0.85 #effectiveness of IPC measures for patients tested positive
alpha <- 0.00978552 #infection progression rate
gammasg <- 0 #testing rate of susceptible patients
gammacg <- 0 #testing rate of colonised patients
gammaig <- 1 / 1 #testing rate of infected patients
gammasx <- 0 #testing rate of susceptible tested positive patients
gammacx <- 0 #testing rate of colonised tested positive patients
gammaix <- 0 #testing rate of infected tested positive patients
gammasy <- 0 #testing rate of susceptible tested negative patients
gammacy <- 0 #testing rate of colonised tested negative patients
gammaiy <- 1 / 1 #testing rate of infected tested negative patients

# Length of stay #

lossg <- 7.34 #susceptible length of stay
loscg <- lossg + 12.21470588 #colonised length of stay
losig <- lossg + 12.21470588 #infected length of stay
losst <- lossg #susceptible testing length of stay
losct <- loscg #colonised testing length of stay
losit <- losig #infected testing length of stay
lossx <- lossg #susceptible tested positive length of stay
loscx <- loscg #colonised tested positive length of stay
losix <- losig #infected tested positive length of stay
lossy <- lossg #susceptible tested negative length of stay
loscy <- loscg #colonised tested negative length of stay
losiy <- losig #infected tested negative length of stay

# Discharges from hospital #

dsg <- 1 / lossg #discharge rate of susceptible patients
dcg <- 1 / loscg #discharge rate of colonised patients
dig <- 1 / losig #discharge rate of infected patients
dst <- 1 / losst #discharge rate of susceptible testing patients
dct <- 1 / losct #discharge rate of colonised testing patients
dit <- 1 / losit #discharge rate of infected testing patients
dsx <- 1 / lossx #discharge rate of susceptible tested positive patients
dcx <- 1 / loscx #discharge rate of colonised tested positive patients
dix <- 1 / losix #discharge rate of infected tested positive patients
dsy <- 1 / lossy #discharge rate of susceptible tested negative patients
dcy <- 1 / loscy #discharge rate of colonised tested negative patients
diy <- 1 / losiy #discharge rate of infected tested negative patients

# Transitions within community #

dp <- 365 #discharge period (days)
rdp <- 30 #recently discharged period (days)
nrdp <- dp - rdp #non-recent discharge period
psi <- 1 / rdp #1/recently discharged period
delta <- 1 / nrdp #1/non-recent discharged period
cropf <- 365 #CRO positive flag period on patient records (days)
omega <- 1 / (cropf - rdp) #1/CRO positive flag period on patient records - recent discharge period
decol <- 220.0517471 #average time to decolonisation (days)
mu <- 1 / decol #decolonisation rate in the community

# Hospital readmission from recently discharged #

rhos <- prob2rate(0.15065221, 30) #readmission rate among recently discharged susceptible patients
rhoc <- prob2rate(0.297, 30) #readmission rate among recently discharged colonised patients
rhop <- prob2rate(0.097462, (365 - 30)) #readmission rate among non-recently discharged patients, i.e. general population
nudx <- 1 #proportion discharged with positive flag screened on readmission
nud <- 1 #proportion discharged without positive flag screened on readmission

# Hospital admission from wider community #

sigmaec <- 0.000532 #prevalence of colonisation upon admission from UK wider community
sigmanec <- 0.01393728 #prevalence of colonisation upon admission from foreign wider community
tau <- 0.0031088 #probability of colonised patient presenting with infection
zeta <- 0.02296281 #proportion of all admissions from abroad
nuec <- 1 #proportion of UK wider community patients screened on admission
nunec <- 1 #proportion of foreign wider community patients screened on admission

# Mortality #

msg <- prob2rate(0.02147048, lossg) #mortality rate susceptible patients
mcg <- msg #mortality rate colonised patients
mig <- prob2rate(0.27702703, losig) #mortality rate infected patients
mst <- msg #mortality rate of susceptible testing patients
mct <- mcg #mortality rate of colonised testing patients
mit <- mig #mortality rate of infected testing patients
msx <- msg #mortality rate of susceptible tested positive patients
mcx <- mcg #mortality rate of colonised tested positive patients
mix <- mig #mortality rate of infected tested positive patients
msy <- msg #mortality rate of susceptible tested negative patients
mcy <- mcg #mortality rate of colonised tested negative patients
miy <- mig #mortality rate of infected tested negative patients
msda <- prob2rate(0.00977582, 30) #mortality rate recently discharged susceptible
mcda <- msda #mortality rate recently discharged unknown colonised
msdb <- prob2rate(0.21525354, 365) #mortality rate non-recently discharged susceptible
mcdb <- msdb #mortality rate non-recently discharged unknown colonised
msdxa <- msda #mortality rate recently discharged susceptible positive flag
mcdxa <- msda #mortality rate recently discharged colonised positive flag
msdxb <- msdb #mortality rate non-recently discharged susceptible positive flag
mcdxb <- msdb #mortality rate non-recently discharged colonised positive flag

# Unit costs #

#bed-day unit costs
cGenWardBed <- 594.90 #general ward bed day cost, per patient per day
cICUBed <- 2202.03 #ICU bed day cost, per patient per day

#screening test unit costs
cCultureTest <- 14.63 #one-off cost of a culture test, per test
cPCRTest <- 56.37 #one-off cost of a PCR test, per test
cScreenTest <- cCultureTest #culture test treated as standard approach

#screening test staff unit opportunity costs
cScreenStaff <- 0 #6.85 #unit cost of staff time opportunity cost, per test
#NOTE: cScreenStaff is set to 0 because screening test unit costs from Manoukian et al. (2022) incorporate all costs, including staff time

#treatment (i.e., abx) unit costs
cInfectionTreat <- 275.37 #infection treatment cost, per patient per day
cInfectionTreatMeet <- 406.33 #infection treatment staff costs (i.e., multidisciplinary meetings re patient treatment), per new infection

#toxicity test unit costs
cToxInit <- 147.99 #initial toxicity test cost, per new infection
cToxOngo <- 21.15 #ongoing toxicity test cost, per patient per day

#IPC costs
cPPEequip <- 26.91 #PPE (gloves and aprons) for contact precautions, per patient per day
cIPCStaff <- 51.56 #staff opportunity cost of IPC measures, per patient per day
cRoomClean <- 25.84 #cost of room cleaning for patients under IPC measures, per patient per day
cStockDisp <- 155.84 #cost of stock disposal after regular discharge, per patient discharge (dead or alive)
cStockDispICU <- 514.99 #cost of stock disposal after ICU discharge, per patient discharge (dead or alive)
cInfectWaste <- 0.46 #cost of infectious waste, per isolating patient per day

#mortality review panel
cMortRevPanel <- 335.96 #cost of mortality review panel, per positive patient death

#outbreak costs
cOutbreakpI <- 0 #35672 #outbreak cost, per infected patient (= 0 OR total 'outbreak' costs in Otter et al. (2017) divided by 18 patients that had CPE+ infection)

# ========================================================== #

#5) Create a function to generate the model parameters for each scenario ####

make_mod_parms <- function(
    ttt,
    theta,
    phi,
    cScreenTest
) {
  stopifnot(length(ttt) == 1, length(theta) == 1, length(phi) == 1)
  stopifnot(ttt > 0, theta >= 0, theta <= 1, phi >= 0, phi <= 1)
  n <- 1 / ttt
  
  list(
    beta0 = beta0, beta1 = beta1,
    Zg = Zg, Zt = Zt, Zx = Zx, Zy = Zy,
    chi = chi, epsilon = epsilon,
    alpha = alpha,
    gammasg = gammasg, gammacg = gammacg, gammaig = gammaig,
    gammasx = gammasx, gammacx = gammacx, gammaix = gammaix,
    gammasy = gammasy, gammacy = gammacy, gammaiy = gammaiy,
    
    dsg = dsg, dcg = dcg, dig = dig,
    dst = dst, dct = dct, dit = dit,
    dsx = dsx, dcx = dcx, dix = dix,
    dsy = dsy, dcy = dcy, diy = diy,
    
    psi = psi, delta = delta, omega = omega, mu = mu,
    
    rhos = rhos, rhoc = rhoc, rhop = rhop,
    nudx = nudx, nud = nud,
    sigmaec = sigmaec, sigmanec = sigmanec,
    tau = tau, zeta = zeta, nuec = nuec, nunec = nunec,
    
    msg = msg, mcg = mcg, mig = mig,
    mst = mst, mct = mct, mit = mit,
    msx = msx, mcx = mcx, mix = mix,
    msy = msy, mcy = mcy, miy = miy,
    msda = msda, mcda = mcda, msdb = msdb, mcdb = mcdb,
    msdxa = msdxa, mcdxa = mcdxa, msdxb = msdxb, mcdxb = mcdxb,
    
    ttt = ttt, n = n, theta = theta, phi = phi,
    
    cGenWardBed = cGenWardBed, cICUBed = cICUBed,
    cCultureTest = cCultureTest, cPCRTest = cPCRTest,
    cScreenTest = cScreenTest, cScreenStaff = cScreenStaff,
    cInfectionTreat = cInfectionTreat, cInfectionTreatMeet = cInfectionTreatMeet,
    cToxInit = cToxInit, cToxOngo = cToxOngo,
    cPPEequip = cPPEequip, cIPCStaff = cIPCStaff, cRoomClean = cRoomClean,
    cStockDisp = cStockDisp, cStockDispICU = cStockDispICU,
    cInfectWaste = cInfectWaste, cMortRevPanel = cMortRevPanel,
    cOutbreakpI = cOutbreakpI
  )
}

# ========================================================== #

# 6) Generate the model parameters for the base case and intervention scenarios ####

mod.parms.i0 <- make_mod_parms(ttt = 2, theta = 0.834, phi = 0.933,  cScreenTest = cCultureTest)
mod.parms.i1 <- make_mod_parms(ttt = 2, theta = 0.960, phi = 0.966,  cScreenTest = cPCRTest)
mod.parms.i2 <- make_mod_parms(ttt = 0.5, theta = 0.960, phi = 0.966, cScreenTest = cCultureTest)

# ========================================================== #

# END OF SCRIPT ####
