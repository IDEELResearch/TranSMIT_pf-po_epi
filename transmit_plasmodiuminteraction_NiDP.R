###Implementing NiDP model for pathogen interactions
### to TranSMIT screening data
### methods from Hamelin et al. PLOS Biology 2019


setwd('~/Documents/P. ovale Epi Bagamoyo')

#
# Delete all data
#
rm(list=ls())

#
# Load the utility routines for model fitting
#
source("NiDP_modelFitting.R")

#
# Load the package for reading xlsx files
#
library(XNomial)

#
# Package for matching rows of matrices
#
library(prodlim)

library(gdata)



###Analysis for Pf-Po antagonism in TranSMIT
#
# data file to use
#
thisFileName <- "transmit_pfpo_counts.csv"

#
# load in the data file
#
allData <- read.csv(thisFileName)

#
# extract the relevant data as required by the fitting function
#
n <- 2
mPathCodes <- allData[1:n]
vCounts <- allData$Count
nCounts <- length(vCounts)

#
# fit the n parameter NiDP model (i.e. without specific clearance)
#
vNiDPMod <- fit_NiDP_nParam(n,mPathCodes,vCounts,nCounts,nFits=10,minP1=-10,maxP1=10)

#
# extract the model fit and parameter information and display to screen
#
print("Fit of n parameter NiDP model (gamma_i = 0 for all i)")
vNiDPMod_Info <- disp_NiDP_ZeroGamma(n,vNiDPMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
vModPred <- numeric(4)
vReorderedCounts <- numeric(4)
for(j in 1:4)
{
  infCode <- vNiDPMod_Info[[2]][j,1:2]
  z <- row.match(infCode,allData[,1:2])
  vModPred[j] <- vNiDPMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,2]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# fit the multinomial model
#
vMultiMod <- fit_Multinomial(n,mPathCodes,vCounts,nCounts,nFits=10,minLogit=-10,maxLogit=10)

#
# extract the model fit and parameter information and display to the screen
#
print("Fit of multinomial model")
vMultiMod_Info <- disp_Multinomial(n,vMultiMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
for(j in 1:4)
{
  infCode <- vMultiMod_Info[[2]][j,1:2]
  z <- row.match(infCode,allData[,1:2])
  vModPred[j] <- vMultiMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,2]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# display the difference in AIC values; note both models have same number of parameters
#
deltaAIC <- 2*(vNiDPMod_Info[[3]]-vMultiMod_Info[[3]])
if(deltaAIC>0) # NiDP is better fit
{
  sComp <- sprintf("Comp. The NiDP model is a better fit than the multinomial model (deltaAIC=%.3f)", deltaAIC)
}else{ # Binomial is better fit
  sComp <- sprintf("Comp. The multinomial model is a better fit than the NiDP model (deltaAIC=%.3f)", deltaAIC)
}
print(sComp)



### Analysis for Pf-Poc-Pow antagonism-synergism
#
# data file to use
#
thisFileName <- "transmit_pfpocpow_counts.csv"

#
# load in the data file
#
allData <- read.csv(thisFileName)


#
# extract the relevant data as required by the fitting function
#
n <- 3
mPathCodes <- allData[1:n]
vCounts <- allData$Count
nCounts <- length(vCounts)

#
# fit the n parameter NiDP model (i.e. without specific clearance)
#
vNiDPMod <- fit_NiDP_nParam(n,mPathCodes,vCounts,nCounts,nFits=10,minP1=-10,maxP1=10)

#
# extract the model fit and parameter information and display to screen
#
print("Fit of n parameter NiDP model (gamma_i = 0 for all i)")
vNiDPMod_Info <- disp_NiDP_ZeroGamma(n,vNiDPMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
vModPred <- numeric(8)
vReorderedCounts <- numeric(8)
for(j in 1:8)
{
  infCode <- vNiDPMod_Info[[2]][j,1:3]
  z <- row.match(infCode,allData[,1:3])
  vModPred[j] <- vNiDPMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,4]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# fit the multinomial model
#
vMultiMod <- fit_Multinomial(n,mPathCodes,vCounts,nCounts,nFits=10,minLogit=-10,maxLogit=10)

#
# extract the model fit and parameter information and display to the screen
#
print("Fit of multinomial model")
vMultiMod_Info <- disp_Multinomial(n,vMultiMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
for(j in 1:8)
{
  infCode <- vMultiMod_Info[[2]][j,1:3]
  z <- row.match(infCode,allData[,1:3])
  vModPred[j] <- vMultiMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,4]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# display the difference in AIC values; note both models have same number of parameters
#
deltaAIC <- 2*(vNiDPMod_Info[[3]]-vMultiMod_Info[[3]])
if(deltaAIC>0) # NiDP is better fit
{
  sComp <- sprintf("Comp. The NiDP model is a better fit than the multinomial model (deltaAIC=%.3f)", deltaAIC)
}else{ # Binomial is better fit
  sComp <- sprintf("Comp. The multinomial model is a better fit than the NiDP model (deltaAIC=%.3f)", deltaAIC)
}
print(sComp)


### Analysis for Pf-Poc-Pow antagonism-synergism using conservative mixed-infection scheme
#
# data file to use
#
thisFileName <- "transmit_pfpocpowconservative_counts.csv"

#
# load in the data file
#
allData <- read.csv(thisFileName)


#
# extract the relevant data as required by the fitting function
#
n <- 3
mPathCodes <- allData[1:n]
vCounts <- allData$Count
nCounts <- length(vCounts)

#
# fit the n parameter NiDP model (i.e. without specific clearance)
#
vNiDPMod <- fit_NiDP_nParam(n,mPathCodes,vCounts,nCounts,nFits=10,minP1=-10,maxP1=10)

#
# extract the model fit and parameter information and display to screen
#
print("Fit of n parameter NiDP model (gamma_i = 0 for all i)")
vNiDPMod_Info <- disp_NiDP_ZeroGamma(n,vNiDPMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
vModPred <- numeric(8)
vReorderedCounts <- numeric(8)
for(j in 1:8)
{
  infCode <- vNiDPMod_Info[[2]][j,1:3]
  z <- row.match(infCode,allData[,1:3])
  vModPred[j] <- vNiDPMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,4]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# fit the multinomial model
#
vMultiMod <- fit_Multinomial(n,mPathCodes,vCounts,nCounts,nFits=10,minLogit=-10,maxLogit=10)

#
# extract the model fit and parameter information and display to the screen
#
print("Fit of multinomial model")
vMultiMod_Info <- disp_Multinomial(n,vMultiMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
for(j in 1:8)
{
  infCode <- vMultiMod_Info[[2]][j,1:3]
  z <- row.match(infCode,allData[,1:3])
  vModPred[j] <- vMultiMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,4]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# display the difference in AIC values; note both models have same number of parameters
#
deltaAIC <- 2*(vNiDPMod_Info[[3]]-vMultiMod_Info[[3]])
if(deltaAIC>0) # NiDP is better fit
{
  sComp <- sprintf("Comp. The NiDP model is a better fit than the multinomial model (deltaAIC=%.3f)", deltaAIC)
}else{ # Binomial is better fit
  sComp <- sprintf("Comp. The multinomial model is a better fit than the NiDP model (deltaAIC=%.3f)", deltaAIC)
}
print(sComp)


### Analysis for Pf-Poc-Pow antagonism-synergism using conservative mixed-infection scheme and randomized unidentified Po infection assignment to mono-infections
#
# data file to use
#
thisFileName <- "transmit_pfpocpow_cons_randomized_counts.csv"

#
# load in the data file
#
allData <- read.csv(thisFileName)


#
# extract the relevant data as required by the fitting function
#
n <- 3
mPathCodes <- allData[1:n]
vCounts <- allData$Count
nCounts <- length(vCounts)

#
# fit the n parameter NiDP model (i.e. without specific clearance)
#
vNiDPMod <- fit_NiDP_nParam(n,mPathCodes,vCounts,nCounts,nFits=10,minP1=-10,maxP1=10)

#
# extract the model fit and parameter information and display to screen
#
print("Fit of n parameter NiDP model (gamma_i = 0 for all i)")
vNiDPMod_Info <- disp_NiDP_ZeroGamma(n,vNiDPMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
vModPred <- numeric(8)
vReorderedCounts <- numeric(8)
for(j in 1:8)
{
  infCode <- vNiDPMod_Info[[2]][j,1:3]
  z <- row.match(infCode,allData[,1:3])
  vModPred[j] <- vNiDPMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,4]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# fit the multinomial model
#
vMultiMod <- fit_Multinomial(n,mPathCodes,vCounts,nCounts,nFits=10,minLogit=-10,maxLogit=10)

#
# extract the model fit and parameter information and display to the screen
#
print("Fit of multinomial model")
vMultiMod_Info <- disp_Multinomial(n,vMultiMod)

#
# do the Monte Carlo goodness of fit test (taking care to ensure rows are ordered correctly)
#
for(j in 1:8)
{
  infCode <- vMultiMod_Info[[2]][j,1:3]
  z <- row.match(infCode,allData[,1:3])
  vModPred[j] <- vMultiMod_Info[[1]][j]
  vReorderedCounts[j] <- allData[z,4]
}
v <- xmonte(vReorderedCounts,vModPred, histobins=F,statName="LLR",ntrials=1000000,detail=0)
sMC <- sprintf("   GoF.  Monte-Carlo p=%.3f", v$pLLR)
print(sMC)

#
# display the difference in AIC values; note both models have same number of parameters
#
deltaAIC <- 2*(vNiDPMod_Info[[3]]-vMultiMod_Info[[3]])
if(deltaAIC>0) # NiDP is better fit
{
  sComp <- sprintf("Comp. The NiDP model is a better fit than the multinomial model (deltaAIC=%.3f)", deltaAIC)
}else{ # Binomial is better fit
  sComp <- sprintf("Comp. The multinomial model is a better fit than the NiDP model (deltaAIC=%.3f)", deltaAIC)
}
print(sComp)



