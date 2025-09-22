################################################################################
# RUNS BAYESIAN QBA FOR THE APPLIED EXAMPLE

# A. Preliminaries
# loads necessary packages and functions

# B. Creates objects necessary to analyse the applied example
# lists model parameters
# lists priors
# sets parameters necessary to run JAGS

# C. Runs Bayesian QBA for the applied example
# creates trace plots
# displays the summary results
################################################################################

################################################################################
# A. Preliminaries
################################################################################

# REMOVE ALL OBJECTS (IF ANY)
rm(list=ls())
# CLOSE OPEN GRAPHICS DEVICES (IF ANY)
graphics.off()

# INSTALL unmconf IF NECESSARY
if (!("unmconf" %in% rownames(installed.packages())))
{
  install.packages("unmconf")
}
# LOAD unmconf
library(unmconf)
# CHECK THE LOADED VERSION IS 1.0.0
if (sessionInfo()$otherPkgs$unmconf$Version!="1.0.0")
  stop("wrong version of unmconf")
# LOAD EDITED unm_glm FUNCTION CALLED MyUnm_glm
source("../../Bayesian function/MyUnm_glm.R")

################################################################################
# B. Create objects necessary to analyse the applied example
################################################################################

# LIST MODEL PARAMETERS
# model parameter names using JAGS notation (obtained from dimnames(as.array(Posteriors))[2]$var)
Params <- c("beta[1]", "beta[alcoholXnever]", "beta[alcoholXformer]", "beta[alcoholX3to4]", "beta[alcoholXatleast5]", "beta[ageC40to59]", "beta[ageC60to79]", "beta[raceCBlack]", "beta[raceCMexicanAmerican]", "beta[raceCHispanic]", "beta[raceCOther]", "beta[educationCHighSchoolGraduate]", "beta[educationCCollege]", "beta[educationCCollegeGraduate]", "beta[PALCmoderate]", "beta[PALChigh]", "beta[yearC2005to2006]", "beta[yearC2007to2008]", "beta[yearC2009to2010]", "beta[yearC2011to2012]", "beta[recall1Cweekday]", "beta[recall2Cweekday]", "beta[smokingCnever]", "beta[smokingCformer]", "beta[incomeC131to299]", "beta[incomeCatleast300]", "beta[misreportingCaccurate]", "beta[misreportingCunder]", "delta[1]", "delta[alcoholXnever]", "delta[alcoholXformer]", "delta[alcoholX3to4]", "delta[alcoholXatleast5]", "eta[energyU]", "gamma[alcoholXnever]", "gamma[alcoholXformer]", "gamma[alcoholX3to4]", "gamma[alcoholXatleast5]", "lambda[energyU]", "lambda[maritalstatusU]", "sigma[waistcircumferenceY]")
# model parameter names using our manuscript's notation
ParamNames <- c("beta0", "betaalcoholXnever", "betaalcoholXformer", "betaalcoholtX3to4", "betaalcoholXatleast5", "betaageC40to59", "betaageC60to79", "betaraceCBlack", "betaraceCMexicanAmerican", "betaraceCHispanic", "betaraceCOther", "betaeducationCHighSchoolGraduate", "betaeducationCCollege", "betaeducationCCollegeGraduate", "betaPALCmoderate", "betaPALChigh", "betayearC2005to2006", "betayearC2007to2008", "betayearC2009to2010", "betayearC2011to2012", "betarecall1Cweekday", "betarecall2Cweekday", "betasmokingCnever", "betasmokingCformer", "betaincomeC131to299", "betaincomeCatleast300", "betamisreportingCaccurate", "betamisreportingCunder", "alpha20", "alpha2alcoholXnever", "alpha2alcoholXformer", "alpha2alcoholX3to4", "alpha2alcoholXatleast5", "etaenergyU", "alpha1alcoholXnever", "alpha1alcoholXformer", "alpha1alcoholX3to4", "alpha1alcoholXatleast5", "betaenergyU", "betamaritalstatusU", "epsilonY")

# LIST PRIORS
# coefficients in Y|X,C,energyU,maritalstatusU
betaenergyU_mean <- 0.006931
betaenergyU_sd <- 0.0002805
betamaritalstatusU_mean <- -1.312
betamaritalstatusU_sd <- 0.3453
# coefficients in U1|X
alpha1alcoholXnever_mean <- -43.47
alpha1alcoholXnever_sd <- 27.21
alpha1alcoholXformer_mean <- 69.72
alpha1alcoholXformer_sd <- 18.63
alpha1alcoholX3to4_mean <- -65.37
alpha1alcoholX3to4_sd <- 18.07
alpha1alcoholXatleast5_mean <- -94.71
alpha1alcoholXatleast5_sd <- 19.64
# residual standard deviation in U1|X
etaenergyU_min <- 372.6
etaenergyU_max <- 716.5
# coefficients in U2|X
alpha20_mean <- -0.9175
alpha20_sd <- 0.04286
alpha2alcoholXnever_mean <- 0.4524
alpha2alcoholXnever_sd <- 0.1081
alpha2alcoholXformer_mean <- -0.07963
alpha2alcoholXformer_sd <- 0.08107
alpha2alcoholX3to4_mean <- 0.3831
alpha2alcoholX3to4_sd <- 0.07405
alpha2alcoholXatleast5_mean <- 0.3785
alpha2alcoholXatleast5_sd <- 0.07958

# SET PARAMETERS NECESSARY TO RUN JAGS
# number of burn-ins for analysis
Burnin <- 2000
# number of total iterations for analysis
Iter <- 20000
# number of MCMC chains
Chains <- 1

################################################################################
# A. Run Bayesian QBA for the applied example
################################################################################

# READ DATA
Data <- read.csv(file="../Data/Data.csv")
# set alcoholX as a factor as "1to2" as baseline
Data$alcoholX <- factor(Data$alcoholX, levels=c("1to2", "never", "former", "3to4", "atleast5"))
Data$alcoholXnever <- as.numeric(Data$alcoholX=="never")
Data$alcoholXformer <- as.numeric(Data$alcoholX=="former")
Data$alcoholX3to4 <- as.numeric(Data$alcoholX=="3to4")
Data$alcoholXatleast5 <- as.numeric(Data$alcoholX=="atleast5")
Data$ageC <- factor(Data$ageC)
Data$raceC <- factor(Data$raceC, levels=c("White", "Black", "MexicanAmerican", "Hispanic", "Other"))
Data$educationC <- factor(Data$educationC, levels=c("<HighSchool", "HighSchoolGraduate", "College", "CollegeGraduate"))
Data$PALC <- Data$physicalactivitylevelC
Data$PALC <- factor(Data$PALC, levels=c("low", "moderate", "high"))
Data$yearC <- factor(Data$yearC)
Data$recall1C <- factor(Data$recall1C, levels=c("weekend", "weekday"))
Data$recall2C <- factor(Data$recall2C, levels=c("weekend", "weekday"))
Data$smokingC <- factor(Data$smokingC, levels=c("current", "never", "former"))
Data$incomeC <- factor(Data$incomeC)
Data$misreportingC <- factor(Data$misreportingC, levels=c("over", "accurate", "under"))

# CREATE A DATA FRAME FOR RUNNING JAGS
JAGSData <- Data
# SET energyU_tilde (part of energyU uncorrelated with C) to NA
JAGSData$energyU[1:nrow(JAGSData)] <- NA
# SET maritalstatusU_tilde (part of maritalstatusU uncorrelated with C) to NA
JAGSData$maritalstatusU[1:nrow(JAGSData)] <- NA

# SET SEED FOR GENERATING INITIAL VALUES
set.seed(20250721)
# GENERATE A SET OF INITIAL VALUES AND SET .RNG.name AND .RNG.seed
Inits5 <- list("beta"=rnorm(28, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(4, mean=0, sd=1), "delta"=rnorm(5, mean=0, sd=1), "tau_waistcircumferenceY"=rgamma(1, shape=0.1, rate=1), "eta_energyU"=runif(1, min=etaenergyU_min, max=etaenergyU_max), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=5)
# RUN MyUnm_glm FOR 2*Iter iterations
Posteriors <- MyUnm_glm(form1=waistcircumferenceY~alcoholX+ageC+raceC+educationC+PALC+yearC+recall1C+recall2C+smokingC+incomeC+misreportingC+energyU+maritalstatusU, form2=energyU~-1+alcoholXnever+alcoholXformer+alcoholX3to4+alcoholXatleast5, form3=maritalstatusU~alcoholXnever+alcoholXformer+alcoholX3to4+alcoholXatleast5, family1=gaussian(), family2=gaussian(), family3=binomial(), data=JAGSData, n.iter=2*Iter+Burnin, n.adapt=0, n.chains=Chains, priors=c("beta[1]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[alcoholXnever]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[alcoholXformer]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[alcoholX3to4]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[alcoholXatleast5]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[ageC40to59]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[ageC60to79]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[raceCBlack]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[raceCMexicanAmerican]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[raceCHispanic]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[raceCOther]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[educationCHighSchoolGraduate]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[educationCCollege]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[educationCCollegeGraduate]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[PALCmoderate]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[PALChigh]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[yearC2005to2006]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[yearC2007to2008]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[yearC2009to2010]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[yearC2011to2012]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[recall1Cweekday]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[recall2Cweekday]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[smokingCnever]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[smokingCformer]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[incomeC131to299]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[incomeCatleast300]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[misreportingCaccurate]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "beta[misreportingCunder]"=paste0("dnorm(", 0, ",", 1/10^2, ")"), "lambda[energyU]"=paste0("dnorm(", betaenergyU_mean, ",", 1/betaenergyU_sd^2, ")"), "lambda[maritalstatusU]"=paste0("dnorm(", betamaritalstatusU_mean, ",", 1/betamaritalstatusU_sd^2, ")"), "gamma[alcoholXnever]"=paste0("dnorm(", alpha1alcoholXnever_mean, ",", 1/alpha1alcoholXnever_sd^2, ")"), "gamma[alcoholXformer]"=paste0("dnorm(", alpha1alcoholXformer_mean, ",", 1/alpha1alcoholXformer_sd^2, ")"), "gamma[alcoholX3to4]"=paste0("dnorm(", alpha1alcoholX3to4_mean, ",", 1/alpha1alcoholX3to4_sd^2, ")"), "gamma[alcoholXatleast5]"=paste0("dnorm(", alpha1alcoholXatleast5_mean, ",", 1/alpha1alcoholXatleast5_sd^2, ")"), "delta[1]"=paste0("dnorm(", alpha20_mean, ",", 1/alpha20_sd^2, ")"), "delta[alcoholXnever]"=paste0("dnorm(", alpha2alcoholXnever_mean, ",", 1/alpha2alcoholXnever_sd^2, ")"), "delta[alcoholXformer]"=paste0("dnorm(", alpha2alcoholXformer_mean, ",", 1/alpha2alcoholXformer_sd^2, ")"), "delta[alcoholX3to4]"=paste0("dnorm(", alpha2alcoholX3to4_mean, ",", 1/alpha2alcoholX3to4_sd^2, ")"), "delta[alcoholXatleast5]"=paste0("dnorm(", alpha2alcoholXatleast5_mean, ",", 1/alpha2alcoholXatleast5_sd^2, ")")), response_nuisance_priors=paste0("tau_waistcircumferenceY~dgamma(", 0.001, ",", 0.001, "); sigma_waistcircumferenceY <- tau_waistcircumferenceY^-0.5"), response_params_to_track="sigma_waistcircumferenceY", confounder1_nuisance_priors=paste0("tau_energyU <- eta_energyU^-2; eta_energyU~dunif(", etaenergyU_min, ",", etaenergyU_max, ")"), confounder1_params_to_track="eta_energyU", inits=list(Inits5))

for (Rep in 1:3)
{
  # CREATE A TRACE PLOT FOR EACH MODEL PARAMETER
  par(mfrow=c(ceiling(length(Params)/3/2),2), mar=c(4,4,0.5,0.5))
  for (param in Params[((Rep-1)*length(Params)/3+1):(Rep*length(Params)/3)])
  {
    plot(as.array(Posteriors)[,param], type="l", xlab="Iterations including burn-ins", ylab=ParamNames[which(Params==param)], ylim=extendrange(c(as.array(Posteriors)[,param])))
  }
}

# DISPLAY SUMMARY RESULTS INCLUDING THE MEDIAN, 2.5TH AND 97.5TH PERCENTILES OF EXPOSURE EFFECT ESTIMATES
print(data.frame("median"=apply(as.array(Posteriors)[(Burnin+1):(Burnin+Iter),], 2, function(x) quantile(x, probs=0.5)), "sd"=apply(as.array(Posteriors)[(Burnin+1):(Burnin+Iter),], 2, sd), "2.5%"=apply(as.array(Posteriors)[(Burnin+1):(Burnin+Iter),], 2, function(x) quantile(x, probs=0.025)), "97.5%"=apply(as.array(Posteriors)[(Burnin+1):(Burnin+Iter),], 2, function(x) quantile(x, probs=0.975)), check.names=F))
