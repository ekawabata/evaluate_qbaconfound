################################################################################
# RUNS SIMULATION STUDY I FOR BAYESIAN QBA
# SCENARIO B: Y CONTINUOUS, X CONTINUOUS, MULTIPLE U=(U1, U2) CONTINUOUS

# A. Preliminaries
# loads necessary packages and functions

# B. Creates objects necessary to perform simulation study I Scenario B
# lists model parameters
# lists true values
# lists priors
# sets parameters necessary to run JAGS

# C. Chooses the numbers of burn-in and monitoring iterations
# generates a chain of 50,000 iterations on Dataset_1.csv
# creates a trace plot and write it to a png file
# generates three chains of 4,000 iterations on 10 datasets from Dataset_1.csv to Dataset_10.csv
# writes chain-specific effective sample sizes to csv files
# writes Rhat to a csv file
# creates a trace plot and write it to a png file

# D. Runs Bayesian QBA on all 500 datasets
# records the summary results in a csv file
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
source("../../../../../Bayesian function/MyUnm_glm.R")

################################################################################
# B. Create objects necessary to perform simulation study I Scenario B
################################################################################

# LIST MODEL PARAMETERS
# model parameter names using JAGS notation (obtained from dimnames(as.array(Posteriors))[2]$var)
Params <- c("beta[1]", "beta[conX]", "beta[binC]", "beta[conC]", "beta[nomC1]", "beta[nomC2]", "lambda[conU1[tilde]]", "lambda[conU2[tilde]]", "gamma[conX]", "delta[conX]", "sigma[conY]", "sigma_conU1[tilde]", "sigma_conU2[tilde]")
# model parameter names using our manuscript's notation
ParamNames <- c("beta0", "betaconX", "betabinC", "betaconC", "betanomC1", "betanomC2", "betaconU1", "betaconU2", "alpha1conX", "alpha2conX", "epsilonY", "eta1", "eta2")

# LIST TRUE VALUES
# coefficients in Y|X,C,U1,U2
beta0 <- 8.534
betaconX <- 0
betabinC <- -0.2531
betaconC <- -0.1783
betanomC1 <- 0.5421
betanomC2 <- 0.5630
betaconU1 <- 0.1299
betaconU2 <- 0.06747
# coefficient in U1|X,C
alpha1conX <- 1.260
# coefficient in U2|X
alpha2conX <- 1.281
# residual standard deviation in Y|X,C,U1,U2
epsilonY <- 4.322
# residual standard deviation in U1|X
eta1 <- 9.394
# residual standard deviation in U2|X
eta2 <- 10.30
# combine to create a data frame
TrueValues <- data.frame('par'=Params, 'true'=c(beta0, betaconX, betabinC, betaconC, betanomC1, betanomC2, betaconU1, betaconU2, alpha1conX, alpha2conX, epsilonY, eta1, eta2))

# LIST VERY INFORMATIVE UNIFORM PRIORS
betaconU1_sd1 <- 0.006628
betaconU2_sd1 <- 0.003442
alpha1conX_sd1 <- 0.06429
eta1_min1 <- 8.405
eta1_max1 <- 10.38
alpha2conX_sd1 <- 0.06536
eta2_min1 <- 9.216
eta2_max1 <- 11.38

# LIST INFORMATIVE UNIFORM PRIORS
betaconU1_sd2 <- 0.01988
betaconU2_sd2 <- 0.01033
alpha1conX_sd2 <- 0.1929
eta1_min2 <- 6.427
eta1_max2 <- 12.36
alpha2conX_sd2 <- 0.1961
eta2_min2 <- 7.047
eta2_max2 <- 13.55

# LIST INFORMATIVE GAMMA PRIORS
inverseeta1squared_shape4 <- 13.19
inverseeta1squared_rate4 <- 1164
inverseeta2squared_shape4 <- 13.19
inverseeta2squared_rate4 <- 1399

# LIST NAMES OF THREE PRIOR SETTINGS
# VeryInformative;Uniform (95% of the sampled values were within plus/minus 10% of the truth; uniform prior on the two residual standard deviations that are bias parameters)
# Informative;Uniform (95% of the sampled values were within plus/minus 30% of the truth; uniform prior on the two residual standard deviations that are bias parameters)
# third prior setting is no longer investigated
# Informative;Gamma (95% of the sampled values were within plus/minus 30% of the truth; gamma prior on the two residual precisions that are bias parameters)
Informativeness <- c("VeryInformative;Uniform", "Informative;Uniform", "", "Informative;Gamma")

# SET PARAMETERS NECESSARY TO RUN JAGS
# number of burn-ins for the one long chain
LongBurnin <- 0
# total iterations for the one long chain
LongIter <- 50000
# number of burn-ins for analysis
Burnin <- 1000
# number of total iterations for analysis
Iter <- 4000
# number of MCMC chains
Chains <- 3
# read a list of random integers for seeds
Seeds <- read.csv(file="../Seeds/RandomIntegers.csv")$seed

################################################################################
# C. Choose the numbers of burn-in and monitoring iterations
################################################################################

for (informativeness in c(2, 4))
{
  # FOR Dataset_1.csv,
  for (dataset in 1)
  {
    # READ SIMULATED DATA
    SimulatedData <- read.csv(file=paste0("../../Data/Dataset_", dataset, ".csv"))
    SimulatedData$nomC <- SimulatedData$cnom
    SimulatedData$binC <- SimulatedData$cbin
    SimulatedData$conC <- SimulatedData$ccon
    SimulatedData$conX <- SimulatedData$xcon
    SimulatedData$conU1 <- SimulatedData$ucon1
    SimulatedData$conU2 <- SimulatedData$ucon2
    SimulatedData$conY <- SimulatedData$ycon
    
    # MEAN-CENTRE COVARIATES
    SimulatedData$conX_cen <- SimulatedData$conX-mean(SimulatedData$conX)
    SimulatedData$binC_cen <- SimulatedData$binC-mean(SimulatedData$binC)
    SimulatedData$conC_cen <- SimulatedData$conC-mean(SimulatedData$conC)
    SimulatedData$nomC <- factor(SimulatedData$nomC, levels=0:2)
    SimulatedData$nomC1 <- as.numeric(SimulatedData$nomC==1)
    SimulatedData$nomC1_cen <- SimulatedData$nomC1-mean(SimulatedData$nomC1)
    SimulatedData$nomC2 <- as.numeric(SimulatedData$nomC==2)
    SimulatedData$nomC2_cen <- SimulatedData$nomC2-mean(SimulatedData$nomC2)
    
    # CREATE A DATA FRAME FOR RUNNING JAGS
    JAGSData <- data.frame(conY=SimulatedData[,"conY"], conX=SimulatedData[,"conX_cen"], binC=SimulatedData[,"binC_cen"], conC=SimulatedData[,"conC_cen"], nomC1=SimulatedData[,"nomC1_cen"], nomC2=SimulatedData[,"nomC2_cen"], conU1_tilde=1:nrow(SimulatedData), conU2_tilde=1:nrow(SimulatedData))
    # SET conU1_tilde (part of conU1 uncorrelated with C) to NA
    JAGSData$conU1_tilde[1:nrow(JAGSData)] <- NA
    # SET conU2_tilde (part of conU2 uncorrelated with C) to NA
    JAGSData$conU2_tilde[1:nrow(JAGSData)] <- NA
    
    # SET SEED FOR GENERATING INITIAL VALUES
    set.seed(2)
    
    if (informativeness==2)
    {
      # GENERATE A SET OF INITIAL VALUES AND SET .RNG.name AND .RNG.seed
      Inits1 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "sigma_conU1_tilde"=runif(1, c(min=eta1_min1, eta1_min2)[informativeness], max=c(min=eta1_max1, eta1_max2)[informativeness]), "sigma_conU2_tilde"=runif(1, c(min=eta2_min1, eta2_min2)[informativeness], max=c(min=eta2_max1, eta2_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=2)
      # RUN MyUnm_glm
      Posteriors <- MyUnm_glm(form1=conY~conX+binC+conC+nomC1+nomC2+conU1_tilde+conU2_tilde, form2=conU1_tilde~conX-1, form3=conU2_tilde~conX-1, family1=gaussian(), family2=gaussian(), family3=gaussian(), data=JAGSData, n.iter=LongIter-LongBurnin, n.adapt=LongBurnin, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU1, ",", 1/c(betaconU1_sd1, betaconU1_sd2)[informativeness]^2, ")"), "lambda[2]"=paste0("dnorm(", betaconU2, ",", 1/c(betaconU2_sd1, betaconU2_sd2)[informativeness]^2, ")"), "gamma[1]"=paste0("dnorm(", alpha1conX, ",", 1/c(alpha1conX_sd1, alpha1conX_sd2)[informativeness]^2, ")"), "delta[1]"=paste0("dnorm(", alpha2conX, ",", 1/c(alpha2conX_sd1, alpha2conX_sd2)[informativeness]^2, ")")), response_nuisance_priors="sigma_{y} <- tau_{y}^-0.5; tau_{y}~dgamma(0.001, 0.001)", confounder1_nuisance_priors=paste0("tau_{u1} <- sigma_{u1}^-2; sigma_{u1}~dunif(", c(eta1_min1, eta1_min2)[informativeness], ",", c(eta1_max1, eta1_max2)[informativeness], ")"), confounder2_nuisance_priors=paste0("tau_{u2} <- sigma_{u2}^-2; sigma_{u2}~dunif(", c(eta2_min1, eta2_min2)[informativeness], ",", c(eta2_max1, eta2_max2)[informativeness], ")"), response_params_to_track="sigma_{y}", confounder1_params_to_track="sigma_{u1}", confounder2_params_to_track="sigma_{u2}", inits=list(Inits1))
    } else {
      Inits1 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "tau_conU1_tilde"=rgamma(1, shape=0.1, rate=1), "tau_conU2_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=3)
      Posteriors <- MyUnm_glm(form1=conY~conX+binC+conC+nomC1+nomC2+conU1_tilde+conU2_tilde, form2=conU1_tilde~conX-1, form3=conU2_tilde~conX-1, family1=gaussian(), family2=gaussian(), family3=gaussian(), data=JAGSData, n.iter=LongIter-LongBurnin, n.adapt=LongBurnin, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU1, ",", 1/betaconU1_sd2^2, ")"), "lambda[2]"=paste0("dnorm(", betaconU2, ",", 1/betaconU2_sd2^2, ")"), "gamma[1]"=paste0("dnorm(", alpha1conX, ",", 1/alpha1conX_sd2^2, ")"), "delta[1]"=paste0("dnorm(", alpha2conX, ",", 1/alpha2conX_sd2^2, ")")), response_nuisance_priors="sigma_{y} <- tau_{y}^-0.5; tau_{y}~dgamma(0.001, 0.001)", confounder1_nuisance_priors=paste0("sigma_{u1} <- tau_{u1}^-0.5; tau_{u1}~dgamma(", inverseeta1squared_shape4, ",", inverseeta1squared_rate4, ")"), confounder2_nuisance_priors=paste0("sigma_{u2} <- tau_{u2}^-0.5; tau_{u2}~dgamma(", inverseeta2squared_shape4, ",", inverseeta2squared_rate4, ")"), response_params_to_track="sigma_{y}", confounder1_params_to_track="sigma_{u1}", confounder2_params_to_track="sigma_{u2}", inits=list(Inits1))
    }
    
    # CREATE A PNG FILE
    png(paste0("../Results/TracePlot", Informativeness[informativeness], "PriorsOneChainDataset", dataset, ".png"), height=10, width=10, units="in", res=300)
    # CREATE A TRACE PLOT FOR EACH MODEL PARAMETER
    par(mfrow=c(ceiling(length(Params)/2),2), mar=c(4,4,0.5,0.5))
    for (param in Params)
    {
      plot(as.array(Posteriors)[,param], type="l", xlab="Iterations including burn-ins", ylab=ParamNames[which(Params==param)], ylim=extendrange(c(Posteriors[,param], TrueValues$true[TrueValues$par==param])))
      # add a red horizontal line of the true value
      abline(h=TrueValues$true[TrueValues$par==param], col=2)
    }
    graphics.off()
  }
  
  # FOR Dataset_1.csv TO Dataset_10.csv,
  for (dataset in 1:10)
  {
    # READ SIMULATED DATA
    SimulatedData <- read.csv(file=paste0("../../Data/Dataset_", dataset, ".csv"))
    SimulatedData$nomC <- SimulatedData$cnom
    SimulatedData$binC <- SimulatedData$cbin
    SimulatedData$conC <- SimulatedData$ccon
    SimulatedData$conX <- SimulatedData$xcon
    SimulatedData$conU1 <- SimulatedData$ucon1
    SimulatedData$conU2 <- SimulatedData$ucon2
    SimulatedData$conY <- SimulatedData$ycon
    
    # MEAN-CENTRE COVARIATES
    SimulatedData$conX_cen <- SimulatedData$conX-mean(SimulatedData$conX)
    SimulatedData$binC_cen <- SimulatedData$binC-mean(SimulatedData$binC)
    SimulatedData$conC_cen <- SimulatedData$conC-mean(SimulatedData$conC)
    SimulatedData$nomC <- factor(SimulatedData$nomC, levels=0:2)
    SimulatedData$nomC1 <- as.numeric(SimulatedData$nomC==1)
    SimulatedData$nomC1_cen <- SimulatedData$nomC1-mean(SimulatedData$nomC1)
    SimulatedData$nomC2 <- as.numeric(SimulatedData$nomC==2)
    SimulatedData$nomC2_cen <- SimulatedData$nomC2-mean(SimulatedData$nomC2)
    
    # CREATE A DATA FRAME FOR RUNNING JAGS
    JAGSData <- data.frame(conY=SimulatedData[,"conY"], conX=SimulatedData[,"conX_cen"], binC=SimulatedData[,"binC_cen"], conC=SimulatedData[,"conC_cen"], nomC1=SimulatedData[,"nomC1_cen"], nomC2=SimulatedData[,"nomC2_cen"], conU1_tilde=1:nrow(SimulatedData), conU2_tilde=1:nrow(SimulatedData))
    # SET conU1_tilde (unmeasured continuous U1) to NA
    JAGSData$conU1_tilde[1:nrow(JAGSData)] <- NA
    # SET conU2_tilde (unmeasured continuous U2) to NA
    JAGSData$conU2_tilde[1:nrow(JAGSData)] <- NA
    
    # SET SEED FOR GENERATING INITIAL VALUES
    set.seed(3)
    
    if (informativeness==2)
    {
      # GENERATE THREE SETS OF INITIAL VALUES AND SET .RNG.name AND .RNG.seed
      Inits2 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "sigma_conU1_tilde"=runif(1, c(min=eta1_min1, eta1_min2)[informativeness], max=c(min=eta1_max1, eta1_max2)[informativeness]), "sigma_conU2_tilde"=runif(1, c(min=eta2_min1, eta2_min2)[informativeness], max=c(min=eta2_max1, eta2_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=4)
      Inits3 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "sigma_conU1_tilde"=runif(1, c(min=eta1_min1, eta1_min2)[informativeness], max=c(min=eta1_max1, eta1_max2)[informativeness]), "sigma_conU2_tilde"=runif(1, c(min=eta2_min1, eta2_min2)[informativeness], max=c(min=eta2_max1, eta2_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=5)
      Inits4 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "sigma_conU1_tilde"=runif(1, c(min=eta1_min1, eta1_min2)[informativeness], max=c(min=eta1_max1, eta1_max2)[informativeness]), "sigma_conU2_tilde"=runif(1, c(min=eta2_min1, eta2_min2)[informativeness], max=c(min=eta2_max1, eta2_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=6)
      # RUN MyUnm_glm
      Posteriors <- MyUnm_glm(form1=conY~conX+binC+conC+nomC1+nomC2+conU1_tilde+conU2_tilde, form2=conU1_tilde~conX-1, form3=conU2_tilde~conX-1, family1=gaussian(), family2=gaussian(), family3=gaussian(), data=JAGSData, n.iter=Iter-Burnin, n.adapt=Burnin, n.chains=Chains, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU1, ",", 1/c(betaconU1_sd1, betaconU1_sd2)[informativeness]^2, ")"), "lambda[2]"=paste0("dnorm(", betaconU2, ",", 1/c(betaconU2_sd1, betaconU2_sd2)[informativeness]^2, ")"), "gamma[1]"=paste0("dnorm(", alpha1conX, ",", 1/c(alpha1conX_sd1, alpha1conX_sd2)[informativeness]^2, ")"), "delta[1]"=paste0("dnorm(", alpha2conX, ",", 1/c(alpha2conX_sd1, alpha2conX_sd2)[informativeness]^2, ")")), response_nuisance_priors="sigma_{y} <- tau_{y}^-0.5; tau_{y}~dgamma(0.001, 0.001)", confounder1_nuisance_priors=paste0("tau_{u1} <- sigma_{u1}^-2; sigma_{u1}~dunif(", c(eta1_min1, eta1_min2)[informativeness], ",", c(eta1_max1, eta1_max2)[informativeness], ")"), confounder2_nuisance_priors=paste0("tau_{u2} <- sigma_{u2}^-2; sigma_{u2}~dunif(", c(eta2_min1, eta2_min2)[informativeness], ",", c(eta2_max1, eta2_max2)[informativeness], ")"), response_params_to_track="sigma_{y}", confounder1_params_to_track="sigma_{u1}", confounder2_params_to_track="sigma_{u2}", inits=list(Inits2, Inits3, Inits4))
    } else {
      Inits2 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "tau_conU1_tilde"=rgamma(1, shape=0.1, rate=1), "tau_conU2_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=4)
      Inits3 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "tau_conU1_tilde"=rgamma(1, shape=0.1, rate=1), "tau_conU2_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=5)
      Inits4 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "tau_conU1_tilde"=rgamma(1, shape=0.1, rate=1), "tau_conU2_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=6)
      Posteriors <- MyUnm_glm(form1=conY~conX+binC+conC+nomC1+nomC2+conU1_tilde+conU2_tilde, form2=conU1_tilde~conX-1, form3=conU2_tilde~conX-1, family1=gaussian(), family2=gaussian(), family3=gaussian(), data=JAGSData, n.iter=Iter, n.adapt=0, n.chains=Chains, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU1, ",", 1/betaconU1_sd2^2, ")"), "lambda[2]"=paste0("dnorm(", betaconU2, ",", 1/betaconU2_sd2^2, ")"), "gamma[1]"=paste0("dnorm(", alpha1conX, ",", 1/alpha1conX_sd2^2, ")"), "delta[1]"=paste0("dnorm(", alpha2conX, ",", 1/alpha2conX_sd2^2, ")")), response_nuisance_priors="sigma_{y} <- tau_{y}^-0.5; tau_{y}~dgamma(0.001, 0.001)", confounder1_nuisance_priors=paste0("sigma_{u1} <- tau_{u1}^-0.5; tau_{u1}~dgamma(", inverseeta1squared_shape4, ",", inverseeta1squared_rate4, ")"), confounder2_nuisance_priors=paste0("sigma_{u2} <- tau_{u2}^-0.5; tau_{u2}~dgamma(", inverseeta2squared_shape4, ",", inverseeta2squared_rate4, ")"), response_params_to_track="sigma_{y}", confounder1_params_to_track="sigma_{u1}", confounder2_params_to_track="sigma_{u2}", inits=list(Inits2, Inits3, Inits4))
    }
    
    # CHECK IF THE ARRAY OF POSTERIOR SAMPLES INCLUDE SAMPLES FROM THE BURN-IN PERIOD
    # IF THE POSTERIOR SAMPLES DO NOT INCLUDE SAMPLES FROM THE BURN-IN PERIOD
    if (summary(Posteriors)$start>1)
    {
      # CREATE A PNG FILE
      png(paste0("../Results/TracePlot", Informativeness[informativeness], "PriorsThreeChainsDataset", dataset, ".png"), height=10, width=10, units="in", res=300)
      # CREATE A TRACE PLOT FOR EACH MODEL PARAMETER
      par(mfrow=c(ceiling(length(Params)/2),2), mar=c(4,4,0.5,0.5))
      for (param in Params)
      {
        plot(as.array(Posteriors)[,param,1], type="l", ylim=extendrange(c(Posteriors[,param,], TrueValues$true[TrueValues$par==param])), xlab=paste0("Iterations after ", Burnin, " burn-ins"), ylab=ParamNames[which(Params==param)])
        for (chain in 2:dim(as.array(Posteriors))[[3]])
        {
          lines(as.array(Posteriors)[,param,chain], col=c(1,3,4)[chain])
        }
        # add a red horizontal line of the true value
        abline(h=TrueValues$true[TrueValues$par==param], col=2)
      }
      graphics.off()
      
      # CALCULATE RHAT
      Rhat <- gelman.diag(Posteriors)[[1]][,2]
      write.csv(t(Rhat), paste0("../Results/Rhat", Informativeness[informativeness], "PriorsDataset", dataset, ".csv"), row.names=F)
      
      # CALCULATE CHAIN-SPECIFIC EFFECTIVE SAMPLE SIZES (ESS)
      for (chain in 1:summary(Posteriors)$nchain)
      {
        ESS <- lapply(Posteriors, effectiveSize)[[chain]][Params]
        write.csv(t(ESS), paste0("../Results/ESS", Informativeness[informativeness], "PriorsDataset", dataset, "Chain", chain, ".csv"), row.names=F)
      }
      
      # IF THE POSTERIOR SAMPLES INCLUDES SAMPLES FROM THE BURN-IN PERIOD
    } else if (summary(Posteriors)$start==1)
    {
      png(paste0("../Results/TracePlot", Informativeness[informativeness], "PriorsThreeChainsDataset", dataset, ".png"), height=10, width=10, units="in", res=300)
      par(mfrow=c(ceiling(length(Params)/2),2), mar=c(4,4,0.5,0.5))
      for (param in Params)
      {
        plot(as.array(Posteriors)[,param,1], type="l", ylim=extendrange(c(Posteriors[,param,], TrueValues$true[TrueValues$par==param])), xlab="Iterations including burn-ins", ylab=ParamNames[which(Params==param)])
        for (chain in 2:dim(as.array(Posteriors))[[3]])
        {
          lines(as.array(Posteriors)[,param,chain], col=c(1,3,4)[chain])
        }
        abline(h=TrueValues$true[TrueValues$par==param], col=2)
      }
      graphics.off()
      
      Rhat <- gelman.diag(mcmc.list(mcmc(as.array(Posteriors)[(Burnin+1):Iter,,1], start=Burnin+1, end=Iter), mcmc(as.array(Posteriors)[(Burnin+1):Iter,,2], start=Burnin+1, end=Iter), mcmc(as.array(Posteriors)[(Burnin+1):Iter,,3], start=Burnin+1, end=Iter)))[[1]][,2]
      write.csv(t(Rhat), paste0("../Results/", Informativeness[informativeness], "PriorsDataset", dataset, "Rhat.csv"), row.names=F)
      
      for (chain in 1:summary(Posteriors)$nchain)
      {
        ESS <- effectiveSize(mcmc(as.array(Posteriors)[(Burnin+1):Iter,,chain], start=Burnin+1, end=Iter))
        write.csv(t(ESS), paste0("../Results/", Informativeness[informativeness], "PriorsDataset", dataset, "ESS", 1+chain, ".csv"), row.names=F)
      }
    }
  }
}

################################################################################
# D. Run Bayesian QBA on all 500 datasets
################################################################################

# FOR Dataset_1.csv TO Dataset_500.csv,
for (dataset in 1:500)
{
  # READ SMULATED DATA
  SimulatedData <- read.csv(file=paste0("../../Data/Dataset_", dataset, ".csv"))
  SimulatedData$nomC <- SimulatedData$cnom
  SimulatedData$binC <- SimulatedData$cbin
  SimulatedData$conC <- SimulatedData$ccon
  SimulatedData$conX <- SimulatedData$xcon
  SimulatedData$conU1 <- SimulatedData$ucon1
  SimulatedData$conU2 <- SimulatedData$ucon2
  SimulatedData$conY <- SimulatedData$ycon
  
  # MEAN-CENTRE COVARIATES
  SimulatedData$conX_cen <- SimulatedData$conX-mean(SimulatedData$conX)
  SimulatedData$binC_cen <- SimulatedData$binC-mean(SimulatedData$binC)
  SimulatedData$conC_cen <- SimulatedData$conC-mean(SimulatedData$conC)
  SimulatedData$nomC <- factor(SimulatedData$nomC, levels=0:2)
  SimulatedData$nomC1 <- as.numeric(SimulatedData$nomC==1)
  SimulatedData$nomC1_cen <- SimulatedData$nomC1-mean(SimulatedData$nomC1)
  SimulatedData$nomC2 <- as.numeric(SimulatedData$nomC==2)
  SimulatedData$nomC2_cen <- SimulatedData$nomC2-mean(SimulatedData$nomC2)
  
  for (informativeness in c(1:2, 4))
  {
    # CREATE A DATA FRAME FOR RUNNING JAGS
    JAGSData <- data.frame(conY=SimulatedData[,"conY"], conX=SimulatedData[,"conX_cen"], binC=SimulatedData[,"binC_cen"], conC=SimulatedData[,"conC_cen"], nomC1=SimulatedData[,"nomC1_cen"], nomC2=SimulatedData[,"nomC2_cen"], conU1_tilde=1:nrow(SimulatedData), conU2_tilde=1:nrow(SimulatedData))
    # SET conU1_tilde (unmeasured continuous U1) to NA
    JAGSData$conU1_tilde[1:nrow(JAGSData)] <- NA
    # SET conU2_tilde (unmeasured continuous U2) to NA
    JAGSData$conU2_tilde[1:nrow(JAGSData)] <- NA
    
    # SET SEED FOR GENERATING INITIAL VALUES
    set.seed(Seeds[1000*(informativeness-1)+2*(dataset-1)+1])
    
    if (informativeness>=1 & informativeness<=2)
    {
      # GENERATE A SET OF INITIAL VALUES AND SET .RNG.name AND .RNG.seed
      Inits5 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "sigma_conU1_tilde"=runif(1, c(min=eta1_min1, eta1_min2)[informativeness], max=c(min=eta1_max1, eta1_max2)[informativeness]), "sigma_conU2_tilde"=runif(1, c(min=eta2_min1, eta2_min2)[informativeness], max=c(min=eta2_max1, eta2_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=Seeds[1000*(informativeness-1)+2*(dataset-1)+2])
      # RECORD CURRENT DATE AND TIME
      StartTime <- Sys.time()
      # RUN MyUnm_glm
      Posteriors <- MyUnm_glm(form1=conY~conX+binC+conC+nomC1+nomC2+conU1_tilde+conU2_tilde, form2=conU1_tilde~conX-1, form3=conU2_tilde~conX-1, family1=gaussian(), family2=gaussian(), family3=gaussian(), data=JAGSData, n.iter=Iter-Burnin, n.adapt=Burnin, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU1, ",", 1/c(betaconU1_sd1, betaconU1_sd2)[informativeness]^2, ")"), "lambda[2]"=paste0("dnorm(", betaconU2, ",", 1/c(betaconU2_sd1, betaconU2_sd2)[informativeness]^2, ")"), "gamma[1]"=paste0("dnorm(", alpha1conX, ",", 1/c(alpha1conX_sd1, alpha1conX_sd2)[informativeness]^2, ")"), "delta[1]"=paste0("dnorm(", alpha2conX, ",", 1/c(alpha2conX_sd1, alpha2conX_sd2)[informativeness]^2, ")")), response_nuisance_priors="sigma_{y} <- tau_{y}^-0.5; tau_{y}~dgamma(0.001, 0.001)", confounder1_nuisance_priors=paste0("tau_{u1} <- sigma_{u1}^-2; sigma_{u1}~dunif(", c(eta1_min1, eta1_min2)[informativeness], ",", c(eta1_max1, eta1_max2)[informativeness], ")"), confounder2_nuisance_priors=paste0("tau_{u2} <- sigma_{u2}^-2; sigma_{u2}~dunif(", c(eta2_min1, eta2_min2)[informativeness], ",", c(eta2_max1, eta2_max2)[informativeness], ")"), response_params_to_track="sigma_{y}", confounder1_params_to_track="sigma_{u1}", confounder2_params_to_track="sigma_{u2}", inits=list(Inits5))
      # RECORD CURRENT DATE AND TIME AGAIN TO CALCULATE RUN TIME FOR MyUnm_glm
      RunTime <- difftime(Sys.time(), StartTime, units='mins')
      
      # RECORD SUMMARY RESULTS INCLUDING THE MEDIAN, 2.5TH AND 97.5TH PERCENTILES OF EXPOSURE EFFECT ESTIMATE AND RUNTIME
      write.table(data.frame(t(c(summary(Posteriors)$statistics[Params,"Mean"], summary(Posteriors)$statistics[Params,"SD"], summary(Posteriors)$statistics[Params,"Naive SE"], summary(Posteriors)$statistics[Params,"Time-series SE"], summary(Posteriors)$quantiles[Params,"2.5%"], summary(Posteriors)$quantiles[Params,"25%"], summary(Posteriors)$quantiles[Params,"50%"], summary(Posteriors)$quantiles[Params,"75%"], summary(Posteriors)$quantiles[Params,"97.5%"], effectiveSize(Posteriors)[Params], as.numeric(unlist(Inits5)[setdiff(names(unlist(Inits5)), ".RNG.name")]), attr(Posteriors, "n.adapt"), attr(Posteriors, "n.adapt")+attr(Posteriors, "n.iter"), RunTime))), file=paste0("../Results/Summary", Informativeness[informativeness], "PriorsDataset", dataset, ".csv"), append=F, sep=',', row.names=F, col.names=c(paste0(Params, "_mean"), paste0(Params, "_sd"), paste0(Params, "_naivese"), paste0(Params, "_timeseriesse"), paste0(Params, "_lci"), paste0(Params, "_lqt"), paste0(Params, "_median"), paste0(Params, "_uqt"), paste0(Params, "_uci"), paste0(Params, "_es"), paste0(setdiff(names(unlist(Inits5)), ".RNG.name"), "_inits"), "Burnins", "TotalIterations", "RunTime"))
    } else {
      Inits5 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(2, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "delta"=rnorm(1, mean=0, sd=1), "tau_conY"=rgamma(1, shape=0.1, rate=1), "tau_conU1_tilde"=rgamma(1, shape=0.1, rate=1), "tau_conU2_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=Seeds[1000*(informativeness-1)+2*(dataset-1)+2])
      StartTime <- Sys.time()
      Posteriors <- MyUnm_glm(form1=conY~conX+binC+conC+nomC1+nomC2+conU1_tilde+conU2_tilde, form2=conU1_tilde~conX-1, form3=conU2_tilde~conX-1, family1=gaussian(), family2=gaussian(), family3=gaussian(), data=JAGSData, n.iter=Iter, n.adapt=0, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU1, ",", 1/betaconU1_sd2^2, ")"), "lambda[2]"=paste0("dnorm(", betaconU2, ",", 1/betaconU2_sd2^2, ")"), "gamma[1]"=paste0("dnorm(", alpha1conX, ",", 1/alpha1conX_sd2^2, ")"), "delta[1]"=paste0("dnorm(", alpha2conX, ",", 1/alpha2conX_sd2^2, ")")), response_nuisance_priors="sigma_{y} <- tau_{y}^-0.5; tau_{y}~dgamma(0.001, 0.001)", confounder1_nuisance_priors=paste0("sigma_{u1} <- tau_{u1}^-0.5; tau_{u1}~dgamma(", inverseeta1squared_shape4, ",", inverseeta1squared_rate4, ")"), confounder2_nuisance_priors=paste0("sigma_{u2} <- tau_{u2}^-0.5; tau_{u2}~dgamma(", inverseeta2squared_shape4, ",", inverseeta2squared_rate4, ")"), response_params_to_track="sigma_{y}", confounder1_params_to_track="sigma_{u1}", confounder2_params_to_track="sigma_{u2}", inits=list(Inits5))
      RunTime <- difftime(Sys.time(), StartTime, units='mins')
      
      write.table(data.frame(t(c(apply(as.array(Posteriors)[(Burnin+1):Iter,Params], 2, mean), apply(as.array(Posteriors)[(Burnin+1):Iter,Params], 2, sd), apply(as.array(Posteriors)[(Burnin+1):Iter,Params], 2, function(x) quantile(x, probs=0.025)), apply(as.array(Posteriors)[(Burnin+1):Iter,Params], 2, function(x) quantile(x, probs=0.25)), apply(as.array(Posteriors)[(Burnin+1):Iter,Params], 2, function(x) quantile(x, probs=0.5)), apply(as.array(Posteriors)[(Burnin+1):Iter,Params], 2, function(x) quantile(x, probs=0.75)), apply(as.array(Posteriors)[(Burnin+1):Iter,Params], 2, function(x) quantile(x, probs=0.975)), as.numeric(unlist(Inits5)[setdiff(names(unlist(Inits5)), ".RNG.name")]), attr(Posteriors, "n.adapt")+Burnin, attr(Posteriors, "n.adapt")+attr(Posteriors, "n.iter"), RunTime))), file=paste0("../Results/Summary", Informativeness[informativeness], "PriorsDataset", dataset, ".csv"), append=F, sep=',', row.names=F, col.names=c(paste0(Params, "_mean"), paste0(Params, "_sd"), paste0(Params, "_lci"), paste0(Params, "_lqt"), paste0(Params, "_median"), paste0(Params, "_uqt"), paste0(Params, "_uci"), paste0(setdiff(names(unlist(Inits5)), ".RNG.name"), "_inits"), "Burnins", "TotalIterations", "RunTime"))
    }
  }
}
