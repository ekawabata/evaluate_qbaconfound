################################################################################
# RUNS SIMULATION STUDY I FOR BAYESIAN QBA
# SCENARIO A: Y BINARY, X BINARY, SINGLE U=(U1) CONTINUOUS

# A. Preliminaries
# loads necessary packages and functions

# B. Creates objects necessary to perform simulation study I Scenario A
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
# B. Create objects necessary to perform simulation study I Scenario A
################################################################################

# LIST MODEL PARAMETERS
# model parameter names using JAGS notation (obtained from dimnames(as.array(Posteriors))[2]$var)
Params <- c("beta[1]", "beta[binX]", "beta[binC]", "beta[conC]", "beta[nomC1]", "beta[nomC2]", "lambda[conU[tilde]]", "gamma[binX]", "sigma_conU[tilde]")
# model parameter names using our manuscript's notation
ParamNames <- c("beta0", "betabinX", "betabinC", "betaconC", "betanomC1", "betanomC2", "betaconU", "alphabinX", "eta")

# LIST TRUE VALUES
# coefficients in Y|X,C,U1,U2
beta0 <- -2.37
betabinX <- 0
betabinC <- -0.3559
betaconC <- -0.009095
betanomC1 <- 0.1753
betanomC2 <- 0.04492
betaconU <- 0.06588
# coefficient in U|X,C
alphabinX <- 6.636
# residual standard deviation in U|X
eta <- 9.523
# combine to create a data frame
TrueValues <- data.frame('par'=Params, 'true'=c(beta0, betabinX, betabinC, betaconC, betanomC1, betanomC2, betaconU, alphabinX, eta))

# LIST VERY INFORMATIVE UNIFORM PRIORS
betaconU_sd1 <- 0.003361
alphabinX_sd1 <- 0.3386
eta_min1 <- 8.521
eta_max1 <- 10.53

# LIST INFORMATIVE UNIFORM PRIORS
betaconU_sd2 <- 0.01008
alphabinX_sd2 <- 1.016
eta_min2 <- 6.516
eta_max2 <- 12.53

# LIST INFORMATIVE GAMMA PRIORS
tauU_shape4 <- 13.19
tauU_rate4 <- 1196

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

for (informativeness in c(2,4))
{
  # FOR Dataset_1.csv,
  for (dataset in 1)
  {
    # READ SMULATED DATA
    SimulatedData <- read.csv(file=paste0("../../Data/Dataset_", dataset, ".csv"))
    SimulatedData$nomC <- SimulatedData$cnom
    SimulatedData$binC <- SimulatedData$cbin
    SimulatedData$conC <- SimulatedData$ccon
    SimulatedData$binX <- SimulatedData$xbin
    SimulatedData$conU <- SimulatedData$ucon1
    SimulatedData$binY <- SimulatedData$ybin
    
    # MEAN-CENTRE COVARIATES
    SimulatedData$binX_cen <- SimulatedData$binX-mean(SimulatedData$binX)
    SimulatedData$binC_cen <- SimulatedData$binC-mean(SimulatedData$binC)
    SimulatedData$conC_cen <- SimulatedData$conC-mean(SimulatedData$conC)
    SimulatedData$nomC <- factor(SimulatedData$nomC, levels=0:2)
    SimulatedData$nomC1 <- as.numeric(SimulatedData$nomC==1)
    SimulatedData$nomC1_cen <- SimulatedData$nomC1-mean(SimulatedData$nomC1)
    SimulatedData$nomC2 <- as.numeric(SimulatedData$nomC==2)
    SimulatedData$nomC2_cen <- SimulatedData$nomC2-mean(SimulatedData$nomC2)
    
    # CREATE A DATA FRAME FOR RUNNING JAGS
    JAGSData <- data.frame(binY=SimulatedData[,"binY"], binX=SimulatedData[,"binX_cen"], binC=SimulatedData[,"binC_cen"], conC=SimulatedData[,"conC_cen"], nomC1=SimulatedData[,"nomC1_cen"], nomC2=SimulatedData[,"nomC2_cen"], conU_tilde=1:nrow(SimulatedData))
    # SET conU_tilde (part of conU uncorrelated with C) to NA
    JAGSData$conU_tilde[1:nrow(JAGSData)] <- NA
    
    # SET SEED FOR GENERATING INITIAL VALUES
    set.seed(1)
    
    if (informativeness==2)
    {
      # GENERATE A SET OF INITIAL VALUES AND SET .RNG.name AND .RNG.seed
      Inits1 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "sigma_conU_tilde"=runif(1, c(min=eta_min1, eta_min2)[informativeness], max=c(min=eta_max1, eta_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=2)
      # RUN MyUnm_glm
      Posteriors <- MyUnm_glm(form1=binY~binX+binC+conC+nomC1+nomC2+conU_tilde, form2=conU_tilde~binX-1, family1=binomial(), family2=gaussian(), data=JAGSData, n.iter=LongIter-LongBurnin, n.adapt=LongBurnin, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU, ",", 1/c(betaconU_sd1, betaconU_sd2)[informativeness]^2, ")"), "gamma[1]"=paste0("dnorm(", alphabinX, ",", 1/c(alphabinX_sd1, alphabinX_sd2)[informativeness]^2, ")")), confounder1_nuisance_priors=paste0("tau_{u1} <- sigma_{u1}^-2; sigma_{u1}~dunif(", c(eta_min1, eta_min2)[informativeness], ",", c(eta_max1, eta_max2)[informativeness], ")"), confounder1_params_to_track="sigma_{u1}", inits=list(Inits1))
    } else {
      Inits1 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "tau_conU_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=2)
      Posteriors <- MyUnm_glm(form1=binY~binX+binC+conC+nomC1+nomC2+conU_tilde, form2=conU_tilde~binX-1, family1=binomial(), family2=gaussian(), data=JAGSData, n.iter=LongIter-LongBurnin, n.adapt=LongBurnin, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU, ",", 1/betaconU_sd2^2, ")"), "gamma[1]"=paste0("dnorm(", alphabinX, ",", 1/alphabinX_sd2^2, ")")), confounder1_nuisance_priors=paste0("sigma_{u1} <- tau_{u1}^-0.5; tau_{u1}~dgamma(", tauU_shape4, ",", tauU_rate4, ")"), confounder1_params_to_track="sigma_{u1}", inits=list(Inits1))
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
    SimulatedData$binX <- SimulatedData$xbin
    SimulatedData$conU <- SimulatedData$ucon1
    SimulatedData$binY <- SimulatedData$ybin
    
    # MEAN-CENTRE COVARIATES
    SimulatedData$binX_cen <- SimulatedData$binX-mean(SimulatedData$binX)
    SimulatedData$binC_cen <- SimulatedData$binC-mean(SimulatedData$binC)
    SimulatedData$conC_cen <- SimulatedData$conC-mean(SimulatedData$conC)
    SimulatedData$nomC <- factor(SimulatedData$nomC, levels=0:2)
    SimulatedData$nomC1 <- as.numeric(SimulatedData$nomC==1)
    SimulatedData$nomC1_cen <- SimulatedData$nomC1-mean(SimulatedData$nomC1)
    SimulatedData$nomC2 <- as.numeric(SimulatedData$nomC==2)
    SimulatedData$nomC2_cen <- SimulatedData$nomC2-mean(SimulatedData$nomC2)
    
    # CREATE A DATA FRAME FOR RUNNING JAGS
    JAGSData <- data.frame(binY=SimulatedData[,"binY"], binX=SimulatedData[,"binX_cen"], binC=SimulatedData[,"binC_cen"], conC=SimulatedData[,"conC_cen"], nomC1=SimulatedData[,"nomC1_cen"], nomC2=SimulatedData[,"nomC2_cen"], conU_tilde=1:nrow(SimulatedData))
    # SET conU_tilde (unmeasured continuous U) to NA
    JAGSData$conU_tilde[1:nrow(JAGSData)] <- NA
    
    # SET SEED FOR GENERATING INITIAL VALUES
    set.seed(3)
    
    if (informativeness==2)
    {
      # GENERATE THREE SETS OF INITIAL VALUES AND SET .RNG.name AND .RNG.seed
      Inits2 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "sigma_conU_tilde"=runif(1, c(min=eta_min1, eta_min2)[informativeness], max=c(min=eta_max1, eta_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=4)
      Inits3 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "sigma_conU_tilde"=runif(1, c(min=eta_min1, eta_min2)[informativeness], max=c(min=eta_max1, eta_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=5)
      Inits4 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "sigma_conU_tilde"=runif(1, c(min=eta_min1, eta_min2)[informativeness], max=c(min=eta_max1, eta_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=6)
      # RUN MyUnm_glm
      Posteriors <- MyUnm_glm(form1=binY~binX+binC+conC+nomC1+nomC2+conU_tilde, form2=conU_tilde~binX-1, family1=binomial(), family2=gaussian(), data=JAGSData, n.iter=Iter-Burnin, n.adapt=Burnin, n.chains=Chains, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU, ",", 1/c(betaconU_sd1, betaconU_sd2)[informativeness]^2, ")"), "gamma[1]"=paste0("dnorm(", alphabinX, ",", 1/c(alphabinX_sd1, alphabinX_sd2)[informativeness]^2, ")")), confounder1_nuisance_priors=paste0("tau_{u1} <- sigma_{u1}^-2; sigma_{u1}~dunif(", c(eta_min1, eta_min2)[informativeness], ",", c(eta_max1, eta_max2)[informativeness], ")"), confounder1_params_to_track="sigma_{u1}", inits=list(Inits2, Inits3, Inits4))
    } else {
      Inits2 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "tau_conU_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=4)
      Inits3 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "tau_conU_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=5)
      Inits4 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "tau_conU_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=6)
      Posteriors <- MyUnm_glm(form1=binY~binX+binC+conC+nomC1+nomC2+conU_tilde, form2=conU_tilde~binX-1, family1=binomial(), family2=gaussian(), data=JAGSData, n.iter=Iter-Burnin, n.adapt=Burnin, n.chains=Chains, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU, ",", 1/betaconU_sd2^2, ")"), "gamma[1]"=paste0("dnorm(", alphabinX, ",", 1/alphabinX_sd2^2, ")")), confounder1_nuisance_priors=paste0("sigma_{u1} <- tau_{u1}^-0.5; tau_{u1}~dgamma(", tauU_shape4, ",", tauU_rate4, ")"), confounder1_params_to_track="sigma_{u1}", inits=list(Inits2, Inits3, Inits4))
    }
    
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
  }
}

################################################################################
# D. Run Bayesian QBA on all 500 datasets
################################################################################

# FOR Dataset_1.csv TO Dataset_500.csv,
for (dataset in 1:500)
{
  # READ SIMULATED DATA
  SimulatedData <- read.csv(file=paste0("../../Data/Dataset_", dataset, ".csv"))
  SimulatedData$nomC <- SimulatedData$cnom
  SimulatedData$binC <- SimulatedData$cbin
  SimulatedData$conC <- SimulatedData$ccon
  SimulatedData$binX <- SimulatedData$xbin
  SimulatedData$conU <- SimulatedData$ucon1
  SimulatedData$binY <- SimulatedData$ybin
  
  # MEAN-CENTRE COVARIATES
  SimulatedData$binX_cen <- SimulatedData$binX-mean(SimulatedData$binX)
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
    JAGSData <- data.frame(binY=SimulatedData[,"binY"], binX=SimulatedData[,"binX_cen"], binC=SimulatedData[,"binC_cen"], conC=SimulatedData[,"conC_cen"], nomC1=SimulatedData[,"nomC1_cen"], nomC2=SimulatedData[,"nomC2_cen"], conU_tilde=1:nrow(SimulatedData))
    # SET conU_tilde (unmeasured continuous U) to NA
    JAGSData$conU_tilde[1:nrow(JAGSData)] <- NA
    
    # SET SEED FOR GENERATING INITIAL VALUES
    set.seed(Seeds[1000*(informativeness-1)+2*(dataset-1)+1])
    
    if (informativeness>=1 & informativeness<=2)
    {
      # GENERATE A SET OF INITIAL VALUES AND SET .RNG.name AND .RNG.seed
      Inits5 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "sigma_conU_tilde"=runif(1, c(min=eta_min1, eta_min2)[informativeness], max=c(min=eta_max1, eta_max2)[informativeness]), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=Seeds[1000*(informativeness-1)+2*(dataset-1)+2])
      # RECORD CURRENT DATE AND TIME
      StartTime <- Sys.time()
      # RUN MyUnm_glm
      Posteriors <- MyUnm_glm(form1=binY~binX+binC+conC+nomC1+nomC2+conU_tilde, form2=conU_tilde~binX-1, family1=binomial(), family2=gaussian(), data=JAGSData, n.iter=Iter-Burnin, n.adapt=Burnin, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU, ",", 1/c(betaconU_sd1, betaconU_sd2)[informativeness]^2, ")"), "gamma[1]"=paste0("dnorm(", alphabinX, ",", 1/c(alphabinX_sd1, alphabinX_sd2)[informativeness]^2, ")")), confounder1_nuisance_priors=paste0("tau_{u1} <- sigma_{u1}^-2; sigma_{u1}~dunif(", c(eta_min1, eta_min2)[informativeness], ",", c(eta_max1, eta_max2)[informativeness], ")"), confounder1_params_to_track="sigma_{u1}", inits=list(Inits5))
      # RECORD CURRENT DATE AND TIME AGAIN TO CALCULATE RUN TIME FOR MyUnm_glm
      RunTime <- difftime(Sys.time(), StartTime, units='mins')
    } else {
      Inits5 <- list("beta"=rnorm(6, mean=0, sd=1), "lambda"=rnorm(1, mean=0, sd=1), "gamma"=rnorm(1, mean=0, sd=1), "tau_conU_tilde"=rgamma(1, shape=0.1, rate=1), ".RNG.name"="base::Wichmann-Hill", ".RNG.seed"=Seeds[1000*(informativeness-1)+2*(dataset-1)+2])
      StartTime <- Sys.time()
      Posteriors <- MyUnm_glm(form1=binY~binX+binC+conC+nomC1+nomC2+conU_tilde, form2=conU_tilde~binX-1, family1=binomial(), family2=gaussian(), data=JAGSData, n.iter=Iter-Burnin, n.adapt=Burnin, n.chains=1, priors=c("beta[1]"="dnorm(0, 0.1)", "beta[2]"="dnorm(0, 0.1)", "beta[3]"="dnorm(0, 0.1)", "beta[4]"="dnorm(0, 0.1)", "beta[5]"="dnorm(0, 0.1)", "beta[6]"="dnorm(0, 0.1)", "lambda[1]"=paste0("dnorm(", betaconU, ",", 1/betaconU_sd2^2, ")"), "gamma[1]"=paste0("dnorm(", alphabinX, ",", 1/alphabinX_sd2^2, ")")), confounder1_nuisance_priors=paste0("sigma_{u1} <- tau_{u1}^-0.5; tau_{u1}~dgamma(", tauU_shape4, ",", tauU_rate4, ")"), confounder1_params_to_track="sigma_{u1}", inits=list(Inits5))
      RunTime <- difftime(Sys.time(), StartTime, units='mins')
    }
    
    # RECORD SUMMARY RESULTS INCLUDING THE MEDIAN, 2.5TH AND 97.5TH PERCENTILES OF EXPOSURE EFFECT ESTIMATE AND RUNTIME
    write.table(data.frame(t(c(summary(Posteriors)$statistics[Params,"Mean"], summary(Posteriors)$statistics[Params,"SD"], summary(Posteriors)$statistics[Params,"Naive SE"], summary(Posteriors)$statistics[Params,"Time-series SE"], summary(Posteriors)$quantiles[Params,"2.5%"], summary(Posteriors)$quantiles[Params,"25%"], summary(Posteriors)$quantiles[Params,"50%"], summary(Posteriors)$quantiles[Params,"75%"], summary(Posteriors)$quantiles[Params,"97.5%"], effectiveSize(Posteriors)[Params], as.numeric(unlist(Inits5)[setdiff(names(unlist(Inits5)), ".RNG.name")]), attr(Posteriors, "n.adapt"), attr(Posteriors, "n.adapt")+attr(Posteriors, "n.iter"), RunTime))), file=paste0("../Results/Summary", Informativeness[informativeness], "PriorsDataset", dataset, ".csv"), append=F, sep=',', row.names=F, col.names=c(paste0(Params, "_mean"), paste0(Params, "_sd"), paste0(Params, "_naivese"), paste0(Params, "_timeseriesse"), paste0(Params, "_lci"), paste0(Params, "_lqt"), paste0(Params, "_median"), paste0(Params, "_uqt"), paste0(Params, "_uci"), paste0(Params, "_es"), paste0(setdiff(names(unlist(Inits5)), ".RNG.name"), "_inits"), "Burnins", "TotalIterations", "RunTime"))
  }
}
