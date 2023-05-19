library(ggplot2) # plotting
library(coda) # test convergence
library(stats)

configureDataSchliep <- function(subdir, sim_data_raw, locList, params, mode = "N", timePts = timePts, locIndices = locIndices) {
  #timePts<-1:params$numGens #time points to "collect"
  #timePts <- seq(from = 1, to = params$numGens, by = 1000/num_samples)
  
  
  # Wrangle data
  locList.asDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
  locDF <- locList.asDF[locIndices,]
  distMx <- dist(locDF, diag = T, upper = T)
  
  Y <- list()
  X <- list()
  for (t_index in 1:length(timePts)) {
    t <- timePts[t_index]
    if (mode == "X") { # configure data differently for extended mode
      # Y: list containing t items where each item is also a list of length n, 
      # each of which is a matrix with S columns and J_(i,j) rows. The values are the ordinal abundance scores for each replicate at location i and time j for species i (column i). 
      Y[[t_index]] <- lapply(X = sim_data_raw[[t]]$y[locIndices], FUN = function(y_s){return(matrix(y_s, ncol = params$numSpecies, nrow = 1))})
    } else {
      #dataSchliep$Y=Y #list containing t items where each item is an n x S matrix of ordinal responses 
      Y[[t_index]] <- matrix(as.double(unlist(sim_data_raw[[t]]$y[locIndices])), nrow = length(locIndices), ncol = params$numSpecies, byrow = TRUE)
    }

    sampledCovs <- list()
    for (covNum in 1:length(sim_data_raw[[t]]$covs)){
        if(paste0("Cov", covNum) %in% params$names_cov){
            sampledCovs[[covNum]] <- sim_data_raw[[t]]$covs[[covNum]][locIndices]
        }
    }
    #cov3<-sim_data_raw[[t]]$covs[[3]]
    X[[t_index]] <- matrix(c(rep(1, times = length(locIndices)), unlist(sampledCovs)), nrow = length(locIndices), ncol = length(params$names_cov) + 1, byrow=FALSE)
  }
  
  dataSchliep = list()
  dataSchliep$t = length(timePts) # number of survey times
  dataSchliep$L = 2 # number of ordinal categories (2 = binary)
  dataSchliep$S = params$numSpecies # number of species (2 or more)
  dataSchliep$p = length(params$names_cov) + 1 # number of covariates in model
  dataSchliep$n = length(locIndices) # number of locations
  dataSchliep$dist = as.matrix(distMx) # distance matrix between locations
  dataSchliep$Y = Y # list containing t items where each item is an n x S matrix of ordinal responses 
  dataSchliep$X = X # list containing t items where each item is an n x p matrix of covariates
  
  if(mode == "X"){
    #J is matrix with n rows and t columns. Each J_(i,j) is the number of replicate observations for plot i and year j. Even if no replicates, this matrix needs to be provided with 1 (when the location is observed that year) or 0 (when the location is not observed that year)
    J = matrix(rep(x = as.integer(1), times = length(locIndices) * length(timePts)), nrow = length(locIndices), ncol = length(timePts))
    dataSchliep$J = J
  }
  
  # dataSchliep is list of the format of HWAdataM.Rdata (from paper supplement)
  #load("../schliep_code/HWAdataMBinary.Rdata")
  # below is from schliep supplement
  return(dataSchliep)
}

#priors=get.priors(HWAdataMBinary,temporalRE=T) 
#pars=get.startingValues(HWAdataMBinary,temporalRE=T) 
#run=Multivariate.Ordinal.Spatial.Model(HWAdataMBinary,pars,priors,temporalRE=T,iters=1000,print.out=100)

testData <- function(data, mode = "N") {
  goodData = TRUE
  # Y: list containing t items where each item is also a list of length n, 
  # each of which is a matrix with S columns and J_(i,j) rows. The values are the ordinal abundance scores for each replicate at location i and time j for species i (column i). 
  
  if(length(data[["Y"]])!=data[["t"]]){
    print("t does not match length of Y")
    goodData = FALSE
  }
  if(length(data[["X"]])!=data[["t"]]){
    print("t does not match length of X")
    goodData = FALSE
  }
  for(time in 1:data[["t"]]){
    if(mode == "X"){
      if(length(unique(unlist(as.vector(data[["Y"]][[time]]))))!=data[["L"]]){
        print(paste("number of ordinal categories (L) doesnt match values in Y at time", time))
        goodData = FALSE
      }
      if(dim(data[["Y"]][[time]][[1]])[1]!=data[["J"]][1,time]){
        print("dimension of Y not matched to J_ij")
        goodData = FALSE
      }
      if(dim(data[["Y"]][[time]][[1]])[2]!=data[["S"]]){
        print("dimension of Y not matched to S")
        goodData = FALSE
      }

      if(length(data[["Y"]][[time]])!=data[["n"]]){
        print("length of Y not matched to n")
        goodData = FALSE
      }
      if(typeof(data[["Y"]][[time]][[1]])!="integer"){
        print("type of Y[[t]] wrong")
        goodData = FALSE
      }
      if(dim(data[["J"]])[1] != data[["n"]]){
        print("J dimension wrong")
        goodData = FALSE
      }
      if(dim(data[["J"]])[2] != data[["t"]]){
        print("J dimension wrong")
        goodData = FALSE
      }
      if(typeof(data[["J"]])!="integer"){
        print("check typeof J" )
        goodData = FALSE
      }
      if(class(data[["J"]])[1]!="matrix"){
        print("check typeof J" )
        goodData = FALSE
      }
    } else {
      if(length(unique(as.vector(data[["Y"]][[time]])))!=data[["L"]]){
        print("number of ordinal categories (L) doesnt match values in Y")
        goodData = FALSE
      }
      if(dim(data[["Y"]][[time]])[1]!=data[["n"]]){
        print("dimension of Y not matched to n")
        goodData = FALSE
      }
      if(dim(data[["Y"]][[time]])[2]!=data[["S"]]){
        print("dimension of Y not matched to S")
        goodData = FALSE
      }
      if(typeof(data[["Y"]][[time]])!="double"){
        print("type of Y[[t]] wrong")
        goodData = FALSE
      }
    }
    if(dim(data[["X"]][[time]])[1]!=data[["n"]]){
      print("dimension of X not matched to n")
      goodData = FALSE
    } 
    if(dim(data[["X"]][[time]])[2]!=data[["p"]]){
      print("dimension of X not matched to p")
      goodData = FALSE
    } 
    if(typeof(data[["X"]][[time]])!="double"){
      print("type of X[[t]] wrong")
      goodData = FALSE
    }
  }
  if(typeof(data[["X"]])!="list"){
    print("X type wrong")
    goodData = FALSE
  }
  if(typeof(data[["Y"]])!="list"){
    print("Y type wrong")
    goodData = FALSE
  }
  if(typeof(data$dist)!="double"){
    print("dist type wrong")
    goodData = FALSE
  }
  if(!(all.equal(data[["dist"]],t(data[["dist"]])))){
    print("distance mx not symmetric")
    goodData = FALSE
  }
  if(dim(data[["dist"]])[1] != data[["n"]]){
    print("dist mx not the right dimensions")
    goodData = FALSE
  }
  if(dim(data[["dist"]])[2] != data[["n"]]){
    print("dist mx not the right dimensions")
    goodData = FALSE
  }
  if(class(data$dist)[1]!="matrix"){
    print("class of dist wrong")
    goodData = FALSE
  }
  
  return(goodData)
}

#############
# plot #
#############
plotPosteriorTable <- function(schliepOutput, ci_percent = 0.95, mode = "N", burn_in = 500) {
  p <- schliepOutput$data$p #number of covariates
  S <- schliepOutput$data$S #number of species
  t <- schliepOutput$data$t #number of time points
  n <- schliepOutput$data$n #number of sites
  iter <- dim(schliepOutput$A)[3]
  
  credibleIntervalsL <- list()
  
  A <- schliepOutput$A[, , burn_in:iter]
  alphamx <- matrix(rep(0, times = S*S), nrow = S, ncol = S)
  for (species1 in 1:S) {
    for (species2 in 1:S) {
      varName = paste0("A",species1,species2)
      mean = mean(A[species1, species2, ])
      quants = quantile(x = A[species1, species2, ], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
      if (species1 == species2) {
        alphamx[species1, species2] <- 0
      } else if (as.numeric(quants[1]) > 0) {
        alphamx[species1, species2] <- 1
        alphamx[species2, species1] <- 1
      } else if (as.numeric(quants[2]) < 0) {
        alphamx[species1, species2] <- -1
        alphamx[species2, species1] <- -1
      } 
      credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
    }
  }
  
  beta <- schliepOutput$beta
  betamx <- matrix(rep(NA, times = S*(p-1)), nrow = S, ncol = (p-1))
  if(mode == "X"){
    beta <- beta[burn_in:iter,]
    i = 1
    for(species in 1:S){
      for(cov in 0:(p-1)){
        varName <- paste0("beta_cov",cov,"_sp",species)
        #### BUG FIXED
        mean = mean(beta[,i])
        quants = quantile(x = beta[,i], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
        if (cov != 0) {
            if (sign(as.numeric(quants[1])) > 0) {
                betamx[species, cov] <- 1
            } else if (sign(as.numeric(quants[2])) < 0) {
                betamx[species, cov] <- -1
            } else {
                betamx[species, cov] <- 0
            }
        }
        credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
        i=i+1
      }
    }
  } else {
    beta <- beta[,,burn_in:iter]
    for(species in 1:S){
      for(cov in 0:(p-1)){
        varName<-paste0("beta_cov",cov,"_sp",species)
        mean = mean(beta[cov+1, species, ])
        quants = quantile(x = beta[cov+1, species, ], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
        if (cov != 0) {
            if (sign(as.numeric(quants[1])) > 0) {
                betamx[species, cov] <- 1
            } else if (sign(as.numeric(quants[2])) < 0) {
                betamx[species, cov] <- -1
            } else {
                betamx[species, cov] <- 0
            }
        }
        credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
      }
    }
  }
  
  if(mode == "X"){
    rho <- schliepOutput$rho[,,burn_in:iter]
    for(species1 in 1:S){
      for(species2 in 1:S){
        varName = paste0("rho",species1,species2)
        mean = mean(rho[species1, species2, ])
        quants = quantile(x = rho[species1, species2, ], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
        if(species1 != species2){
            if(as.numeric(quants[1]) > 0) {
                alphamx[species2,species1] <- 1
            } else if (as.numeric(quants[2]) < 0) {
                alphamx[species2,species1] <- -1
            }
        }
        credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
      }
    }
  }
  
  credibleIntervals <- data.frame(t(data.frame(credibleIntervalsL)))
  names(credibleIntervals) <- c("mean", "lowBound", "highBound")
  credibleIntervals$varName = row.names(credibleIntervals)

  returnL <- list(credibleIntervals=credibleIntervals, betamx = betamx, alphamx = alphamx)
  return(returnL)
}


######## count mistakes ##########
numMistakes <- function(plotPosteriorTableOut, simParms, trial) {
    
    # runL <- rep(NA, times = numRuns * numTrials)
    # trialL <- rep(NA, times = numRuns * numTrials)
    num_correctInferencesL <- rep(NA, times = numRuns * numTrials) # true pos
    num_incorrect_alphaL <- rep(NA, times = numRuns * numTrials) # species interaction incorrectInference
    num_incorrect_betaL <- rep(NA, times = numRuns * numTrials) # covariate interaction incorrectInference
    num_incorrectInferencesL <- rep(NA, times = numRuns * numTrials) # false pos (or wrong direction)
    num_actualEffectsL <- rep(NA, times = numRuns * numTrials) # total actual effects
    num_possibleEffectsL <- rep(NA, times = numRuns * numTrials) # n^2 + n*p -n (minus n because of the diagonal of alpha)
    num_missedEffects_alphaL <- rep(NA, times = numRuns * numTrials)
    num_missedEffects_betaL <- rep(NA, times = numRuns * numTrials)
    num_missedEffectsL <- rep(NA, times = numRuns * numTrials)# false negatives

    # simParms <- readRDS(paste0(dataDir, "/params.Rdata"))
    # figure out which beta column to remove based on just being an intercept variable
    covVars <- simParms$covVars
    const_covs <- unlist(lapply(covVars, FUN = function(covVar) {covVar$type != "constant"}))

    actualBeta <- sign(simParms$beta)
    actualBeta <- actualBeta[, const_covs]
    actualAlpha <- sign(simParms$alpha)

    # check if the inference for this trial finished and store info
    #finished_INLA_tr <- file.exists(paste0(inla_dir, "trial", trial, "/inferenceRes.Rdata"))
    #finished_INLA_trL[i] <- finished_INLA_tr

    #if (finished_INLA_tr) {
    #TODO vvvvvv
    #inferredParms <- list() #readRDS(paste0(inla_dir, "trial", trial, "/inferenceRes.Rdata"))

    betaInferred <- plotPosteriorTableOut$betamx
    alphaInferred <- plotPosteriorTableOut$alphamx
    
    # Note: inferredParms$betaInferred * actualBeta == 1 if and only if both are 1 or both are -1
    num_correct <- sum(betaInferred * actualBeta == 1) + sum(inferredParms$alphaInferred * actualAlpha == 1)

    # Note: here we want to add together inferrence in the wrong direction with saying there is an effect when it's actually 0
    # type 1 error -- inferring an effect where there is none or wrong direction of effect
    count_incorrectT1_beta <- sum(betaInferred * actualBeta == -1) + sum(actualBeta == 0 & betaInferred != 0)
    count_incorrectT1_alpha <- sum(inferredParms$alphaInferred * actualAlpha == -1) + sum(actualAlpha == 0 & inferredParms$alphaInferred != 0)
    count_incorrectT1 <- count_incorrectT1_beta + count_incorrectT1_alpha
    # type 2 error
    # count missed effects (number of times that actual effect is 1 or -1 and inferred effect is 0)
    num_missedEffects_alpha <- sum(actualAlpha != 0 & inferredParms$alphaInferred == 0)
    num_missedEffects_beta <- sum(actualBeta != 0 & betaInferred == 0) 
    num_missedEffects <- num_missedEffects_alpha + num_missedEffects_beta

    # add to running lists
    num_incorrect_alphaL[i] <- count_incorrectT1_alpha
    num_incorrect_betaL[i] <- count_incorrectT1_beta
    num_correctInferencesL[i] <- num_correct
    num_incorrectInferencesL[i] <- count_incorrectT1
    num_missedEffects_alphaL[i] <- num_missedEffects_alpha
    num_missedEffects_betaL[i] <- num_missedEffects_beta
    num_missedEffectsL[i] <- num_missedEffects
    
    # count total number of actual effects in the model
    num_actualEffects <- sum(abs(actualAlpha), abs(actualBeta))

    num_actualEffectsL[i] <- num_actualEffects
    num_possibleEffectsL[i] <- dim(actualBeta)[1] * dim(actualBeta)[2] + dim(actualAlpha)[1]^2 - dim(actualAlpha)[1]

    df <- data.frame(sim_run = runL, 
                trial = trialL, 
                finished_INLA = finished_INLA_trL,
                INLA_runtime = INLA_timeL,
                num_incorrect_alpha = num_incorrect_alphaL,
                num_incorrect_beta = num_incorrect_betaL,
                num_correctInferences = num_correctInferencesL, 
                num_incorrectInferences = num_incorrectInferencesL, 
                num_actualEffects = num_actualEffectsL,
                num_missedEffects_alpha = num_missedEffects_alphaL,
                num_missedEffects_beta = num_missedEffects_betaL,
                num_missedEffectsL = num_missedEffectsL,
                num_possibleEffectsL = num_possibleEffectsL)

    # print(df)

    write.csv(df, paste0(dataDir, "/SchliepRes/trial", trial, "/schliep_mistakes.csv"))
}



testConvergence <- function(schliepRes1, schliepRes2, iter, mode){

    # burn_in = 10000 # this isn't currently doing anything -- autoburnin or something  --look into this
    #iter <- 50000
    #mode <- "X"

    #trial <- 1
    #schliepRes1<-readRDS(paste0(subdir, "trial",trial,"X/","SchliepRes", trial, "Burn", burn, "_t10.Rdata"))
    #res1name="HWA_resultsX1_50k.Rdata"
    #res1name <- args[1]
    #schliepRes1 <- readRDS(paste0("/home/fiona_callahan/HWA_schliep_reanalyze/results/", res1name))

    #trial <- 2
    #res2name <- args[2]
    #res2name="HWA_resultsX2_50k_disp.Rdata"
    #schliepRes2<-readRDS(paste0(subdir, "trial",trial,"X/","SchliepRes", trial, "Burn", burn, "_t10.Rdata"))
    #schliepRes2 <- readRDS(paste0("/home/fiona_callahan/HWA_schliep_reanalyze/results/", res2name))

    #print(res1name)
    #print(res2name)

    S <- schliepRes1$data$S
    p <- schliepRes1$data$p
    t <- schliepRes1$data$t

    variableL <- c()
    gelmanStatL <- c()
    effSizeL <- c()

    # beta
    i <- 1
    for (species in 1:S) {
        for (cov in 0:(p - 1)) {
            varName <- paste0("beta_cov", cov, "_sp", species)
            varValues1 <- schliepRes1$beta[, i]
            varValues2 <- schliepRes2$beta[, i]
            varMCMC <- mcmc.list(mcmc(varValues1), mcmc(varValues2))
            # want this to be around 1 (<1.1 is apparently a good rule of thumb)
            # the two chains should have overdispersed starting points-- need to look into this!!
            # autoburnin causes the first half of the chain to be thrown out -- makes a huge difference
            gelman <- gelman.diag(varMCMC, confidence = 0.95, transform = FALSE, autoburnin = TRUE,
                    multivariate = TRUE)
            effective_size <- effectiveSize(varMCMC)
            #put values into lists
            variableL <- c(variableL, varName)
            gelmanStatL <- c(gelmanStatL, gelman$psrf[1, 1])
            effSizeL <- c(effSizeL, effective_size)
            i <- i + 1
        }
    } ## all of these look fine for 10000 iters/autoburnin with some randomness in starting points EXCEPT the intercepts
    ## why is the estimate of the intercept so different than the others -- did I do something wrong? Is there a reason this should fail to converge

    # add alpha
    for (time in 1:t){
        for (species in 1:S){
            varName <- paste0("alpha_time", time, "_sp", species)
            varValues1 <- schliepRes1$alpha[time, species, ]
            varValues2 <- schliepRes2$alpha[time, species, ]
            varMCMC <- mcmc.list(mcmc(varValues1), mcmc(varValues2))
            # want this to be around 1 (<1.1 is apparently a good rule of thumb)
            # the two chains should have overdispersed starting points-- need to look into this!!
            # autoburnin causes the first half of the chain to be thrown out -- makes a huge difference
            gelman <- gelman.diag(varMCMC, confidence = 0.95, transform = FALSE, autoburnin = TRUE,
                    multivariate = TRUE)
            effective_size <- effectiveSize(varMCMC)
            #values in lists
            variableL <- c(variableL, varName)
            gelmanStatL <- c(gelmanStatL, gelman$psrf[1, 1])
            effSizeL <- c(effSizeL, effective_size)
        }
    }

    # add A
    for (species1 in 1:S){
        for (species2 in 1:S){
            varName <- paste0("A_", species1, species2)
            varValues1 <- schliepRes1$A[species1, species2, ]
            varValues2 <- schliepRes2$A[species1, species2, ]
            varMCMC <- mcmc.list(mcmc(varValues1), mcmc(varValues2))
            gelman <- gelman.diag(varMCMC, confidence = 0.95, transform = FALSE, autoburnin = TRUE,
                    multivariate = TRUE)
            effective_size <- effectiveSize(varMCMC)
            #values in lists
            variableL <- c(variableL, varName)
            gelmanStatL <- c(gelmanStatL, gelman$psrf[1, 1])
            effSizeL <- c(effSizeL, effective_size)
        }
    }


    # add phi 
    for (s1 in 1:S){
        varName <- paste0("phi_", s1)
        varValues1 <- schliepRes1$phi[, s1]
        varValues2 <- schliepRes2$phi[, s1]
        varMCMC <- mcmc.list(mcmc(varValues1), mcmc(varValues2))
        gelman <- gelman.diag(varMCMC, confidence = 0.95, transform = FALSE, autoburnin = TRUE,
                    multivariate = TRUE)
        effective_size <- effectiveSize(varMCMC)
        #values in lists
        variableL <- c(variableL, varName)
        gelmanStatL <- c(gelmanStatL, gelman$psrf[1, 1])
        effSizeL <- c(effSizeL, effective_size)
    }

    # and rho
    if (mode == "X") {
        for (species1 in 1:S){
            for (species2 in 1:S){
                varName <- paste0("rho_", species1, species2)
                varValues1 <- schliepRes1$rho[species1, species2, ]
                varValues2 <- schliepRes2$rho[species1, species2, ]
                varMCMC <- mcmc.list(mcmc(varValues1), mcmc(varValues2))
                gelman <- gelman.diag(varMCMC, confidence = 0.95, transform = FALSE, autoburnin = TRUE,
                            multivariate = TRUE)
                effective_size <- effectiveSize(varMCMC)
                #values in lists
                variableL <- c(variableL, varName)
                gelmanStatL <- c(gelmanStatL, gelman$psrf[1, 1])
                effSizeL <- c(effSizeL, effective_size)
            }
        }
    }

    # make df
    allVar_stats <- data.frame(Variable = variableL, GelmanRubin = gelmanStatL, EffectiveSize = effSizeL)
    
    return(allVar_stats)
}