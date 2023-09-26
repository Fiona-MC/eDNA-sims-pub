# takes folder with the same locations for many simulations and scrambles them so that spatial autocorrelation is eliminated
# samples just a few time points from late in the sims (eg halfway through and at the end) to avoid autocorrelation
library(tidyr)

data_dir <- "/space/s1/fiona_callahan/multiSim_rw3/"
# note rw3 is 100 space pts and 100 sims
#outDir <- paste0(data_dir, "scrambled_4/")
#dir.create(outDir)
numSims <- 100
nSpace <- 100
n_sim_out <- min(numSims, nSpace)
nTime <- 1
numGens <- 1000
abd <- FALSE

# initialize
new_sitetab_sampledL <- lapply(X = 1:n_sim_out, FUN = function(x) {
    df <- data.frame(matrix(nrow = n_sim_out * nTime, ncol = 14))
    names(df) <- c("labID", "Age", "Lat", "Long", "Cov1", "Cov2", "Cov3", "Cov4", "Sp1", "Sp2", "Sp3", "Sp1_abd", "Sp2_abd", "Sp3_abd")
    return(df)})

for(run in 1:n_sim_out) {
    sitetab <- read.csv(paste0(data_dir, "randomRun", run, "/sitetab_longer.csv"))
    outDir <- paste0(data_dir, "randomRun", run, "/")
    sitetab <- sitetab[, -1]
    ages <- seq(1, numGens - 1, round(numGens / nTime)) # note age is gens before present
    sitetab <- sitetab[sitetab$Age %in% ages, ]
    sitetab2 <- pivot_wider(sitetab, names_from = Species, values_from = c(Presence, Abundance))

    # params <- readRDS(paste0(dataDir, "randomRun", run, "/params.Rdata"))
    for(loc in 1:n_sim_out) {
        for(age in ages) {
            # put location in data frame i
            destination_index <- ((run + n_sim_out - loc) %% n_sim_out) + 1
            dest_row <- which(is.na(new_sitetab_sampledL[[destination_index]]))[1]
            new_sitetab_sampledL[[destination_index]][dest_row, ] <- sitetab2[sitetab2$Age == age & sitetab2$labID == loc, ]
        }
    }
}   

# save new sitetabs
for (i in 1:n_sim_out) {
    write.csv(new_sitetab_sampledL[[i]], file = paste0(data_dir, "randomRun", i, "/sitetab_scrambled.csv"))
}





#############################################
# logistic regression and counting mistakes
library(stats)

runs <- 1:n_sim_out
numTrials <- 1
trials <- 1:1

runL <- rep(NA, times = n_sim_out * numTrials)
trialL <- rep(NA, times = n_sim_out * numTrials)
num_correctInferencesL <- rep(NA, times = n_sim_out * numTrials) # true pos
num_incorrect_alphaL <- rep(NA, times = n_sim_out * numTrials) # species interaction incorrectInference
num_incorrect_betaL <- rep(NA, times = n_sim_out * numTrials) # covariate interaction incorrectInference
num_incorrectInferencesL <- rep(NA, times = n_sim_out * numTrials) # false pos (or wrong direction)
num_actualEffectsL <- rep(NA, times = n_sim_out * numTrials) # total actual effects
num_possibleEffectsL <- rep(NA, times = n_sim_out * numTrials) # n^2 + n*p -n (minus n because of the diagonal of alpha)
num_missedEffects_alphaL <- rep(NA, times = n_sim_out * numTrials)
num_missedEffects_betaL <- rep(NA, times = n_sim_out * numTrials)
num_missedEffectsL <- rep(NA, times = n_sim_out * numTrials)# false negatives

# this is a hacky way to get the number of parms
simParms <- readRDS(paste0(data_dir, "randomRun", 1, "/params.Rdata"))
parmVals <- unlist(simParms)
# take out functions from parmVals
parmVals <- parmVals[lapply(parmVals, FUN = class) != "function"]
parmNames <- names(parmVals)

# initialize place to store parm values
parmsDF <- data.frame(matrix(data = NA, nrow = n_sim_out * numTrials, ncol = length(parmNames)))
names(parmsDF) <- parmNames
i <- 1
inferredParmsL <- list()
for (run in runs) {
    for (trial in trials) { # basically ignore the trials thing -- I think this is deterministic so trials should be irrelevant
    outDir <- paste0(data_dir, "randomRun", run, "/")
        if(file.exists(paste0(outDir, "sitetab_scrambled.csv"))) { # this is for the runs that were deleted
            # store run and trial info
            runL[i] <- run
            trialL[i] <- trial

            ############### LOAD ACTUAL PARAMS ######################
            simParms <- readRDS(paste0(data_dir, "randomRun", run, "/params.Rdata"))
            
            parmVals <- unlist(simParms)
            # take out functions from parmVals
            parmVals <- parmVals[lapply(parmVals, FUN = class) != "function"]
            parmNames <- names(parmVals)
            parmsDF[i, ] <- parmVals

            # figure out which beta column to remove based on just being an intercept variable
            covVars <- simParms$covVars
            const_covs <- unlist(lapply(covVars, FUN = function(covVar) {covVar$type != "constant"}))

            actualBeta <- sign(simParms$beta)
            actualBeta <- actualBeta[, const_covs]
            actualAlpha <- sign(simParms$alpha)

            ############### DO LOGISTIC REGRESSION ######################
            # load sitetab
            sim_sitetab_sampled <- read.csv(paste0(outDir, "sitetab_scrambled.csv"), header = TRUE)
            
            if (abd) {
                # sp1
                model_sp1 <- glm(Sp1 ~ 1 + Cov1 + Cov2 + Cov3 + Sp2_abd + Sp3_abd, 
                                family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_sp1_summary <- data.frame(summary(model_sp1)$coefficients)
                # sp2
                model_sp2 <- glm(Sp2 ~ 1 + Cov1 + Cov2 + Cov3 + Sp1_abd + Sp3_abd, 
                                family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_sp2_summary <- data.frame(summary(model_sp2)$coefficients)
                # sp3
                model_sp3 <- glm(Sp3 ~ 1 + Cov1 + Cov2 + Cov3 + Sp1_abd + Sp2_abd, 
                                family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_sp3_summary <- data.frame(summary(model_sp3)$coefficients)

                ############### GET INFERRED ALPHA AND BETA ######################
                B11 <- (model_sp1_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov1", "Estimate"])
                B12 <- (model_sp1_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov2", "Estimate"])
                B13 <- (model_sp1_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov3", "Estimate"])

                B21 <- (model_sp2_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov1", "Estimate"])
                B22 <- (model_sp2_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov2", "Estimate"])
                B23 <- (model_sp2_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov3", "Estimate"])

                B31 <- (model_sp3_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov1", "Estimate"])
                B32 <- (model_sp3_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov2", "Estimate"])
                B33 <- (model_sp3_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov3", "Estimate"])
                
                A12 <- (model_sp1_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp2_abd", "Estimate"])
                A13 <- (model_sp1_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp3_abd", "Estimate"])

                A21 <- (model_sp2_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp1_abd", "Estimate"])
                A23 <- (model_sp2_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp3_abd", "Estimate"])

                A31 <- (model_sp3_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp1_abd", "Estimate"])
                A32 <- (model_sp3_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp2_abd", "Estimate"])
            } else {
                # sp1
                model_sp1 <- glm(Sp1 ~ 1 + Cov1 + Cov2 + Cov3 + Sp2 + Sp3, 
                                family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_sp1_summary <- data.frame(summary(model_sp1)$coefficients)
                # sp2
                model_sp2 <- glm(Sp2 ~ 1 + Cov1 + Cov2 + Cov3 + Sp1 + Sp3, 
                                family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_sp2_summary <- data.frame(summary(model_sp2)$coefficients)
                # sp3
                model_sp3 <- glm(Sp3 ~ 1 + Cov1 + Cov2 + Cov3 + Sp1 + Sp2, 
                                family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_sp3_summary <- data.frame(summary(model_sp3)$coefficients)

                ############### GET INFERRED ALPHA AND BETA ######################
                B11 <- (model_sp1_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov1", "Estimate"])
                B12 <- (model_sp1_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov2", "Estimate"])
                B13 <- (model_sp1_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov3", "Estimate"])

                B21 <- (model_sp2_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov1", "Estimate"])
                B22 <- (model_sp2_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov2", "Estimate"])
                B23 <- (model_sp2_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov3", "Estimate"])

                B31 <- (model_sp3_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov1", "Estimate"])
                B32 <- (model_sp3_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov2", "Estimate"])
                B33 <- (model_sp3_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov3", "Estimate"])
                
                A12 <- (model_sp1_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp2", "Estimate"])
                A13 <- (model_sp1_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp3", "Estimate"])

                A21 <- (model_sp2_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp1", "Estimate"])
                A23 <- (model_sp2_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp3", "Estimate"])

                A31 <- (model_sp3_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp1", "Estimate"])
                A32 <- (model_sp3_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp2", "Estimate"])
            }

            betaInferred <- matrix(c(B11, B12, B13, B21, B22, B23, B31, B32, B33), nrow = 3, ncol = 3, byrow = TRUE)
            alphaInferred <- matrix(c(0, A12, A13, A21, 0, A23, A31, A32, 0), nrow = 3, ncol = 3, byrow = TRUE)

            inferredParmsL[[i]] <- list(alphaInferred = alphaInferred, betaInferred = betaInferred)

            ############### COMPARE INFERRECE RESULTS TO ACTUAL ALPHA AND BETA ######################

            # Note: inferredParms$betaInferred * actualBeta == 1 if and only if both are 1 or both are -1
            num_correct <- sum(betaInferred * actualBeta == 1) + sum(alphaInferred * actualAlpha == 1)

            # Note: here we want to add together inferrence in the wrong direction with saying there is an effect when it's actually 0
            # type 1 error -- inferring an effect where there is none or wrong direction of effect
            count_incorrectT1_beta <- sum(betaInferred * actualBeta == -1) + sum(actualBeta == 0 & betaInferred != 0)
            count_incorrectT1_alpha <- sum(alphaInferred * actualAlpha == -1) + sum(actualAlpha == 0 & alphaInferred != 0)
            count_incorrectT1 <- count_incorrectT1_beta + count_incorrectT1_alpha
            # type 2 error
            # count missed effects (number of times that actual effect is 1 or -1 and inferred effect is 0)
            num_missedEffects_alpha <- sum(actualAlpha != 0 & alphaInferred == 0)
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
        }
        i <- i + 1
    }
}

df <- data.frame(sim_run = runL, 
            trial = trialL, 
            num_incorrect_alpha = num_incorrect_alphaL,
            num_incorrect_beta = num_incorrect_betaL,
            num_correctInferences = num_correctInferencesL, 
            num_incorrectInferences = num_incorrectInferencesL, 
            num_actualEffects = num_actualEffectsL,
            num_missedEffects_alpha = num_missedEffects_alphaL,
            num_missedEffects_beta = num_missedEffects_betaL,
            num_missedEffectsL = num_missedEffectsL,
            num_possibleEffectsL = num_possibleEffectsL)

#print(df)

parmsDF$sim_run <- runL
parmsDF$trial <- trialL
#print(parmsDF)

fulldf <- merge(x = df, y = parmsDF, by = c("sim_run", "trial"))

#print(fulldf)
write.csv(fulldf, paste0(data_dir, "/logistic_mistakes.csv"))









##############################################
# specific parms
nSpecies <- 3
nCov <- 4

runL <- rep(NA, times = n_sim_out * numTrials)
trialL <- rep(NA, times = n_sim_out * numTrials)
alphaInferredL <- matrix(nrow = n_sim_out * numTrials, ncol = nSpecies * nSpecies)
betaInferredL <- matrix(nrow = n_sim_out * numTrials, ncol = nSpecies * (nCov - 1))

i <- 1
simParms <- readRDS(paste0(data_dir, "randomRun", 1, "/params.Rdata"))
for (run in 1:n_sim_out) {
    for (trial in 1:numTrials) {
        # store run and trial info
        runL[i] <- run
        trialL[i] <- trial
    
        # figure out which beta column to remove based on just being an intercept variable
        covVars <- simParms$covVars
        const_covs <- unlist(lapply(covVars, FUN = function(covVar) {covVar$type != "constant"}))

        actualBeta <- sign(simParms$beta)
        actualBeta <- actualBeta[, const_covs]
        actualAlpha <- sign(simParms$alpha)
 
        inferredParms <- inferredParmsL[[i]]

        betaInferred <- inferredParms$betaInferred
        alphaInferred <- inferredParms$alphaInferred
        
        betaInferredL[i, ] <- as.vector(betaInferred)
        alphaInferredL[i, ] <- as.vector(alphaInferred)
        
        i <- i + 1
    }
}

alphaInferredDF <- as.data.frame(alphaInferredL)
names(alphaInferredDF) <- c("alpha11", "alpha21", "alpha31", 
                            "alpha12", "alpha22", "alpha32",
                            "alpha13", "alpha23", "alpha33")
betaInferredDF <- as.data.frame(betaInferredL)
names(betaInferredDF) <- c("beta11", "beta21", "beta31", 
                            "beta12", "beta22", "beta32",
                            "beta13", "beta23", "beta33")

df <- data.frame(sim_run = runL, 
            trial = trialL, 
            alphaInferredDF,
            betaInferredDF)

summary(df)

average_alphaInferred <- matrix(data = apply(X = alphaInferredDF[!is.na(alphaInferredDF$alpha11), ], MARGIN = 2, FUN = mean), 
                                nrow = nSpecies, ncol = nSpecies)

average_betaInferred <- matrix(data = apply(X = betaInferredDF[!is.na(betaInferredDF$beta11), ], MARGIN = 2, FUN = mean), 
                                nrow = nSpecies, ncol = nCov - 1)

average_abs_alphaInferred <- matrix(data = apply(X = alphaInferredDF[!is.na(alphaInferredDF$alpha11), ], MARGIN = 2, 
                                FUN = function(x) {mean(abs(x))}), 
                                nrow = nSpecies, ncol = nSpecies)

average_abs_betaInferred <- matrix(data = apply(X = betaInferredDF[!is.na(betaInferredDF$beta11), ], MARGIN = 2, 
                                FUN = function(x) {mean(abs(x))}), 
                                nrow = nSpecies, ncol = nCov - 1)

average_alphaInferred
average_betaInferred

average_abs_alphaInferred
average_abs_betaInferred



##############################################
# specific parms for INLA
nSpecies <- 3
nCov <- 4
numTrials <- 2

runL <- rep(NA, times = n_sim_out * numTrials)
trialL <- rep(NA, times = n_sim_out * numTrials)
alphaInferredL <- matrix(nrow = n_sim_out * numTrials, ncol = nSpecies * nSpecies)
betaInferredL <- matrix(nrow = n_sim_out * numTrials, ncol = nSpecies * (nCov - 1))

i <- 1
simParms <- readRDS(paste0(data_dir, "randomRun", 1, "/params.Rdata"))
for (run in 1:n_sim_out) {
    for (trial in 1:numTrials) {
        # store run and trial info
        runL[i] <- run
        trialL[i] <- trial
    
        # figure out which beta column to remove based on just being an intercept variable
        covVars <- simParms$covVars
        const_covs <- unlist(lapply(covVars, FUN = function(covVar) {covVar$type != "constant"}))

        actualBeta <- sign(simParms$beta)
        actualBeta <- actualBeta[, const_covs]
        actualAlpha <- sign(simParms$alpha)

        if(file.exists(paste0(data_dir, "randomRun", run, "/INLA_res_paper/trial", trial, "/inferenceRes.Rdata"))){
            inferredParms <- readRDS(paste0(data_dir, "randomRun", run, "/INLA_res_paper/trial", trial, "/inferenceRes.Rdata"))

            betaInferred <- inferredParms$betaInferred
            betaInferred <- betaInferred[, const_covs]

            alphaInferred <- inferredParms$alphaInferred
            
            betaInferredL[i, ] <- as.vector(betaInferred)
            alphaInferredL[i, ] <- as.vector(alphaInferred)
        } else {
            betaInferredL[i, ] <- NA
            alphaInferredL[i, ] <- NA
        }
        i <- i + 1
    }
}

alphaInferredDF <- as.data.frame(alphaInferredL)
names(alphaInferredDF) <- c("alpha11", "alpha21", "alpha31", 
                            "alpha12", "alpha22", "alpha32",
                            "alpha13", "alpha23", "alpha33")
betaInferredDF <- as.data.frame(betaInferredL)
names(betaInferredDF) <- c("beta11", "beta21", "beta31", 
                            "beta12", "beta22", "beta32",
                            "beta13", "beta23", "beta33")

df <- data.frame(sim_run = runL, 
            trial = trialL, 
            alphaInferredDF,
            betaInferredDF)

summary(df)

average_alphaInferred <- matrix(data = apply(X = alphaInferredDF[!is.na(alphaInferredDF$alpha11), ], MARGIN = 2, FUN = mean), 
                                nrow = nSpecies, ncol = nSpecies)

average_betaInferred <- matrix(data = apply(X = betaInferredDF[!is.na(betaInferredDF$beta11), ], MARGIN = 2, FUN = mean), 
                                nrow = nSpecies, ncol = nCov - 1)

average_abs_alphaInferred <- matrix(data = apply(X = alphaInferredDF[!is.na(alphaInferredDF$alpha11), ], MARGIN = 2, 
                                FUN = function(x) {mean(abs(x))}), 
                                nrow = nSpecies, ncol = nSpecies)

average_abs_betaInferred <- matrix(data = apply(X = betaInferredDF[!is.na(betaInferredDF$beta11), ], MARGIN = 2, 
                                FUN = function(x) {mean(abs(x))}), 
                                nrow = nSpecies, ncol = nCov - 1)

round(average_alphaInferred, digits = 2)
round(average_betaInferred, digits = 2)

round(average_abs_alphaInferred, digits = 2)
round(average_abs_betaInferred, digits = 2)
