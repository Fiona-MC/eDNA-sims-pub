#!/usr/bin/env Rscript

# to run
# Rscript count_INLAmistakes.R /space/s1/fiona_callahan/multiSim2/randomRun1/ /space/s1/fiona_callahan/multiSim2/randomRun1/INLA_res_faster/

#data_dir <- "/space/s1/fiona_callahan/multiSim2/randomRun1/"
#inla_dir <- "/space/s1/fiona_callahan/multiSim2/randomRun1/INLA_res_faster/"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input folders need to be supplied", call. = FALSE)
} 

# TODO! fix this so that it won't count the constant covariate
# but be careful because this is currently running and I don't wanna ruin it

#data_dir <- "/home/fiona_callahan/simData/testing/run1/"
#save_dir <- "/home/fiona_callahan/simData/testing/INLAres/"
data_dir <- args[1]
inla_dir <- args[2]

numRuns <- 1 # runs per folder (1 is correct)
numTrials <- 1

runL <- rep(NA, times = numRuns * numTrials)
trialL <- rep(NA, times = numRuns * numTrials)
num_correctInferencesL <- rep(NA, times = numRuns * numTrials) # true pos
num_incorrect_alphaL <- rep(NA, times = numRuns * numTrials) # species interaction incorrectInference
num_incorrect_betaL <- rep(NA, times = numRuns * numTrials) # covariate interaction incorrectInference
num_incorrectInferencesL <- rep(NA, times = numRuns * numTrials) # false pos (or wrong direction)
num_actualEffectsL <- rep(NA, times = numRuns * numTrials) # total actual effects
num_possibleEffectsL <- rep(NA, times = numRuns * numTrials) # n^2 + n*p -n (minus n because of the diagonal of alpha)
num_missedEffects_alphaL <- rep(NA, times = numRuns * numTrials)
num_missedEffects_betaL <- rep(NA, times = numRuns * numTrials)
num_missedEffectsL <- rep(NA, times = numRuns * numTrials)# false negatives
INLA_timeL <- rep(NA, times = numRuns * numTrials)
finished_INLA_trL <- rep(NA, times = numRuns * numTrials) # T or F whether this trial finished the INLA part

i <- 1
for (run in 1:numRuns) {
    for (trial in 1:numTrials) {
        # store run and trial info
        runL[i] <- run
        trialL[i] <- trial

        simParms <- readRDS(paste0(data_dir, "/params.Rdata"))
        # figure out which beta column to remove based on just being an intercept variable
        covVars <- simParms$covVars
        const_covs <- unlist(lapply(covVars, FUN = function(covVar) {covVar$type != "constant"}))

        actualBeta <- sign(simParms$beta)
        actualBeta <- actualBeta[, const_covs]
        actualAlpha <- sign(simParms$alpha)

        # check if the inference for this trial finished and store info
        finished_INLA_tr <- file.exists(paste0(inla_dir, "trial", trial, "/inferenceRes.Rdata"))
        finished_INLA_trL[i] <- finished_INLA_tr

        if (finished_INLA_tr) {
          inferredParms <- readRDS(paste0(inla_dir, "trial", trial, "/inferenceRes.Rdata"))

          betaInferred <- inferredParms$betaInferred
          betaInferred <- betaInferred[, const_covs]
        
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
          # note this "if" is just because I changed this during a run -- can be taken out later (or not)
          if ("time" %in% names(inferredParms)) {
            INLA_timeL[i] <- inferredParms$time
          } 
        } 
        
        # count total number of actual effects in the model
        num_actualEffects <- sum(abs(actualAlpha), abs(actualBeta))

        num_actualEffectsL[i] <- num_actualEffects
        num_possibleEffectsL[i] <- dim(actualBeta)[1] * dim(actualBeta)[2] + dim(actualAlpha)[1]^2 - dim(actualAlpha)[1]

        i <- i + 1
    }
}

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

write.csv(df, paste0(inla_dir, "INLA_mistakes.csv"))
