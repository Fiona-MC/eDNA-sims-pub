#!/usr/bin/env Rscript
# to run
# Rscript count_Specificmistakes.R /space/s1/fiona_callahan/multiSim_ParmSet2/

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input folders need to be supplied", call. = FALSE)
} 

sim_dir <- "/space/s1/fiona_callahan/multiSim_rw3/"
sim_dir <- args[1]
INLAname <- "paper"
#INLAname <- "faster"

numRuns <- 100
numTrials <- 2
nSpecies <- 3
nCov <- 4

runL <- rep(NA, times = numRuns * numTrials)
trialL <- rep(NA, times = numRuns * numTrials)
alphaInferredL <- matrix(nrow = numRuns * numTrials, ncol = nSpecies * nSpecies)
betaInferredL <- matrix(nrow = numRuns * numTrials, ncol = nSpecies * nCov)

i <- 1
for (run in 1:numRuns) {
    for (trial in 1:numTrials) {
        data_dir <- paste0(sim_dir, "randomRun", run, "/")
        inla_dir <- paste0(sim_dir, "randomRun", run, "/INLA_res_", INLAname, "/")
        # store run and trial info
        runL[i] <- run
        trialL[i] <- trial
        if(dir.exists(paste0(sim_dir, "randomRun", run, "/"))) {
          simParms <- readRDS(paste0(data_dir, "params.Rdata"))
          # figure out which beta column to remove based on just being an intercept variable
          covVars <- simParms$covVars
          const_covs <- unlist(lapply(covVars, FUN = function(covVar) {covVar$type != "constant"}))

          actualBeta <- sign(simParms$beta)
          actualBeta <- actualBeta[, const_covs]
          actualAlpha <- sign(simParms$alpha)

          # check if the inference for this trial finished and store info
          finished_INLA_tr <- file.exists(paste0(inla_dir, "trial", trial, "/inferenceRes.Rdata"))
          #finished_INLA_trL[i] <- finished_INLA_tr

          if (finished_INLA_tr) {
            inferredParms <- readRDS(paste0(inla_dir, "trial", trial, "/inferenceRes.Rdata"))

            betaInferred <- inferredParms$betaInferred
            alphaInferred <- inferredParms$alphaInferred
            
            betaInferredL[i, ] <- as.vector(betaInferred)
            alphaInferredL[i, ] <- as.vector(alphaInferred)
          } 
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
                            "beta13", "beta23", "beta33",
                            "beta14", "beta24", "beta34")

df <- data.frame(sim_run = runL, 
            trial = trialL, 
            alphaInferredDF,
            betaInferredDF)

summary(df)

average_alphaInferred <- matrix(data = apply(X = alphaInferredDF[!is.na(alphaInferredDF$alpha11), ], MARGIN = 2, FUN = mean), 
                                nrow = nSpecies, ncol = nSpecies)

average_betaInferred <- matrix(data = apply(X = betaInferredDF[!is.na(betaInferredDF$beta11), ], MARGIN = 2, FUN = mean), 
                                nrow = nSpecies, ncol = nCov)

average_abs_alphaInferred <- matrix(data = apply(X = alphaInferredDF[!is.na(alphaInferredDF$alpha11), ], MARGIN = 2, 
                                FUN = function(x) {mean(abs(x))}), 
                                nrow = nSpecies, ncol = nSpecies)

average_abs_betaInferred <- matrix(data = apply(X = betaInferredDF[!is.na(betaInferredDF$beta11), ], MARGIN = 2, 
                                FUN = function(x) {mean(abs(x))}), 
                                nrow = nSpecies, ncol = nCov)

average_alphaInferred
average_betaInferred

hist(alphaInferredDF$alpha13)

# print(df)
write.csv(df, paste0(sim_dir, "specific_InferredParms.csv"))
