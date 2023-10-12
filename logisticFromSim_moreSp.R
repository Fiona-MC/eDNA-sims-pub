# from sim output (multisim folder with lots of sims), run a logistic model on sitetab_sampled
# outputs a csv of the number of logistic regression inference mistakes per sim

# NOTE: five covariates and three species is hard coded here -- need to change this to generalize

library(stats)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input folder needs to be supplied", call. = FALSE)
} 

data_dir <- "/space/s1/fiona_callahan/multiSim_5sp_random/"
numRuns <- 1000

data_dir <- args[1]
numRuns <- as.numeric(args[2])

# to run
# Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim.R /space/s1/fiona_callahan/multiSim11_testSet/ 200

runs <- 1:numRuns
numTrials <- 1
trials <- 1:1

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

# this is a hacky way to get the number of parms
simParms <- readRDS(paste0(data_dir, "randomRun", 1, "/params.Rdata"))
parmVals <- unlist(simParms)
# take out functions from parmVals
parmVals <- parmVals[lapply(parmVals, FUN = class) != "function"]
parmNames <- names(parmVals)
# initialize place to store parm values
parmsDF <- data.frame(matrix(data = NA, nrow = numRuns * numTrials, ncol = length(parmNames)))
names(parmsDF) <- parmNames

i <- 1
for (run in runs) {
    for (trial in trials) { # basically ignore the trials thing -- I think this is deterministic so trials should be irrelevant
        if(file.exists(paste0(data_dir, "randomRun", run))) { # this is for the runs that were deleted
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
            sim_sitetab_sampled <- read.csv(file = paste0(data_dir, "randomRun", run, "/sim_sitetab_sampled.csv"), header = TRUE)
            
            # sp1
            model_sp1 <- glm(Sp1 ~ 1 + Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Sp2 + Sp3 + Sp4 + Sp5, 
                            family = binomial(link = 'logit'), data = sim_sitetab_sampled)
            model_sp1_summary <- data.frame(summary(model_sp1)$coefficients)
            # sp2
            model_sp2 <- glm(Sp2 ~ 1 + Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Sp1 + Sp3 + Sp4 + Sp5, 
                            family = binomial(link = 'logit'), data = sim_sitetab_sampled)
            model_sp2_summary <- data.frame(summary(model_sp2)$coefficients)
            # sp3
            model_sp3 <- glm(Sp3 ~ 1 + Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Sp1 + Sp2 + Sp4 + Sp5, 
                            family = binomial(link = 'logit'), data = sim_sitetab_sampled)
            model_sp3_summary <- data.frame(summary(model_sp3)$coefficients)
            # sp4
            model_sp4 <- glm(Sp3 ~ 1 + Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Sp1 + Sp2 + Sp3 + Sp5, 
                            family = binomial(link = 'logit'), data = sim_sitetab_sampled)
            model_sp4_summary <- data.frame(summary(model_sp3)$coefficients)
            # sp5
            model_sp5 <- glm(Sp3 ~ 1 + Cov1 + Cov2 + Cov3 + Cov4 + Cov5 + Sp1 + Sp2 + Sp3 + Sp4, 
                            family = binomial(link = 'logit'), data = sim_sitetab_sampled)
            model_sp5_summary <- data.frame(summary(model_sp3)$coefficients)

            # TODO this isn't done -- gotta store these and decide whether they're "right"
            ############### GET INFERRED ALPHA AND BETA ######################
            B11 <- (model_sp1_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov1", "Estimate"])
            B12 <- (model_sp1_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov2", "Estimate"])
            B13 <- (model_sp1_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov3", "Estimate"])
            B14 <- (model_sp1_summary["Cov4", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov4", "Estimate"])
            B15 <- (model_sp1_summary["Cov5", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Cov5", "Estimate"])

            B21 <- (model_sp2_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov1", "Estimate"])
            B22 <- (model_sp2_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov2", "Estimate"])
            B23 <- (model_sp2_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov3", "Estimate"])
            B24 <- (model_sp2_summary["Cov4", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov4", "Estimate"])
            B25 <- (model_sp2_summary["Cov5", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Cov5", "Estimate"])


            B31 <- (model_sp3_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov1", "Estimate"])
            B32 <- (model_sp3_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov2", "Estimate"])
            B33 <- (model_sp3_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov3", "Estimate"])
            B34 <- (model_sp3_summary["Cov4", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov4", "Estimate"])
            B35 <- (model_sp3_summary["Cov5", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Cov5", "Estimate"])

            B41 <- (model_sp4_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Cov1", "Estimate"])
            B42 <- (model_sp4_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Cov2", "Estimate"])
            B43 <- (model_sp4_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Cov3", "Estimate"])
            B44 <- (model_sp4_summary["Cov4", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Cov4", "Estimate"])
            B45 <- (model_sp4_summary["Cov5", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Cov5", "Estimate"])

            B51 <- (model_sp5_summary["Cov1", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Cov1", "Estimate"])
            B52 <- (model_sp5_summary["Cov2", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Cov2", "Estimate"])
            B53 <- (model_sp5_summary["Cov3", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Cov3", "Estimate"])
            B54 <- (model_sp5_summary["Cov4", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Cov4", "Estimate"])
            B55 <- (model_sp5_summary["Cov5", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Cov5", "Estimate"])

            betaInferred <- matrix(c(B11, B12, B13, B14, B15, 
                                    B21, B22, B23, B24, B25, 
                                    B31, B32, B33, B34, B35,
                                    B41, B42, B43, B44, B45,
                                    B51, B52, B53, B54, B55), nrow = 5, ncol = 5, byrow = TRUE)
            
            A12 <- (model_sp1_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp2", "Estimate"])
            A13 <- (model_sp1_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp3", "Estimate"])
            A14 <- (model_sp1_summary["Sp4", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp4", "Estimate"])
            A15 <- (model_sp1_summary["Sp5", "Pr...z.."] < 0.05) * sign(model_sp1_summary["Sp5", "Estimate"])
            
            A21 <- (model_sp2_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp1", "Estimate"])
            A23 <- (model_sp2_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp3", "Estimate"])
            A24 <- (model_sp2_summary["Sp4", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp4", "Estimate"])
            A25 <- (model_sp2_summary["Sp5", "Pr...z.."] < 0.05) * sign(model_sp2_summary["Sp5", "Estimate"])

            A31 <- (model_sp3_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp1", "Estimate"])
            A32 <- (model_sp3_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp2", "Estimate"])
            A34 <- (model_sp3_summary["Sp4", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp4", "Estimate"])
            A35 <- (model_sp3_summary["Sp5", "Pr...z.."] < 0.05) * sign(model_sp3_summary["Sp5", "Estimate"])

            A41 <- (model_sp4_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Sp1", "Estimate"])
            A42 <- (model_sp4_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Sp2", "Estimate"])
            A43 <- (model_sp4_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Sp3", "Estimate"])
            A45 <- (model_sp4_summary["Sp5", "Pr...z.."] < 0.05) * sign(model_sp4_summary["Sp5", "Estimate"])

            A51 <- (model_sp5_summary["Sp1", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Sp1", "Estimate"])
            A52 <- (model_sp5_summary["Sp2", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Sp2", "Estimate"])
            A53 <- (model_sp5_summary["Sp3", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Sp3", "Estimate"])
            A54 <- (model_sp5_summary["Sp4", "Pr...z.."] < 0.05) * sign(model_sp5_summary["Sp4", "Estimate"])
            alphaInferred <- matrix(c(0, A12, A13, A14, A15,
                                    A21, 0, A23, A24, A25,
                                    A31, A32, 0, A34, A35,
                                    A41, A42, A43, 0, A45,
                                    A51, A52, A53, A54, 0), nrow = 5, ncol = 5, byrow = TRUE)



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
