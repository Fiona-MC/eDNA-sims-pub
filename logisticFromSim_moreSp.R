# from sim output (multisim folder with lots of sims), run a logistic model on sitetab_sampled
# outputs a csv of the number of logistic regression inference mistakes per sim

# NOTE: five covariates and five species is hard coded here -- need to change this to generalize

library(stats)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input folder needs to be supplied", call. = FALSE)
} 

data_dir <- "/space/s1/fiona_callahan/multiSim_10sp_dep/"
numRuns <- 20
dumb <- TRUE

#data_dir <- args[1]
#numRuns <- as.numeric(args[2])

# to run
# Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim_moreSp.R /space/s1/fiona_callahan/multiSim_5sp_random/ 1000

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

avg_alphInferred <- matrix(0, nrow = simParms$numSpecies, ncol = simParms$numSpecies)
avg_betInferred <- matrix(0, nrow = simParms$numSpecies, ncol = simParms$numCovs - 1)
nCompleteB <- 0
nCompleteA <- 0

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
            if (dumb) {
                sim_sitetab_sampled <- read.csv(file = paste0(data_dir, "randomRun", run, "/sitetab_dumb.csv"), header = TRUE)
            } else {
                sim_sitetab_sampled <- read.csv(file = paste0(data_dir, "randomRun", run, "/sim_sitetab_sampled.csv"), header = TRUE)
            }
            
            sp_glm_L <- list()
            for (speciesName in simParms$names_species) {
                # glm
                fmla <- as.formula(paste(speciesName, "~ ", 
                                    paste(c(simParms$names_species[simParms$names_species != speciesName]), collapse = "+"), "+",
                                    paste(c(simParms$names_cov), collapse = "+")))
                model <- glm(fmla, family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_summary <- data.frame(summary(model)$coefficients)
                sp_glm_L[[speciesName]] <- model_summary
            }
            
            ############### GET INFERRED ALPHA AND BETA ######################
            betaInferred <- matrix(NA, nrow = simParms$numSpecies, ncol = (simParms$numCovs - 1))
            for (spNum in 1:simParms$numSpecies) {
                for (covNum in 1:(simParms$numCovs - 1)) {
                    speciesName <- paste0("Sp", spNum)
                    covName <- paste0("Cov", covNum)
                    model_summary <- sp_glm_L[[speciesName]]
                    betaInferred[spNum, covNum] <- (model_summary[covName, "Pr...z.."] < 0.05) * sign(model_summary[covName, "Estimate"])
                }
            }
            
            if(!is.na(sum(betaInferred))) {
                avg_betInferred <- avg_betInferred + betaInferred
                nCompleteB <- nCompleteB + 1
            }


            alphaInferred <- matrix(NA, nrow = simParms$numSpecies, ncol = simParms$numSpecies)
            for (spNum1 in 1:simParms$numSpecies) {
                for (spNum2 in 1:simParms$numSpecies) {
                    speciesName1 <- paste0("Sp", spNum1)
                    speciesName2 <- paste0("Sp", spNum2)
                    model_summary <- sp_glm_L[[speciesName1]]
                    if (spNum1 == spNum2) {
                        alphaInferred[spNum1, spNum2] <- 0
                    } else {
                        alphaInferred[spNum1, spNum2] <- (model_summary[speciesName2, "Pr...z.."] < 0.05) * 
                                                        sign(model_summary[speciesName2, "Estimate"])
                    }
                }
            }

            if(!is.na(sum(alphaInferred))) {
                avg_alphInferred <- avg_alphInferred + alphaInferred
                nCompleteA <- nCompleteA + 1
            }

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

# done w loop

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

avg_alphInferred <- avg_alphInferred / nCompleteA
avg_betInferred <- avg_betInferred / nCompleteB

round(avg_alphInferred, 4)
round(abs(avg_alphInferred), 4)

parmsDF$sim_run <- runL
parmsDF$trial <- trialL
#print(parmsDF)

fulldf <- merge(x = df, y = parmsDF, by = c("sim_run", "trial"))

#print(fulldf)
if(dumb) {
    write.csv(fulldf, paste0(data_dir, "/logistic_mistakes_dumb.csv"))
} else {
    write.csv(fulldf, paste0(data_dir, "/logistic_mistakes.csv"))
}
