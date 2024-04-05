# this is edited to work for filtered species names
# should actually work for both though

library(stats)
library(igraph)
library(stringr)

data_dir <- "/space/s1/fiona_callahan/multiSim_50sp/"
numRuns <- 100
covs <- FALSE
logi <- (as.numeric("1") == 0)
sitetab_name <- "sim_sitetab_sampled1000_filtered100.csv"
outName <- "logistic_mistakes_sampled1000_noCov_100runs_filtered100"
#sitetab_name <- "sim_sitetab_sampled1000.csv"
#outName <- "logistic_mistakes_sampled1000_noCov_100runs"

randomSim <- str_detect(data_dir, "random")

# to run
# Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim_moreSp.R /space/s1/fiona_callahan/multiSim_5sp_random/ 1000

runs <- 1:numRuns
numTrials <- 1
trials <- 1:1

cutoffs <- c(0, 1e-128, 1e-64, 1e-32, 1e-16, 1e-8, 1e-4, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15,
             0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.9999, 0.99999, 0.999999, 0.9999999, 0.99999999, 1)
#cutoffs <- c(1e-128)

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
timeL <- rep(NA, times = numRuns * numTrials)
finished_trL <- rep(NA, times = numRuns * numTrials) # T or F whether this trial finished the INLA part
num_correct_clusterL <- rep(NA, times = numRuns * numTrials)
num_incorrect_clusterL <- rep(NA, times = numRuns * numTrials)
num_missed_clusterL <- rep(NA, times = numRuns * numTrials)
alpha_direction_mistakesL <- rep(NA, times = numRuns * numTrials)
alpha_incorrect_undirectedL <- rep(NA, times = numRuns * numTrials)
alpha_correct_undirectedL <- rep(NA, times = numRuns * numTrials)
TP_ignoreSign <- rep(NA, times = numRuns * numTrials)
FP_ignoreSign <- rep(NA, times = numRuns * numTrials)
TN_ignoreSign <- rep(NA, times = numRuns * numTrials)
FN_ignoreSign <- rep(NA, times = numRuns * numTrials)
TP_cluster <- rep(NA, times = numRuns * numTrials)
FP_cluster <- rep(NA, times = numRuns * numTrials)
TN_cluster <- rep(NA, times = numRuns * numTrials)
FN_cluster <- rep(NA, times = numRuns * numTrials)

# this is a hacky way to get the number of parms
simParms <- readRDS(paste0(data_dir, "randomRun", 1, "/params.Rdata"))

logisticRes <- list()
for (run in runs) {
    for (trial in trials) { # basically ignore the trials thing -- I think this is deterministic so trials should be irrelevant
        if (file.exists(paste0(data_dir, "randomRun", run))) { # this is for the runs that were deleted
            # LOAD ACTUAL PARAMS 
            simParms <- readRDS(paste0(data_dir, "randomRun", run, "/params.Rdata"))

            ############### DO LOGISTIC REGRESSION ######################
            # load sitetab
            sim_sitetab_sampled <- read.csv(file = paste0(data_dir, "randomRun", run, "/", sitetab_name), header = TRUE)
            
            sp_glm_L <- list()
            speciesNames <- names(sim_sitetab_sampled)[grep("Sp", names(sim_sitetab_sampled))]
            numSpecies <- length(speciesNames)
            for (speciesName in speciesNames) {
                if (covs) {
                    # glm
                    fmla <- as.formula(paste(speciesName, "~ ", 
                                        paste(c(speciesNames[speciesNames != speciesName]), collapse = "+"), "+",
                                        paste(c(simParms$names_cov), collapse = "+")))
                } else {
                    fmla <- as.formula(paste(speciesName, "~ ", 
                                    paste(c(speciesNames[speciesNames != speciesName]), collapse = "+")))
                }
                model <- glm(fmla, family = binomial(link = 'logit'), data = sim_sitetab_sampled)
                model_summary <- data.frame(summary(model)$coefficients)
                sp_glm_L[[speciesName]] <- model_summary
            }
            logisticRes[[run]] <- sp_glm_L
        }
    }
}

z.values_pos <- c()
z.values_neg <- c()
z.values_0 <- c()


for (cutoff in cutoffs) {
    i <- 1

    # TODO get df of actual and predicted values to test if I am making the confusion mx wrong
    actual_predicted_vals <- data.frame()
    alpha_actualL <- c()
    alpha_predictedL <- c()
    beta_actualL <- c()
    beta_predictedL <- c()

    # initialize place to store parm values
    parmsDF <- data.frame()

    for (run in runs) {
        for (trial in trials) { # basically ignore the trials thing -- I think this is deterministic so trials should be irrelevant
            if (file.exists(paste0(data_dir, "randomRun", run))) { # this is for the runs that were deleted
                sim_sitetab_sampled <- read.csv(file = paste0(data_dir, "randomRun", run, "/", sitetab_name), header = TRUE)
                speciesNames <- names(sim_sitetab_sampled)[grep("Sp", names(sim_sitetab_sampled))]
                numSpecies <- length(speciesNames)
                sp_glm_L <- logisticRes[[run]]
                # store run and trial info
                runL[i] <- run
                trialL[i] <- trial

                ############### LOAD ACTUAL PARAMS ######################
                simParms <- readRDS(paste0(data_dir, "randomRun", run, "/params.Rdata"))
                
                parmVals <- unlist(simParms)
                # take out functions from parmVals
                parmVals <- parmVals[lapply(parmVals, FUN = class) != "function"]
                # take out parms that are not in all of the runs
                if (run != runs[1]) {
                    parmsDF <- parmsDF[, !duplicated(names(parmsDF))]
                    parmVals <- parmVals[!duplicated(names(parmVals))]
                    parmsDF <- parmsDF[, names(parmsDF) %in% names(parmVals)]
                    parmVals <- parmVals[names(parmVals) %in% names(parmsDF)]
                    if (!(all(names(parmsDF) == names(parmVals)))) {
                        stop("rbind in logisticFromSim_moreSp is failing because the columns don't match")
                    }
                }
                parmsDF <- rbind(parmsDF, parmVals)
                colnames(parmsDF) <- names(parmVals)

                # figure out which beta column to remove based on just being an intercept variable
                covVars <- simParms$covVars
                const_covs <- unlist(lapply(covVars, FUN = function(covVar) {covVar$type != "constant"}))

                # Extract numeric parts using regular expressions
                numeric_parts <- gsub("\\D", "", speciesNames)
                # Convert to numeric
                numeric_species <- as.numeric(numeric_parts)

                actualBeta <- sign(simParms$beta)
                actualBeta <- actualBeta[, const_covs]
                actualBeta <- actualBeta[numeric_species, ]

                actualAlpha <- sign(simParms$alpha)
                actualAlpha <- actualAlpha[numeric_species, numeric_species]

                #simParmsFiltered <- readRDS(paste0(data_dir, "randomRun", run, "/paramsFiltered.Rdata"))
                #actualBeta <- simParmsFiltered$filteredBeta
                #actualAlpha <- simParmsFiltered$filteredAlpha

                ############### GET INFERRED ALPHA AND BETA ######################
                betaInferred <- matrix(NA, nrow = numSpecies, ncol = (simParms$numCovs - 1))
                if (covs) {
                    for (spNum in 1:numSpecies) {
                        for (covNum in 1:(simParms$numCovs - 1)) {
                            speciesName <- speciesNames[spNum]
                            covName <- paste0("Cov", covNum)
                            model_summary <- sp_glm_L[[speciesName]]
                            betaInferred[spNum, covNum] <- (model_summary[covName, "Pr...z.."] < cutoff) * sign(model_summary[covName, "Estimate"])
                            if (cutoff == 0) {
                                if (actualBeta[spNum, covNum] == 0) {
                                    z.values_0 <- c(z.values_0, model_summary[covName, "z.value"])
                                } else if (actualBeta[spNum, covNum] < 0) {
                                    z.values_neg <- c(z.values_neg, model_summary[covName, "z.value"])
                                } else {
                                    z.values_pos <- c(z.values_pos, model_summary[covName, "z.value"])
                                }
                            }
                        }
                    }
                    if(!is.na(sum(betaInferred)) && !randomSim) {
                        avg_betInferred <- avg_betInferred + betaInferred
                        nCompleteB <- nCompleteB + 1
                    }   
                }

                alphaInferred <- matrix(NA, nrow = numSpecies, ncol = numSpecies)
                for (spNum1 in 1:numSpecies) {
                    for (spNum2 in 1:numSpecies) {
                        speciesName1 <- speciesNames[spNum1]
                        speciesName2 <- speciesNames[spNum2]
                        model_summary <- sp_glm_L[[speciesName1]]
                        if (spNum1 == spNum2) {
                            alphaInferred[spNum1, spNum2] <- 0
                        } else {
                            alphaInferred[spNum1, spNum2] <- (model_summary[speciesName2, "Pr...z.."] < cutoff) * 
                                                            sign(model_summary[speciesName2, "Estimate"])
                            if (cutoff == 0) {
                                if (actualAlpha[spNum1, spNum2] == 0) {
                                    z.values_0 <- c(z.values_0, model_summary[speciesName2, "z.value"])
                                } else if (actualAlpha[spNum1, spNum2] < 0) {
                                    z.values_neg <- c(z.values_neg, model_summary[speciesName2, "z.value"])
                                } else {
                                    z.values_pos <- c(z.values_pos, model_summary[speciesName2, "z.value"])
                                }
                            }
                        }
                    }
                }


                ############### COMPARE INFERRECE RESULTS TO ACTUAL ALPHA AND BETA ######################
                inferredParms <- list()
                inferredParms$alphaInferred <- alphaInferred
                inferredParms$betaInferred <- betaInferred
                inferredParms$cutoff <- cutoff

                alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
                inferredAlphaG <- graph_from_adjacency_matrix(inferredParms$alphaInferred != 0, mode = "max")
                connected_alpha_actual <- (distances(alphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                            (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
                connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                            (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
                # Note: inferredParms$betaInferred * actualBeta == 1 if and only if both are 1 or both are -1
                num_correct_alpha <- sum(inferredParms$alphaInferred * actualAlpha == 1, na.rm = TRUE)
                if (covs) {
                    num_correct_beta <- sum(inferredParms$betaInferred * actualBeta == 1, na.rm = TRUE)
                    num_correct <- num_correct_alpha + num_correct_beta
                } else {
                    num_correct <- num_correct_alpha
                }

                num_correct_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == 1, na.rm = TRUE)

                # Note: here we want to add together inferrence in the wrong direction with saying there is an effect when it's actually 0
                # type 1 error -- inferring an effect where there is none or wrong direction of effect
                count_incorrectT1_alpha <- sum(inferredParms$alphaInferred * actualAlpha == -1, na.rm = TRUE) + 
                                                sum(actualAlpha == 0 & inferredParms$alphaInferred != 0, na.rm = TRUE)
                alpha_direction_mistakes <- sum(inferredParms$alphaInferred * actualAlpha == -1, na.rm = TRUE)

                if (covs) {
                    count_incorrectT1_beta <- sum(inferredParms$betaInferred * actualBeta == -1, na.rm = TRUE) + 
                                                sum(actualBeta == 0 & inferredParms$betaInferred != 0, na.rm = TRUE)
                    count_incorrectT1 <- count_incorrectT1_beta + count_incorrectT1_alpha
                } else {
                    count_incorrectT1 <- count_incorrectT1_alpha
                }

                count_incorrect_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == -1, na.rm = TRUE) + 
                                                sum(connected_alpha_actual == 0 & connected_alpha_inferred != 0, na.rm = TRUE)

                if (covs) {
                TP_cluster[i] <- sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1, na.rm = TRUE) + 
                                    sum(abs(betaInferred) * abs(actualBeta) == 1, na.rm = TRUE)
                FP_cluster[i] <- sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0), na.rm = TRUE) + 
                                    sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0), na.rm = TRUE)
                FN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1), na.rm = TRUE) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1), na.rm = TRUE)
                TN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0), na.rm = TRUE) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0), na.rm = TRUE)
                } else {
                    TP_cluster[i] <- sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1, na.rm = TRUE) 
                    FP_cluster[i] <- sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0), na.rm = TRUE) 
                    FN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1), na.rm = TRUE) 
                    TN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0), na.rm = TRUE) 
                }
                
                # type 2 error
                # count missed effects (number of times that actual effect is 1 or -1 and inferred effect is 0)
                num_missedEffects_alpha <- sum(actualAlpha != 0 & inferredParms$alphaInferred == 0, na.rm = TRUE)

                if (covs) {
                    num_missedEffects_beta <- sum(actualBeta != 0 & inferredParms$betaInferred == 0, na.rm = TRUE)
                }

                if (covs) {
                    num_missedEffects <- num_missedEffects_alpha + num_missedEffects_beta
                } else {
                    num_missedEffects <- num_missedEffects_alpha 
                }

                num_missedEffects_cluster <- sum(connected_alpha_actual != 0 & connected_alpha_inferred == 0, na.rm = TRUE)

                # undirected meaning that a-->b iff b-->a in "actual alpha". This also ignores the sign.
                undirected_alpha_actual <- distances(alphaG, v = 1:numSpecies, to = 1:numSpecies) == 1
                undirected_alpha_inferred <- distances(inferredAlphaG, v = 1:numSpecies, to = 1:numSpecies) == 1

                alpha_incorrect_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == -1, na.rm = TRUE) + 
                                                sum(undirected_alpha_actual == 0 & undirected_alpha_inferred != 0, na.rm = TRUE)
                alpha_correct_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == 1, na.rm = TRUE)

                # add to running lists
                if (covs) {
                    num_incorrect_betaL[i] <- count_incorrectT1_beta
                    num_missedEffects_betaL[i] <- num_missedEffects_beta
                } 

                if (covs) {
                    TP_ignoreSign[i] <- sum(abs(alphaInferred) * abs(actualAlpha) == 1, na.rm = TRUE) + 
                                        sum(abs(betaInferred) * abs(actualBeta) == 1, na.rm = TRUE)
                    FP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0), na.rm = TRUE) + 
                                        sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0), na.rm = TRUE)
                    FN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1), na.rm = TRUE) + 
                                        sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1), na.rm = TRUE)
                    TN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0), na.rm = TRUE) + 
                                        sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0), na.rm = TRUE)
                } else {
                    TP_ignoreSign[i] <- sum(abs(alphaInferred) * abs(actualAlpha) == 1, na.rm = TRUE) 
                    FP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0), na.rm = TRUE) 
                    FN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1), na.rm = TRUE) 
                    TN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0), na.rm = TRUE) 
                }

                num_incorrect_alphaL[i] <- count_incorrectT1_alpha
                num_correctInferencesL[i] <- num_correct
                num_incorrectInferencesL[i] <- count_incorrectT1
                num_missedEffects_alphaL[i] <- num_missedEffects_alpha
                num_missedEffectsL[i] <- num_missedEffects
                alpha_direction_mistakesL[i] <- alpha_direction_mistakes
                num_correct_clusterL[i] <- num_correct_cluster
                num_incorrect_clusterL[i] <- count_incorrect_cluster
                num_missed_clusterL[i] <- num_missedEffects_cluster
                alpha_incorrect_undirectedL[i] <- alpha_incorrect_undirected
                alpha_correct_undirectedL[i] <- alpha_correct_undirected

                if ("time" %in% names(inferredParms)) {
                    timeL[i] <- inferredParms$time
                }
                    
                
            }

                # count total number of actual effects in the model
                num_actualEffects <- sum(abs(actualAlpha))

                num_actualEffectsL[i] <- num_actualEffects
                num_possibleEffectsL[i] <- dim(actualAlpha)[1]^2 - dim(actualAlpha)[1]
                i <- i + 1
        }
    }

    # done w loop

    df <- data.frame(sim_run = runL, 
                trial = trialL, 
                finished = finished_trL,
                runtime = timeL,
                num_incorrect_alpha = num_incorrect_alphaL,
                num_incorrect_beta = num_incorrect_betaL,
                num_correctInferences = num_correctInferencesL, 
                num_incorrectInferences = num_incorrectInferencesL, 
                num_actualEffects = num_actualEffectsL,
                num_missedEffects_alpha = num_missedEffects_alphaL,
                num_missedEffects_beta = num_missedEffects_betaL,
                num_missedEffectsL = num_missedEffectsL,
                num_possibleEffectsL = num_possibleEffectsL,
                num_correct_cluster = num_correct_clusterL,
                num_incorrect_cluster = num_incorrect_clusterL,
                num_missed_cluster = num_missed_clusterL,
                alpha_direction_mistakes = alpha_direction_mistakesL,
                alpha_incorrect_undirected = alpha_incorrect_undirectedL,
                alpha_correct_undirected = alpha_correct_undirectedL,
                TP_ignoreSign = TP_ignoreSign,
                FP_ignoreSign = FP_ignoreSign,
                TN_ignoreSign = TN_ignoreSign,
                FN_ignoreSign = FN_ignoreSign,
                TP_cluster = TP_cluster,
                FP_cluster = FP_cluster,
                TN_cluster = TN_cluster,
                FN_cluster = FN_cluster)

    #print(df)

    parmsDF$sim_run <- runL
    parmsDF$trial <- trialL
    #print(parmsDF)

    fulldf <- merge(x = df, y = parmsDF, by = c("sim_run", "trial"))

    #print(fulldf)
    write.csv(fulldf, paste0(data_dir, "/", outName, "_cutoff", cutoff, "_", numRuns, "sims.csv"))
     
}

plotZvals <- TRUE
if (plotZvals) {
    library(ggplot2)
    z.values_0 <- z.values_0[!is.na(z.values_0)]
    z.values_pos <- z.values_pos[!is.na(z.values_pos)]
    z.values_neg <- z.values_neg[!is.na(z.values_neg)]
    z_vals_df <- data.frame(actual_beta = as.factor(c(rep(0, times = length(z.values_0)), 
                                        rep(1, times = length(z.values_pos)), 
                                        rep(-1, times = length(z.values_neg)))), 
                        zValues = c(z.values_0, z.values_pos, z.values_neg))

    plot1 <- ggplot(z_vals_df, aes(x = zValues, fill = actual_beta, color = actual_beta)) +
        geom_density(alpha = 0.5) +  # Add transparency for better visibility of overlapping densities
        coord_cartesian(xlim = c(-4, 4)) +  # Set x-axis limits
        labs(title = paste0("Density Plot of zValues: ", outName),
            x = "zValues",
            y = "Density") +
        scale_fill_manual(values = rainbow(length(unique(z_vals_df$actual_beta)))) +  # Adjust colors
        scale_color_manual(values = rainbow(length(unique(z_vals_df$actual_beta)))) +  # Adjust colors
        theme_minimal()  # Optional: Change the theme
    
    ggsave(filename = paste0(data_dir, "zValsDistribution_", outName, "_", numRuns, "sims.png"), plot = plot1)

}