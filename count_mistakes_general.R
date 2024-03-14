#!/usr/bin/env Rscript
library(igraph)
# to run
# Rscript countEcoCopulaMistakes.R /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun1/ /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun1/ecoCopula_res/

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input folders need to be supplied", call. = FALSE)
} 

#data_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun1/"
#outdir <- "/space/s1/fiona_callahan/multiSim_100/randomRun1/ecoCopula_res_cov500/"
#outdir <- "/space/s1/fiona_callahan/multiSim_100/randomRun2/INLA_res_faster/"
#covs <- TRUE # is beta inferred?

#echo ${folder}/ ${folder}/${resDirName}/ ${covs} ${cutoff}
#/space/s1/fiona_callahan/multiSim_100/randomRun1/ /space/s1/fiona_callahan/multiSim_100/randomRun1/ecoCopula_res_cov500/ 1
#data_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun1/"
#outdir <- "/space/s1/fiona_callahan/multiSim_100/randomRun1/ecoCopula_res_cov500/"
#covs <- (as.numeric(0) == 1)
#cutoff <- as.numeric(NA)

data_dir <- args[1]
outdir <- args[2]
covs <- (as.numeric(args[3]) == 1)
cutoff <- as.numeric(args[4]) # MUST BE LAST ARGUMENT because of if(!is.na(cutoff))

numRuns <- 1 # runs per folder (1 is correct)
numTrials <- 1
cluster <- FALSE #just controls cluster plotting


runL <- rep(NA, times = numRuns * numTrials)
trialL <- rep(NA, times = numRuns * numTrials)
num_correctInferencesL <- rep(NA, times = numRuns * numTrials) # true pos
num_incorrect_alphaL <- rep(NA, times = numRuns * numTrials) # species interaction incorrectInference
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
num_incorrect_betaL <- rep(NA, times = numRuns * numTrials)
TP_ignoreSign <- rep(NA, times = numRuns * numTrials)
FP_ignoreSign <- rep(NA, times = numRuns * numTrials)
TN_ignoreSign <- rep(NA, times = numRuns * numTrials)
FN_ignoreSign <- rep(NA, times = numRuns * numTrials)
TP_cluster <- rep(NA, times = numRuns * numTrials)
FP_cluster <- rep(NA, times = numRuns * numTrials)
TN_cluster <- rep(NA, times = numRuns * numTrials)
FN_cluster <- rep(NA, times = numRuns * numTrials)
TP_sign <- rep(NA, times = numRuns * numTrials)
FP_sign <- rep(NA, times = numRuns * numTrials)
TN_sign <- rep(NA, times = numRuns * numTrials)
FN_sign <- rep(NA, times = numRuns * numTrials)

i <- 1
for (run in 1:numRuns) {
    for (trial in 1:numTrials) {
        # store run and trial info
        runL[i] <- run
        trialL[i] <- trial

        simParms <- readRDS(paste0(data_dir, "/params.Rdata"))
        actualAlpha <- sign(simParms$alpha)

        if (covs) {
            actualBeta <- sign(simParms$beta)
            covTypes <- unlist(lapply(simParms$covVars, FUN = function(x) {return(x[["type"]])}))
            if (dim(actualBeta)[2] == length(covTypes)) {
                actualBeta <- actualBeta[, covTypes != "constant"]
            }
        }

        if (!is.na(cutoff)) {
            # check if the inference for this trial finished and store info
            finished_tr <- file.exists(paste0(outdir, "trial", trial, "/inferenceRes_cutoff", cutoff, ".Rdata"))
        } else {
            finished_tr <- file.exists(paste0(outdir, "trial", trial, "/inferenceRes.Rdata"))
        }
        finished_trL[i] <- finished_tr

        if (finished_tr) {
            if (!is.na(cutoff)) {
                inferredParms <- readRDS(paste0(outdir, "trial", trial, "/inferenceRes_cutoff", cutoff, ".Rdata"))
            } else {
                inferredParms <- readRDS(paste0(outdir, "trial", trial, "/inferenceRes.Rdata"))
            }
            
            if (covs) { # take out constant covariate
                covTypes <- unlist(lapply(simParms$covVars, FUN = function(x) {return(x[["type"]])}))
                if (dim(inferredParms$betaInferred)[2] == length(covTypes)) {
                    inferredParms$betaInferred <- inferredParms$betaInferred[, covTypes != "constant"]
                }
            }

            try(if (sum(diag(inferredParms$alphaInferred)) != 0) {stop("Elements of inferred alpha diagonal are not 0")})
            alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
            inferredAlphaG <- graph_from_adjacency_matrix(inferredParms$alphaInferred != 0, mode = "max")
            connected_alpha_actual <- (distances(alphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) != Inf) * 
                                        (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
            connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) != Inf) * 
                                        (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
            # Note: inferredParms$betaInferred * actualBeta == 1 if and only if both are 1 or both are -1
            num_correct_alpha <- sum(inferredParms$alphaInferred * actualAlpha == 1)
            if (covs) {
                num_correct_beta <- sum(inferredParms$betaInferred * actualBeta == 1)
                num_correct <- num_correct_alpha + num_correct_beta
            } else {
                num_correct <- num_correct_alpha
            }

            num_correct_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == 1)

            # Note: here we want to add together inferrence in the wrong direction with saying there is an effect when it's actually 0
            # type 1 error -- inferring an effect where there is none or wrong direction of effect
            count_incorrectT1_alpha <- sum(inferredParms$alphaInferred * actualAlpha == -1) + sum(actualAlpha == 0 & inferredParms$alphaInferred != 0)
            alpha_direction_mistakes <- sum(inferredParms$alphaInferred * actualAlpha == -1)

            if (covs) {
                count_incorrectT1_beta <- sum(inferredParms$betaInferred * actualBeta == -1) + sum(actualBeta == 0 & inferredParms$betaInferred != 0)
                count_incorrectT1 <- count_incorrectT1_beta + count_incorrectT1_alpha
            } else {
                count_incorrectT1 <- count_incorrectT1_alpha
            }

            count_incorrect_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == -1) + 
                                            sum(connected_alpha_actual == 0 & connected_alpha_inferred != 0)

            # type 2 error
            # count missed effects (number of times that actual effect is 1 or -1 and inferred effect is 0)
            num_missedEffects_alpha <- sum(actualAlpha != 0 & inferredParms$alphaInferred == 0)

            if (covs) {
                num_missedEffects_beta <- sum(actualBeta != 0 & inferredParms$betaInferred == 0)
            }

            if (covs) {
                num_missedEffects <- num_missedEffects_alpha + num_missedEffects_beta
            } else {
                num_missedEffects <- num_missedEffects_alpha 
            }

            num_missedEffects_cluster <- sum(connected_alpha_actual != 0 & connected_alpha_inferred == 0)

            alphaInferred <- inferredParms$alphaInferred
            betaInferred <- inferredParms$betaInferred
            if (covs) {
                TP_cluster[i] <- sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1) + 
                                    sum(abs(betaInferred) * abs(actualBeta) == 1)
                FP_cluster[i] <- sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0)) + 
                                    sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0))
                FN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1))
                TN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0)) -
                                    simParms$numSpecies # subtract diagonal
                try(if (TP_cluster[i] + FP_cluster[i] + FN_cluster[i] + TN_cluster[i] != 
                        (simParms$numSpecies^2 - simParms$numSpecies + simParms$numSpecies * (simParms$numCovs - 1))) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx in cluster")})

            } else {
                TP_cluster[i] <- sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1) 
                FP_cluster[i] <- sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0)) 
                FN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1)) 
                TN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0)) -
                                    simParms$numSpecies # subtract diagonal
                try(if (TP_cluster[i] + FP_cluster[i] + FN_cluster[i] + TN_cluster[i] != simParms$numSpecies^2 - simParms$numSpecies) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx in cluster")})
            }

            # undirected meaning that a-->b iff b-->a in "actual alpha". This also ignores the sign.
            undirected_alpha_actual <- distances(alphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) == 1
            undirected_alpha_inferred <- distances(inferredAlphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) == 1

            alpha_incorrect_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == -1) + 
                                            sum(undirected_alpha_actual == 0 & undirected_alpha_inferred != 0)
            alpha_correct_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == 1)

            # add to running lists
            if (covs) {
                num_incorrect_betaL[i] <- count_incorrectT1_beta
                num_missedEffects_betaL[i] <- num_missedEffects_beta
            } 
            
            if (covs) {
                TP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 1)) + 
                                    sum((abs(betaInferred) == 1) & (abs(actualBeta) == 1))
                FP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0)) + 
                                    sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0))
                FN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1))
                TN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0)) -
                                    simParms$numSpecies # subtract diagonal
                try(if (TN_ignoreSign[i] + FN_ignoreSign[i] + FP_ignoreSign[i] + TP_ignoreSign[i] != 
                        (simParms$numSpecies^2 - simParms$numSpecies + simParms$numSpecies * (simParms$numCovs - 1))) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx")})

            } else {
                TP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 1))
                FP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0)) 
                FN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1)) 
                TN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0)) -
                                        simParms$numSpecies # subtract diagonal
                try(if (TP_ignoreSign[i] + FP_ignoreSign[i] + FN_ignoreSign[i] + TN_ignoreSign[i] != simParms$numSpecies^2 - simParms$numSpecies) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx")})
            }

            if (covs) {
                TP_sign[i] <- sum(alphaInferred == 1 & actualAlpha == 1) + 
                                    sum(betaInferred == 1 & actualBeta == 1) +
                                    sum(alphaInferred == -1 & actualAlpha == -1) + 
                                    sum(betaInferred == -1 & actualBeta == -1)
                FP_sign[i] <- sum(alphaInferred != 0 & actualAlpha == 0) + 
                                    sum(betaInferred != 0 & actualBeta == 0) +
                                    sum(alphaInferred == 1 & actualAlpha == -1) + 
                                    sum(betaInferred == 1 & actualBeta == -1) +
                                    sum(alphaInferred == -1 & actualAlpha == 1) + 
                                    sum(betaInferred == -1 & actualBeta == 1) 
                FN_sign[i] <- sum(alphaInferred == 0 & actualAlpha != 0) + 
                                    sum(betaInferred == 0 & actualBeta != 0)
                TN_sign[i] <- sum(alphaInferred == 0 & actualAlpha == 0) + 
                                    sum(betaInferred == 0 & actualBeta == 0) -
                                    simParms$numSpecies # subtract diagonal
                try(if (TP_sign[i] + FP_sign[i] + FN_sign[i] + TN_sign[i] != 
                        (simParms$numSpecies^2 - simParms$numSpecies + simParms$numSpecies * (simParms$numCovs - 1))) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx")})

            } else {
                TP_sign[i] <- sum(alphaInferred != 0 & actualAlpha == alphaInferred)
                FP_sign[i] <- sum(alphaInferred != 0 & actualAlpha == 0) +
                                sum(alphaInferred == 1 & actualAlpha == -1) + 
                                sum(alphaInferred == -1 & actualAlpha == 1) 
                FN_sign[i] <- sum(alphaInferred == 0 & actualAlpha != 0)
                TN_sign[i] <- sum(alphaInferred == 0 & actualAlpha == 0) -
                                    simParms$numSpecies # subtract diagonal
                try(if (TP_sign[i] + FP_sign[i] + FN_sign[i] + TN_sign[i] != simParms$numSpecies^2 - simParms$numSpecies) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx")})
            }

            # add to running lists
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
        if (covs) {
            num_actualEffects <- sum(abs(actualAlpha)) + sum(abs(actualBeta))
        } else {
            num_actualEffects <- sum(abs(actualAlpha))
        }
        num_actualEffectsL[i] <- num_actualEffects

        if (covs) {
            num_possibleEffects <- dim(actualAlpha)[1]^2 - dim(actualAlpha)[1] + (dim(actualBeta)[1] * dim(actualBeta)[2])
        } else {
            num_possibleEffects <- dim(actualAlpha)[1]^2 - dim(actualAlpha)[1]
        }
        num_possibleEffectsL[i] <- num_possibleEffects

        if (cluster) {
          layout <- layout_with_mds(graph = alphaG)
          plot(alphaG, main = "Actual Network", vertex.size = 0)
          png(filename = paste0(outdir, "actualNetworkPlot.png"), height = 800, width = 800, units = "px")
          dev.off()
          layout <- layout_with_mds(graph = inferredAlphaG)
          plot(inferredAlphaG, main = paste0("Inferred Network: ", outdir), layout = layout, vertex.size = 0)
          png(filename = paste0(outdir, "inferredNetworkPlot.png"), height = 800, width = 800, units = "px")
          dev.off()
        }
        i <- i + 1
    }
}

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
            FN_cluster = FN_cluster,
            TP_sign = TP_sign,
            FP_sign = FP_sign,
            TN_sign = TN_sign,
            FN_sign = FN_sign)

# print(df)

if (!is.na(cutoff)) {
    write.csv(df, paste0(outdir, "mistakes", cutoff, ".csv"))
} else {
    write.csv(df, paste0(outdir, "mistakes.csv"))
}