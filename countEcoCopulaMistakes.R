#!/usr/bin/env Rscript
library(igraph)
# to run
# Rscript countEcoCopulaMistakes.R /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun1/ /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun1/ecoCopula_res/

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input folders need to be supplied", call. = FALSE)
} 

# TODO! fix this so that it won't count the constant covariate
# but be careful because this is currently running and I don't wanna ruin it
data_dir <- "/space/s1/fiona_callahan/multiSim_manySp_testing2/randomRun3/"
eco_dir <- "/space/s1/fiona_callahan/multiSim_manySp_testing2/randomRun3/spiecEasi_res_mb/"
data_dir <- args[1]
eco_dir <- args[2]

numRuns <- 1 # runs per folder (1 is correct)
numTrials <- 1
cluster <- TRUE

runL <- rep(NA, times = numRuns * numTrials)
trialL <- rep(NA, times = numRuns * numTrials)
num_correctInferencesL <- rep(NA, times = numRuns * numTrials) # true pos
num_incorrect_alphaL <- rep(NA, times = numRuns * numTrials) # species interaction incorrectInference
num_incorrectInferencesL <- rep(NA, times = numRuns * numTrials) # false pos (or wrong direction)
num_actualEffectsL <- rep(NA, times = numRuns * numTrials) # total actual effects
num_possibleEffectsL <- rep(NA, times = numRuns * numTrials) # n^2 + n*p -n (minus n because of the diagonal of alpha)
num_missedEffects_alphaL <- rep(NA, times = numRuns * numTrials)
num_missedEffectsL <- rep(NA, times = numRuns * numTrials)# false negatives
timeL <- rep(NA, times = numRuns * numTrials)
finished_trL <- rep(NA, times = numRuns * numTrials) # T or F whether this trial finished the INLA part
num_correct_clusterL <- rep(NA, times = numRuns * numTrials)
num_incorrect_clusterL <- rep(NA, times = numRuns * numTrials)
num_missed_clusterL <- rep(NA, times = numRuns * numTrials)

i <- 1
for (run in 1:numRuns) {
    for (trial in 1:numTrials) {
        # store run and trial info
        runL[i] <- run
        trialL[i] <- trial

        simParms <- readRDS(paste0(data_dir, "/params.Rdata"))
        actualAlpha <- sign(simParms$alpha)

        # check if the inference for this trial finished and store info
        finished_tr <- file.exists(paste0(eco_dir, "trial", trial, "/inferenceRes.Rdata"))
        finished_trL[i] <- finished_tr

        if (finished_tr) {
          inferredParms <- readRDS(paste0(eco_dir, "trial", trial, "/inferenceRes.Rdata"))
          alphaG <- graph_from_adjacency_matrix(actualAlpha)
          inferredAlphaG <- graph_from_adjacency_matrix(inferredParms$alphaInferred)
          connected_alpha_actual <- (distances(alphaG, v = 1:100, to = 1:100) != Inf) * 
                                    (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
          connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:100, to = 1:100) != Inf) * 
                                    (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
          # Note: inferredParms$betaInferred * actualBeta == 1 if and only if both are 1 or both are -1
          num_correct <- sum(inferredParms$alphaInferred * actualAlpha == 1)
          num_correct_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == 1)

          # Note: here we want to add together inferrence in the wrong direction with saying there is an effect when it's actually 0
          # type 1 error -- inferring an effect where there is none or wrong direction of effect
          count_incorrectT1_alpha <- sum(inferredParms$alphaInferred * actualAlpha == -1) + sum(actualAlpha == 0 & inferredParms$alphaInferred != 0)
          count_incorrectT1 <- count_incorrectT1_alpha

          count_incorrect_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == -1) + 
                                          sum(connected_alpha_actual == 0 & connected_alpha_inferred != 0)

          # type 2 error
          # count missed effects (number of times that actual effect is 1 or -1 and inferred effect is 0)
          num_missedEffects_alpha <- sum(actualAlpha != 0 & inferredParms$alphaInferred == 0)
          num_missedEffects <- num_missedEffects_alpha 

          num_missedEffects_cluster <- sum(connected_alpha_actual != 0 & connected_alpha_inferred == 0)

          # add to running lists
          num_incorrect_alphaL[i] <- count_incorrectT1_alpha
          num_correctInferencesL[i] <- num_correct
          num_incorrectInferencesL[i] <- count_incorrectT1
          num_missedEffects_alphaL[i] <- num_missedEffects_alpha
          num_missedEffectsL[i] <- num_missedEffects

          num_correct_clusterL[i] <- num_correct_cluster
          num_incorrect_clusterL[i] <- count_incorrect_cluster
          num_missed_clusterL[i] <- num_missedEffects_cluster

          # note this "if" is just because I changed this during a run -- can be taken out later (or not)
          if ("time" %in% names(inferredParms)) {
            timeL[i] <- inferredParms$time
          } 
        } 
        
        # count total number of actual effects in the model
        num_actualEffects <- sum(abs(actualAlpha))

        num_actualEffectsL[i] <- num_actualEffects
        num_possibleEffectsL[i] <- dim(actualAlpha)[1]^2 - dim(actualAlpha)[1]

        if (cluster) {
          layout <- layout_with_mds(graph = alphaG)
          plot(alphaG, main = "Actual Network", vertex.size = 0)
          png(filename = paste0(eco_dir, "actualNetworkPlot.png"), height = 800, width = 800, units = "px")
          dev.off()
          layout <- layout_with_mds(graph = inferredAlphaG)
          plot(inferredAlphaG, main = paste0("Inferred Network: ", eco_dir), layout = layout, vertex.size = 0)
          png(filename = paste0(eco_dir, "inferredNetworkPlot.png"), height = 800, width = 800, units = "px")
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
            num_correctInferences = num_correctInferencesL, 
            num_incorrectInferences = num_incorrectInferencesL, 
            num_actualEffects = num_actualEffectsL,
            num_missedEffects_alpha = num_missedEffects_alphaL,
            num_missedEffectsL = num_missedEffectsL,
            num_possibleEffectsL = num_possibleEffectsL,
            num_correct_cluster = num_correct_clusterL,
            num_incorrect_cluster = num_incorrect_clusterL,
            num_missed_cluster = num_missed_clusterL)

# print(df)

write.csv(df, paste0(eco_dir, "mistakes.csv"))
