# this is for after running runSchliepSimAnalysis.sh on a bunch of sims
# note: the making of the schliep_runNums.csv is in analyze_multiSim.R

dataDir <- "/space/s1/fiona_callahan/multiSim11/"

schliep_runNums <- read.csv(paste0(dataDir, "schliep_runNums2.csv"), header = FALSE)

source("/home/fiona_callahan/eDNA_sims_code/Schliep_sims_fcns.R")

###############
############## convergence res
numTrials <- 2
for(run in schliep_runNums[, 1]) {
    subdir <- paste0(dataDir, "randomRun", run, "/")
    schliepParms <- readRDS(paste0(subdir, "SchliepRes/trial", 1, "/schliepParms.Rdata"))

    if (numTrials > 1) {
        # takes a long time to load the 50k iterations ones
        if (file.exists(paste0(subdir, "SchliepRes/trial", 1, "/SchliepRes.Rdata")) && 
                        file.exists(paste0(subdir, "SchliepRes/trial", 2, "/SchliepRes.Rdata"))) {
            schliepRes1 <- readRDS(paste0(subdir, "SchliepRes/trial", 1, "/SchliepRes.Rdata"))
            schliepRes2 <- readRDS(paste0(subdir, "SchliepRes/trial", 2, "/SchliepRes.Rdata"))
            allVar_stats <- testConvergence(schliepRes1, schliepRes2, iter = schliepParms$numIters, mode = schliepParms$mode)
            write.csv(allVar_stats, paste0(subdir, "SchliepRes/convergence_diagnostics.csv"))
        }
    }
}
################

#look at convergence results
for(run in schliep_runNums[, 1]){
    if(file.exists(paste0(dataDir, "randomRun", run, "/SchliepRes/convergence_diagnostics.csv"))) {
        print(run)
        covergenceRes <- read.csv(paste0(dataDir, "randomRun", run, "/SchliepRes/convergence_diagnostics.csv"))
        print(covergenceRes)
        print()
    }
}

schliep_inferenceRes <- read.csv(paste0(dataDir, "schliep_inferenceRes_gathered.csv"), header = TRUE)
inla_multiSimRes <- read.csv(paste0(dataDir, "infResGathered.csv"), header = TRUE)

schliep_inferenceRes$totalMistakes <- schliep_inferenceRes$num_incorrectInferences + schliep_inferenceRes$num_missedEffects
inla_multiSimRes$totalMistakes <- inla_multiSimRes$num_incorrectInferences + inla_multiSimRes$num_missedEffectsL

schliepResMerged <- merge(schliep_inferenceRes, inla_multiSimRes, by = "r", all.x = TRUE, all.y = FALSE, suffixes = c("_JSDM", "_INLA"))

plot(schliepResMerged$totalMistakes_JSDM, schliepResMerged$totalMistakes_INLA, xlab = "totalMistakes_JSDM", ylab = "totalMistakes_INLA")
