# this is for after running runSchliepSimAnalysis.sh on a bunch of sims
# note: the making of the schliep_runNums.csv is in analyze_multiSim.R

dataDir <- "/space/s1/fiona_callahan/multiSim11/"

schliep_runNums <- read.csv(paste0(dataDir, "schliep_runNums_50k.csv"), header = FALSE)

############## convergence res -- this is in schliep_examineRes.R so it's taken out here
#numTrials <- 2
#source("/home/fiona_callahan/eDNA_sims_code/Schliep_sims_fcns.R")
#for(run in schliep_runNums[, 1]) {
#    subdir <- paste0(dataDir, "randomRun", run, "/")
#    schliepParms <- readRDS(paste0(subdir, "SchliepRes/trial", 1, "/schliepParms.Rdata"))

#    if (numTrials > 1) {
        # takes a long time to load the 50k iterations ones
#        if (file.exists(paste0(subdir, "SchliepRes/trial", 1, "/SchliepRes.Rdata")) && 
#                        file.exists(paste0(subdir, "SchliepRes/trial", 2, "/SchliepRes.Rdata"))) {
#            schliepRes1 <- readRDS(paste0(subdir, "SchliepRes/trial", 1, "/SchliepRes.Rdata"))
#            schliepRes2 <- readRDS(paste0(subdir, "SchliepRes/trial", 2, "/SchliepRes.Rdata"))
#            allVar_stats <- testConvergence(schliepRes1, schliepRes2, iter = schliepParms$numIters, mode = schliepParms$mode)
#            write.csv(allVar_stats, paste0(subdir, "SchliepRes/convergence_diagnostics.csv"))
#        }
#    }
#}
################

# look at convergence results
for(run in schliep_runNums[, 1]) {
    print(run)
    if(file.exists(paste0(dataDir, "randomRun", run, "/SchliepRes/convergence_diagnostics.csv"))) {
        covergenceRes <- read.csv(paste0(dataDir, "randomRun", run, "/SchliepRes/convergence_diagnostics.csv"))
        schliepParms <- readRDS(paste0(dataDir, "randomRun", run, "/SchliepRes/trial1/schliepParms.Rdata"))
        INLA_mistakes <- read.csv(paste0(dataDir, "randomRun", run, "/INLA_res_faster/INLA_mistakes.csv"))
        print(paste("Num iterations:", schliepParms$numIters))
        print(paste("INLA mistakes:", INLA_mistakes$num_incorrectInferences + INLA_mistakes$num_missedEffectsL))
        print(head(covergenceRes, n = 20))
    } else {
        print("No convergence results")
    }
}

for(run in schliep_runNums[, 1]) {
    print(run)
    print(paste("Num Samples Total:", params$num_samples_space,  params$num_samples_time))

    if(file.exists(paste0(dataDir, "randomRun", run, "/SchliepRes/trial1/SchliepRes.Rdata"))) {
        print("Trial1 file exists")
        params <- readRDS(paste0(dataDir, "randomRun", run, "/params.Rdata"))
    } 
    if(file.exists(paste0(dataDir, "randomRun", run, "/SchliepRes/trial2/SchliepRes.Rdata"))) {
        print("Trial2 file exists")
    }
}

schliep_inferenceRes <- read.csv(paste0(dataDir, "schliep_inferenceRes_gathered.csv"), header = TRUE)
inla_multiSimRes <- read.csv(paste0(dataDir, "infResGathered.csv"), header = TRUE)

schliep_inferenceRes$totalMistakes <- schliep_inferenceRes$num_incorrectInferences + schliep_inferenceRes$num_missedEffects
inla_multiSimRes$totalMistakes <- inla_multiSimRes$num_incorrectInferences + inla_multiSimRes$num_missedEffectsL

schliepResMerged <- merge(schliep_inferenceRes, inla_multiSimRes, by = "r", all.x = TRUE, all.y = FALSE, suffixes = c("_JSDM", "_INLA"))

plot(schliepResMerged$totalMistakes_JSDM, schliepResMerged$totalMistakes_INLA, xlab = "totalMistakes_JSDM", ylab = "totalMistakes_INLA")
