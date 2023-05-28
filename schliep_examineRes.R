# take folder with a bunch of multisim results and a file called  
# schliep_runNums.csv which indicates the runs that you've done schliep_sims.R on
# makes a file with the number of mistakes for each of these runs and the parmameters for those runs for analysis

args <- commandArgs(trailingOnly = TRUE)

source("/home/fiona_callahan/eDNA_sims_code/Schliep_sims_fcns.R")

# Rscript /home/fiona_callahan/eDNA_sims_code/schliep_examineRes.R /space/s1/fiona_callahan/multiSim11/ 2

dataDir <- "/space/s1/fiona_callahan/multiSim11/"
numTrials <- 2

dataDir <- args[1]
numTrials <- as.numeric(args[2])

runNums <- read.csv(paste0(dataDir, "schliep_runNums.csv"))[, 1]
trials <- 1:numTrials

# initialize results df
res_df <- data.frame()
runsL <- c()
trialsL <- c()
for (run in runNums) {
    if (file.exists(paste0(dataDir, "randomRun", run, "/SchliepRes/trial1/SchliepRes.Rdata")) && file.exists(paste0(dataDir, "randomRun", run, "/SchliepRes/trial2/SchliepRes.Rdata"))) {
        print(paste("Running run number", run))
        subdir <- paste0(dataDir, "randomRun", run, "/")
        # get sim parms
        simParms <- readRDS(paste0(subdir, "params.Rdata"))
        parmVals <- unlist(simParms)
        # take out functions from parmVals
        parmVals <- parmVals[lapply(parmVals, FUN = class) != "function"]
        for (trial in trials) {
            schliepParms <- readRDS(paste0(subdir, "SchliepRes/trial",trial,"/schliepParms.Rdata"))
            # get results from schliep
            schliepOutput <- readRDS(paste0(subdir, "SchliepRes/trial",trial, "/SchliepRes.Rdata"))
            plotPosteriorTableOut <- plotPosteriorTable(schliepOutput, ci_percent = 0.95, mode = schliepParms$mode, burn_in = schliepParms$burn)
            thisMistakes <- numMistakes(plotPosteriorTableOut = plotPosteriorTableOut, simParms = simParms, trial = trial, simRun = run)
            write.csv(thisMistakes, paste0(subdir, "SchliepRes/trial", trial, "/schliep_mistakes.csv"))
            thisRow <- data.frame(c(parmVals, thisMistakes[1,]))
            res_df <- merge(res_df, thisRow, all = TRUE)
            runsL <- c(runsL, run)
            trialsL <- c(trialsL, trial)
        } # end trials loop
        # convergence testing
        if (numTrials > 1) {
            schliepRes1 <- readRDS(paste0(subdir, "SchliepRes/trial",1, "/SchliepRes.Rdata"))
            schliepRes2 <- readRDS(paste0(subdir, "SchliepRes/trial",2, "/SchliepRes.Rdata"))
            allVar_stats <- testConvergence(schliepRes1, schliepRes2, iter = schliepParms$numIters, mode = schliepParms$mode)
            write.csv(allVar_stats, paste0(subdir, "SchliepRes/convergence_diagnostics.csv"))
        }
    }
} # end runs loop
names(res_df) <- c(names(parmVals), names(thisMistakes))
res_df$RunNum <- runsL
res_df$SchiepTrial <- trialsL

write.csv(res_df, file = paste0(dataDir, "schliep_inferenceRes_gathered.csv"))