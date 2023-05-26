# take folder with a bunch of multisim results and a file called  
# schliep_runNums.csv which indicates the runs that you've done schliep_sims.R on
# makes a file with the number of mistakes for each of these runs and the parmameters for those runs for analysis

args <- commandArgs(trailingOnly = TRUE)

source("/home/fiona_callahan/eDNA_sims_code/Schliep_sims_fcns.R")

# Rscript /home/fiona_callahan/eDNA_sims_code/schliep_examineRes.R /space/s1/fiona_callahan/multiSim11/ 2

dataDir <- args[1]
numTrials <- as.numeric(args[2])

dataDir <- "/space/s1/fiona_callahan/multiSim11/"
numTrials <- 2

runNums <- read.csv(paste0(dataDir, "schliep_runNums.csv"))[, 1]
trials <- 1:numTrials

# initialize results df
res_df <- data.frame()
for (run in runNums) {
    subdir <- paste0(dataDir, "randomRun", run, "/")
    schliepParms <- readRDS(paste0(subdir, "SchliepRes/trial",trial,"/schliepParms.Rdata"))
    # get sim parms
    simParms <- readRDS(paste0(subdir, "params.Rdata"))
    parmVals <- unlist(simParms)
    # take out functions from parmVals
    parmVals <- parmVals[lapply(parmVals, FUN = class) != "function"]
    for (trial in trials) {
        # get results from schliep
        schliepOutput <- readRDS(paste0(subdir, "SchliepRes/trial",trial, "/SchliepRes.Rdata"))
        plotPosteriorTableOut <- plotPosteriorTable(schliepOutput, ci_percent = 0.95, mode = schliepParms$mode, burn_in = schliepParms$burn)
        thisMistakes <- numMistakes(plotPosteriorTableOut = plotPosteriorTableOut, simParms = simParms, trial = trial, simRum = simRun)
        thisRow <- c(parmVals, thisMistakes[1,])
        res_df <- rbind(res_df, thisRow)
    } # end trials loop
    # convergence testing
    if (numTrials > 1) {
        schliepRes1 <- readRDS(paste0(subdir, "SchliepRes/trial",1, "/SchliepRes.Rdata"))
        schliepRes2 <- readRDS(paste0(subdir, "SchliepRes/trial",2, "/SchliepRes.Rdata"))
        allVar_stats <- testConvergence(schliepRes1, schliepRes2, iter = schliepParms$numIters, mode = schliepParms$mode)
        write.csv(allVar_stats, paste0(subdir, "SchliepRes/convergence_diagnostics.csv"))
    }
} # end runs loop
names(res_df) <- c(names(parmVals), names(thisMistakes))

write.csv(res_df, file = paste0(dataDir, "schliep_inferenceRes_gathered.csv"))