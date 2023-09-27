# put all of the parameter values and results into a table from a bunch of inla runs
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("input folders need to be supplied", call. = FALSE)
} 

#Rscript /home/fiona_callahan/eDNA_sims_code/gather_inferenceRes.R /space/s1/fiona_callahan/multiSim11/ 1000 1

thisDir <- "/space/s1/fiona_callahan/multiSim11/"
numRuns <- 1000
numTrials <- 1

thisDir="/space/s1/fiona_callahan/multiSim_rw3/"
#numRuns=100
numRuns=1
numTrials=1

thisDir <- args[1]
numRuns <- as.numeric(args[2])
numTrials <- as.numeric(args[3])
scramble <- (as.numeric(args[4]) == 1)
if (scramble) {
    print("scramble")
}
INLAname <- "paper" # name of folder with INLA results is INLA_res_[INLAname]
#INLAname <- "faster"

runs <- 1:numRuns
trials <- 1:numTrials

infResDF <- data.frame()
for (run in runs) {
    # note there is already a check for this in the count_INLAmistakes, 
    # so this is just in case I wanna run this on a whole folder that I havent finished yet
    # load parameter values
    if(file.exists(paste0(thisDir, "randomRun", run))) { # for runs that were sorted into the unrealistic folder
        params <- readRDS(paste0(thisDir, "randomRun", run, "/params.Rdata"))
        parmVals <- unlist(params)
        # take out functions from parmVals
        parmVals <- parmVals[lapply(parmVals, FUN = class) != "function"]
        parmNames <- names(parmVals)
    }
    sitetabName <- paste0(thisDir, "randomRun", run, "/sim_sitetab_sampled.csv")
    if (scramble) {
        sitetabName <- paste0(thisDir, "randomRun", run, "/sitetab_scrambled.csv")
    }
    if (file.exists(paste0(thisDir, "randomRun", run, "/INLA_res_", INLAname, "/INLA_mistakes.csv")) 
                && file.exists(sitetabName)
                && file.info(paste0(thisDir, "randomRun", run, "/INLA_res_", INLAname, "/INLA_mistakes.csv"))$size > 0) {
        # load number of mistakes
        INLA_mistakes <- read.csv(paste0(thisDir, "randomRun", run, "/INLA_res_", INLAname, "/INLA_mistakes.csv"), header = TRUE)
        sim_sitetab_sampled <- read.csv(sitetabName, header = TRUE)
        # get percent presence per species
        percent_presence <- list()
        for(species in params$names_species) {
            # get #presences/#totalDataPts
            percent_presence[[paste0(species, "PercentPresence")]] <- sum(sim_sitetab_sampled[, species]) / length(sim_sitetab_sampled[, species])
        }
        percent_presence <- unlist(percent_presence)
        for (trial in trials) {
            # load inference results (i think this is not necessary)
            # INLAres <- readRDS(paste0(thisDir, "randomRun", run, "/INLA_res/trial", trial, "/inferenceRes.Rdata"))
            INLA_mistakes_tr <- INLA_mistakes[INLA_mistakes$trial == trial, ]
            # put data in DF
            thisRow <- c(run, trial, parmVals, percent_presence, INLA_mistakes_tr)
            infResDF <- rbind(infResDF, thisRow)
            colnames(infResDF) <- names(thisRow)
        }
    }
}
saveRDS(parmNames, paste0(thisDir, "parmNames.Rdata"))

# this will fail if the INLA_mistakes file was never made but I think that's ok
colnames(infResDF) <- c("RunNum", "INLA_trial", parmNames, names(percent_presence), names(INLA_mistakes)) 

write.csv(infResDF, paste0(thisDir, "infResGathered.csv"))
# note -- beta1, beta2,..., beta9 (and same for alphas) are read off by column
