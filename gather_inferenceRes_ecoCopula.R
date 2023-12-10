# put all of the parameter values and results into a table from a bunch of runs
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("input folders need to be supplied", call. = FALSE)
} 

#Rscript /home/fiona_callahan/eDNA_sims_code/gather_inferenceRes_ecoCopula.R /space/s1/fiona_callahan/multiSim_5sp_testing/ 5 1

thisDir="/space/s1/fiona_callahan/multiSim_manySp_testing2/"
numRuns=10
numTrials=1
resFolderName <- "spiecEasi_res_mb"

thisDir <- args[1]
numRuns <- as.numeric(args[2])
numTrials <- as.numeric(args[3])
resFolderName <- args[4]
cutoff <- as.numeric(args[5])

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
    #if (scramble) {
    #    sitetabName <- paste0(thisDir, "randomRun", run, "/sitetab_scrambled.csv")
    #}
    if (!is.na(cutoff)) {
        mistakes_file <- paste0(thisDir, "randomRun", run, "/", resFolderName, "mistakes", cutoff, ".csv")
    } else {
        mistakes_file <- paste0(thisDir, "randomRun", run, "/", resFolderName, "/mistakes.csv")
    }

    if (file.exists(mistakes_file) 
                && file.exists(sitetabName)
                && file.info(mistakes_file)$size > 0) {
        # load number of mistakes
        mistakes <- read.csv(mistakes_file, header = TRUE)
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
            mistakes_tr <- mistakes[mistakes$trial == trial, ]
            # put data in DF
            thisRow <- c(run, trial, parmVals, percent_presence, mistakes_tr)
            infResDF <- rbind(infResDF, thisRow)
            colnames(infResDF) <- names(thisRow)
        }
    }
}
#saveRDS(parmNames, paste0(thisDir, "parmNames.Rdata"))

# this will fail if the mistakes file was never made but I think that's ok
colnames(infResDF) <- c("RunNum", "inference_trial", parmNames, names(percent_presence), names(mistakes)) 

if (!is.na(cutoff)) {
    write.csv(infResDF, paste0(thisDir, resFolderName, "_infResGathered_cutoff", cutoff, ".csv"))
} else {
    write.csv(infResDF, paste0(thisDir, resFolderName, "_infResGathered.csv"))
}
# note -- beta1, beta2,..., beta9 (and same for alphas) are read off by column
