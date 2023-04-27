library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("input folder needs to be supplied", call. = FALSE)
} 

# Rscript /home/fiona_callahan/filter_sims.R /space/s1/fiona_callahan/multiSim5/
#thisDir <- "/space/s1/fiona_callahan/multiSim2/"
thisDir <- args[1]
numRuns <- 1000

badRuns <- c()
for(run in 1:numRuns) {
    keepRun = TRUE
    sim_sitetab_sampled <- read.csv(paste0(thisDir, "randomRun", run, "/sim_sitetab_sampled.csv"), header = TRUE)
    params <- readRDS(paste0(thisDir, "randomRun", run, "/params.Rdata"))
    percent_presence <- list()
    for(species in params$names_species) {
        # get #presences/#totalDataPts
        percent_presence[[paste0(species, "PercentPresence")]] <- sum(sim_sitetab_sampled[, species]) / length(sim_sitetab_sampled[, species])
    }
    # get the last time index at which each species was alive (in any location)
    sp_presence_time <- sim_sitetab_sampled %>%
        select(Age, params$names_species) %>%
        group_by(Age) %>%
        summarise(across(everything(), list(max)))
    last_presence_index <- apply(X = sp_presence_time, MARGIN = 2, FUN = function(col){tail(which(col == 1), 1)})
    extinct_time <- sapply(last_presence_index, FUN = function(i){sp_presence_time$Age[i]})[-1]
    # conditions for deleting this run
    if(length(extinct_time) == 0) {
        print(paste0("extinct for all time for run ", run))
        keepRun <- FALSE
    } else if (min(unlist(extinct_time)) < params$num_gens / 2) {
        print(paste0("extinct too early for run ", run))
        keepRun <- FALSE
    }
    if(min(unlist(percent_presence)) < 0.05) {
        print(paste0("at least one species %detection under 5% for run ", run))
        keepRun <- FALSE
    }
    if (max(unlist(percent_presence)) > 0.95) {
        print(paste0("at least one species %detection over 95% for run ", run))
        keepRun <- FALSE
    }
    if (!keepRun) {
        badRuns <- c(badRuns, run)
    }
}

write.csv(data.frame(badRuns = badRuns), file = paste0(thisDir, "unrealistic_runNums.csv"))
