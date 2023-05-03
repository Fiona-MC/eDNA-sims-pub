library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("input folder needs to be supplied", call. = FALSE)
} 

# Rscript /home/fiona_callahan/filter_sims.R /space/s1/fiona_callahan/multiSim5/
#thisDir <- "/space/s1/fiona_callahan/multiSimTest/"
thisDir <- args[1]
numRuns <- 1000


badRuns <- c()

runL <- c()
speciesL <- c()
percent_presenceL <- c()
extinct_timeL <- c()
K_cov_avgL <- c()
K_sp_avgL <- c()

for(run in 1:numRuns) {
    keepRun = TRUE
    sim_sitetab_sampled <- read.csv(paste0(thisDir, "randomRun", run, "/sim_sitetab_sampled.csv"), header = TRUE)
    params <- readRDS(paste0(thisDir, "randomRun", run, "/params.Rdata"))
    sim_data <- readRDS(paste0(subdir, "sim_data.Rdata"))

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
    
    # balance of species effect and cov effect
    k_cov_sum <- rep(0, times = params$numSpecies)
    k_sp_sum <- rep(0, times = params$numSpecies)
    total_n <- 0
    for(t in seq_along(sim_data)){
        for(s in seq_along(sim_data[[t]]$K_sp)){
            total_n <- total_n + 1
            k_cov_sum <- k_cov_sum + sim_data[[t]]$K_cov[[s]] # vector of length numSpecies
            k_sp_sum <- k_sp_sum + sim_data[[t]]$K_sp[[s]] # vector of length numSpecies
        }
    }
    K_cov_avg <- k_cov_sum / total_n
    K_sp_avg <- k_sp_sum / total_n

    runL <- c(runL, rep(run, times = params$numSpecies))
    speciesL <- c(speciesL, params$names_species)
    percent_presenceL <- c(percent_presenceL, unlist(percent_presence))
    extinct_timeL <- c(extinct_timeL ,unlist(extinct_time))
    K_cov_avgL <- c(K_cov_avgL, K_cov_avg)
    K_sp_avgL <- c(K_sp_avgL, K_sp_avg)

    # conditions for deleting this run
    if(length(extinct_time) == 0) {
        thisReason <- "extinct for all time"
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    } else if (min(unlist(extinct_time)) < params$num_gens / 2) {
        thisReason <- "extinct too early"
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    }
    if(min(unlist(percent_presence)) < 0.05) {
        thisReason <- "at least one species %detection under 5%"
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    }
    if (max(unlist(percent_presence)) > 0.95) {
        thisReason <- "at least one species %detection over 95%"
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    }
    if (max(abs(K_sp_avg/K_cov_avg)) < 0.1) {
        # note if there are species with no covariate effects this will never happen
        thisReason <- "max avg (species effect / covariate effect) < 0.1 "
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    }
    if (min(abs(K_sp_avg/K_cov_avg)) > 0.9) { 
        # note if there are species with no species effects this will never happen
        thisReason <- "min avg (species effect / covariate effect) > 0.9"
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    }

    if (mean(abs(K_sp_avg/K_cov_avg)) < 0.05) {
        # note if there are species with no covariate effects this will never happen
        thisReason <- "mean avg (species effect / covariate effect) < 0.05 "
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    }
    if (mean(abs(K_sp_avg/K_cov_avg)) > 0.95) { 
        # note if there are species with no species effects this will never happen
        thisReason <- "mean avg (species effect / covariate effect) > 0.95"
        print(paste0(thisReason, " for run ", run))
        keepRun <- FALSE
    }

    # add bad runs and reasons to running list
    if (!keepRun) {
        badRuns <- c(badRuns, run)
        reasons <- c(reasons, thisReason)
    }
}

stats_df <- data.frame(run = runL, 
                        species = speciesL,
                        percent_presence = percent_presenceL,
                        extinct_time = extinct_timeL,
                        K_cov_avg = K_cov_avgL,
                        K_sp_avg = K_sp_avgL)

write.csv(data.frame(badRuns = badRuns, reason = reasons), file = paste0(thisDir, "unrealistic_runNums.csv"))
write.csv(stats_df, file = paste0(thisDir, "filtering_stats.csv"))
