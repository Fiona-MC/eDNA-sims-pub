# this is for making plots and examining parms for a specific simulation
# like if you simulated a lot with diff parms and then want to examine specific ones 
# input is the output from INLA_simAnalysis_faster.R

#to run: Rscript examine_sim.R /space/s1/fiona_callahan/multiSim2/randomRun1/
source("/home/fiona_callahan/eDNA_sims_code/multiSimFunctions.R")
library(dplyr)

#args <- commandArgs(trailingOnly = TRUE)
thisDir <- "/space/s1/fiona_callahan/multiSim11/"
#thisDir <- args[1]
#runNums <- c(447, 254, 487, 169, 426, 143, 102, 385, 373, 335) # these are the ones that failed in multiSim 3
#runNums <- seq(from = 10, to = 1000, by = 10) # this is doing it for every 10th one--this now happens in multiSim_Mar2023.R anyways
#runNums <- c(175, 178, 203, 448, 668, 958, 988) # this is the ones that did VERY well in multiSim7
#runNums<-c(1:10)
runNums <- read.csv(paste0(thisDir, "schliep_runNums.csv"))[,1]
runNums <- read.csv(paste0(thisDir, "best_multiSim11_runNums.csv"))[,1]
#locPlots <- c(14, 15, 23, 36, 55, 61)
locPlots <- c(19, 8, 6, 37, 64, 46)

for (run in runNums) {
    subdir <- paste0(thisDir, "randomRun", run, "/")
    sim_data <- readRDS(paste0(subdir, "sim_data.Rdata"))
    params <- readRDS(paste0(subdir, "params.Rdata"))
    locList <- readRDS(paste0(subdir, "locList.Rdata"))

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

    print(run)
    print("K_cov_avg")
    print(K_cov_avg)
    print("K_sp_avg")
    print(K_sp_avg)

    # make plots
    makePlots(sim_data = sim_data, params = params, locList = locList, locPlots = locPlots, outDir = subdir)

    # look at how many presences per species were actually passed to INLA
    sim_sitetab_sampled <- read.csv(paste0(subdir, "sim_sitetab_sampled.csv"), header = TRUE)
    dim(sim_sitetab_sampled)
    
    sim_sitetab_sampled %>%
        summarise_if(is.numeric, mean, na.rm = TRUE)

    # record average numbers of each species over time
}
