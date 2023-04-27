# this is for making plots and examining parms for a specific simulation
# like if you simulated a lot with diff parms and then want to examine specific ones 
# input is the output from INLA_simAnalysis_faster.R

#to run: Rscript examine_sim.R /space/s1/fiona_callahan/multiSim2/randomRun1/
source("/home/fiona_callahan/multiSimFunctions.R")
library(dplyr)

#args <- commandArgs(trailingOnly = TRUE)

#subdir <- args[1]
#runNums <- c(447, 254, 487, 169, 426, 143, 102, 385, 373, 335)
runNums <- seq(from = 10, to = 1000, by = 10)

#locPlots <- c(14, 15, 23, 36, 55, 61)
locPlots <- c(19, 8, 6, 37, 64, 46)

for (run in runNums){
    subdir <- paste0("/space/s1/fiona_callahan/multiSim3/randomRun", run, "/")
    sim_data <- readRDS(paste0(subdir, "sim_data.Rdata"))
    params <- readRDS(paste0(subdir, "params.Rdata"))
    locList <- readRDS(paste0(subdir, "locList.Rdata"))

    # make plots
    makePlots(sim_data = sim_data, params = params, locList = locList, locPlots = locPlots, outDir = subdir)

    # look at how many presences per species were actually passed to INLA
    sim_sitetab_sampled <- read.csv(paste0(subdir, "sim_sitetab_sampled.csv"), header = TRUE)
    dim(sim_sitetab_sampled)
    
    sim_sitetab_sampled %>%
        summarise_if(is.numeric, mean, na.rm = TRUE)

    # record average numbers of each species over time
}
