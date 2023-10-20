# NOTE -- ALL THIS CODE HAS NOW BEEN TRANSPLANTED INTO multiSim_Mar2023.R to happen if readAbdMode == T






# takes the results of multiSim_Mar2023.R and samples the abundances of animals into read abundances
# NOTE -- ALL THIS CODE HAS NOW BEEN TRANSPLANTED INTO multiSim_Mar2023.R to happen if readAbdMode == T
# then normalizes the read abundances (?)
library(ggplot2)
source("/home/fiona_callahan/eDNA_sims_code/multiSimFunctions.R")

# Rscript sample_readAbd.R /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun4/ ${random} ${plot}
# Rscript sample_readAbd.R /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun4/ 1 1

args <- commandArgs(trailingOnly = TRUE)

data_dir <- "/space/s1/fiona_callahan/multiSim_5sp_testing/randomRun4/"
data_dir <- args[1]

random <- TRUE
plot <- TRUE
random <- (as.integer(args[2]) == 1)
plot <- (as.integer(args[3]) == 1)

# load data
sim_data_raw <- readRDS(paste0(data_dir, "sim_data.Rdata"))
locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))

# indivSampleRate is the poission rate of sampling individuals to deposit DNA 
# readSampleRate is the number of reads per individual that deposits DNA (avg) 
# TODO (I think probably I want this to be NB not poission but for now it's poisson)
if (random) {
    indivSampleProb <- runif(n = 1, 0.001, 0.02)
    readSampleRate <- runif(n = 1, 0.5, 100)
    corr_mode <- sample(c("independent")) # this is the only one implemented thus far
} else {
    indivSampleProb <- 0.01
    readSampleRate <- 1
    corr_mode <- "independent"
}

abdParms <- list(indivSampleProb = indivSampleProb, readSampleRate = readSampleRate, corr_mode = corr_mode)
saveRDS(abdParms, paste0(data_dir, "abdParms.Rdata"))

for(t in 1:params$num_gens) {
    for(s in seq_along(locList)) {
        thisN <- sim_data_raw[[t]]$N[[s]] # vector of N for all species
        lambdaDep <- abdParms$indivSampleRate * thisN
        N_depDNA <- rbinom(n = params$numSpecies, size = as.integer(thisN), p = abdParms$indivSampleProb)
        lambdaReads <- abdParms$readSampleRate * N_depDNA 
        Nreads <- rpois(n = params$numSpecies, lambda = lambdaReads)
        sim_data_raw[[t]]$readAbd[[s]] <- Nreads
    }
}

#timePts<-1:params$num_gens #time points to "collect"
timePts <- seq(from = 1, to = params$num_gens, by = 1)

# Wrangle data
sitetab_data <- numeric(length(timePts) * length(params$locList) * params$numSpecies * (7 + params$numCovs))
i <- 1
for (t in timePts){
for (s_index in seq_along(locList)){
    for (k in 1:params$numSpecies){
    sitetab_data[i] <- s_index
    sitetab_data[i + 1] <- params$num_gens - t # age
    sitetab_data[i + 2] <- locList[s_index][[1]][1] #x_location
    sitetab_data[i + 3] <- locList[s_index][[1]][2] #y_location
    for (cov in 1:params$numCovs) { # add covariate values 
        sitetab_data[i + 3 + cov] <- sim_data_raw[[t]]$covs[[cov]][s_index]
    }
    sitetab_data[i + 4 + params$numCovs] <- k
    sitetab_data[i + 5 + params$numCovs] <- sim_data_raw[[t]]$readAbd[[s_index]][k] # add presence absence for this species
    sitetab_data[i + 6 + params$numCovs] <- sim_data_raw[[t]]$N[[s_index]][k] #add abundances here
    i <- i + 4 + params$numCovs + 3
    }
}
}

sim_sitetab_longer <- as.data.frame(matrix(data = sitetab_data, nrow = length(timePts) * length(locList) * params$numSpecies, 
                                ncol = 7 + params$numCovs, byrow = TRUE))

namesCov <- sapply(1:params$numCov, FUN = function(cov) {paste0("Cov", cov)})
names(sim_sitetab_longer) <- c("labID", "Age", "Lat", "Long", namesCov, "Species", "ReadAbd", "Abundance")

write.csv(sim_sitetab_longer, paste0(data_dir, "sitetab_longer_readAbd.csv"))

# sample from full sitetab
time_pts_sample <- seq(from = 1, to = params$num_gens, by = round(1000 / params$num_samples_time))
locations_sample <- sort(sample(x = seq_along(locList), size = params$num_samples_space))

# Wrangle data into format for analysis
# the 4 is for the spatial index, generation, x-location, and y-location
sitetab_data_sampled <- numeric(length(time_pts_sample) * length(locations_sample) *
                (4 + params$numCovs + params$numSpecies))
i <- 1
for (t in time_pts_sample) {
    for (s_index in locations_sample) { 
        sitetab_data_sampled[i] <- s_index
        sitetab_data_sampled[i + 1] <- params$num_gens - t # age
        sitetab_data_sampled[i + 2] <- locList[s_index][[1]][1]
        sitetab_data_sampled[i + 3] <- locList[s_index][[1]][2]
        for (cov in 1:params$numCovs){ # add covariate values to sitetab
            # ADD COV MEASUREMENT NOISE HERE
            sitetab_data_sampled[i + 3 + cov] <- sim_data_raw[[t]]$covs[[cov]][s_index] + rnorm(n = 1, mean = 0, sd = params$covMeasureNoise_sd)
        }
        for (k in 1:params$numSpecies){ # add read abd absence to sitetab
            sitetab_data_sampled[i + 3 + params$numCovs + k] <- sim_data_raw[[t]]$readAbd[[s_index]][k]
        }
        i <- i + 4 + params$numCovs + params$numSpecies
    }
}

sim_sitetab_sampled <- as.data.frame(matrix(data = sitetab_data_sampled, nrow = length(time_pts_sample) * length(locations_sample), 
                                ncol = 4 + params$numCovs + params$numSpecies, byrow = TRUE)) 
names(sim_sitetab_sampled) <- c("labID", "Age", "Lat", "Long", params$names_all_cov, params$names_species)

if (plot) {
    makePlots(sim_data = sim_data_raw, params = params, locList = locList, outDir = data_dir, sitetab_sampled = sim_sitetab_sampled)
}
saveRDS(sim_data_raw, paste0(data_dir, "sim_data_abd.Rdata"))
write.csv(sim_sitetab_sampled, paste0(data_dir, "sim_sitetab_readAbd_sampled.csv"))


if (plot == TRUE) {
    p <- plotDetectProb(abdParms = abdParms)
    ggsave(p, file = paste0(data_dir, "detProbPlot_abd.png"), width = 7, height = 7)
}
