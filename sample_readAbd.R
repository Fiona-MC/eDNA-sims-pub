# takes the results of multiSim_Mar2023.R and samples the abundances of animals into read abundances
# then normalizes the read abundances (?)
library(ggplot2)
data_dir <- "/space/s1/fiona_callahan/multiSim11a/randomRun30/"

# load data
sim_data_raw <- readRDS(paste0(data_dir, "sim_data.Rdata"))
locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))

# indivSampleRate is the poission rate of sampling individuals to deposit DNA 
# readSampleRate is the number of reads per individual that deposits DNA (avg) 
# TODO (I think probably I want this to be NB not poission but for now it's poisson)
abdParms <- list(indivSampleProb = 0.01, readSampleRate = 1, corr_mode = "independent")

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

saveRDS(sim_data_raw, paste0(data_dir, "sim_data_abd.Rdata"))

#plot P(#reads>0) vs #indiv
N_max <- 1000
N_set <- seq(1, N_max, 10)
P_detect <- rep(NA, times = length(N_set))
for (i in seq_along(N_set)) {
    N <- N_set[i]
    z_set <- seq(0, N, 1)
    P_detect[i] <- sum((1 - dpois(x = 0, lambda = abdParms$readSampleRate * z_set)) * dbinom(x = z_set, size = N, p = abdParms$indivSampleProb))
}
plot(x = N_set, y = P_detect)
