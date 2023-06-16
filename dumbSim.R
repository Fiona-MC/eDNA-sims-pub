#library(Rlab)
library(ggplot2)
library(reshape)
library(gridExtra)
#install.packages(c('gapminder','gganimate','gifski'))
library(gapminder)
#library(gganimate)
source("/home/fiona_callahan/eDNA_sims_code/multiSimFunctions.R")

# Note to self about a bug you're going to create someday: 
# count_INLAmistakes is going to fail if there is a constant covariate that is not the last one in the list

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("input folder and runstart to run end need to be supplied", call. = FALSE)
} 
# Rscript /home/fiona_callahan/eDNA_sims_code/dumbSim.R /space/s1/fiona_callahan/dumbSim3/ 1 100
#thisdir<-"/space/s1/fiona_callahan/dumbSim/"

thisdir <- args[1]
dir.create(thisdir)

runstart <- args[2]
runend <- args[3]
runs <- runstart:runend

random <- FALSE
parmSet <- 2

betaParm <- 2

for (run in runs) {

  subdir <- paste0(thisdir, "randomRun", run, "/")
  dir.create(subdir)
  params <- getParms(random = random, parmSet = parmSet)
  params$beta <- matrix(data = c(1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1) * betaParm, nrow = params$numSpecies, ncol = params$numCovs, byrow = TRUE)
  saveRDS(params, paste0(subdir, "params.Rdata"))

  num_locations <- params$x_split * params$y_split

  #Get locations
  locList <- get_locations(xdim = params$xdim, ydim = params$ydim, mode = "grid", xSplit = params$x_split, ySplit = params$y_split)
  # get distMx
  locDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
  distMx <- as.matrix(dist(locDF, diag = TRUE, upper = TRUE))
  
  #Run dumb sim
  ######################################################################
  sim_data <- list()
  for (t in 1:params$num_gens) {
    sim_data[[t]] <- list()
    thisY <- list()
    
    # covariates are either constant or iid normal(0,1)
    covL_norm <- list()
    for (covNum in 1:params$numCovs) {
      if (params$covVars[[covNum]]$type != "constant") { # dont normalize if constant covariate
        covL_norm[[covNum]] <- rnorm(n = num_locations, mean = 0, sd = 1)
      } else {
        covL_norm[[covNum]] <- rep(1, times = num_locations)
      }
    }
    
    covList <- list()
    for (s_index in seq_along(locList)){
      covList[[s_index]] <- unlist(lapply(X = covL_norm, FUN = function(covNorm) {covNorm[s_index]}))
    }

    for (s in 1:num_locations) { # start location loop
        linear <- params$beta %*% covList[[s]]
        p_i <-  exp(linear) / (1 + exp(linear))
        y_i <- rbinom(size = 1, n = params$numSpecies, prob = p_i)
      
        thisY[[s]] <- y_i
    } # end location loop
    sim_data[[t]]$y <- thisY
    sim_data[[t]]$covs <- covL_norm
  } # end time loop
  

  #timePts<-1:params$num_gens #time points to "collect"
  timePts <- seq(from = 1, to = params$num_gens, by = 1)
  
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
              # ADD NO COV MEASUREMENT NOISE HERE
              sitetab_data_sampled[i + 3 + cov] <- sim_data[[t]]$covs[[cov]][s_index] #+ rnorm(n = 1, mean = 0, sd = params$covMeasureNoise_sd)
          }
          for (k in 1:params$numSpecies){ # add species presence absence to sitetab
              sitetab_data_sampled[i + 3 + params$numCovs + k] <- sim_data[[t]]$y[[s_index]][k]
          }
          i <- i + 4 + params$numCovs + params$numSpecies
      }
  }

  sim_sitetab_sampled <- as.data.frame(matrix(data = sitetab_data_sampled, nrow = length(time_pts_sample) * length(locations_sample), 
                                    ncol = 4 + params$numCovs + params$numSpecies, byrow = TRUE)) 
  names(sim_sitetab_sampled) <- c("labID", "Age", "Lat", "Long", params$names_all_cov, params$names_species)

  saveRDS(sim_data, paste0(subdir, "sim_data.Rdata"))
  saveRDS(locList, paste0(subdir, "locList.Rdata"))
  saveRDS(covList, paste0(subdir, "covList.Rdata"))
  write.csv(sim_sitetab_sampled, paste0(subdir, "sim_sitetab_sampled.csv"))
}
