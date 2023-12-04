#setwd("/home/fiona_callahan/simData/FALSE_POS/independent_FPR/")
#library(Rlab)
library(ggplot2)
library(reshape)
library(gridExtra)
#install.packages(c('gapminder','gganimate','gifski'))
library(gapminder)
#install.packages("data.table")
library(data.table)
#library(gganimate)
source("/home/fiona_callahan/eDNA_sims_code/multiSimFunctions.R")
source("/home/fiona_callahan/eDNA_sims_code/multiSimParms.R")

# Note to self about a bug you're going to create someday: 
# count_INLAmistakes is going to fail if there is a constant covariate that is not the last one in the list

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("input folder and runstart to run end need to be supplied", call. = FALSE)
} 
# Rscript /home/fiona_callahan/eDNA_sims_code/multiSim_Mar2023.R /space/s1/fiona_callahan/multiSim_100/ 1 100
# thisdir<-"/space/s1/fiona_callahan/multiSim5/"
random <- FALSE
#parmSet <- "indep" # indep means that all alphas will be 0
parmSet <- 1 # this is default
# amount of generations to prime the sim before starting to record -- (not implemented TODO)
burn <- 100

spNumMode <- "many"
#nSp <- 2
nSp <- 10
readAbdMode <- TRUE

thisdir <- "/space/s1/fiona_callahan/multiSim_100/"
thisdir <- args[1]
dir.create(thisdir)

runstart <- args[2]
runend <- args[3]
runs <- runstart:runend

# controls whether sim_sitetab_longer.csv is created with all time points or just the sampled table
longer <- FALSE 

print(Sys.time())
for (run in runs) {
  if(run %% 10 == 1) {
    makePlotsTF <- TRUE
  } else {
    makePlotsTF <- FALSE
  }
  print(run)
  subdir <- paste0(thisdir, "randomRun", run, "/")
  dir.create(subdir)

  if(spNumMode == 5) {
    params <- getParms_5(random = random, parmSet = parmSet)
  } else if(spNumMode == "many") {
    params <- getParms_many(random = random, parmSet = parmSet, numSpecies = nSp)
  } else {
    params <- getParms(random = random, parmSet = parmSet)
  }
  saveRDS(params, paste0(subdir, "params.Rdata"))
  
  num_locations <- params$x_split * params$y_split
  
  #Get locations
  locList <- get_locations(xdim = params$xdim, ydim = params$ydim, mode = "grid", xSplit = params$x_split, ySplit = params$y_split)
  # get distMx
  locDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
  distMx <- as.matrix(dist(locDF, diag = TRUE, upper = TRUE))
  
  #Initialize abundances with all 10's
  N0 <- list()
  for (s_index in seq_along(locList)) {
    N0[[s_index]] <- rep(x = params$N0_0, times = params$numSpecies)
  }
  
  #Run simulation
  #5 seconds for 64 locations and 1000 time steps
  ######################################################################
  #Run sim
  ######################################################################
  #
  sim_data <- list()
  t_minus_one <- list()
  for (t in 1:params$num_gens){
    #if(t %% 100 == 1){
    #  print("t=")
    #  print(t)
    #}
    sim_data[[t]] <- list()
    thisK <- list()
    thisK_cov <- list()
    thisK_sp <- list()
    thisN <- list()
    thisY <- list()
    if(readAbdMode) {
      thisReadAbd <- list()
    }
    
    FP_count <- 0
    
    # TODO where should I implement adding squares of covariates?
    # TODO think more about this normalization scheme -- normalizing at each time point across space
    # get covariates (normalize)
    covL_norm <- list()
    for (covNum in 1:params$numCovs) {
      thisF1 <- params$covVars[[covNum]]$coef1func
      thisF2 <- params$covVars[[covNum]]$coef2func
      thisType <- params$covVars[[covNum]]$type
      if (thisType == "randomWalk") {
        if (t == 1) {
          cov <- rnorm(n = length(locList), mean = 0, sd = 1)
          t_minus_one[[covNum]] <- cov
        } else {
          cov <- unlist(lapply(X = seq_along(locList), FUN = 
                           function(loc_index) {return(get_cov(locList[[loc_index]], type = thisType, 
                                                        coef1func = thisF1, 
                                                        coef2func = thisF2, 
                                                        t = t, t_minus_1 = t_minus_one[[covNum]][loc_index], params = params))}))
          t_minus_one[[covNum]] <- cov
        }
      } else if (thisType == "randomNormal") {
        cov <- rnorm(n = length(locList), mean = 0, sd = 1)
      } else {
      cov <- unlist(lapply(X = locList, FUN = 
                           function(loc) {return(get_cov(loc, type = thisType, 
                                                        coef1func = thisF1, 
                                                        coef2func = thisF2, 
                                                        t = t, params = params))}))
      }
      if (params$covVars[[covNum]]$type != "constant") { # dont normalize if constant covariate
        covL_norm[[covNum]] <- (cov - mean(cov)) / sd(cov)
        # add noise (NOTE THIS IS PROCESS NOISE NOT MEASUREMENT NOISE)
        covL_norm[[covNum]] <- covL_norm[[covNum]] + rnorm(n = num_locations, mean = 0, sd = params$covNoise_sd)
      } else {
        covL_norm[[covNum]] <- cov
      }
    }
    
    covList <- list()
    for (s_index in seq_along(locList)){
      covList[[s_index]] <- unlist(lapply(X = covL_norm, FUN = function(covNorm) {covNorm[s_index]}))
    }
    
    # start space loop
    for (s in seq_along(locList)){
      if (t == 1) {
        lastN <- N0[[s]]
      } else {
        lastN <- sim_data[[t - 1]]$N[[s]]
      }
      # calculate carrying capacities -- divide into covariate and species effects for later analysis
      thisK_cov[[s]] <- params$beta %*% covList[[s]]
      thisK_sp[[s]] <- params$alpha %*% abundanceEffectTransform(lastN, type = params$abundanceEffectType)
      thisK[[s]] <- thisK_cov[[s]] + thisK_sp[[s]]
      thisK[[s]] <- apply(thisK[[s]], FUN = function(x) {max(x, 0.001)}, MARGIN = 1) # truncate at 0-ish


      # calculate new abundances
      # calculate noise
      dWt <- rnorm(n = params$numSpecies, mean = 0, sd = 1)
      
      # get numNeighbors (total number of animals in neighboring locations)
      # neighbors<-get_neighbors(location = s, xSplit = x_split, ySplit = y_split)
      distFromHere <- distMx[s, ] # get row in distance mx that corresponds to this spatial location
      isNeighbor <- distFromHere > 0 & distFromHere <= params$radius
      neighbors <- c(1:num_locations)[isNeighbor]
      numNeighbors <- rep(0, times = params$numSpecies)
      for (neighbor in neighbors){
        if (t == 1) {
          numNeighbors <- numNeighbors + N0[[neighbor]] # these are vectors of length numSpecies
        } else {
          numNeighbors <- numNeighbors + sim_data[[t - 1]]$N[[neighbor]]
        }
      }
      
      M <- rpois(n = params$numSpecies, lambda = params$mig_rates * numNeighbors) 
      # note: mig_rates*numNeighbors is a vector -- these rates are different per species
      lastN <- lastN + M
      dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      thisN[[s]] <- lastN + dN
      # truncate at 0 -- no negative populations
      thisN[[s]] <- apply(matrix(thisN[[s]]), FUN = function(x) {max(x, 0)}, MARGIN = 1)  # nolint: brace_linter.
      
      if(readAbdMode) {
        #lambdaDep <- params$indivSampleRate * thisN
        N_depDNA <- rbinom(n = params$numSpecies, size = as.integer(thisN[[s]]), p = params$indivSampleProb)
        lambdaReads <- params$readSampleRate * N_depDNA 
        Nreads <- rpois(n = params$numSpecies, lambda = lambdaReads)
        thisReadAbd[[s]] <- Nreads
        thisY[[s]] <- as.integer(Nreads >= params$readThreshold)
      } else {
        # calculate new response data
        detect_prob <- (thisN[[s]]^params$det_prob_exp) / (params$det_prob_add + thisN[[s]]^params$det_prob_exp)
        thisY[[s]] <- rbinom(n = params$numSpecies, size = 1, prob = detect_prob)
        
        # add false detections (positives)
        # TODO how to add option of whether fpr changes with location and time?
        # fpr_mode = "constant" # fpr is constant 
        # sample an fpr for this time and location
        fpr_mode <- params$fpr$mode
        if (fpr_mode == "none") {
          # no false positives (ie rate is 0)
          fpr <- rep(0, params$numSpecies)
        } else if (fpr_mode == "constant") {
          # fpr is a set constant  -- this is basically just for experimenting
          fpr <- rep(params$fpr$constant_fpr, params$numSpecies)
        } else if (fpr_mode == "independent") {
          # draw separate fpr for each species
          fpr <- rbeta(n = params$numSpecies, shape1 = params$fpr$alpha_fpr, shape2 = params$fpr$beta_fpr)
        } else if (fpr_mode == "dependent_sp") {
          # draw one fpr for this location and time and apply to all species
          fpr <- rep(rbeta(n = 1, shape1 = params$fpr$alpha_fpr, shape2 = params$fpr$beta_fpr), params$numSpecies)
        } else {
          print("your fpr_mode is not implemented")
        }
        # apply this fpr to all species (bernoulli trials)
        for (species in 1:params$numSpecies){
          if (thisY[[s]][species] == 0) { #if species is not detected (thisY[[s]] == 0) then set to 1 with prob fpr
            FP <- rbinom(n = 1, size = 1, prob = fpr[species]) 
            thisY[[s]][species] <- FP
            if (FP) {
              FP_count <- FP_count + 1
            } 
          }
        }
      }
    } # end location loop
    sim_data[[t]]$K <- thisK
    sim_data[[t]]$K_cov <- thisK_cov
    sim_data[[t]]$K_sp <- thisK_sp
    sim_data[[t]]$N <- thisN
    sim_data[[t]]$y <- thisY
    sim_data[[t]]$covs <- covL_norm
    sim_data[[t]]$FP_count <- FP_count
    if(readAbdMode) {
      sim_data[[t]]$readAbd <- thisReadAbd
    }
  } # end time loop

  #timePts<-1:params$num_gens #time points to "collect"
  timePts <- seq(from = 1, to = params$num_gens, by = 1)
  if (longer) {
    # Wrangle data
    #TODO -- this is taking WAY too long with more species lol
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
              sitetab_data[i + 3 + cov] <- sim_data[[t]]$covs[[cov]][s_index]
          }
          sitetab_data[i + 4 + params$numCovs] <- k
          sitetab_data[i + 5 + params$numCovs] <- sim_data[[t]]$y[[s_index]][k] # add presence absence for this species
          sitetab_data[i + 6 + params$numCovs] <- sim_data[[t]]$N[[s_index]][k] #add abundances here
          i <- i + 4 + params$numCovs + 3
        }
      }
    }
    
    sim_sitetab_longer <- as.data.frame(matrix(data = sitetab_data, nrow = length(timePts) * length(locList) * params$numSpecies, 
                                      ncol = 7 + params$numCovs, byrow = TRUE))

    namesCov <- sapply(1:params$numCov, FUN = function(cov) {paste0("Cov", cov)})
    names(sim_sitetab_longer) <- c("labID", "Age", "Lat", "Long", namesCov, "Species", "Presence", "Abundance")
    
    # TODO this is taking a stupid long time
    fwrite(sim_sitetab_longer, paste0(subdir, "sitetab_longer.csv"))
    #write.csv(sim_sitetab_longer, paste0(subdir, "sitetab_longer.csv"))
  }
  
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
              sitetab_data_sampled[i + 3 + cov] <- sim_data[[t]]$covs[[cov]][s_index] + rnorm(n = 1, mean = 0, sd = params$covMeasureNoise_sd)
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

  if (makePlotsTF) {
    if(params$numSpecies < 10) {
      makePlots(sim_data = sim_data, params = params, locList = locList, outDir = subdir, 
              sitetab_sampled = sim_sitetab_sampled, spList = 1:params$numSpecies)
    } else {
      makePlots(sim_data = sim_data, params = params, locList = locList, outDir = subdir, 
              sitetab_sampled = sim_sitetab_sampled, spList = c(1, 2, 3))
    }
  }
  saveRDS(sim_data, paste0(subdir, "sim_data.Rdata"))
  saveRDS(locList, paste0(subdir, "locList.Rdata"))
  saveRDS(covList, paste0(subdir, "covList.Rdata"))
  fwrite(sim_sitetab_sampled, paste0(subdir, "sim_sitetab_sampled.csv"))
  #write.csv(sim_sitetab_sampled, paste0(subdir, "sim_sitetab_sampled.csv"))

  if(readAbdMode) {
    if (longer) {
      # Wrangle data for abd
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
              sitetab_data[i + 3 + cov] <- sim_data[[t]]$covs[[cov]][s_index]
          }
          sitetab_data[i + 4 + params$numCovs] <- k
          sitetab_data[i + 5 + params$numCovs] <- sim_data[[t]]$readAbd[[s_index]][k] # add presence absence for this species
          sitetab_data[i + 6 + params$numCovs] <- sim_data[[t]]$N[[s_index]][k] #add abundances here
          i <- i + 4 + params$numCovs + 3
          }
      }
      }

      sim_sitetab_longer <- as.data.frame(matrix(data = sitetab_data, nrow = length(timePts) * length(locList) * params$numSpecies, 
                                      ncol = 7 + params$numCovs, byrow = TRUE))

      namesCov <- sapply(1:params$numCov, FUN = function(cov) {paste0("Cov", cov)})
      names(sim_sitetab_longer) <- c("labID", "Age", "Lat", "Long", namesCov, "Species", "ReadAbd", "Abundance")

      fwrite(sim_sitetab_longer, paste0(subdir, "sitetab_longer_readAbd.csv"))
      #write.csv(sim_sitetab_longer, paste0(subdir, "sitetab_longer_readAbd.csv"))
    }

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
                sitetab_data_sampled[i + 3 + cov] <- sim_data[[t]]$covs[[cov]][s_index] + rnorm(n = 1, mean = 0, sd = params$covMeasureNoise_sd)
            }
            for (k in 1:params$numSpecies){ # add read abd absence to sitetab
                sitetab_data_sampled[i + 3 + params$numCovs + k] <- sim_data[[t]]$readAbd[[s_index]][k]
            }
            i <- i + 4 + params$numCovs + params$numSpecies
        }
    }

    sim_sitetab_sampled <- as.data.frame(matrix(data = sitetab_data_sampled, nrow = length(time_pts_sample) * length(locations_sample), 
                                    ncol = 4 + params$numCovs + params$numSpecies, byrow = TRUE)) 
    names(sim_sitetab_sampled) <- c("labID", "Age", "Lat", "Long", params$names_all_cov, params$names_species)

    saveRDS(sim_data, paste0(subdir, "sim_data_abd.Rdata"))
    fwrite(sim_sitetab_sampled, paste0(subdir, "sim_sitetab_readAbd_sampled.csv"))
    #write.csv(sim_sitetab_sampled, paste0(subdir, "sim_sitetab_readAbd_sampled.csv"))


    if (makePlotsTF == TRUE) {
        p <- plotDetectProb(abdParms = params)
        ggsave(p, file = paste0(subdir, "detProbPlot_abd.png"), width = 7, height = 7)
    }
  }
  print(Sys.time())
}
