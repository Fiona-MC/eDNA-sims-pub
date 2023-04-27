#!/usr/bin/env Rscript

# to run
# Rscript INLA_simAnalysis_faster.R /home/fiona_callahan/simData/testing/run1/ /home/fiona_callahan/simData/testing/INLAres/

source("/home/fiona_callahan/inla_code/INLA_ST_func.R")
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input and output files need to be supplied", call. = FALSE)
} 

#data_dir <- "/home/fiona_callahan/simData/testing/run1/"
#save_dir <- "/home/fiona_callahan/simData/testing/INLAres/"
data_dir <- args[1]
save_dir <- args[2]
dir.create(save_dir)

proj <- FALSE
numTrials <- 1

# load data
sim_data_raw <- readRDS(paste0(data_dir, "sim_data.Rdata"))
sitetab <- read.csv(paste0(data_dir, "sim_sitetab_sampled.csv"))
#locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))

# note this does not include constant covs because those are just intercepts
names_cov <- params$names_cov 
names_species <- params$names_species

# run analysis three times on the same data
for (trial in 1:numTrials) {
  subdir <- paste0(save_dir, "trial", trial, "/")
  dir.create(subdir)
  # Select what you'll be working with
  k <- 11 # k is number of time knots
  #k <- 4
  # Prepare spatial and temporal meshes
  # alltidx <- seq(1,k)
  #reducedtidx <- seq(2,k,2)
  tknots <- seq(min(sitetab$Age), max(sitetab$Age), length = k)
  #reducedtknots <- tknots[seq(2,k,2)]
  mesh_t <- inla.mesh.1d(tknots)

  xytab <- cbind(sitetab$Long, sitetab$Lat)
  #xytab <- LonLatToAzi(xytab)
  n <- dim(sitetab)[1]
  
  bound <- inla.nonconvex.hull(as.matrix(xytab))
  #TODO what are these params??
  #mesh_s <- inla.mesh.2d(boundary = bound, max.edge = c(200, 400), max.n=200, cutoff = 500)
  mesh_s <- inla.mesh.2d(boundary = bound, max.edge = c(10, 25))#, max.edge = c(2, 4), max.n=200, cutoff = 500)
  spde_s <- inla.spde2.matern(mesh_s)
  
  #reslist1<-reslist
  # Define lists to store results
  # ADD MODELS
  res_all <- list()
  
  for (response in names_species){ 
    # Run with other animals as covariates
    start_time <- Sys.time()
    all_names_cov <- c(names_cov, names_species)
    all_names_cov <- all_names_cov[all_names_cov != response]
    res_all[[response]] <- #withTimeout({
           RunInlaBin(response, sitetab, xytab, mesh_s, mesh_t, 1, 
                              namescov = all_names_cov, normcov = 0, project = proj)
            #}, timeout = 2 * 3600, onTimeout = "warning")
    total_time <- (as.numeric(Sys.time()) - as.numeric(start_time)) / 3600 # time elapsed in hours
  }
  
  sim_lists <- list()
  sim_lists[["Cov+Animal"]] <- res_all
  
  #save results for later
  saveRDS(sim_lists, paste0(subdir, "sim_lists-iid-", trial, ".Rdata"))
  run_data <- list()
  run_data$params <- params
  run_data$sim_data_raw <- sim_data_raw
  run_data$sitetab <- sitetab
  run_data$time <- total_time

  saveRDS(run_data, paste0(subdir, "rundata-iid-", trial, ".Rdata"))
  # CPO 
  modelcands <- names(sim_lists)

  # get inferred alpha and beta (just 1's -1's and 0's )
  inferenceRes <- list()
  inferenceRes$time <- total_time
  inferenceRes$alphaInferred <- matrix(data = 0, nrow = length(names_species), ncol = length(names_species))
  inferenceRes$betaInferred <- matrix(data = 0, nrow = length(names_species), ncol = params$numCovs)
  for (sp_index in seq_along(names_species)) {
    this_best_model <- "Cov+Animal"
    best_model_res <- sim_lists[[this_best_model]][[names_species[sp_index]]]
    significant_pos_vars <- row.names(best_model_res$summary.fixed)[best_model_res$summary.fixed$`0.025quant` > 0]
    significant_neg_vars <- row.names(best_model_res$summary.fixed)[best_model_res$summary.fixed$`0.975quant` < 0]
    for (sp_index2 in seq_along(names_species)) { # put 1's and -1's into alpha and beta inferred based on 95% CI
      if (names_species[sp_index2] %in% significant_pos_vars) {inferenceRes$alphaInferred[sp_index, sp_index2] <- 1}
      if (names_species[sp_index2] %in% significant_neg_vars) {inferenceRes$alphaInferred[sp_index, sp_index2] <- -1}
    }
    for (cov_index in seq_along(names_cov)) {
      if (names_cov[cov_index] %in% significant_pos_vars) {inferenceRes$betaInferred[sp_index, cov_index] <- 1}
      if (names_cov[cov_index] %in% significant_neg_vars) {inferenceRes$betaInferred[sp_index, cov_index] <- -1}
    }
  }

  saveRDS(inferenceRes, paste0(subdir, "inferenceRes.Rdata"))
}
