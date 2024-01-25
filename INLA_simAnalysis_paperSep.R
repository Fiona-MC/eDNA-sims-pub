#!/usr/bin/env Rscript

# to run
# Rscript INLA_simAnalysis_paperSep.R /space/s1/fiona_callahan/multiSim_rw3/randomRun1/ /space/s1/fiona_callahan/multiSim_rw3/randomRun1/INLAres/

# NOTE this is the original unedited version!!!
source("/home/fiona_callahan/Arctic_eDNA_2021/script/INLA_ST_functions.R")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input and output files need to be supplied", call. = FALSE)
} 

#data_dir <- "/home/fiona_callahan/simData/testing/run1/"
#save_dir <- "/home/fiona_callahan/simData/testing/INLAres/"
#data_dir <- "/space/s1/fiona_callahan/multiSim_rw3/randomRun1/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_rw3/randomRun1/INLAres/"
#sitetabName <- "sim_sitetab_sampled.csv"
#${folder}/ ${folder}/INLA_res_${INLA_type}/#
data_dir <- args[1]
save_dir <- args[2]
#dir.create(save_dir)
sitetabName <- args[3]
modelParms <- args[4] # "none" "cov" "sp" "spCov"


plot <- FALSE
proj <- FALSE
numTrials <- 1

# load data
# sim_data_raw <- readRDS(paste0(data_dir, "sim_data.Rdata"))
sitetab <- read.csv(paste0(data_dir, sitetabName))
#locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))


names_cov <- params$names_cov
names_species <- params$names_species

# run analysis three times on the same data
for (trial in 1:numTrials){
  subdir <- paste0(save_dir, "trial", trial, "/")
  dir.create(subdir)
  # Select what you'll be working with
  k <- 11 # k is number of time knots
  #k <- 4
  # Prepare spatial and temporal meshes
  # alltidx <- seq(1,k)
  #reducedtidx <- seq(2,k,2)
  tknots <- seq(min(sitetab$Age), max(sitetab$Age), length = k)
  #edge case, this is probably the wrong way to handle this
  #(fixes bug if only one time pt)
  if(min(sitetab$Age) == max(sitetab$Age)) {
    tknots <- seq(min(sitetab$Age) - 1, max(sitetab$Age), length = k)
  }
  #reducedtknots <- tknots[seq(2,k,2)]
  mesh_t <- inla.mesh.1d(tknots)

  xytab <- cbind(sitetab$Long, sitetab$Lat)
  #xytab <- LonLatToAzi(xytab)
  n <- dim(sitetab)[1]
  
  bound <- inla.nonconvex.hull(as.matrix(xytab))
  #TODO what are these params??
  #mesh_s <- inla.mesh.2d(boundary = bound, max.edge = c(200, 400), max.n=200, cutoff = 500)
  mesh_s <- inla.mesh.2d(boundary = bound, max.edge = c(10, 25))#, max.edge = c(2, 4), max.n=200, cutoff = 500)
  plot(mesh_s)
  spde_s <- inla.spde2.matern(mesh_s)
  
  
  
  #goodanimalnames<-c("Sp1", "Sp2")#, "Sp3")#, "Sp4", "Sp5")
  # this takes a long time even with only 10 time points
  # Model with no covariates
  reslist <- list()
  if (modelParms == "none") {
    all_names_cov <- c()
    norm_cov <- 0
  } else if (modelParms == "sp") {
    all_names_cov <- names_species
    norm_cov <- 2
  } else if (modelParms == "cov") {
    all_names_cov <- names_cov
    norm_cov <- 2
  } else if (modelParms == "spCov") {
    all_names_cov <- c(names_cov, names_species)
    norm_cov <- 2
  } else {
    stop("modelParms parameter is not a valid value")
  }

  for (response in names_species) { # NOTE -- this takes a few minutes
    names_cov_temp <- all_names_cov[all_names_cov != response]
    # Run INLA - Binomial model with number of trials = 1
    reslist[[response]] <- RunInlaBin(response, sitetab, xytab, mesh_s, mesh_t, 1, 
                                namescov = names_cov_temp, normcov = norm_cov, project = proj)
  }
  
  #save results for later
  saveRDS(reslist, paste0(subdir, "resList_", modelParms, ".Rdata"))

  run_data <- list()
  run_data$params <- params
  run_data$modelParms <- modelParms
  # run_data$sim_data_raw <- sim_data_raw
  run_data$sitetab <- sitetab

  saveRDS(run_data, paste0(subdir, "rundata", trial, ".Rdata"))
}
