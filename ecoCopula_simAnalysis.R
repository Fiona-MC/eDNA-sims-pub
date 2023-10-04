#!/usr/bin/env Rscript
library(ecoCopula)

# to run
# Rscript ecoCopula_simAnalysis.R /space/s1/fiona_callahan/multiSim11a/randomRun10/ /space/s1/fiona_callahan/multiSim11a/randomRun10/ecoCopula_res/ 0

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input and output files need to be supplied", call. = FALSE)
} 

#data_dir <- "/home/fiona_callahan/simData/testing/run1/"
#save_dir <- "/home/fiona_callahan/simData/testing/INLAres/"
data_dir <- "/space/s1/fiona_callahan/multiSim11a/randomRun30/"
save_dir <- "/space/s1/fiona_callahan/multiSim11a/randomRun30/ecoCopula_res/"
#${folder}/ ${folder}/INLA_res_${INLA_type}/#
data_dir <- args[1]
save_dir <- args[2]
dir.create(save_dir)
scramble <- (as.numeric(args[3]) == 1)
scramble <- FALSE

plot <- FALSE
numTrials <- 2

# load data
sim_data_raw <- readRDS(paste0(data_dir, "sim_data.Rdata"))
sitetab <- read.csv(paste0(data_dir, "sim_sitetab_sampled.csv"))
#locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))
if (scramble == TRUE) {
  sitetab <- read.csv(paste0(data_dir, "sitetab_scrambled.csv"))
}

names_cov <- params$names_cov
names_species <- params$names_species

# run analysis three times on the same data
for (trial in 1:numTrials){
  # example from documentation
  # https://cran.r-project.org/web/packages/ecoCopula/vignettes/the_basics.html

  fit0 <- stackedsdm(sitetab[100:200, names_species], formula_X = ~ Cov1 + Cov2 + Cov3, 
                    data = sitetab[100:200, names_cov], family = "binomial", ncores = 1) 
  cgr_sim <- cgr(fit0)
  plot(cgr_sim, pad = 1)
 


  # get inferred alpha and beta (just 1's -1's and 0's )
  inferenceRes <- list()
  inferenceRes$alphaInferred <- matrix(data = 0, nrow = length(names_species), ncol = length(names_species))
  inferenceRes$betaInferred <- matrix(data = 0, nrow = length(names_species), ncol = params$numCovs)
  for (sp_index in seq_along(names_species)) {
    this_best_model <- modelcands[which.min(waictab[, sp_index])]
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
