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
#data_dir <- "/space/s1/fiona_callahan/multiSim_5sp/randomRun22/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_5sp/randomRun22/ecoCopula_res/"
#${folder}/ ${folder}/INLA_res_${INLA_type}/#
data_dir <- args[1]
save_dir <- args[2]
dir.create(save_dir)
scramble <- FALSE
scramble <- (as.numeric(args[3]) == 1)


plot <- FALSE
numTrials <- 1
prec <- FALSE

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


for (trial in 1:numTrials) {
  # example from documentation
  # https://cran.r-project.org/web/packages/ecoCopula/vignettes/the_basics.html
  subdir <- paste0(save_dir, "trial", trial, "/")
  dir.create(subdir)
  fit0 <- stackedsdm(sitetab[, names_species], formula_X = ~ Cov1 + Cov2 + Cov3 + Cov4 + Cov5, 
                    data = sitetab[, names_cov], family = "binomial", ncores = 1) 
  cgr_sim <- cgr(fit0)
  #plot(cgr_sim, pad = 1)
 
  if(prec == TRUE) {
    inferred_mx <- cgr_sim$best_graph$prec
  } else {
    inferred_mx <- cgr_sim$best_graph$cor
  }

  # get inferred alpha and beta (just 1's -1's and 0's )
  inferenceRes <- list()
  inferenceRes$alphaInferred <- sign(inferred_mx)
  inferenceRes$betaInferred <- NA
  # put 0s on diagonal
  inferenceRes$alphaInferred <- inferenceRes$alphaInferred * as.integer(diag(nrow = params$numSpecies, ncol = params$numSpecies) == 0)

  saveRDS(inferenceRes, paste0(subdir, "inferenceRes.Rdata"))

}
