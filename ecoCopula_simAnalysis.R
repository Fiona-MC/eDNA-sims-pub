#!/usr/bin/env Rscript
library(ecoCopula)

# to run
# Rscript ecoCopula_simAnalysis.R /space/s1/fiona_callahan/multiSim11a/randomRun10/ /space/s1/fiona_callahan/multiSim11a/randomRun10/ecoCopula_res/ 0

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input and output files need to be supplied", call. = FALSE)
} 

data_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun2/"
save_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun2/ecoCopula_res_noCov/"

data_dir <- args[1]
save_dir <- args[2]
dir.create(save_dir)
scramble <- FALSE
scramble <- (as.numeric(args[3]) == 1)
sitetab_name <- "sim_sitetab_sampled.csv"
sitetab_name <- args[4]


plot <- FALSE
numTrials <- 1
prec <- TRUE

#sitetab <- read.csv(paste0(data_dir, "sitetab_dumb.csv"))
#sitetab <- read.csv(paste0(data_dir, "sitetab_dumb_dir.csv"))
sitetab <- read.csv(paste0(data_dir, sitetab_name))
#locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))
if (scramble == TRUE) {
  sitetab <- read.csv(paste0(data_dir, sitetab_name))
}

names_cov <- params$names_cov
names_species <- params$names_species
cov <- FALSE


for (trial in 1:numTrials) {
  # example from documentation
  # https://cran.r-project.org/web/packages/ecoCopula/vignettes/the_basics.html
  subdir <- paste0(save_dir, "trial", trial, "/")
  dir.create(subdir)
  #formula <- as.formula(paste(" ~ ", paste(names_cov, collapse = "+")))
  if(cov) {
    formula <- as.formula(paste(" ~ ", paste(names_cov, collapse = "+")))
    fit0 <- stackedsdm(sitetab[, names_species], formula_X = formula, 
                      data = sitetab[, names_cov], family = "binomial", ncores = 1) 
  } else {
    fit0 <- stackedsdm(sitetab[, names_species], formula_X = ~1, 
                          data = data.frame(intercept = rep(1, times = dim(sitetab)[1])), family = "binomial", ncores = 1) 
  }
  # vvv changing these within reason did not reduce the number of inferred edges by much
  #cgr_sim <- cgr(fit0, n.lambda = 100, n.samp = 1000, method = "BIC")
  cgr_sim <- cgr(fit0)

  saveRDS(cgr_sim, file = paste0(subdir, "EC_res.Rdata"))

  #plot(cgr_sim, pad = 1)
 
  if (prec == TRUE) {
    inferred_mx <- -1 * cgr_sim$best_graph$prec
  } else {
    inferred_mx <- cgr_sim$best_graph$cov
  }

  # get inferred alpha and beta (just 1's -1's and 0's )
  inferenceRes <- list()
  inferenceRes$alphaInferred <- sign(inferred_mx)
  inferenceRes$betaInferred <- NA
  # put 0s on diagonal
  inferenceRes$alphaInferred <- inferenceRes$alphaInferred * as.integer(diag(nrow = params$numSpecies, ncol = params$numSpecies) == 0)

  saveRDS(inferenceRes, paste0(subdir, "inferenceRes.Rdata"))
}
