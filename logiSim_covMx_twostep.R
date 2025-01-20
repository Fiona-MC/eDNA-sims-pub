library(mvtnorm)
#install.packages("highfrequency")
library(data.table)
library(highfrequency)

# this is the code that was used for revisions 2 and 3

source("/home/fiona_callahan/eDNA_sims_code/revisions_code/MVN_poisson_twoStepSim.r")

args <- commandArgs(trailingOnly = TRUE)
# Rscript ./revisions_code/logiSim_covMx_twostep.R /space/s1/fiona_callahan/sim_paper_stuff/multiSim_100sp_revision_test3/randomRun1/ 100

data_dir <- "/space/s1/fiona_callahan/sim_paper_stuff/multiSim_100sp_revision_test/randomRun1/"
nspecies <- 100

data_dir <- args[1]
nspecies <- as.numeric(args[2])
no_covs <- FALSE

# simulate reads from covariance matrix
samples <- 25000
ncov <- 4 # 4 is from wilkinson

if (nspecies == 10) {
    block_sizes <- c(4)
} else if (nspecies == 100) {
    block_sizes <- c(5, 5, 5, 5, 5)
} else {
    block_sizes <- c(floor(sqrt(nspecies)))
}

ninteract <- ceiling(sqrt(nspecies))
corrMx <- sigma_wilkinson(nspecies, block_sizes = block_sizes)
# simulate and output the covariates and covariate effects, resulting in means per species
sim_means_L <- mu_wilkinson(samples, nspecies, ncov = ncov)
mvnMean <- sim_means_L$mu
beta <- sim_means_L$beta
covs <- sim_means_L$covs

if (no_covs) {
    beta <-  beta * 0
    mvnMean <- mvnMean * 0
}

actualAlpha <- sign(corrMx)
diag(actualAlpha) <- 0

params <- list()
params$alpha <- actualAlpha
params$names_species <- sapply(1:nspecies, FUN = function(x) {paste0("Sp", x)})
params$numSpecies <- nspecies
params$numCovs <- ncov + 1
params$beta <- t(beta)
params$names_cov <- sapply(1:params$numCovs, FUN = function(x) {paste0("Cov", x)})
params$covVars <- list()

for(i in 1:(params$numCovs - 1)) {
    params$covVars[[i]] <- list()
    params$covVars[[i]]$type <- "rnorm"
}
params$covVars[[params$numCovs]] <- list()
params$covVars[[params$numCovs]]$type <- "constant"

saveRDS(params, file = paste0(data_dir, "params.Rdata"))
saveRDS(list(covarMx = corrMx, beta = beta), file = paste0(data_dir, "rev_covar_mx.Rdata"))

factor <- 100 # scale up the MVN samples so that number of reads has a higher mean

pois_samples <- generate_poisson_twoStep(samples, mvnMean = mvnMean, corrMx, factor)
reads <- pois_samples$poisson_samples
presAbs <- (pois_samples$mvn_samples > 0) * 1 # probit samples

# df of read counts
sitetab_abd <- data.frame(site = 1:samples, reads)
names_species <- sapply(1:nspecies, FUN = function(i) {paste0("Sp", i)})
names(sitetab_abd) <- c("site", names_species)

for (i in 1:(ncov + 1)) {
    sitetab_abd[[paste0("Cov", i)]] <- covs[, i]
}

sitetab_abd$Age <- sample(x = 1:10000, size = samples, replace = TRUE)
sitetab_abd$Lat <- runif(n = samples, min = 0, max = 100)
sitetab_abd$Long <- runif(n = samples, min = 0, max = 100)

fwrite(x = sitetab_abd, file = paste0(data_dir, "sitetab_abd_logi.csv"))

# df of pres abs
sitetab <- data.frame(site = 1:samples, presAbs)
names(sitetab) <- c("site", names_species)

for (i in 1:(ncov + 1)) {
    sitetab[[paste0("Cov", i)]] <- covs[, i]
}

sitetab$Age <- sitetab_abd$Age
sitetab$Lat <- sitetab_abd$Lat
sitetab$Long <- sitetab_abd$Long

fwrite(x = sitetab, file = paste0(data_dir, "sitetab_logi.csv"))