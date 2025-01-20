# 1/17/25 obsidian
# make sure the poisson values have the correct means (after drawing x and $\beta$)
library(truncnorm)
library(mvtnorm)
#install.packages("highfrequency")
library(data.table)
library(highfrequency)

# this is the code that was used for revisions 2 and 3

source("/home/fiona_callahan/eDNA_sims_code/revisions_code/MVN_poisson_twoStepSim.r")

# Rscript ./revisions_code/logiSim_covMx_twostep.R /space/s1/fiona_callahan/sim_paper_stuff/multiSim_100sp_revision_test3/randomRun1/ 100

data_dir <- "/space/s1/fiona_callahan/sim_paper_stuff/multiSim_10sp_revision_test/randomRun1/"
nspecies <- 10
no_covs <- FALSE

ntrials <- 1000
# simulate reads from covariance matrix
samples <- 10000
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

for (trial in 1:ntrials) {
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

    fwrite(x = sitetab_abd, file = paste0(data_dir, "sitetab_abd_logi", trial, ".csv"))

}






########################################################################

params <- readRDS(paste0(data_dir, "params.Rdata"))
sitetab1 <- read.csv(paste0(data_dir, "sitetab_abd_logi", 1, ".csv"))

beta <- t(params$beta)
sampled_X <- t(as.matrix(sitetab1[, c("Cov1", "Cov2", "Cov3", "Cov4", "Cov5")]))

z_means <- t(sampled_X) %*% beta 

# are the means the same as X * \beta
all(abs(z_means - mvnMean) < 10^8)

# Calculate the mean of the truncated normal
truncated_means <- matrix(etruncnorm(a = 0, mean = z_means, sd = 1), 
                          nrow = nrow(z_means), 
                          ncol = ncol(z_means))

poisson_expectations <- truncated_means * 100

## get average of actual draws
# ntrials <- 100
running_sum <- matrix(0, nrow = 10000, ncol = 10)
n_nonzero <- matrix(0, nrow = 10000, ncol = 10)
for (trial in 1:ntrials) {
    sitetab <- read.csv(paste0(data_dir, "sitetab_abd_logi", trial, ".csv"))
    names_species <- sapply(1:10, FUN = function(x) {paste0("Sp", x)})
    actual_values <- as.matrix(sitetab[, names_species])

    running_sum <- running_sum + actual_values
    n_nonzero <- n_nonzero + (actual_values != 0)
}

avg_actual <- running_sum / n_nonzero

mean(abs(poisson_expectations - avg_actual), na.rm = TRUE)

mean(poisson_expectations - avg_actual, na.rm = TRUE)

mean(poisson_expectations)
mean(avg_actual, na.rm = TRUE)
