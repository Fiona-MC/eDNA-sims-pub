library(mvtnorm)
#install.packages("highfrequency")
library(data.table)
library(highfrequency)

args <- commandArgs(trailingOnly = TRUE)
# Rscript dumbSim_covMx.R /space/s1/fiona_callahan/multiSim_10sp_indep/randomRun1/

# bash code to do this for all subfolders within this folder
#sim_dir="/space/s1/fiona_callahan/multiSim_10sp_dep"
#for folder in ${sim_dir}/randomRun*; do (Rscript dumbSim_covMx.R ${folder}/) done

data_dir <- "/space/s1/fiona_callahan/multiSim_10sp_dep/randomRun1/"
data_dir <- args[1]
params <- readRDS(paste0(data_dir, "params.Rdata"))

actualAlpha <- params$alpha
prec <- FALSE
sim_covs <- TRUE

# simulate reads from covariance matrix
samples <- params$num_samples_time * params$num_samples_space

sigma <- 1 # variance of read counts [per species]
rho <- 0.45 # covariance of species with interaction 

corrMx <- matrix(0, nrow = params$numSpecies, ncol = params$numSpecies)
for (i in 1:params$numSpecies) {
    for (j in 1:params$numSpecies) {
        if (actualAlpha[i, j] < 0) {
            corrMx[i, j] <- -rho
            corrMx[j, i] <- -rho
        } else if (actualAlpha[i, j] > 0) {
            corrMx[i, j] <- rho
            corrMx[j, i] <- rho
        } 
    }
}
corrMx <- corrMx + diag(x = 1, nrow = params$numSpecies, ncol = params$numSpecies)
precMx <- corrMx

eig <- eigen(precMx)
min(eig$values)

if (prec) {
    covarMx <- solve(precMx)
} else {
    covarMx <- corrMx
}
#covarMx <- makePsd(corrMx, method = "covariance")

sum(covarMx == 0)
sum(covarMx != 0)
sum(abs(covarMx) < 0.001)
sum(abs(covarMx) > 0.001)

sum(actualAlpha != 0)
sum(precMx != 0)
#96 * 2 + 100

# this is not going to work because covarMx is not positive semidefinite -- 
# need to think of anohter way to get this mx
reads <- rmvnorm(n = samples, mean = rep(0, times = params$numSpecies), sigma = covarMx)
if (prec) {
    saveRDS(covarMx, file = paste0(data_dir, "dumb_covar_mx_prec.Rdata"))
} else {
   saveRDS(covarMx, file = paste0(data_dir, "dumb_covar_mx.Rdata"))
}
# truncate at 0
reads <- reads * (reads > 0)
# turn into reasonable numbers 
reads <- reads * 100
reads <- round(reads)

#sim random normal covariates


# df of read counts
if (sim_covs) {
    covVals <- matrix(data = rnorm(n = samples * (params$numCovs - 1)), nrow = samples, ncol = (params$numCovs - 1))
    sitetab_abd <- data.frame(site = 1:samples, reads, covVals)
    names(sitetab_abd) <- c("site", params$names_species, params$names_cov)
} else {
    sitetab_abd <- data.frame(site = 1:samples, reads)
    names(sitetab_abd) <- c("site", params$names_species)
}

if (prec) {
    fwrite(x = sitetab_abd, file = paste0(data_dir, "sitetab_abd_dumb_prec.csv"))
} else {
    fwrite(x = sitetab_abd, file = paste0(data_dir, "sitetab_abd_dumb.csv"))
}

# truncate to presence absence 
presAbs <- (reads > params$readThreshold) * 1
if (sim_covs) {
    sitetab <- data.frame(site = 1:samples, presAbs, covVals)
    names(sitetab) <- c("site", params$names_species, params$names_cov)
} else {
    sitetab <- data.frame(site = 1:samples, presAbs)
    names(sitetab) <- c("site", params$names_species)
}

if (prec) {
    fwrite(x = sitetab, file = paste0(data_dir, "sitetab_dumb_prec.csv"))
} else {
    fwrite(x = sitetab, file = paste0(data_dir, "sitetab_dumb.csv"))
}
