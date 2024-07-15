library(mvtnorm)
#install.packages("highfrequency")
library(data.table)
library(highfrequency)

args <- commandArgs(trailingOnly = TRUE)
# Rscript logiSim_covMx_higherCorr.R /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples/randomRun1/


# bash code to do this for all subfolders within this folder
#sim_dir="/space/s1/fiona_callahan/multiSim_100sp_logi0.5"
#mkdir $sim_dir
#for folder in ${sim_dir}/randomRun{1..100}; do (Rscript logiSim_covMx_higherCorr.r ${folder}/) done

save_dir <- "/space/s1/fiona_callahan/multiSim_100sp_logi_0.5/randomRun1/"
save_dir <- args[1]
numSpecies <- 100
absRho <- 0.5

dir.create(save_dir)

actualAlpha <- matrix(data = 0, nrow = numSpecies, ncol = numSpecies)
corrMx <- diag(numSpecies)

for (pair in 1:(numSpecies / 2)) {
    rho <- sample(c(absRho, -absRho), size = 1)
    corrMx[pair, (numSpecies / 2) + pair] <- rho
    corrMx[(numSpecies / 2) + pair, pair] <- rho

    actualAlpha[pair, (numSpecies / 2) + pair] <- sign(rho)
    actualAlpha[(numSpecies / 2) + pair, pair] <- sign(rho)
}

prec <- FALSE
sim_covs <- TRUE

params <- list()
params$alpha <- actualAlpha
params$names_species <- sapply(1:numSpecies, FUN = function(x) {paste0("Sp", x)})
params$numSpecies <- numSpecies
params$numCovs <- numSpecies + 1
params$beta <- matrix(0, nrow = numSpecies, ncol = params$numCovs)
params$names_cov <- sapply(1:params$numCovs, FUN = function(x) {paste0("Cov", x)})
params$covVars <- list()

for(i in 1:(params$numCovs - 1)) {
    params$covVars[[i]] <- list()
    params$covVars[[i]]$type <- "rnorm"
}
params$covVars[[params$numCovs]] <- list()
params$covVars[[params$numCovs]]$type <- "constant"

saveRDS(params, file = paste0(save_dir, "params.Rdata"))

# simulate reads from covariance matrix
samples <- 25000

eig <- eigen(corrMx)
min(eig$values)

if (prec) {
    covarMx <- - solve(corrMx)
} else {
    covarMx <- corrMx
}
#covarMx <- makePsd(corrMx, method = "covariance")

#sum(covarMx == 0)
#sum(covarMx != 0)
#sum(abs(covarMx) < 0.001)
#sum(abs(covarMx) > 0.001)

##sum(actualAlpha != 0)
#sum(precMx != 0)
#96 * 2 + 100

# this is not going to work because covarMx is not positive semidefinite -- 
# need to think of anohter way to get this mx
reads <- rmvnorm(n = samples, mean = rep(0, times = numSpecies), sigma = covarMx)
if (prec) {
    saveRDS(covarMx, file = paste0(save_dir, "logi_covar_mx_prec.Rdata"))
} else {
   saveRDS(covarMx, file = paste0(save_dir, "logi_covar_mx.Rdata"))
}
# truncate at 0
reads <- reads * (reads > 0)
# turn into reasonable numbers 
reads <- reads * 100
reads <- round(reads)

#sim random normal covariates
sitetab_abd <- data.frame(site = 1:samples, reads)
names(sitetab_abd) <- c("site", sapply(1:numSpecies, FUN = function(x) {paste0("Sp", x)}))

sitetab_abd$Age <- sample(x = 1:50000, size = samples, replace = TRUE)
sitetab_abd$Lat <- runif(n = samples, min = 0, max = 100)
sitetab_abd$Long <- runif(n = samples, min = 0, max = 100)

# truncate to presence absence 
presAbs <- (reads > 5) * 1

sitetab <- data.frame(site = 1:samples, presAbs)
names(sitetab) <- c("site", sapply(1:numSpecies, FUN = function(x) {paste0("Sp", x)}))

sitetab$Age <- sitetab_abd$Age
sitetab$Lat <- sitetab_abd$Lat
sitetab$Long <- sitetab_abd$Long

for (i in 1:(params$numCovs - 1)) {
    thisCov <- rnorm(n = samples)
    sitetab_abd[[paste0("Cov", i)]] <- thisCov
    sitetab[[paste0("Cov", i)]] <- thisCov
}
sitetab_abd[[paste0("Cov", params$numCovs)]] <- rep(x = 1, times = samples)
sitetab[[paste0("Cov", params$numCovs)]] <- rep(x = 1, times = samples)

if (prec) {
    fwrite(x = sitetab_abd, file = paste0(save_dir, "sitetab_abd_logi_prec.csv"))
} else {
    fwrite(x = sitetab_abd, file = paste0(save_dir, "sitetab_abd_logi.csv"))
}

if (prec) {
    fwrite(x = sitetab, file = paste0(save_dir, "sitetab_logi_prec.csv"))
} else {
    fwrite(x = sitetab, file = paste0(save_dir, "sitetab_logi.csv"))
}
