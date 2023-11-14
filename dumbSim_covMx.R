library(mvtnorm)
#install.packages("highfrequency")
library(data.table)
library(highfrequency)

data_dir <- "/space/s1/fiona_callahan/multiSim_2sp_test/randomRun1/"
params <- readRDS(paste0(data_dir, "params.Rdata"))

actualAlpha <- params$alpha

# simulate reads from covariance matrix
samples <- params$num_samples_time * params$num_samples_space

sigma <- 1 # variance of read counts [per species]
rho <- 0.8 # covariance of species with interaction 

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

covarMx <- solve(precMx)
covarMx <- corrMx
#covarMx <- makePsd(corrMx, method = "covariance")

sum(covarMx == 0)
sum(covarMx != 0)
sum(abs(covarMx) < 0.001)
sum(abs(covarMx) > 0.001)

sum(actualAlpha != 0)
sum(precMx != 0)
96 * 2 + 100

# this is not going to work because covarMx is not positive semidefinite -- 
# need to think of anohter way to get this mx
reads <- rmvnorm(n = samples, mean = rep(0, times = params$numSpecies), sigma = covarMx)
# truncate at 0
reads <- reads * (reads > 0)
# turn into reasonable numbers 
reads <- reads * 100
reads <- round(reads)

# df of read counts
sitetab_abd <- data.frame(site = 1:samples, reads)
names(sitetab_abd) <- c("site", params$names_species)

fwrite(x = sitetab_abd, file = paste0(data_dir, "sitetab_abd_dumb_dir.csv"))

# truncate to presence absence 
presAbs <- (reads > params$readThreshold) * 1
sitetab <- data.frame(site = 1:samples, presAbs)
names(sitetab) <- c("site", params$names_species)

fwrite(x = sitetab, file = paste0(data_dir, "sitetab_dumb_dir.csv"))
