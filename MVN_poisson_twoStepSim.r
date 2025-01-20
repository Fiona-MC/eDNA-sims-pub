library(MASS) # For generating multivariate normal random variables
library(ggplot2)

generate_poisson_twoStep <- function(n, mvnMean, target_covar, factor) {
  # n: Number of samples
  # mvnMean: mean of multivariate normal -- if less than 0, there will be fewer 0 counts.
  # target_covar: Desired covariance matrix
  
  if (min(eigen(target_covar)$values) <= 0) {
    stop("Target covariance matrix must be positive definite.")
  }

  if (dim(mvnMean)[1] != n) {
    print("dim(mvnMean):")
    print(dim(mvnMean))
    stop("mvnMean should be n by numSpecies")
  }

  if (dim(mvnMean)[2] != dim(target_covar)[1]) {
    print("dim(mvnMean):")
    print(dim(mvnMean))
    print("dim(target_covar)")
    print(dim(target_covar))
    stop("mvnMean should be n by numSpecies and dim(target_covar) should be numSpecies by numSpecies")
  }

  if (dim(target_covar)[1] != dim(target_covar)[2]) {
    print("dim(target_covar)")
    print(dim(target_covar))
    stop("dim(target_covar) should be numSpecies by numSpecies")
  }
  
  # Generate MVN variables
  num_species <- dim(target_covar)[1]

  mvn_samples <- mvnMean + mvrnorm(n, mu = rep(0, times = num_species), Sigma = target_covar) # Multivariate normal
  mvn_truncated <- mvn_samples * (mvn_samples > 0) * factor
  
  # Use truncated MVN as poisson means
  poisson_samples <- sapply(1:num_species, function(i) {
    rpois(n = dim(mvn_truncated)[1], mvn_truncated[, i])
  })

  return(list(mvn_samples = mvn_samples, mvn_truncated = mvn_truncated, poisson_samples = poisson_samples))
}

sigma_wilkinson <- function(nspecies, block_sizes) {
  # implementing the simulation model from Wilkinson (in supplement 5)
  # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13106
  stopifnot(sum(block_sizes) <= nspecies)
  
  sigma_final <- matrix(0, nrow = nspecies, ncol = nspecies)
  blockstart <- 1
  for (blocksize in block_sizes) {
    factors <- 3
    lambda <- matrix(rnorm(n = blocksize * factors), nrow = blocksize, ncol = factors)

    # low rank correlation matrix
    sigma <- (lambda %*% t(lambda)) + diag(1, nrow = blocksize)

    # rescale to correlation mx
    D_inv <- diag(1 / sqrt(diag(sigma)))
    sigma_corr <- D_inv %*% sigma %*% D_inv

    # add interacting species
    block_location <- seq(blockstart, blockstart + blocksize - 1)
    sigma_final[block_location, block_location] <- sigma_corr
    blockstart <- blockstart + blocksize # start for the next block is right after this block ends
  }

  diag(sigma_final) <- 1

  return(sigma_final)
}

mu_wilkinson <- function(n, nspecies, ncov = 4) {
  # standard normal covariates (plus intercept)
  covs <- matrix(c(rnorm(n = n * ncov), rep(1, times = n)), 
                  nrow = n, ncol = ncov + 1, byrow = FALSE)
  beta <- matrix(rnorm(n = (ncov + 1) * nspecies), nrow = ncov + 1, ncol = nspecies)
  mu <- covs %*% beta

  return(list(covs = covs, beta = beta, mu = mu))
}

testRun <- function(nspecies) {
  # simulate reads from covariance matrix
  n <- 25000
  ncov <- 4 # 4 is from wilkinson

  ninteract <- ceiling(sqrt(nspecies))
  corrMx <- sigma_wilkinson(nspecies, ninteract = ninteract)
  sim_means_L <- mu_wilkinson(n, nspecies, ncov = ncov)
  mvnMean <- sim_means_L$mu
  #mvnMean <- matrix(rep(0, times = n * nspecies), nrow = n, ncol = nspecies)
  samples <- generate_poisson_twoStep(n, mvnMean, corrMx, factor = 100)
  poisson_samples <- samples$poisson_samples
  # Check empirical covariance
  #sample_cov <- cov(poisson_samples)
  sample_cor <- cor(poisson_samples)
  sample_cor <- cor(samples$mvn_truncated)
  sample_cor <- cor(samples$mvn_samples)
  #print("Sample Covariance Matrix:")
  #print(sample_cov)
  #print("Sample Correlations:")
  #print(sample_cor)

  cor_interact <- sample_cor[1:ninteract, 1:ninteract]
  #hist(cor_interact[upper.tri(cor_interact)])

  cor_noInteract <- sample_cor[(ninteract + 1):nspecies, (ninteract + 1):nspecies]
  #hist(cor_noInteract[upper.tri(cor_noInteract)])

  target_cor_interact <- corrMx[1:ninteract, 1:ninteract]
  target_cor_noInteract <- corrMx[(ninteract + 1):nspecies, (ninteract + 1):nspecies]

  return(list(target_cor_interact = target_cor_interact, target_cor_noInteract = target_cor_noInteract, 
            cor_interact = cor_interact, cor_noInteract = cor_noInteract))
}

if (FALSE) { # don't run on source()
noInteract <- c()
interact <- c()
target_noInteract <- c()
target_interact <- c()
for (i in 1:10) {
  test <- testRun(100)
  noInteract <- c(noInteract, test$cor_noInteract[upper.tri(test$cor_noInteract)])
  interact <- c(interact, test$cor_interact[upper.tri(test$cor_interact)])

  target_noInteract <- c(target_noInteract, test$target_cor_noInteract[upper.tri(test$target_cor_noInteract)])
  target_interact <- c(target_interact, test$target_cor_interact[upper.tri(test$target_cor_interact)])
}

group <- c(rep("empirical_interact", times = length(interact)),
          rep("empirical_independent", times = length(noInteract)),
          rep("target_interact", times = length(target_interact)),
          rep("target_independent", times = length(target_noInteract)))

correlation <- c(interact, noInteract, target_interact, target_noInteract)

empirical_v_target <- data.frame(group, correlation)

ggplot(empirical_v_target, aes(x = correlation, fill = group)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~group, scales = "free_y") +
  theme_minimal() +
  labs(title = "Histograms of Correlations by Group",
       x = "Correlation",
       y = "Frequency") +
  theme(legend.position = "none")
}