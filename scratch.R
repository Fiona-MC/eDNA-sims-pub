library(MASS)
library(mvtnorm)



    
    # Generate covariance matrix
    cov_matrix <- matrix(NA, nrow = nrow(locDF), ncol = nrow(locDF))
    for (locindex1 in seq_len(nrow(locDF))) {
      for (locindex2 in seq_len(nrow(locDF))) {
        loc1 <- as.numeric(locDF[locindex1, ])
        loc2 <- as.numeric(locDF[locindex2, ])
        cov_matrix[locindex1, locindex2] <- spatial_covariance(loc1, loc2, covar_scale_space = covar_scale_space)
      }
    }

mvnorm_inverse_covariance <- function(n_samples = 1, mean, inv_covariance) {
  # Check if the inverse covariance matrix is symmetric positive definite
  if (!all(eigen(inv_covariance)$values > 0)) {
    stop("Inverse covariance matrix must be symmetric positive definite.")
  }

  # Cholesky decomposition of the inverse covariance matrix
  chol_decomp <- chol(inv_covariance)

  # Generate standard normal samples
  standard_normal_samples <- matrix(rnorm(length(mean) * n_samples), nrow = n_samples)

  # Transform the samples to obtain multivariate normal samples
  samples <- t(chol_decomp %*% t(standard_normal_samples)) + rep(mean, each = n_samples)

  return(samples)
}

MASS.time <- sapply(seq(500, 2500, 500), function(x) {
  system.time(MASS::mvrnorm(100, rnorm(x), diag(x)))[3]  
})

mvtnorm.time <- sapply(seq(500, 2500, 500), function(x) {
  system.time(mvtnorm::rmvnorm(100, rnorm(x), diag(x)))[3]  
})

inv_cov_time <- sapply(seq(500, 2500, 500), function(x) {
  system.time(draw_mvnorm_inverse_covariance(100, rnorm(x), diag(x)))[3]  
})

MASS.time
mvtnorm.time
inv_cov_time