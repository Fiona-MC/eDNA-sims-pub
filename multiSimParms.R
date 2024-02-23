library(ggplot2)
library(reshape)
library(gridExtra)
library(jsonlite)


getParms_json <- function(parmFile = "parameters.json") { #TODO this is not implemented at all
  params <- fromJSON(parmFile)
  for(parm in params) {
    # check if random
    # sample if random
  }
  # do a bunch of the transformaation and function stuff here
  return(params)
}

#any number of species should work
getParms_many <- function(random = FALSE, parmSet = 1, numSpecies = 100, parmSetCov = "indep") { 
    if (random == TRUE) {
        num_samples_time <- sample(500, size = 1) # sample times per location
        num_samples_space <- sample(50, size = 1) # sample locations per time
        #radius <- 16
        radius <- runif(n = 1, min = 12, max = 50) # for neighborhoods: a little more the 
        # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
        #mean_fpr <- runif(n = 1, min = 0, max = 0.2)
        #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

        #fpr_mode <- sample(x = c("none", "independent", "dependent_sp", "constant"), size = 1)
        #fpr_mode <- "independent"
        covNoise_sd <- runif(n = 1, min = 0, max = 0.1) # process noise
        covMeasureNoise_sd <- runif(n = 1, min = 0, max = 0.2) # measurement noise

        # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
        # intrinsic growth rate
        r <- runif(n = 1, min = 0.01, max = 1) #todo think on reasonableness here
        # noise param for popn growth 
        sigma <- runif(n = 1, min = 0, max = 0.2)

        #N_50 <- sample(x = 1:1000, size = 1) # number of animals for which detection prob is 50%

        #species effect
        c2 <- runif(n = 1, min = 100, max = 500)

        mean_mig_rate <- runif(n = 1, min = 0, max = 0.1) # poisson rate per individual per time

        indivSampleProb <- runif(n = 1, 0.001, 0.02)
        readSampleRate <- runif(n = 1, 0.5, 100)
        corr_mode <- sample(c("independent")) # this is the only one implemented thus far
        numCovs <- sample(1:numSpecies, size = 1)

        covVars <- list()
        for (i in 1:numCovs) {
          period <- runif(1, min = 100, max = 10000)
          covVars[[i]] <- list(type = "spatialRandomField", period = period)
        }
        covar_scale_space <- runif(n = 1, min = 1, max = 100)

        covVars[[numCovs + 1]] <- list(type = "constant")

        numCovs <- length(covVars)
        names_cov <- sapply(1:numCovs, FUN = function(cov) {
                        if (covVars[[cov]]$type != "constant") {
                          return(paste0("Cov", cov))}}) 
        names_cov <- unlist(names_cov)
        names_all_cov <- unlist(sapply(1:numCovs, FUN = function(cov) {paste0("Cov", cov)}))
        names_species <- sapply(1:numSpecies, FUN = function(sp) {paste0("Sp", sp)}) # nolint: brace_linter.

        alpha <- matrix(0, nrow = numSpecies, ncol = numSpecies)
        beta <- matrix(0, nrow = numSpecies, ncol = numCovs)

        if (parmSet != "indep" && numSpecies > 5) { # if the species are independent, there should be 0 interactions
          # choose species interactions for block 1
          for(ii in sample(x = 1:(numSpecies %/% 2), size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
            for(jj in sample(x = 1:(numSpecies %/% 2), size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
              alpha[ii, jj] <- sample(x = c(-1, 1), size = 1)
            }
          }

          # choose species interactions for block 2
          for(ii in sample(x = (numSpecies %/% 2 + 1):numSpecies, size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
            for(jj in sample(x = (numSpecies %/% 2 + 1):numSpecies, size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
              alpha[ii, jj] <- sample(x = c(-1, 1), size = 1)
            }
          }
        }

        if (parmSet != "indep" && numSpecies == 3) {
          alpha <- matrix(c(0, -1, 0, -1, 0, 0, 0, 0, 0), byrow = TRUE)
        }

        if (parmSetCov == "indep") { # if the species are independent, there should be 0 interactions
          # set environment interactions as independent
          for(iii in 1:numSpecies) {
            beta[iii, iii] <- 1 # one independent covariate per species
            alpha[iii, iii] <- 0
            beta[iii, numSpecies + 1] <- 1 # this is the intercept variable
          }
        } else if (parmSetCov == "random") {
          #nrows is num species, ncol is number of covs
          # choose species interactions for block 1
          for(ii in sample(x = 1:(numSpecies %/% 2), size = as.integer(sqrt(numSpecies %/% 4)), replace = TRUE)) {
            for(jj in sample(x = 1:(numCovs %/% 2), size = as.integer(sqrt(numCovs %/% 4)), replace = TRUE)) {
              beta[ii, jj] <- runif(min = -1, max = 1, n = 1)
            }
          }

          # choose species interactions for block 2
          for(ii in sample(x = (numSpecies %/% 2 + 1):numSpecies, size = as.integer(sqrt(numSpecies %/% 4)), replace = TRUE)) {
            for(jj in sample(x = (numCovs %/% 2 + 1):numCovs, size = as.integer(sqrt(numSpecies %/% 4)), replace = TRUE)) {
              beta[ii, jj] <- runif(min = -1, max = 1, n = 1)
            }
          }

          for(iiii in 1:numCovs) {
            beta[iiii, iiii] <- 1 # one independent covariate per species
            beta[iiii, numCovs] <- 1 # this is the intercept variable
          }

          for(iii in 1:numSpecies) {
            alpha[iii, iii] <- 0
          }
        }

        #alpha and beta numbers respectively (for scaling purposes)
        #cov effect
        c3 <- 200

        alpha <- alpha * c2
        beta <- beta * c3
        
        #mig_rates=rep(0.1, times=numSpecies) # OR vvv sample this from gamma per species
        mig_rates <- rgamma(n = numSpecies, shape = mean_mig_rate / 0.01, scale = 0.01) # mean = 0.1, mode 0.09, variance 0.001
        readThreshold <- sample(x = 1:50, size = 1)

    } else { 
        num_samples_time <- 500 # sample times per location
        num_samples_space <- 50 # sample locations per time
        radius <- 16 # for neighborhoods: a little more the 
        # hypotenuse of 11.1 unit grid, to get neighbors be right for this setup
        #mean_fpr <- 0.05
        #fpr_mode <- "none" # "dependent_sp" "constant" "none"
        covNoise_sd <- 0
        #intrinsic growth rate
        r <- 0.05
        #noise param for popn growth 
        sigma <- 0.005

        #N_50 <- 50 # number of animals for which detection prob is 50%

        mean_mig_rate <- 0.01

        c2 <- 300 #changed to be varied -- sp effect
        covMeasureNoise_sd <- 0

        indivSampleProb <- 0.01
        readSampleRate <- 1
        corr_mode <- "independent"

        covVars <- list()
        for (i in 1:numSpecies) {
          covVars[[i]] <- list(type = "randomWalk")
        }
        covVars[[numSpecies + 1]] <- list(type = "constant")

        covar_scale_space <- 1

        numCovs <- length(covVars)
        names_cov <- sapply(1:numCovs, FUN = function(cov) {
                        if (covVars[[cov]]$type != "constant") {
                          return(paste0("Cov", cov))}}) 
        names_cov <- unlist(names_cov)
        names_all_cov <- unlist(sapply(1:numCovs, FUN = function(cov) {paste0("Cov", cov)}))
        names_species <- sapply(1:numSpecies, FUN = function(sp) {paste0("Sp", sp)}) # nolint: brace_linter.

        alpha <- matrix(0, nrow = numSpecies, ncol = numSpecies)
        beta <- matrix(0, nrow = numSpecies, ncol = numCovs)

        if (parmSet != "indep" && numSpecies > 5) { # if the species are independent, there should be 0 interactions
          # choose species interactions for block 1
          for(ii in sample(x = 1:(numSpecies %/% 2), size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
            for(jj in sample(x = 1:(numSpecies %/% 2), size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
              alpha[ii, jj] <- sample(x = c(-1, 1), size = 1)
            }
          }

          # choose species interactions for block 2
          for(ii in sample(x = (numSpecies %/% 2 + 1):numSpecies, size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
            for(jj in sample(x = (numSpecies %/% 2 + 1):numSpecies, size = as.integer(sqrt(numSpecies %/% 2)), replace = TRUE)) {
              alpha[ii, jj] <- sample(x = c(-1, 1), size = 1)
            }
          }
        }

        if (parmSet != "indep" && numSpecies == 3) {
          alpha <- matrix(c(0, -1, 0, -1, 0, 0, 0, 0, 0), byrow = TRUE)
        }

        if (parmSetCov == "indep") { # if the species are independent, there should be 0 interactions
          # set environment interactions as independent
          for(iii in 1:numSpecies) {
            beta[iii, iii] <- 1 # one independent covariate per species
            alpha[iii, iii] <- 0
            beta[iii, numSpecies + 1] <- 1 # this is the intercept variable
          }
        } else if (parmSetCov == "random") {
          #nrows is num species, ncol is number of covs
          # choose species interactions for block 1
          for(ii in sample(x = 1:(numSpecies %/% 2), size = as.integer(sqrt(numSpecies %/% 4)), replace = TRUE)) {
            for(jj in sample(x = 1:(numCovs %/% 2), size = as.integer(sqrt(numCovs %/% 4)), replace = TRUE)) {
              beta[ii, jj] <- runif(min = -1, max = 1, n = 1)
            }
          }

          # choose species interactions for block 2
          for(ii in sample(x = (numSpecies %/% 2 + 1):numSpecies, size = as.integer(sqrt(numSpecies %/% 4)), replace = TRUE)) {
            for(jj in sample(x = (numCovs %/% 2 + 1):numCovs, size = as.integer(sqrt(numSpecies %/% 4)), replace = TRUE)) {
              beta[ii, jj] <- runif(min = -1, max = 1, n = 1)
            }
          }

          for(iii in 1:numSpecies) {
            beta[iii, iii] <- 1 # one independent covariate per species
            alpha[iii, iii] <- 0
            beta[iii, numSpecies + 1] <- 1 # this is the intercept variable
          }
        }

        #alpha and beta numbers respectively (for scaling purposes)
        #cov effect
        c3 <- 200

        alpha <- alpha * c2
        beta <- beta * c3
        
        #mig_rates=rep(0.1, times=numSpecies) # OR vvv sample this from gamma per species
        mig_rates <- rgamma(n = numSpecies, shape = mean_mig_rate / 0.01, scale = 0.01) # mean = 0.1, mode 0.09, variance 0.001
        readThreshold <- 5
        
    }

    num_gens <- 10000
    N0_0 <- 10
    x_split <- 10
    y_split <- 10
    #TODO change to accommodate actual spatial locations -- this is mostly done I think -- 
    # might need something to figure out the lat long to dist mx
    xdim1 <- 100
    ydim1 <- 100
    location_mode <- "grid"
    
    #constant_fpr <- 0.05
    #beta_fpr <- 50 # mean = alpha/(alpha+beta) = 5/(100+5) 
    #alpha_fpr <- (-beta_fpr * mean_fpr) / (mean_fpr - 1) # moment matching
    #fpr_mode = "dependent_sp" # "independent" "dependent_sp" "constant" "none"

    #detection prob param

    #det_prob_exp <- 2 # fix this at 2
    #det_prob_add <- N_50^det_prob_exp # vary this 

    abundanceEffectType <- "arctangent"

    # parameters from old versions that dont exist in abundance version anymore
    det_prob_exp <- NA
    det_prob_add <- NA
    alpha_fpr <- NA
    beta_fpr <- NA
    constant_fpr <- NA
    N_50 <- NA
    mean_fpr <- NA
    fpr_mode <- NA

    params <- list(xdim = xdim1, ydim = ydim1, location_mode = location_mode, x_split = x_split, y_split = y_split, 
               num_gens = num_gens, mig_rates = mig_rates, mean_mig_rate = mean_mig_rate, 
               numCovs = numCovs, numSpecies = numSpecies, r = r, sigma = sigma, 
               det_prob_exp = det_prob_exp, det_prob_add = det_prob_add, N_50 = N_50, c2 = c2, c3 = c3, covVars = covVars, 
               abundanceEffectType = abundanceEffectType, beta = beta, alpha = alpha, fpr = list(mean_fpr = mean_fpr, mode = fpr_mode, 
               alpha_fpr = alpha_fpr, beta_fpr = beta_fpr, constant_fpr = constant_fpr), names_cov = names_cov, names_all_cov = names_all_cov,
               names_species = names_species, N0_0 = N0_0, covNoise_sd = covNoise_sd, covMeasureNoise_sd = covMeasureNoise_sd, 
               num_samples_time = num_samples_time, num_samples_space = num_samples_space, radius = radius, indivSampleProb = indivSampleProb, 
               readSampleRate = readSampleRate, corr_mode = corr_mode, readThreshold = readThreshold, covar_scale_space = covar_scale_space)
    return(params)
}


















getParms_5 <- function(random = FALSE, parmSet = 1) { # only for 5 species
    if (random == TRUE) {
        num_samples_time <- sample(x = 2:50, size = 1) # sample times per location
        num_samples_space <- sample(2:50, size = 1) # sample locations per time
        #radius <- 16
        radius <- runif(n = 1, min = 12, max = 50) # for neighborhoods: a little more the 
        # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
        mean_fpr <- runif(n = 1, min = 0, max = 0.2)
        #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

        fpr_mode <- sample(x = c("none", "independent", "dependent_sp", "constant"), size = 1)
        #fpr_mode <- "independent"
        covNoise_sd <- runif(n = 1, min = 0, max = 0.1) # process noise
        covMeasureNoise_sd <- runif(n = 1, min = 0, max = 0.2) # measurement noise

        # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
        # intrinsic growth rate
        r <- runif(n = 1, min = 0.01, max = 1) #todo think on reasonableness here
        # noise param for popn growth 
        sigma <- runif(n = 1, min = 0, max = 0.2)

        N_50 <- sample(x = 1:1000, size = 1) # number of animals for which detection prob is 50%

        #species effect
        c2 <- runif(n = 1, min = 100, max = 200)

        mean_mig_rate <- runif(n = 1, min = 0, max = 0.1) # poisson rate per individual per time

        N0_0 <- 10
        x_split <- 8
        y_split <- 8
    } else { 
        num_samples_time <- 50 # sample times per location
        num_samples_space <- 50 # sample locations per time
        radius <- 16 # for neighborhoods: a little more the 
        # hypotenuse of 11.1 unit grid, to get neighbors be right for this setup
        mean_fpr <- 0.05
        fpr_mode <- "none" # "dependent_sp" "constant" "none"
        covNoise_sd <- 0
        #intrinsic growth rate
        r <- 0.05
        #noise param for popn growth 
        sigma <- 0.005

        N_50 <- 50 # number of animals for which detection prob is 50%

        mean_mig_rate <- 0.01

        N0_0 <- 10
        x_split <- 8
        y_split <- 8

        c2 <- 300 #changed to be varied -- sp effect
        covMeasureNoise_sd <- 0
    }

    #TODO change to accommodate actual spatial locations -- this is mostly done I think -- 
    # might need something to figure out the lat long to dist mx
    xdim1 <- 100
    ydim1 <- 100
    location_mode <- "grid"
    
    constant_fpr <- 0.05
    beta_fpr <- 50 # mean = alpha/(alpha+beta) = 5/(100+5) 
    alpha_fpr <- (-beta_fpr * mean_fpr) / (mean_fpr - 1) # moment matching
    #fpr_mode = "dependent_sp" # "independent" "dependent_sp" "constant" "none"

    num_gens <- 1000
    #detection prob param

    det_prob_exp <- 2 # fix this at 2
    det_prob_add <- N_50^det_prob_exp # vary this 



    abundanceEffectType <- "arctangent"

    covVars <- list()
    covVars[[1]] <- list(type = "randomWalk")
    covVars[[2]] <- list(type = "randomWalk")
    covVars[[3]] <- list(type = "randomWalk")
    covVars[[4]] <- list(type = "randomWalk")
    covVars[[5]] <- list(type = "randomWalk")
    covVars[[6]] <- list(type = "constant")


    numCovs <- length(covVars)
    names_cov <- sapply(1:numCovs, FUN = function(cov) {
                    if (covVars[[cov]]$type != "constant") {
                      return(paste0("Cov", cov))}}) 
    names_cov <- unlist(names_cov)
    names_all_cov <- unlist(sapply(1:numCovs, FUN = function(cov) {paste0("Cov", cov)}))
    numSpecies <- 5
    names_species <- sapply(1:numSpecies, FUN = function(sp) {paste0("Sp", sp)}) # nolint: brace_linter.

    #TODO
    alpha <- matrix(0, nrow = numSpecies, ncol = numSpecies)
    beta <- matrix(0, nrow = numSpecies, ncol = numCovs)

    alpha[1, 2] <- -1
    alpha[2, 1] <- -1
    alpha[5, 3] <- -1
    alpha[5, 4] <- -1

    beta[1, 1] <- 1
    beta[2, 2]  <- 1
    beta[3, 3]  <- 1
    beta[4, 4]  <- 1
    beta[5, 5]  <- 1
    beta[1, 6] <- 1
    beta[2, 6] <- 1
    beta[3, 6] <- 1
    beta[4, 6] <- 1
    beta[5, 6] <- 1

    #alpha and beta numbers respectively (for scaling purposes)
    #cov effect
    c3 <- 200

    alpha <- alpha * c2
    beta <- beta * c3
    
    #mig_rates=rep(0.1, times=numSpecies) # OR vvv sample this from gamma per species
    mig_rates <- rgamma(n = numSpecies, shape = mean_mig_rate / 0.01, scale = 0.01) # mean = 0.1, mode 0.09, variance 0.001

    params <- list(xdim = xdim1, ydim = ydim1, location_mode = location_mode, x_split = x_split, y_split = y_split, 
               num_gens = num_gens, mig_rates = mig_rates, mean_mig_rate = mean_mig_rate, 
               numCovs = numCovs, numSpecies = numSpecies, r = r, sigma = sigma, 
               det_prob_exp = det_prob_exp, det_prob_add = det_prob_add, N_50 = N_50, c2 = c2, c3 = c3, covVars = covVars, 
               abundanceEffectType = abundanceEffectType, beta = beta, alpha = alpha, fpr = list(mean_fpr = mean_fpr, mode = fpr_mode, 
               alpha_fpr = alpha_fpr, beta_fpr = beta_fpr, constant_fpr = constant_fpr), names_cov = names_cov, names_all_cov = names_all_cov,
               names_species = names_species, N0_0 = N0_0, covNoise_sd = covNoise_sd, covMeasureNoise_sd = covMeasureNoise_sd, 
               num_samples_time = num_samples_time, num_samples_space = num_samples_space, radius = radius)
    return(params)
}











getParms <- function(random = TRUE, parmSet = 1) { # only for 3 species
    if (random == TRUE) {
        num_samples_time <- sample(x = 5:500, size = 1) # sample times per location
        num_samples_space <- sample(5:64, size = 1) # sample locations per time
        #radius <- 16
        radius <- runif(n = 1, min = 12, max = 50) # for neighborhoods: a little more the 
        # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
        mean_fpr <- runif(n = 1, min = 0, max = 0.2)
        #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

        fpr_mode <- sample(x = c("none", "independent", "dependent_sp", "constant"), size = 1)
        #fpr_mode <- "independent"
        covNoise_sd <- runif(n = 1, min = 0, max = 0.1) # process noise
        covMeasureNoise_sd <- runif(n = 1, min = 0, max = 0.2) # measurement noise

        # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
        # intrinsic growth rate
        r <- runif(n = 1, min = 0.01, max = 1) #todo think on reasonableness here
        # noise param for popn growth 
        sigma <- runif(n = 1, min = 0, max = 0.5)

        N_50 <- sample(x = 1:1000, size = 1) # number of animals for which detection prob is 50%

        #species effect
        c2 <- runif(n = 1, min = 100, max = 200)

        #a12 <- sample(x = c(-1, 0, 1), size = 1)
        #a13 <- sample(x = c(-1, 0, 1), size = 1)
        #a21 <- sample(x = c(-1, 0, 1), size = 1)
        #a23 <- sample(x = c(-1, 0, 1), size = 1)
        #a31 <- sample(x = c(-1, 0, 1), size = 1)
        #a32 <- sample(x = c(-1, 0, 1), size = 1)
        #b11 <- sample(x = c(-1, 0, 1), size = 1)
        #b12 <- sample(x = c(-1, 0, 1), size = 1)
        #b13 <- sample(x = c(-1, 0, 1), size = 1)
        #b14 <- 1
        #b21 <- sample(x = c(-1, 0, 1), size = 1)
        #b22 <- sample(x = c(-1, 0, 1), size = 1)
        #b23 <- sample(x = c(-1, 0, 1), size = 1)
        #b24 <- 1
        #b31 <- sample(x = c(-1, 0, 1), size = 1)
        #b32 <- sample(x = c(-1, 0, 1), size = 1)
        #b33 <- sample(x = c(-1, 0, 1), size = 1)
        #b34 <- 1

        mean_mig_rate <- runif(n = 1, min = 0, max = 0.1) # poisson rate per individual per time

        N0_0 <- 10
        x_split <- 8
        y_split <- 8
    } else if (parmSet == 1) {
      num_samples_time <- 15 # sample times per location
      num_samples_space <- 20 # sample locations per time
      radius <- 16
      # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
      mean_fpr <- 0
      #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

      fpr_mode <- "none"
      #fpr_mode <- "independent"
      covNoise_sd <- 0.01 # process noise
      covMeasureNoise_sd <- 0

      # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      # intrinsic growth rate
      r <- 1 #todo think on reasonableness here
      # noise param for popn growth 
      sigma <- 0.001

      N_50 <- 10 # number of animals for which detection prob is 50%

      #species effect
      c2 <- 150

      mean_mig_rate <- 0.001 # poisson rate per individual per time

      N0_0 <- 10
      x_split <- 8
      y_split <- 8
    } else if(parmSet == 2) {
      # more samples, higher sigma, lower species effect (compared to 1) 
      # best so far
      num_samples_time <- 15 # sample times per location
      num_samples_space <- 30 # sample locations per time
      radius <- 16
      # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
      mean_fpr <- 0
      #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

      fpr_mode <- "none"
      #fpr_mode <- "independent"
      covNoise_sd <- 0.05 # process noise
      covMeasureNoise_sd <- 0

      # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      # intrinsic growth rate
      r <- 1 #todo think on reasonableness here
      # noise param for popn growth 
      sigma <- 0.001

      N_50 <- 10 # number of animals for which detection prob is 50%

      #species effect
      c2 <- 100

      mean_mig_rate <- 0.01 # poisson rate per individual per time

      N0_0 <- 10
      x_split <- 8
      y_split <- 8
      } else if (parmSet == 3) {
      # higer species effect (compared to 2) --does worse than 2 on logistic
      num_samples_time <- 15 # sample times per location
      num_samples_space <- 30 # sample locations per time
      radius <- 16
      # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
      mean_fpr <- 0
      #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

      fpr_mode <- "none"
      #fpr_mode <- "independent"
      covNoise_sd <- 0.05 # process noise
      covMeasureNoise_sd <- 0

      # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      # intrinsic growth rate
      r <- 1 #todo think on reasonableness here
      # noise param for popn growth 
      sigma <- 0.001

      N_50 <- 10 # number of animals for which detection prob is 50%

      #species effect
      c2 <- 200

      mean_mig_rate <- 0.01 # poisson rate per individual per time

      N0_0 <- 10
      x_split <- 8
      y_split <- 8
    } else if (parmSet == "optim") {
      # higer species effect (compared to 2) --does worse than 2 on logistic
      num_samples_time <- 294 # sample times per location
      num_samples_space <- 32 # sample locations per time
      radius <- 16
      # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
      mean_fpr <- 0
      #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

      fpr_mode <- "none"
      #fpr_mode <- "independent"
      covNoise_sd <- 0.05 # process noise
      covMeasureNoise_sd <- 0.14

      # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      # intrinsic growth rate
      r <- 0.58 #todo think on reasonableness here
      # noise param for popn growth 
      sigma <- 0.33

      N_50 <- 198 # number of animals for which detection prob is 50%

      #species effect
      c2 <- 50

      mean_mig_rate <- 0.05 # poisson rate per individual per time

      N0_0 <- 10
      x_split <- 8
      y_split <- 8
    } else if(parmSet == 4) {
      # higer species effect (compared to 2) --does worse than 2 on logistic
      num_samples_time <- 25 # sample times per location
      num_samples_space <- 10 # sample locations per time
      radius <- 16
      # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
      mean_fpr <- 0
      #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

      fpr_mode <- "none"
      #fpr_mode <- "independent"
      covNoise_sd <- 0.05 # process noise
      covMeasureNoise_sd <- 0

      # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      # intrinsic growth rate
      r <- .8 #todo think on reasonableness here
      # noise param for popn growth 
      sigma <- 0

      N_50 <- 50 # number of animals for which detection prob is 50%

      #species effect
      c2 <- 50

      mean_mig_rate <- 0.05 # poisson rate per individual per time
      N0_0 <- 10
      x_split <- 8
      y_split <- 8
    } else if(parmSet == 5) {
      # low migration rate, high radius (compared to 2) 
      # 
      num_samples_time <- 15 # sample times per location
      num_samples_space <- 30 # sample locations per time
      radius <- 100
      # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
      mean_fpr <- 0
      #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

      fpr_mode <- "none"
      #fpr_mode <- "independent"
      covNoise_sd <- 0.05 # process noise
      covMeasureNoise_sd <- 0

      # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      # intrinsic growth rate
      r <- 1 #todo think on reasonableness here
      # noise param for popn growth 
      sigma <- 0.001

      N_50 <- 10 # number of animals for which detection prob is 50%

      #species effect
      c2 <- 100

      mean_mig_rate <- 0.001 # poisson rate per individual per time

      N0_0 <- 10
      x_split <- 8
      y_split <- 8
    } else if(parmSet == "rw") {
      # parmSet5 plus random walk covs and no cov noise and r lower, sigma lower, and migration rate lower, c2 higher
      num_samples_time <- 15 # sample times per location
      num_samples_space <- 30 # sample locations per time
      radius <- 100
      # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
      mean_fpr <- 0
      #fpr_mode <- sample(x = c("independent", "dependent_sp"), size = 1)

      fpr_mode <- "none"
      #fpr_mode <- "independent"
      covNoise_sd <- 0 # process noise
      covMeasureNoise_sd <- 0

      # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
      # intrinsic growth rate
      r <- .2 #todo think on reasonableness here
      # noise param for popn growth 
      sigma <- 0.001

      N_50 <- 10 # number of animals for which detection prob is 50%

      #species effect
      c2 <- 200

      mean_mig_rate <- 0.0001 # poisson rate per individual per time

      N0_0 <- 10
      x_split <- 100
      y_split <- 10
    } else { 
        num_samples_time <- 10 # sample times per location
        num_samples_space <- 10 # sample locations per time
        radius <- 16 # for neighborhoods: a little more the 
        # hypotenuse of 11.1 unit grid, to get neighbors be right for this setup
        mean_fpr <- 0.05
        fpr_mode <- "none" # "dependent_sp" "constant" "none"
        covNoise_sd <- 0
        #intrinsic growth rate
        r <- 0.05
        #noise param for popn growth 
        sigma <- 0.005

        N_50 <- 50 # number of animals for which detection prob is 50%

        mean_mig_rate <- 0.01

        N0_0 <- 10
        x_split <- 8
        y_split <- 8
    }

    a12 <- -1
    a13 <- 0
    a21 <- -1
    a23 <- 0
    a31 <- 0
    a32 <- 0

    b11 <- 1
    b12 <- 0 
    b13 <- 0
    b14 <- 1
    b21 <- 0
    b22 <- 1
    b23 <- 0
    b24 <- 1
    b31 <- 0
    b32 <- 0
    b33 <- 1
    b34 <- 1

    #TODO change to accommodate actual spatial locations -- this is mostly done I think -- 
    # might need something to figure out the lat long to dist mx
    xdim1 <- 100
    ydim1 <- 100
    location_mode <- "grid"
    
    constant_fpr <- 0.05
    beta_fpr <- 50 # mean = alpha/(alpha+beta) = 5/(100+5) 
    alpha_fpr <- (-beta_fpr * mean_fpr) / (mean_fpr - 1) # moment matching
    #fpr_mode = "dependent_sp" # "independent" "dependent_sp" "constant" "none"

    num_gens <- 1000
    #detection prob param

    det_prob_exp <- 2 # fix this at 2
    det_prob_add <- N_50^det_prob_exp # vary this 

    #alpha and beta numbers respectively (for scaling purposes)
    # species effect
    #c2 <- 50 #changed to be varied
    #cov effect
    c3 <- 200

    abundanceEffectType <- "arctangent"

    covVars <- list()
    covVars[[1]] <- list(type = "peaksTimesLinear", coef1func = function(t) {return(sin(t * pi / 87))}, coef2func = function(t) {return(-0.5 * t)})
    covVars[[2]] <- list(type = "singleTrough")
    covVars[[3]] <- list(type = "ripples")
    covVars[[4]] <- list(type = "constant")

    if (parmSet == "rw") {
      covVars[[1]] <- list(type = "randomWalk")
      covVars[[2]] <- list(type = "randomWalk")
      covVars[[3]] <- list(type = "randomWalk")
      covVars[[4]] <- list(type = "constant")
    }

    numCovs <- length(covVars)
    names_cov <- sapply(1:numCovs, FUN = function(cov) {
                    if (covVars[[cov]]$type != "constant") {
                      return(paste0("Cov", cov))}}) 
    names_cov <- unlist(names_cov)
    names_all_cov <- unlist(sapply(1:numCovs, FUN = function(cov) {paste0("Cov", cov)}))
    numSpecies <- 3
    names_species <- sapply(1:numSpecies, FUN = function(sp) {paste0("Sp", sp)}) # nolint: brace_linter.
    alpha <- matrix(c2 * c(0, a12, a13, a21, 0, a23, a31, a32, 0), nrow = numSpecies, ncol = numSpecies, byrow = TRUE)
    beta <- matrix(data = c3 * c(b11, b12, b13, b14, b21, b22, b23, b24, b31, b32, b33, b34), nrow = numSpecies, ncol = numCovs, byrow = TRUE)
    
    #mig_rates=rep(0.1, times=numSpecies) # OR vvv sample this from gamma per species
    mig_rates <- rgamma(n = numSpecies, shape = mean_mig_rate / 0.01, scale = 0.01) # mean = 0.1, mode 0.09, variance 0.001

    params <- list(xdim = xdim1, ydim = ydim1, location_mode = location_mode, x_split = x_split, y_split = y_split, 
               num_gens = num_gens, mig_rates = mig_rates, mean_mig_rate = mean_mig_rate, 
               numCovs = numCovs, numSpecies = numSpecies, r = r, sigma = sigma, 
               det_prob_exp = det_prob_exp, det_prob_add = det_prob_add, N_50 = N_50, c2 = c2, c3 = c3, covVars = covVars, 
               abundanceEffectType = abundanceEffectType, beta = beta, alpha = alpha, fpr = list(mean_fpr = mean_fpr, mode = fpr_mode, 
               alpha_fpr = alpha_fpr, beta_fpr = beta_fpr, constant_fpr = constant_fpr), names_cov = names_cov, names_all_cov = names_all_cov,
               names_species = names_species, N0_0 = N0_0, covNoise_sd = covNoise_sd, covMeasureNoise_sd = covMeasureNoise_sd, 
               num_samples_time = num_samples_time, num_samples_space = num_samples_space, radius = radius)
    return(params)
}