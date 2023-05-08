library(ggplot2)
library(reshape)
library(gridExtra)

getParms <- function(random = TRUE) {
    if (random == TRUE) {
        num_samples_time <- sample(x = 5:500, size = 1) # sample times per location
        num_samples_space <- sample(5:64, size = 1) # sample locations per time
        radius <- 16
        # radius <- runif(n = 1, min = 12, max = 50) # for neighborhoods: a little more the 
        # hypotenuse of 11.1 unit grid (16), to get neighbors be "right" for this setup
        mean_fpr <- runif(n = 1, min = 0, max = 0.2)
        fpr_mode <- sample(x = c("none", "dependent_sp", "constant"), size = 1)
        covNoise_sd <- runif(n = 1, min = 0, max = 0.1) # process noise
        covMeasureNoise_sd <- runif(n = 1, min = 0, max = 0.2) # measurement noise

        # this is the code from the sim: dN <- params$r * lastN * (1 - lastN / thisK[[s]]) + params$sigma * lastN * dWt
        # intrinsic growth rate
        r <- runif(n = 1, min = 0.01, max = 1) #todo think on reasonableness here
        # noise param for popn growth 
        sigma <- runif(n = 1, min = 0, max = 0.5)

        N_50 <- sample(x = 1:1000, size = 1) # number of animals for which detection prob is 50%

        #species effect
        c2 <- runif(n = 1, min = 25, max = 150)

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
        sigma <- 0.05

        N_50 <- 50 # number of animals for which detection prob is 50%

        mean_mig_rate <- 0.1

        N0_0 <- 10
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
    x_split <- 8
    y_split <- 8
    
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

# Functions
get_locations <- function(xdim = xdim1, ydim = ydim1, mode = "grid", xSplit = x_split, ySplit = y_split) {
  if (mode == "grid") {
    # need to figure out the data type that's gonna be good for locations?
    # why is r?
    thisX <- 0
    thisY <- 0
    locList <- list()
    i <- 1
    for (x in 1:xSplit){
      thisX <- thisX + xdim / (xSplit + 1)
      for (y in 1:ySplit){
        thisY <- thisY + ydim / (ySplit + 1)
        locList[[i]] <- c(thisX, thisY)
        i <- i + 1
      }
      thisY <- 0
    }
  }else{
    print("mode not implemented yet in get_locations") 
  }
  return(locList)
}

#TODO make it so that you can just read covs in from other place
get_cov <- function(location, type = "peaks", coef1func = NA, coef2func = NA, t = NA, params = params) {
  x <- location[1]
  y <- location[2]
  if (type == "peaks") {
    c1 <- params$xdim1 / 4
    c2 <- params$ydim1 / 4
    temp1 <- (x - 2 * c1) / c1
    temp2 <- (y - 2 * c2) / c2
    return(temp1 * temp2 * exp(-(temp1^2 + temp2^2)))
  } else if (type == "linear") {
    coef1 <- coef1func(t)
    coef2 <- coef2func(t)
    return(coef1 * x + coef2 * y)
  } else if (type == "peaksTimesLinear") {
    coef1 <- coef1func(t)
    coef2 <- coef2func(t)
    c1 <- params$xdim / 4
    c2 <- params$ydim / 4
    temp1 <- (x - 2 * c1) / c1
    temp2 <- (y - 2 * c2) / c2
    return((temp1 * temp2 * exp(-(temp1^2 + temp2^2))) * (coef1 * x + coef2 * y))
  } else if (type == "constant") {
    return(1)
  } else if (type == "singleTrough") {
    temp1 <- x - (t / params$num_gens) * (params$xdim / 2)
    temp2 <- y - (t / params$num_gens) * (params$ydim / 2)
    return(sqrt(temp1^2 + temp2^2))
  }else if (type == "ripples") {
    # x minus some number between 0 and 1 that changes with time (cycles every ~400)
    temp1 <- x / params$xdim - abs(sin((1 / 123) * t))  
    temp2 <- y / params$ydim - abs(cos((1 / 123) * t))
    return(sin(temp1 * 3 * temp2 * 3)) # the times 3 is just to get the function to behave
  } else {
    print("type invalid in function get_cov")
  }
}

abundanceEffectTransform <- function(mx, type = "none") {
  if (type == "sqrt") {
    return(sign(mx) * sqrt(abs(mx)))
  }else if (type == "arctangent") {
    return(atan(mx))
  }else {
    return(mx)
  }
}




makePlots <- function(sim_data = sim_data, params = params, locList = locList, locPlots = c(19, 8, 6, 37, 64, 46), outDir) {
    detProbP <- ggplot(data.frame(x = 0:500, y = unlist(lapply(X = as.list(0:500), 
                                  FUN = function(x) { 
                                    (x^params$det_prob_exp) / (params$det_prob_add + x^params$det_prob_exp)}))),  # nolint: object_usage_linter.
                                    aes(x, y)) + # nolint: object_usage_linter.
                                    geom_point() +
                                    labs(x = "Abundance -- N(s,t)", y = "Detection Probability")
    
    ggsave(detProbP, file = paste0(outDir, "detProbPlot.png"), width = 7, height = 7)
    ######################################################################
    #wrangle data
    ######################################################################
    # this creates a list of length of number of locations and in each index there is a DF with colums time, species1, species2, species3
    species_abundance <- list() 
    for (s_index in seq_along(locList)) {
      # make data frame for that location
      species_abundance[[s_index]] <- data.frame(time = 1:params$num_gens)
      for (k in 1:params$numSpecies) {
        spAbund <- lapply(sim_data, FUN = function(x) {return(x$N[[s_index]][k])}) # nolint: brace_linter.
        species_abundance[[s_index]][[paste0("species", k)]] <- unlist(spAbund)
      }
    }
    
    #wrangle data
    carrying_capacities <- list()
    for (s_index in seq_along(locList)){
      # make data frame for that location
      carrying_capacities[[s_index]] <- data.frame(time = 1:params$num_gens)
      for (k in 1:params$numSpecies){
        Kcapacity <- lapply(sim_data, FUN = function(x) {return(x$K[[s_index]][k])})
        carrying_capacities[[s_index]][[paste0("species", k)]] <- unlist(Kcapacity)
      }
    }
    
    #plot species against each other
    s_indexL <- locPlots
    abd_plotList <- list()
    cc_plotList <- list()
    for (i in seq_along(s_indexL)){
      s_index <- s_indexL[i]
      # species_abundance[[s_index]]
      df <- melt(species_abundance[[s_index]][species_abundance[[s_index]]$time %in% 100:1000, ], id.vars = "time")
      names(df) <- c("time", "species", "abundance")
      abd_plot <- ggplot(df, aes(x = time, y = abundance, color = species)) + # nolint: object_usage_linter.
        geom_line() +
        ggtitle(paste("Species abundance for location", s_index))
      
      carrying_capacities[[s_index]]
      df <- melt(carrying_capacities[[s_index]][carrying_capacities[[s_index]]$time %in% 100:1000, ], id.vars = "time")
      names(df) <- c("time", "species", "carrying_capacity")
      cc_plot <- ggplot(df, aes(x = time, y = carrying_capacity, color = species)) + # nolint: object_usage_linter.
        geom_line() +
        ggtitle(paste("Carrying capacities for location", s_index))
      abd_plotList[[i]] <- abd_plot
      cc_plotList[[i]] <- cc_plot
    }
    
    abd_plots <- arrangeGrob(abd_plotList[[1]], abd_plotList[[2]], abd_plotList[[3]], abd_plotList[[4]],
                            abd_plotList[[5]], abd_plotList[[6]], nrow = 2)
    ccplots <- arrangeGrob(cc_plotList[[1]], cc_plotList[[2]], cc_plotList[[3]], cc_plotList[[4]], 
                            cc_plotList[[5]], cc_plotList[[6]], nrow = 2)
    ggsave(abd_plots, filename = paste0(outDir, "abdPlots.png"), width = 16, height = 8)
    ggsave(ccplots, filename = paste0(outDir, "ccPlots.png"), width = 16, height = 8)

    #Plot average popn over space
    #data wrangling
    x_coords <- lapply(locList, FUN = function(x) {return(x[1])})
    y_coords <- lapply(locList, FUN = function(x) {return(x[2])}) 
    
    sp_average <- lapply(species_abundance, FUN = function(sp_abd) {
      return(apply(X = sp_abd, MARGIN = 2, FUN = mean))
    })
    
    abd_av_df <- data.frame(s_index = seq_along(locList), x = unlist(x_coords), y = unlist(y_coords))
    for (k in 1:params$numSpecies){
      #just cause the first column isnt real
      sp_i <- k + 1
      sp_av <- lapply(sp_average, FUN = function(spav) {return(spav[sp_i])})
      abd_av_df[[paste0("species", k, "AverageAbundance")]] <- unlist(sp_av)
    }
    #abd_av_df
    
    aesNamesString <- lapply(1:params$numSpecies, FUN = function(sp) {paste0("species", sp, "AverageAbundance")})
    #plot average abundance of each species across space
    avgAbdPlotsL <- lapply(aesNamesString, FUN = function(varName) {
      p <- ggplot(abd_av_df, aes_string(x = "x", y = "y", label = "s_index", color = varName)) + # nolint: object_usage_linter.
        geom_text() +
        scale_color_gradient(low = "blue", high = "green")
      #return(p)
    })
    
    avgAbdPlot <- arrangeGrob(grobs = avgAbdPlotsL)
    ggsave(avgAbdPlot, file = paste0(outDir, "avgAbdPlot.png"), height = 7, width = 7)
    
    
    #plot env over space
    
    #data wrangling
    #x_coords <- lapply(locList, FUN = function(x) {return(x[1])})
    #y_coords <- lapply(locList, FUN = function(x) {return(x[2])})
    #env_df <- data.frame(s_index = seq_along(locList), x = unlist(x_coords), 
    #                  y = unlist(y_coords)) 
    #for (cov in 1:params$numCovs){
      # TODO this covariate thing is jsut for one time point and that's misleading if they change thru time
    #  env_df[[paste0("cov", cov, "_norm")]] <- covL_norm[[cov]]
    #}
    #for (k in 1:params$numSpecies){
    #  envSp <- sapply(seq_along(locList), FUN = function(s_index) {return(params$beta[k, ] %*% covList[[s_index]])})
    #  env_df[[paste0("species", k, "Environment")]] <- envSp
    #}
    
    #plot environments
    #aesNamesString <- lapply(1:params$numSpecies, FUN = function(sp) {paste0("species", sp, "Environment")})
    #plot environment for each species over space
    #envPlotsL <- lapply(aesNamesString, FUN = function(varName) {
    #  p <- ggplot(env_df, aes_string(x = "x", y = "y", label = "s_index", color = varName)) +
    #    geom_text() +
    #    scale_color_gradient(low = "blue", high = "red")
    #  return(p)
    #})
    #envPlots <- grid.arrange(grobs = envPlotsL)
    
    #Plot covariates for single time pt
    #aesNamesString <- lapply(1:params$numCovs, FUN = function(cov) {paste0("cov", cov, "_norm")})
    #covPlotsL <- lapply(aesNamesString, FUN = function(varName) {
    #  p <- ggplot(env_df, aes_string(x = "x", y = "y", label = "s_index", color = varName)) +
    #    geom_text() +
    #    scale_color_gradient(low = "blue", high = "red")
    #  return(p)
    #})
    # this is only for one time point --misleading if change thru time
    #covPlots <- grid.arrange(grobs = covPlotsL)
}