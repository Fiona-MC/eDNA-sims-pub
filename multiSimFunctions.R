library(ggplot2)
library(reshape)
library(gridExtra)

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
get_cov <- function(location, type = "peaks", coef1func = NA, coef2func = NA, t = NA, t_minus_1 = NA, params = params) {
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
  } else if (type == "ripples") {
    # x minus some number between 0 and 1 that changes with time (cycles every ~400)
    temp1 <- x / params$xdim - abs(sin((1 / 123) * t))  
    temp2 <- y / params$ydim - abs(cos((1 / 123) * t))
    return(sin(temp1 * 3 * temp2 * 3)) # the times 3 is just to get the function to behave
  } else if (type == "randomWalk") { # not sure if I'm going to use this one--requires a little rearranging on the other end
    return(t_minus_1 + rnorm(n = 1, mean = 0, sd = 0.01))
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




makePlots <- function(sim_data = sim_data, params = params, locList = locList, locPlots = c(19, 8, 6, 37, 64, 46), outDir, sitetab_sampled) {
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
      p <- ggplot(abd_av_df, aes(x = x, y = y, label = s_index, color = !!sym(varName))) + # nolint: object_usage_linter.
        geom_text() +
        scale_color_gradient(low = "blue", high = "green")
      #return(p)
    })
    
    avgAbdPlot <- arrangeGrob(grobs = avgAbdPlotsL)
    ggsave(avgAbdPlot, file = paste0(outDir, "avgAbdPlot.png"), height = 7, width = 7)
    
    # plot covs over time for spatial locations 1 through 8 (first row in current setup)
    covPlotsL <- list()
    for(cov in 1:params$numCovs){
      covPlotsL[[cov]] <- ggplot(data = sitetab_sampled[sitetab_sampled$labID <= 8, ], aes(x = Age, y = Cov1, group = labID, color = labID)) + # nolint: object_usage_linter, line_length_linter.
          geom_point() +
          geom_line()
    }
    covplots <- arrangeGrob(grobs = covPlotsL)
    ggsave(covplots, file = paste0(outDir, "covPlots.png"), height = 8, width = 10)


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

#plot P(#reads>0) vs number of indiv
plotDetectProb <- function(abdParms) {
    N_max <- 1000
    N_set <- seq(1, N_max, 10)
    P_detect <- rep(NA, times = length(N_set))
    for (i in seq_along(N_set)) {
        N <- N_set[i]
        z_set <- seq(0, N, 1)
        P_detect[i] <- sum((1 - dpois(x = 0, lambda = abdParms$readSampleRate * z_set)) * dbinom(x = z_set, size = N, p = abdParms$indivSampleProb))
    }
    p <- ggplot(data.frame(x = N_set, y = P_detect),  # nolint: object_usage_linter.
                                    aes(x, y)) + # nolint: object_usage_linter.
                                    geom_point() +
                                    labs(x = "Abundance -- N(s,t)", y = "Detection Probability")
    return(p)
}