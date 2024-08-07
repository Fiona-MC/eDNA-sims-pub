# plot covariance vs distance across time and space from sim data
library(ggplot2)
library(igraph)
#install.packages("ppcor")
library(ppcor)
library(data.table)

data_dir <- "/space/s1/fiona_callahan/multiSim_10sp/randomRun1/"

# Note: the sim_data thing is only needed for the first two parts with the spatiotemporal part
sim_data <- readRDS(paste0(data_dir, "sim_data.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))

locList <- readRDS(paste0(data_dir, "locList.Rdata"))
# get distMx
locDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
distMx <- as.matrix(dist(locDF, diag = TRUE, upper = TRUE))

spTemp <- FALSE
if (spTemp) {
# for each distance, get the set of pairs that 
distances <- unique(distMx[1, ])
df_L <- list()
for (d in round(distances, 2)) { 
    locPairs <- which(round(distMx, 2) == round(d, 2), arr.ind = TRUE)
    times <- seq(100, params$num_gens, by = 10)
    # initialize for results
    y_df <- data.frame(locInd1 = rep(NA, length(times) * length(locPairs[, 1])),
                        locInd2 = rep(NA, length(times) * length(locPairs[, 1])),
                        time = rep(NA, length(times) * length(locPairs[, 1])),
                        ysp1_1 = rep(NA, length(times) * length(locPairs[, 1])),
                        ysp2_1 = rep(NA, length(times) * length(locPairs[, 1])),
                        ysp3_1 = rep(NA, length(times) * length(locPairs[, 1])),
                        ysp1_2 = rep(NA, length(times) * length(locPairs[, 1])),
                        ysp2_2 = rep(NA, length(times) * length(locPairs[, 1])),
                        ysp3_2 = rep(NA, length(times) * length(locPairs[, 1])))
    df_i <- 1
    for (t in times) { # I think temporal autocorrelation will inflate the spatial autocorrelation on average maybe??
        for (i in seq_along(locPairs[, 1])) {
            y1 <- sim_data[[t]]$y[[locPairs[i, 1]]] # these are vectors of pres abs for all species
            y2 <- sim_data[[t]]$y[[locPairs[i, 2]]] 
            # store data
            y_df$locInd1[df_i] <- locPairs[i, 1]
            y_df$locInd2[df_i] <- locPairs[i, 2]
            y_df$time[df_i] <- t
            y_df$ysp1_1[df_i] <- y1[1]
            y_df$ysp2_1[df_i] <- y1[2]
            y_df$ysp3_1[df_i] <- y1[3]
            y_df$ysp1_2[df_i] <- y2[1]
            y_df$ysp2_2[df_i] <- y2[2]
            y_df$ysp3_2[df_i] <- y2[3]
            # interate index
            df_i <- df_i + 1
        }
    }
    df_L[[as.character(d)]] <- y_df
}

# Same thing but for time pairs
time_splits <- seq(1, 800, 10)
df_L_time <- list()
for (time_spl in time_splits) { 
    timePair1 <- seq(100, 1000 - time_spl, by = 5)
    timePair2 <- timePair1 + time_spl
    locs <- c(32)
    # initialize for results
    y_df_time <- data.frame(time1 = rep(NA, length(timePair1) * length(locs)),
                        time2 = rep(NA, length(timePair1) * length(locs)),
                        location = rep(NA, length(timePair1) * length(locs)),
                        ysp1_1 = rep(NA, length(timePair1) * length(locs)),
                        ysp2_1 = rep(NA, length(timePair1) * length(locs)),
                        ysp3_1 = rep(NA, length(timePair1) * length(locs)),
                        ysp1_2 = rep(NA, length(timePair1) * length(locs)),
                        ysp2_2 = rep(NA, length(timePair1) * length(locs)),
                        ysp3_2 = rep(NA, length(timePair1) * length(locs)))
    df_i <- 1
    for (l in locs) { 
        for (i in seq_along(timePair1)) {
            y1 <- sim_data[[timePair1[i]]]$y[[l]] # these are vectors of pres abs for all species
            y2 <- sim_data[[timePair2[i]]]$y[[l]] 
            # store data
            y_df_time$time1[df_i] <- timePair1[i]
            y_df_time$time2[df_i] <- timePair2[i]
            y_df_time$location[df_i] <- l
            y_df_time$ysp1_1[df_i] <- y1[1]
            y_df_time$ysp2_1[df_i] <- y1[2]
            y_df_time$ysp3_1[df_i] <- y1[3]
            y_df_time$ysp1_2[df_i] <- y2[1]
            y_df_time$ysp2_2[df_i] <- y2[2]
            y_df_time$ysp3_2[df_i] <- y2[3]
            # interate index
            df_i <- df_i + 1
        }
    }
    df_L_time[[as.character(time_spl)]] <- y_df_time
}

# species covariance in space
covars <- c()
for (d in round(distances, 2)) {
    covars <- c(covars, cov(df_L[[as.character(d)]]$ysp3_1, df_L[[as.character(d)]]$ysp3_2))
}
covarData <- data.frame(distance = distances, covariance = covars)

ggplot(covarData, aes(x = distance, y = covariance)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = "blue") +
    geom_vline(xintercept = 0, color = "blue") +
    labs(y = "Covariance of y for species 3") +
    theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16)) 


# species covariance time
covars <- c()
for (time_spl in as.character(time_splits)) {
    covars <- c(covars, cov(df_L_time[[as.character(time_spl)]]$ysp3_1, df_L_time[[as.character(time_spl)]]$ysp3_2))
}
covarData <- data.frame(time_split = time_splits, covariance = covars)

ggplot(covarData, aes(x = time_splits, y = covariance)) + 
    geom_point() +
    geom_hline(yintercept = 0, color = "blue") +
    geom_vline(xintercept = 0, color = "blue") +
    labs(y = "Covariance of y for species 3") +
    theme(axis.text.x = element_text(size = 12),  
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16)) 
}








#######################################################
################## covar by interaction ###############
######################################################

getCovarData <- function(data_dirs, logi = FALSE, prec = FALSE, abd = TRUE, cluster = FALSE, covariates = FALSE, actualMx = FALSE) {
    plotData_all <- data.frame()
    for(data_dir in data_dirs) {
        params <- readRDS(paste0(data_dir, "params.Rdata"))
        # question: does pairwise covariance (between species) correlate to actual connections in the large networks
        if (logi) {
            if (abd) {
                sitetab_path <- paste0(data_dir, "sitetab_abd_logi.csv")
                if (prec) {
                    sitetab_path <- paste0(data_dir, "sitetab_abd_logi_prec.csv")
                }
            } else {
                sitetab_path <- paste0(data_dir, "sitetab_logi.csv")
                if (prec) {
                    sitetab_path <- paste0(data_dir, "sitetab_logi_prec.csv")
                }
            }
        } else {
            sitetab_path <- paste0(data_dir, "sim_sitetab_sampled10000.csv")
        }

        if (!logi && abd) {
            # what about if we are using abundances? vvv
            sitetab_path <- paste0(data_dir, "sim_sitetab_readAbd_sampled10000.csv")
        }

        sitetab <- read.csv(sitetab_path)
        print(sitetab_path)

        if (logi) {
            if (prec) {
                actual_covar_mx <- readRDS(paste0(data_dir, "logi_covar_mx_prec.Rdata"))
            } else {
                actual_covar_mx <- readRDS(paste0(data_dir, "logi_covar_mx.Rdata"))
            }
        }

        corr_mx <- matrix(NA, nrow = params$numSpecies, ncol = params$numSpecies)
        covar_mx <- matrix(NA, nrow = params$numSpecies, ncol = params$numSpecies)
        pcor_mx <- matrix(NA, nrow = params$numSpecies, ncol = params$numSpecies)
        for (sp1 in 1:params$numSpecies) {
            for (sp2 in 1:params$numSpecies) {
                covar_mx[sp1, sp2] <- cor(sitetab[, paste0("Sp", sp1)], sitetab[, paste0("Sp", sp2)])
                corr_mx[sp1, sp2] <- cor(sitetab[, paste0("Sp", sp1)], sitetab[, paste0("Sp", sp2)])
                if (covariates) {
                    # what about after controlling for covariates? 
                    pcors <- pcor(sitetab[, c(paste0("Sp", sp1), paste0("Sp", sp2), params$names_cov)])
                    pcorsEst <- pcors$estimate
                    pcor_mx[sp1, sp2] <- pcorsEst[1, 2]
                }
            }
        }

        if (covariates) {
            # what about after controlling for covariates
            covar_mx <- pcor_mx
        }

        actual_alpha <- sign(params$alpha)
        
        if (cluster) {
            # what about if actual alpha is the clustering version? vvv
            alphaG <- graph_from_adjacency_matrix(actual_alpha, mode = "undirected")
            actual_alpha <- (distances(alphaG, v = 1:params$numSpecies, to = 1:params$numSpecies) != Inf) * 
                                    (diag(nrow = dim(actual_alpha)[1], ncol = dim(actual_alpha)[1]) == 0)
        }

        if (!actualMx) {
            covar_pos_interact <- as.numeric(covar_mx * (actual_alpha > 0))
            covar_pos_interact <- covar_pos_interact[covar_pos_interact != 0]

            covar_neg_interact <- as.numeric(covar_mx * (actual_alpha < 0))
            covar_neg_interact <- covar_neg_interact[covar_neg_interact != 0]

            covar_no_interact <- as.numeric((covar_mx + 0.00001) * (actual_alpha == 0))
            # remove diagonal 
            covar_no_interact <- covar_no_interact * (diag(nrow = params$numSpecies) == 0)
            covar_no_interact <- covar_no_interact[covar_no_interact != 0]
            plotData <- data.frame(actual_interaction = c(rep("positive", times = length(covar_pos_interact)), 
                                                    rep("negative", times = length(covar_neg_interact)), 
                                                    rep("none", times = length(covar_no_interact))),
                                covariance = c(covar_pos_interact, covar_neg_interact, covar_no_interact))
            plotData_all <- rbind(plotData_all, plotData)
        }

        if (actualMx) {
            covar_mx <- actual_covar_mx
            covar_pos_interact <- as.numeric(covar_mx * (actual_alpha > 0))
            covar_pos_interact <- covar_pos_interact[covar_pos_interact != 0]

            covar_neg_interact <- as.numeric(covar_mx * (actual_alpha < 0))
            covar_neg_interact <- covar_neg_interact[covar_neg_interact != 0]

            # hacky way to get the ones that are 0 in alpha but not make it all 0's if covariance is actually 0
            covar_no_interact <- as.numeric((covar_mx + 0.00001) * (actual_alpha == 0))
            # remove diagonal 
            covar_no_interact <- covar_no_interact * (diag(nrow = params$numSpecies) == 0)
            covar_no_interact <- covar_no_interact[covar_no_interact != 0]
            plotData <- data.frame(actual_interaction = c(rep("positive", times = length(covar_pos_interact)), 
                                                    rep("negative", times = length(covar_neg_interact)), 
                                                    rep("none", times = length(covar_no_interact))),
                                covariance = c(covar_pos_interact, covar_neg_interact, covar_no_interact))
            plotData_all <- rbind(plotData_all, plotData)
        }
    }
    
    return(plotData_all)
}


data_dir <- "/space/s1/fiona_callahan/multiSim_100sp/"
abd = TRUE
data_dirs <- sapply(X = 1:100, FUN = function(num) {paste0(data_dir, "randomRun", num, "/")})
plotData <- getCovarData(data_dirs, logi = FALSE, prec = FALSE, abd = abd, cluster = FALSE, covariates = FALSE, actualMx = FALSE)
# what if we do absolute value of covariance?
#plotData <- data.frame(actual_interaction = c(rep("positive", times = length(covar_pos_interact)), 
#                                            rep("negative", times = length(covar_neg_interact)), 
#                                            rep("none", times = length(covar_no_interact))),
#                        covariance = abs(c(covar_pos_interact, covar_neg_interact, covar_no_interact)))

ggplot(data = plotData[plotData$covariance < 1, ], aes(x = actual_interaction, y = covariance)) +
            geom_boxplot() #+ 
            #ggtitle(paste0(data_dir))

plot1 <- ggplot(data = plotData[abs(plotData$covariance) < 1000, ], aes(x = covariance, fill = actual_interaction)) +
            geom_histogram(position = "dodge", stat = "density") #+
            #ggtitle(paste0(data_dir))
plot1
if (abd) {
    ggsave(plot1, filename = paste0(data_dir, "empiricalCovariance_smooth_abd.png"))
} else {
    ggsave(plot1, filename = paste0(data_dir, "empiricalCovariance_smooth_presAbs.png"))
}


plot2 <- ggplot(data = plotData[abs(plotData$covariance) < 100, ], aes(x = covariance, fill = actual_interaction)) +
            geom_histogram(position = "dodge") #+
            #ggtitle(paste0(data_dir))
plot2

if (abd) {
    ggsave(plot2, filename = paste0(data_dir, "empiricalCovariance_abd.png"))
} else {
    ggsave(plot2, filename = paste0(data_dir, "empiricalCovariance_presAbs.png"))
}

