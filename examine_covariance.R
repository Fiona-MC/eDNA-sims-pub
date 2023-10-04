# plot covariance vs distance across time and space from sim data
library(ggplot2)

data_dir <- "/space/s1/fiona_callahan/multiSim11/randomRun1/"

sim_data <- readRDS(paste0(data_dir, "sim_data.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))

locList <- readRDS(paste0(data_dir, "locList.Rdata"))
# get distMx
locDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
distMx <- as.matrix(dist(locDF, diag = TRUE, upper = TRUE))

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

# species covariance
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
