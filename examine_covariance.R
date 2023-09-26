# plot covariance vs distance across time and space from sim data

data_dir <- "/space/s1/fiona_callahan/multiSim_ParmSet5/randomRun1/"

sim_data <- readRDS(paste0(data_dir, "sim_data.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))

locList <- readRDS(paste0(data_dir, "locList.Rdata"))
# get distMx
locDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
distMx <- as.matrix(dist(locDF, diag = TRUE, upper = TRUE))
  

getLocationPairs <- function(d, locList, distMx) {
    #TODO
    pairL <- list()

    return(pairL) # this is the list of pairs of location indices that are exactly d apart
}

# for each distance, get the set of pairs that 
distances <- unique(distMx[1, ])
avgCov <- rep(NA, times = length(distances))
for (d in distances) { 
    for (t in 100:params$num_gens) { # I think temporal autocorrelation will inflate the spatial autocorrelation on average maybe??
        locPairs <- getLocationPairs(d, locList, distMx)
        for (locs in locPairs) {
            y1 <- sim_data[[t]]$y[[locs[1]]] # these are vectors of pres abs for all species
            y2 <- sim_data[[t]]$y[[locs[2]]] 
            # 
        }
        #not done
    }

}
