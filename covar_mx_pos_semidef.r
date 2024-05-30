setwd("/home/fiona_callahan/eDNA_sims_code")
source("./multiSimFunctions.R")
source("./multiSimParms.R")


for (i in 1:100) {
    print(i)
    #Get locations
    locList <- get_locations(xdim = 100, ydim = 100, mode = "grid", xSplit = 10, ySplit = 10)
    # get distMx
    locDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
    distMx <- as.matrix(dist(locDF, diag = TRUE, upper = TRUE))

    covar_scale_space <- runif(n = 1, min = 1, max = 40)
    print(covar_scale_space)

    covar_matrix <- matrix(NA, nrow = nrow(locDF), ncol = nrow(locDF))
    for (locindex1 in seq_len(nrow(locDF))) {
    for (locindex2 in seq_len(nrow(locDF))) {
        loc1 <- as.numeric(locDF[locindex1, ])
        loc2 <- as.numeric(locDF[locindex2, ])
        covar_matrix[locindex1, locindex2] <- spatial_covariance(loc1, loc2, covar_scale_space = covar_scale_space)
    }
    }
    #print(eigen(covar_matrix)$values)

    if (any(is.complex(eigen(covar_matrix)$values))) {
    print(covar_scale_space)
    print("covariance matrix has complex eigenvals and therefore is not pos semidef (line 94 in multiSim_Mar2023.R)", call.=FALSE)
    }

    if (!(all(eigen(covar_matrix)$values >= 0))) {
    print(covar_scale_space)
    print("covariance matrix has negative values and therefore is not positive semidefinite (line 94 in multiSim_Mar2023.R)", call.=FALSE)
    }

    inv_covar_mx <- solve(covar_matrix)

    if (!isSymmetric(inv_covar_mx)) {
        print("inv_covar_mx not symmetric")
        #print(inv_covar_mx)
    }
    
    if (any(is.complex(eigen(inv_covar_mx)$values))) {
    print(covar_scale_space)
    print("covariance matrix has complex eigenvals and therefore is not pos semidef (line 94 in multiSim_Mar2023.R)", call.=FALSE)
    }

    if (!(all(eigen(inv_covar_mx)$values >= 0))) {
    print(covar_scale_space)
    print("covariance matrix has negative values and therefore is not positive semidefinite (line 94 in multiSim_Mar2023.R)", call.=FALSE)
    }

}
