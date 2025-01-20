library(ggplot2)
library(tidyr)
library(stringr)
library(SpiecEasi)
library(igraph)

numRuns <- 100
numSamples <- 100
logi <- FALSE
filtered <- FALSE

dirName <- c("multiSim_10sp_random")
#dirName <- c("multiSim_test2x10sp")
multiSimRes <- data.frame()
#resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
#resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
#              "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")

#ls /space/s1/fiona_callahan/sim_paper_stuff/multiSim_100
logistic_cutoffs <- c(0.05)
#log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _cov_2runs_cutoff", x, ".csv")})
#log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _noCov_2runs_cutoff", x, ".csv")})

# "/space/s1/fiona_callahan/multiSim_10sp/spiecEasi_res_sampled1000_glasso_filtered100_infResGathered_100sims.csv"

if (logi) {
  log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                              paste0("logistic_mistakes_sampled", numSamples, "_cov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
  log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                              paste0("logistic_mistakes_sampled", numSamples, "_noCov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
} else {
  if (filtered) {
    log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                          paste0("logistic_mistakes_sampled", numSamples, "_cov_", numRuns, "runs_filtered_cutoff", x, "_", numRuns, "sims.csv")})
    log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                          paste0("logistic_mistakes_sampled", numSamples, "_noCov_", numRuns, "runs_filtered_cutoff", x, "_", numRuns, "sims.csv")})

  } else {
    log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                              paste0("logistic_mistakes_sampled", numSamples, "_cov_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
    log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                              paste0("logistic_mistakes_sampled", numSamples, "_noCov_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})

  }

  #if (numRuns == 100 && numSamples == 100 && dirName == c("multiSim_10sp")) {
  #  log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
  #                            paste0("logistic_mistakes_sampled", numSamples, "_cov_", numRuns, "runs_cutoff", x, ".csv")})
  #  log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
  #                            paste0("logistic_mistakes_sampled", numSamples, "_noCov_", numRuns, "runs_cutoff", x, ".csv")})
  #}
}

#inla_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_infResGathered_cutoff", x, ".csv")})

#inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paperSep_infResGathered_cutoff", x, ".csv")})

inla_cutoffs <- c(0.05)
if (logi) {
  inla_resnames_cov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_cov_logi_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
  inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_noCov_logi_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
} else {
  if(filtered) {
    inla_resnames_cov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_cov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
    inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_noCov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
  } else {
    inla_resnames_cov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_cov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
    inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_noCov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
  }
}

#inla_cutoffs <- c(0, 0.0000001, .3, .5, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_sampled500_infResGathered_cutoff", x, ".csv")})

#inla_cutoffs <- c(0, 0.0000001, .3, .5, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames_500 <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_sampled500_infResGathered_cutoff", x, ".csv")})
if(logi) {
  seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb_logi")
  seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso_logi")
  seName3 <- paste0("spiecEasi_res_sampled", numSamples, "_sparcc_logi")
} else {
  if (filtered) {
    seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb_filtered")
    seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso_filtered")
    seName3 <- paste0("spiecEasi_res_sampled", numSamples, "_sparcc_filtered") 
  } else {
    seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb")
    seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso")
    seName3 <- paste0("spiecEasi_res_sampled", numSamples, "_sparcc")
  }
}

if (logi) {
  ecName1 <- paste0("ecoCopula_res_sampled", numSamples, "_noCov_logi")
  ecName2 <- paste0("ecoCopula_res_sampled", numSamples, "_cov_logi")
} else {
  if(filtered) {
    ecName1 <- paste0("ecoCopula_res_sampled", numSamples, "_noCov_filtered")
    ecName2 <- paste0("ecoCopula_res_sampled", numSamples, "_cov_filtered")
  } else {
    ecName1 <- paste0("ecoCopula_res_sampled", numSamples, "_noCov")
    ecName2 <- paste0("ecoCopula_res_sampled", numSamples, "_cov")
  }
}

if (logi) {
  inla1 <- paste0("INLA_res_paperSep_sampled", numSamples, "_logi_infResGathered_", numRuns, "sims.csv")
  jags1 <- paste0("JAGS_infResGathered_sampled", numSamples, "_logi_", numRuns, "sims.csv")
} else {
  if (filtered) {
    inla1 <- paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_infResGathered_", numRuns, "sims.csv")
    jags1 <- NA
  } else {
    inla1 <- paste0("INLA_res_paperSep_sampled", numSamples, "_infResGathered_", numRuns, "sims.csv")
    jags1 <- paste0("JAGS_infResGathered_sampled", numSamples, "_", numRuns, "sims.csv")
  }
}

resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName3, "_infResGathered_", numRuns, "sims.csv"),
              inla1, 
              inla_resnames_noCov,
              inla_resnames_cov,
              log_resnames_cov,
              log_resnames_noCov,
              jags1)

#if (filtered) {
#  seName1 <- "spiecEasi_res_sampled1000_mb_filtered100"
#  seName2 <- "spiecEasi_res_sampled1000_glasso_filtered100"
#  inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
#                        paste0("INLA_res_paperSep_sampled", numSamples, "_noCov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
#  resNames <- c(log_resnames_noCov, 
#              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
#              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"),
#              inla_resnames_noCov)
#}

# check if runs exist
thisDir <-  paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/")
for (i in seq_along(resNames)) {
  if(!file.exists(paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/", resNames[i]))) {
    print("This file does not exist:")
    print(paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/", resNames[i]))
  }
}

#for (i in seq_along(resNames)) {
#  if(file.exists(paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/", resNames[i]))) {
#    print("This file exists:")
#    print(paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/", resNames[i]))
##  }
#}

#resNames <- c(inla_resnames,
 #             inla_resnames_500)

#logistic_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
#log_resnames <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_logi_cutoff", x, ".csv")})
#resNames <- log_resnames

file.exists(paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/", paste0("JAGS_infResGathered_sampled", numSamples, "_", numRuns, "sims.csv")))

thisDir <-  paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/")
# load results into list
multiSimResL <- list()
for (i in seq_along(resNames)) {
  if(file.exists(paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/", resNames[i]))) {
    thisMultiSimRes <- read.csv(paste0("/space/s1/fiona_callahan/sim_paper_stuff/", dirName, "/", resNames[i]), header = TRUE)
    thisMultiSimRes$totalMistakes <- thisMultiSimRes$num_incorrectInferences + thisMultiSimRes$num_missedEffectsL
    thisMultiSimRes$fp_fp_tp <- thisMultiSimRes$num_incorrectInferences / 
                                (thisMultiSimRes$num_correctInferences + thisMultiSimRes$num_incorrectInferences)
    #thisMultiSimRes$fp_fp_tp_cluster <- thisMultiSimRes$num_incorrect_cluster / 
    #                            (thisMultiSimRes$num_correct_cluster + thisMultiSimRes$num_incorrect_cluster)
    actual_mean_fpr <- rep(NA, times = dim(thisMultiSimRes)[1])
    if("fpr.mode" %in% names(thisMultiSimRes)) {
      for (row in seq_len(dim(thisMultiSimRes)[1])) { #1:dim(multiSimRes)[1]
        if(is.na(thisMultiSimRes$fpr.mode[row])) { # TODO this is an error with the jags one 
          thisFPR <- NA
        } else {
          if(thisMultiSimRes$fpr.mode[row] == "constant") {
            thisFPR <- thisMultiSimRes$fpr.constant_fpr[row]
          } else if (thisMultiSimRes$fpr.mode[row] == "none") {
            thisFPR <- 0
          } else if (thisMultiSimRes$fpr.mode[row] == "dependent_sp") {
            thisFPR <- thisMultiSimRes$fpr.mean_fpr[row]
          } else if (thisMultiSimRes$fpr.mode[row] == "independent") {
            thisFPR <- thisMultiSimRes$fpr.mean_fpr[row]
          }
        }
          actual_mean_fpr[row] <- thisFPR
      }
    } 
    thisMultiSimRes$actual_mean_fpr <- actual_mean_fpr
    if (str_detect(resNames[i], "JAGS")) {
      for (cutoffIndex in unique(thisMultiSimRes$cutoff_val)) {
        multiSimResL[[paste0(resNames[i], "_cutoffIndex", cutoffIndex)]] <- thisMultiSimRes[thisMultiSimRes$cutoff_val == cutoffIndex, ]
      }
    } else {
      multiSimResL[[resNames[i]]] <- thisMultiSimRes
    }
  }

}

get_falseDiscovery <- function(data = multiSimRes, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
      TP <- data$TP_cluster
      FP <- data$FP_cluster
    } else if (mode == "ignore_sign") {
      TP <- data$TP_ignoreSign
      FP <- data$FP_ignoreSign 
    } else {      
      TP <- data$TP_sign
      FP <- data$FP_sign 
    }
    FDR <- FP / (FP + TP)
    if (return_components) {
      return(list(FDR = FDR, FP = FP, TP = TP))
    } else {
      return(FDR)
    }
}


# method | threshold | avg_TPR | avg_FPR 
file <- rep(NA, times = length(multiSimResL))
methods <- rep(NA, times = length(multiSimResL))
modelSelect <- rep(NA, times = length(multiSimResL))
numSamplesL <- rep(NA, times = length(multiSimResL))
numSims <- rep(NA, times = length(multiSimResL))
numSp <- rep(NA, times = length(multiSimResL))
thresholds <- rep(NA, times = length(multiSimResL))
avg_FDR <- rep(NA, times = length(multiSimResL))

for (i in seq_along(multiSimResL)) {
    multiSimRes <- multiSimResL[[i]]
    file[i] <- names(multiSimResL)[i]
    if (str_detect(names(multiSimResL)[i], "logistic")) {
      if (str_detect(names(multiSimResL)[i], "noCov")) {
        methods[i] <- "logistic_noCov"
        modelSelect[i] <- FALSE
      } else {
        methods[i] <- "logistic"
        modelSelect[i] <- FALSE
      }
    } else if (str_detect(names(multiSimResL)[i], "INLA")) {
      if (str_detect(names(multiSimResL)[i], "noCov")) {
        thisName <- "INLA_noCov"
      } else {
         thisName <- "INLA_cov"
      }
      if (str_detect(names(multiSimResL)[i], "cutoff")) {
        modelSelect[i] <- FALSE
      } else {
        modelSelect[i] <- TRUE
      }
      methods[i] <- thisName
    } else if (str_detect(names(multiSimResL)[i], "ecoCopula")) {
      if (str_detect(names(multiSimResL)[i], "noCov")) {
        methods[i] <- "ecoCopula_noCov"
      } else {
         methods[i] <- "ecoCopula_cov"
      }
      modelSelect[i] <- TRUE
    } else if (str_detect(names(multiSimResL)[i], "spiecEasi")) {
      if (str_detect(names(multiSimResL)[i], "mb")) {
        methods[i] <- "spiecEasi_mb"
      } else if (str_detect(names(multiSimResL)[i], "glasso")) {
        methods[i] <- "spiecEasi_glasso"
      } else {
        methods[i] <- "spiecEasi_sparcc"
      }
      modelSelect[i] <- TRUE
    } else if (str_detect(names(multiSimResL)[i], "JAGS")) {
      methods[i] <- "JAGS"
      modelSelect[i] <- FALSE
    }

    if (str_detect(names(multiSimResL)[i], "sampled")) {
      thisNumSamples <- str_extract(names(multiSimResL)[i], "sampled(\\d+)")
      thisNumSamples <- str_extract(thisNumSamples, "\\d+")
      numSamplesL[i] <- thisNumSamples
    } 

    if (str_detect(names(multiSimResL)[i], "sims")) {
      thisNumSims <- str_extract(names(multiSimResL)[i], "(\\d+)sims")
      thisNumSims <- str_extract(thisNumSims, "\\d+")
      numSims[i] <- thisNumSims
    } 

    if (str_detect(names(multiSimResL)[i], "sp")) {
      thisNumSp <- str_extract(names(multiSimResL)[i], "(\\d+)sp")
      thisNumSp <- str_extract(thisNumSp, "\\d+")
      numSp[i] <- thisNumSp
    } 

    #if (str_detect(names(multiSimResL)[i], "cutoff")) {
    #  thresholds[i] <- str_split(c(names(multiSimResL)[i]), pattern = ".")
    #}

    FDR_temp <- get_falseDiscovery(data = multiSimRes, mode = "ignore_sign")
    avg_FDR[i] <- mean(FDR_temp, na.rm = TRUE)
}


res_summary <- data.frame(file = file,
                            method = methods,
                            modelSelect = modelSelect,
                            numSamples = numSamplesL,
                            numSims = numSims,
                            numSpecies = numSp, 
                            threshold = thresholds,
                            avg_FDR = avg_FDR)

res_summary[,c(2,4,5,8)]

