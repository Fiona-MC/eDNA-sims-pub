library(ggplot2)
library(tidyr)
library(stringr)
library(SpiecEasi)
library(igraph)

cluster <- FALSE
ratio_of_avg <- FALSE #do we compute the average of the ratio or ratio of averages
numRuns <- 100
numSamples <- 10000
logi <- FALSE
saveRes <- FALSE
covMode <- "noCount" # "all" "noCov" "cov" "covNoCount" "noCount"

if (logi) {
  filtered <- FALSE
} else {
  filtered <- TRUE
}

dirName <- c("multiSim_100sp_random_moreSamples")
#dirName <- c("multiSim_test2x10sp")
multiSimRes <- data.frame()
#resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
#resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
#              "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")

#ls /space/s1/fiona_callahan/multiSim_100
logistic_cutoffs <- c(0, 1e-128, 1e-64, 1e-32, 1e-16, 1e-8, 1e-4, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15,
                      0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.9999, 0.99999, 0.999999, 0.9999999, 0.99999999, 1, 1.01)

if (logi) {
  log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                  paste0("logistic_mistakes_sampled", numSamples, "_cov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
  log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                                paste0("logistic_mistakes_sampled", numSamples, "_noCov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
  log_resnames_covNoCount <- sapply(X = logistic_cutoffs, FUN = function(x) {
                      paste0("logistic_mistakes_sampled", numSamples, "_covNoCount_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})     
} else {
  if (filtered) {
    log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                          paste0("logistic_mistakes_sampled", numSamples, "_cov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
    log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                          paste0("logistic_mistakes_sampled", numSamples, "_noCov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
    log_resnames_covNoCount <- sapply(X = logistic_cutoffs, FUN = function(x) {
                      paste0("logistic_mistakes_sampled", numSamples, "_covNoCount_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})                          
  } else {
    log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                              paste0("logistic_mistakes_sampled", numSamples, "_cov_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
    log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                              paste0("logistic_mistakes_sampled", numSamples, "_noCov_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})

  }
}

linear_cutoffs <- c(0, 1e-128, 1e-64, 1e-32, 1e-16, 1e-8, 1e-4, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15,
                      0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.9999, 0.99999, 0.999999, 0.9999999, 0.99999999, 1, 1.01)

if (logi) {
  lin_resnames_cov <- sapply(X = linear_cutoffs, FUN = function(x) {
                  paste0("linearReg_mistakes_sampled", numSamples, "_cov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
  lin_resnames_noCov <- sapply(X = linear_cutoffs, FUN = function(x) {
                                paste0("linearReg_mistakes_sampled", numSamples, "_noCov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
  lin_resnames_covNoCount <- sapply(X = linear_cutoffs, FUN = function(x) {
                      paste0("linearReg_mistakes_sampled", numSamples, "_covNoCount_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})     
} else {
  if (filtered) {
    lin_resnames_cov <- sapply(X = linear_cutoffs, FUN = function(x) {
                          paste0("linearReg_mistakes_sampled", numSamples, "_cov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
    lin_resnames_noCov <- sapply(X = linear_cutoffs, FUN = function(x) {
                          paste0("linearReg_mistakes_sampled", numSamples, "_noCov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
    lin_resnames_covNoCount <- sapply(X = linear_cutoffs, FUN = function(x) {
                      paste0("linearReg_mistakes_sampled", numSamples, "_covNoCount_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})                          
  } else {
    lin_resnames_cov <- sapply(X = linear_cutoffs, FUN = function(x) {
                              paste0("linearReg_mistakes_sampled", numSamples, "_cov_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
    lin_resnames_noCov <- sapply(X = linear_cutoffs, FUN = function(x) {
                              paste0("linearReg_mistakes_sampled", numSamples, "_noCov_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
  }
}

#inla_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_infResGathered_cutoff", x, ".csv")})

#inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paperSep_infResGathered_cutoff", x, ".csv")})

inla_cutoffs <- c(0, 1, 0.0000001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, .3, .5)
if (logi) {
  inla_resnames_cov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_logi_cov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
  inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_logi_noCov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
  inla_resnames_covNoCount <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_logi_covNoCount_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
} else {
  if(filtered) {
    inla_resnames_cov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_cov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
    inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_noCov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
    inla_resnames_covNoCount <- sapply(X = inla_cutoffs, FUN = function(x) {
                        paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_covNoCount_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
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
} else {
  if (filtered) {
    seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb_filtered")
    seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso_filtered")
  } else {
    seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb")
    seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso")
  }
}

sparcc_cutoffs <- c(0.00001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 0.005, 0.04, 0.06, 0.08, 0.125, 0.15, 0.175, 0.25)
if (logi) {
  sparcc_resnames <- sapply(X = sparcc_cutoffs, FUN = function(x) {
                        paste0("sparcc_res_sampled", numSamples, "_logi_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
} else {
  if(filtered) {
    sparcc_resnames <- sapply(X = sparcc_cutoffs, FUN = function(x) {
                        paste0("sparcc_res_sampled", numSamples, "_filtered_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
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
    jags1 <- paste0("JAGS_infResGathered_sampled", numSamples, "_filtered_", numRuns, "sims.csv")
  } else {
    inla1 <- paste0("INLA_res_paperSep_sampled", numSamples, "_infResGathered_", numRuns, "sims.csv")
    jags1 <- paste0("JAGS_infResGathered_sampled", numSamples, "_", numRuns, "sims.csv")
  }
}

if(covMode == "all") { 
  se_include <- TRUE
  resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              sparcc_resnames,
              inla1, 
              inla_resnames_noCov,
              inla_resnames_cov,
              inla_resnames_covNoCount,
              log_resnames_cov,
              log_resnames_noCov,
              log_resnames_covNoCount,
              lin_resnames_cov,
              lin_resnames_noCov,
              lin_resnames_covNoCount,
              jags1)
} else if (covMode == "cov") {
  se_include <- FALSE
  resNames <- c(inla1, 
                inla_resnames_cov,
                log_resnames_cov,
                lin_resnames_cov)
} else if (covMode == "noCov") {
  se_include <- TRUE
  resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              sparcc_resnames,
              inla_resnames_noCov,
              log_resnames_noCov,
              lin_resNames_noCov)
} else if (covMode == "covNoCount") {
  se_include <- FALSE
  resNames <- c(paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              inla_resnames_covNoCount,
              log_resnames_covNoCount,
              lin_resnames_covNoCount,
              jags1)
} else if (covMode == "noCount") {
  se_include <- TRUE
  resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              sparcc_resnames,
              inla1, 
              inla_resnames_noCov,
              inla_resnames_covNoCount,
              log_resnames_noCov,
              log_resnames_covNoCount,
              lin_resnames_noCov,
              lin_resnames_covNoCount,
              jags1)
} else {
  print("covMode not a valid option")
}
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
thisDir <-  paste0("/space/s1/fiona_callahan/", dirName, "/")
for (i in seq_along(resNames)) {
  if(!file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))) {
    print("This file does not exist:")
    print(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))
  }
}

#for (i in seq_along(resNames)) {
#  if(file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))) {
#    print("This file exists:")
#    print(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))
##  }
#}

#resNames <- c(inla_resnames,
 #             inla_resnames_500)

#logistic_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
#log_resnames <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_logi_cutoff", x, ".csv")})
#resNames <- log_resnames


thisDir <-  paste0("/space/s1/fiona_callahan/", dirName, "/")
# load results into list
multiSimResL <- list()
for (i in seq_along(resNames)) {
  if(file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))) {
    thisMultiSimRes <- read.csv(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]), header = TRUE)
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
# default ones
#multiSimRes <- multiSimResL$logistic_mistakes_sampled50_noCov_100runs_cutoff0.9999_100sims.csv
#multiSimRes <- multiSimResL$spiecEasi_res_mb_infResGathered.csv

# if you need to check these later, look at obsidian ROC curve note
# obsidian://open?vault=simulations&file=ROC%20curves
get_TPR <- function(data = multiSimRes, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
        TP <- data$TP_cluster 
        FN <- data$FN_cluster
    } else if (mode == "ignore_sign") {
        TP <- data$TP_ignoreSign
        FN <- data$FN_ignoreSign
    } else {
        # tp/(tp+fn)
        TP <- data$TP_sign
        FN <- data$FN_sign
    } 
    TPR <- TP / (TP + FN)
    if (return_components) {
      return(list(TPR = TPR, TP = TP, FN = FN))
    } else {
      return(TPR)
    }
}

get_FPR <- function(data = multiSimRes, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
      FP <- data$FP_cluster
      TN <- data$TN_cluster
    } else if (mode == "ignore_sign") {
      FP <- data$FP_ignoreSign
      TN <- data$TN_ignoreSign
    } else {
      FP <- data$FP_sign
      TN <- data$TN_sign
    } 
    FPR <- FP / (FP + TN)
    if (return_components) {
      return(list(FPR = FPR, FP = FP, TN = TN))
    } else {
      return(FPR)
    }
}

get_precision <- function(data = multiSimRes, mode = "ignore_sign", return_components = FALSE) {
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
    prec <- TP / (TP + FP)
    if (return_components) {
      return(list(prec = prec, TP = TP, FP = FP))
    } else {
      return(prec)
    }
}

get_recall <- function(data = multiSimRes, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
      TP <- data$TP_cluster
      FN <- data$FN_cluster
    } else if (mode == "ignore_sign") {
      TP <- data$TP_ignoreSign
      FN <- data$FN_ignoreSign 
    } else {      
      TP <- data$TP_sign
      FN <- data$FN_sign 
    }
    recall <- TP / (TP + FN)
    if (return_components) {
      return(list(recall = recall, TP = TP, FN = FN))
    } else {
      return(recall)
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
avg_TPR <- rep(NA, times = length(multiSimResL))
avg_FPR <- rep(NA, times = length(multiSimResL))
TPR_sd <- rep(NA, times = length(multiSimResL))
FPR_sd <- rep(NA, times = length(multiSimResL))
avg_precision <- rep(NA, times = length(multiSimResL))
avg_recall <- rep(NA, times = length(multiSimResL))

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
    } else if(str_detect(names(multiSimResL)[i], "sparcc")) {
      methods[i] <- "sparcc"
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
    if (cluster) {
      if(ratio_of_avg) {
        # TPR <- TP / (TP + FN)
        avg_TPR[i] <- mean(get_TPR(data = multiSimRes, mode = "cluster", return_components = TRUE)$TP, na.rm = TRUE) / 
                      (mean(get_TPR(data = multiSimRes, mode = "cluster", return_components = TRUE)$TP, na.rm = TRUE) + 
                            mean(get_TPR(data = multiSimRes, mode = "cluster", return_components = TRUE)$FN, na.rm = TRUE))
        # FPR <- FP / (FP + TN)
        avg_FPR[i] <- mean(get_FPR(data = multiSimRes, mode = "cluster", return_components = TRUE)$FP, na.rm = TRUE) / 
                      (mean(get_FPR(data = multiSimRes, mode = "cluster", return_components = TRUE)$FP, na.rm = TRUE) + 
                            mean(get_FPR(data = multiSimRes, mode = "cluster", return_components = TRUE)$TN, na.rm = TRUE))
        # prec <- TP / (TP + FP)
        avg_precision[i] <-  mean(get_precision(data = multiSimRes, mode = "cluster", return_components = TRUE)$TP, na.rm = TRUE) / 
                      (mean(get_precision(data = multiSimRes, mode = "cluster", return_components = TRUE)$TP, na.rm = TRUE) + 
                            mean(get_precision(data = multiSimRes, mode = "cluster", return_components = TRUE)$FP, na.rm = TRUE))
        # recall <- TP / (TP + FN)
        avg_recall[i] <-  mean(get_recall(data = multiSimRes, mode = "cluster", return_components = TRUE)$TP, na.rm = TRUE) / 
                      (mean(get_recall(data = multiSimRes, mode = "cluster", return_components = TRUE)$TP, na.rm = TRUE) + 
                            mean(get_recall(data = multiSimRes, mode = "cluster", return_components = TRUE)$FN, na.rm = TRUE)) 
      } else {
        TPR_temp <- get_TPR(data = multiSimRes, mode = "cluster")
        avg_TPR[i] <- mean(TPR_temp, na.rm = TRUE)
        TPR_sd[i] <- sqrt(var(TPR_temp, na.rm = TRUE) / length(TPR_temp))
        FPR_temp <- get_FPR(data = multiSimRes, mode = "cluster")
        avg_FPR[i] <- mean(FPR_temp, na.rm = TRUE)
        FPR_sd[i] <- sqrt(var(FPR_temp, na.rm = TRUE) / length(FPR_temp))
        avg_precision[i] <- mean(get_precision(data = multiSimRes, mode = "cluster"), na.rm = TRUE)
        avg_recall[i] <- mean(get_recall(data = multiSimRes, mode = "cluster"), na.rm = TRUE)
      }
    } else {
      if(ratio_of_avg) {
        # TPR <- TP / (TP + FN)
        avg_TPR[i] <- mean(get_TPR(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$TP, na.rm = TRUE) / 
                      (mean(get_TPR(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$TP, na.rm = TRUE) + 
                            mean(get_TPR(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$FN, na.rm = TRUE))
        # FPR <- FP / (FP + TN)
        avg_FPR[i] <- mean(get_FPR(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$FP, na.rm = TRUE) / 
                      (mean(get_FPR(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$FP, na.rm = TRUE) + 
                            mean(get_FPR(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$TN, na.rm = TRUE))
        # prec <- TP / (TP + FP)
        avg_precision[i] <-  mean(get_precision(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$TP, na.rm = TRUE) / 
                      (mean(get_precision(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$TP, na.rm = TRUE) + 
                            mean(get_precision(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$FP, na.rm = TRUE))
        # recall <- TP / (TP + FN)
        avg_recall[i] <-  mean(get_recall(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$TP, na.rm = TRUE) / 
                      (mean(get_recall(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$TP, na.rm = TRUE) + 
                            mean(get_recall(data = multiSimRes, mode = "ignore_sign", return_components = TRUE)$FN, na.rm = TRUE)) 
      } else {
        TPR_temp <- get_TPR(data = multiSimRes, mode = "ignore_sign")
        avg_TPR[i] <- mean(TPR_temp, na.rm = TRUE)
        TPR_sd[i] <- sqrt(var(TPR_temp, na.rm = TRUE) / length(TPR_temp))
        FPR_temp <- get_FPR(data = multiSimRes, mode = "ignore_sign")
        avg_FPR[i] <- mean(FPR_temp, na.rm = TRUE)
        FPR_sd[i] <- sqrt(var(FPR_temp, na.rm = TRUE) / length(FPR_temp))
        avg_precision[i] <- mean(get_precision(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
        avg_recall[i] <- mean(get_recall(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
      }
    }
}

if(se_include) {
multiples <- FALSE
if(multiples) {
if (!cluster) {
  # for spiec-easi, just put separate lines for a couple of individual sims rather than average
  for (run in sample(1:numRuns, size = 5)) {
    subdir <- paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/", seName2, "/trial1/")
    se <- readRDS(paste0(subdir, "se_rawRes.Rdata"))
    if (filtered) {
      params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/paramsFiltered", numSamples, ".Rdata"))
      actualAlpha <- params$filteredAlpha
    } else {
      params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
      actualAlpha <- params$alpha
    }

    alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
    # theta = true_interactions
    se.roc <- huge::huge.roc(se$est$path, theta = alphaG, verbose = FALSE)
    avg_TPR <- c(avg_TPR, se.roc$tp)
    avg_FPR <- c(avg_FPR, se.roc$fp)

    TPR_sd <- c(TPR_sd, rep(NA, times = length(se.roc$fp))) 
    FPR_sd <- c(FPR_sd, rep(NA, times = length(se.roc$fp)))

    modelSelect <- c(modelSelect, rep(FALSE, times = length(se.roc$fp)))
    numSamplesL <- c(numSamplesL, rep(NA, times = length(se.roc$fp)))
    numSims <- c(numSims, rep(NA, times = length(se.roc$fp)))
    numSp <- c(numSp, rep(NA, times = length(se.roc$fp)))

    avg_precision <- c(avg_precision, rep(NA, times = length(se.roc$fp)))
    avg_recall <- c(avg_recall, rep(NA, times = length(se.roc$fp)))
    methods <- c(methods, rep(paste0("SpiecEasi", run), times = length(se.roc$fp)))
    file <- c(file, rep(NA, times = length(se.roc$fp)))
    thresholds <- c(thresholds, rep(NA, times = length(se.roc$fp)))
  } 
}
}


if (!cluster) { 
  if(file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", paste0(seName1, "_infResGathered_", numRuns, "sims.csv")))) {

    # for spiecEasi, get average for mb across a lot of lambdas
    ######################################
    TPRs <- matrix(nrow = 100, ncol = numRuns) # rows are lambda values, columns are runs
    FPRs <- matrix(nrow = 100, ncol = numRuns)
    for (run in 1:numRuns) {
      subdir <- paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/", seName1, "/trial1/")
      if (file.exists(paste0(subdir, "se_rawRes.Rdata"))) {
        se <- readRDS(paste0(subdir, "se_rawRes.Rdata"))
        if (filtered) {
          params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/paramsFiltered", numSamples, ".Rdata"))
          actualAlpha <- params$filteredAlpha
        } else {
          params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
          actualAlpha <- params$alpha
        }
        # direction of alpha
        #alphaG <- graph_from_adjacency_matrix(actualAlpha < 0, mode = "max")
        alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
        # theta = true_interactions
        se.roc <- huge::huge.roc(se$est$path, theta = alphaG, verbose = FALSE)
        TPRs[, run] <- se.roc$tp
        FPRs[, run] <- se.roc$fp
      }
    } 
    
    # average over ordered lambdas -- NOT SAME LAMBDAS
    this_avg_tpr <- apply(TPRs, MARGIN = 1, FUN = function(x) {mean(x, na.rm = TRUE)})
    this_avg_fpr <- apply(FPRs, MARGIN = 1, FUN = function(x) {mean(x, na.rm = TRUE)})
    this_tpr_sd <- sqrt(apply(TPRs, MARGIN = 1, FUN = var) / apply(TPRs, MARGIN = 1, FUN = length))
    this_fpr_sd <- sqrt(apply(FPRs, MARGIN = 1, FUN = var) / apply(FPRs, MARGIN = 1, FUN = length))
    avg_TPR <- c(avg_TPR, this_avg_tpr)
    avg_FPR <- c(avg_FPR, this_avg_fpr)
    TPR_sd <- c(TPR_sd, this_tpr_sd) 
    FPR_sd <- c(FPR_sd, this_fpr_sd)
    avg_precision <- c(avg_precision, rep(NA, times = length(se.roc$fp)))
    avg_recall <- c(avg_recall, rep(NA, times = length(se.roc$fp)))
    methods <- c(methods, rep(paste0("SpiecEasi_mb"), times = length(se.roc$fp)))

    modelSelect <- c(modelSelect, rep(FALSE, times = length(se.roc$fp)))
    numSamplesL <- c(numSamplesL, rep(NA, times = length(se.roc$fp)))
    numSims <- c(numSims, rep(NA, times = length(se.roc$fp)))
    numSp <- c(numSp, rep(NA, times = length(se.roc$fp)))

    file <- c(file, rep(NA, times = length(se.roc$fp)))
    thresholds <- c(thresholds, rep(NA, times = length(se.roc$fp)))
  }

  if(file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", paste0(seName2, "_infResGathered_", numRuns, "sims.csv")))) {

    # for spiecEasi, get average for glasso across a lot of lambdas
    ######################################
    TPRs <- matrix(nrow = 100, ncol = numRuns) # rows are lambda values, columns are runs
    FPRs <- matrix(nrow = 100, ncol = numRuns)
    for (run in 1:numRuns) {
      subdir <- paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/", seName2, "/trial1/")
      if (file.exists(paste0(subdir, "se_rawRes.Rdata"))) {
        se <- readRDS(paste0(subdir, "se_rawRes.Rdata"))
        if (filtered) {
          params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/paramsFiltered", numSamples, ".Rdata"))
          actualAlpha <- params$filteredAlpha
        } else {
          params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
          actualAlpha <- params$alpha
        }
        # direction
        #alphaG <- graph_from_adjacency_matrix(actualAlpha < 0, mode = "max")
        alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
        # theta = true_interactions
        se.roc <- huge::huge.roc(se$est$path, theta = alphaG, verbose = FALSE)
        TPRs[, run] <- se.roc$tp
        FPRs[, run] <- se.roc$fp
      }
    } 
    # average over ordered lambdas -- NOT SAME LAMBDAS
    this_avg_tpr <- apply(TPRs, MARGIN = 1, FUN = function(x) {mean(x, na.rm = TRUE)})
    this_avg_fpr <- apply(FPRs, MARGIN = 1, FUN = function(x) {mean(x, na.rm = TRUE)})
    this_tpr_sd <- sqrt(apply(TPRs, MARGIN = 1, FUN = var) / apply(TPRs, MARGIN = 1, FUN = length))
    this_fpr_sd <- sqrt(apply(FPRs, MARGIN = 1, FUN = var) / apply(FPRs, MARGIN = 1, FUN = length))
    avg_TPR <- c(avg_TPR, this_avg_tpr)
    avg_FPR <- c(avg_FPR, this_avg_fpr)
    TPR_sd <- c(TPR_sd, this_tpr_sd)
    FPR_sd <- c(FPR_sd, this_fpr_sd)
    avg_precision <- c(avg_precision, rep(NA, times = length(se.roc$fp)))
    avg_recall <- c(avg_recall, rep(NA, times = length(se.roc$fp)))
    methods <- c(methods, rep(paste0("SpiecEasi_glasso"), times = length(se.roc$fp)))

    modelSelect <- c(modelSelect, rep(FALSE, times = length(se.roc$fp)))
    numSamplesL <- c(numSamplesL, rep(NA, times = length(se.roc$fp)))
    numSims <- c(numSims, rep(NA, times = length(se.roc$fp)))
    numSp <- c(numSp, rep(NA, times = length(se.roc$fp)))

    file <- c(file, rep(NA, times = length(se.roc$fp)))
    thresholds <- c(thresholds, rep(NA, times = length(se.roc$fp)))
  }
}
}

# this is not working yet for ec
ecAdjust <- FALSE
if (ecAdjust) {
  TPRs <- matrix(nrow = 100, ncol = numRuns) # rows are lambda values, columns are runs
  FPRs <- matrix(nrow = 100, ncol = numRuns)
  for (run in 1:numRuns) {
    subdir <- paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/ecoCopula_res_noCov/trial1/")
    ec <- readRDS(paste0(subdir, "EC_res.Rdata"))
    params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
    actualAlpha <- params$alpha
    alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
    # TODO
    #ec.roc <- 
    #TPRs[, run] <- 
    #FPRs[, run] <- 
  } 
  this_avg_tpr <- apply(TPRs, MARGIN = 1, FUN = mean)
  this_avg_fpr <- apply(FPRs, MARGIN = 1, FUN = mean)
  avg_TPR <- c(avg_TPR, this_avg_tpr)
  avg_FPR <- c(avg_FPR, this_avg_fpr)
  avg_precision <- c(avg_precision, rep(NA, times = dim(TPRs)[1]))
  avg_recall <- c(avg_recall, rep(NA, times = dim(TPRs)[1]))
  methods <- c(methods, rep(paste0("SpiecEasi_avg"), times = dim(TPRs)[1]))
  file <- c(file, rep(NA, times = dim(TPRs)[1]))
  thresholds <- c(thresholds, rep(NA, times = dim(TPRs)[1]))
}

ROC_data <- data.frame(file = file, method = as.factor(methods), threshold = thresholds, modelSelect = modelSelect,
                      numSamples = numSamplesL, numSims = numSims, numSp = numSp, 
                      avg_FPR = as.numeric(avg_FPR), avg_TPR = as.numeric(avg_TPR), 
                      FPR_sd = as.numeric(FPR_sd), TPR_sd = as.numeric(TPR_sd),
                      avg_precision = as.numeric(avg_precision), avg_recall = as.numeric(avg_recall))

####################
#thisDir <-  paste0("/space/s1/fiona_callahan/multiSim_50sp/")
#ROC_data <- read.csv(paste0(thisDir, "ROC_data_cluster.csv"))
#cluster <- FALSE
#if (cluster) {
#  ROC_data <- read.csv(paste0(thisDir, "ROC_data_cluster.csv"))
#} else {
#  ROC_data <- read.csv(paste0(thisDir, "ROC_data.csv"))
#}
####################

ROC_plot <- ggplot(ROC_data, aes(x = avg_FPR, y = avg_TPR, color = method, group = method)) +
  scale_shape_manual(values = 1:12) +
  geom_point(size = 5, aes(shape = modelSelect)) +
  stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
  labs(x = "FPR", y = "TPR") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +   # Add y=x line (no skill)
  #geom_errorbar(aes(ymin = pmax(0, avg_TPR - TPR_sd), ymax = pmin(1, avg_TPR + TPR_sd)), width = 0.03) +  # Add TPR error bars
  #geom_errorbarh(aes(xmin = pmax(0, avg_FPR - FPR_sd), xmax = pmin(1, avg_FPR + FPR_sd)), height = 0.03) +  # Add FPR error bars
  lims(x = c(0, 1), y = c(0, 1)) + # Set x and y-axis limits
  scale_color_manual(values = c("SpiecEasi_glasso" = "red", 
                                "SpiecEasi_mb" = "deeppink3", 
                                "sparcc" = "darkred",
                                "INLA_noCov" = "blue", 
                                "INLA_cov" = "deepskyblue", 
                                "ecoCopula_cov" = "green", 
                                "ecoCopula_noCov" = "darkgreen", 
                                "logistic_noCov" = "yellow",
                                "logistic_cov" = "brown1",
                                "logistic" = "darkorange",
                                "JAGS" = "darkorchid")) + # Manual color scale for method
  theme(text = element_text(size = 24))  # Set the base size for all text elements

#ROC_data[ROC_data$method == "SpiecEasi_glasso",]

PR_plot <- ggplot(ROC_data, aes(x = avg_recall, y = avg_precision, color = method, shape = method)) +
  scale_shape_manual(values = 1:12) +
  geom_point(size = 3) +
  stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
  labs(title = paste("P-R: cluster = ", cluster), x = "Recall = TP/(TP+FN)", y = "Precision = TP/(TP+FP)") +
  geom_abline(intercept = 0.5, slope = 0, linetype = "dashed", color = "gray") +   # Add y=0.5 line (no skill)
  lims(x = c(0, 1), y = c(0, 1)) + # Set x and y-axis limits
  theme(text = element_text(size = 24))  # Set the base size for all text elements

if (saveRes) {
  if (cluster) {
    if (ratio_of_avg) {
      saveName <- paste0("_cluster_ratioOfAv_", numSamples, "samples_", numRuns, "runs")
    } else {
      saveName <- paste0("_cluster_", numSamples, "samples_", numRuns, "runs")
    }
  } else {
    if (ratio_of_avg) {
      saveName <- paste0("_ratioOfAv_", numSamples, "samples_", numRuns, "runs")
    } else {
      saveName <- paste0("_", numSamples, "samples_", numRuns, "runs")
    }
  }

  if (logi) {
    saveName <- paste0(saveName, "_logi")
  }


  ggsave(ROC_plot, filename = paste0(thisDir, "ROC_plot_", saveName, ".png"))
  write.csv(ROC_data, file = paste0(thisDir, "ROC_data_", saveName, ".csv"))
  ggsave(PR_plot, filename = paste0(thisDir, "PR_plot_", saveName, ".png"))
}

ROC_plot