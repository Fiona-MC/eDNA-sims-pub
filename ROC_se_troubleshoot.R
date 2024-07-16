library(ggplot2)
library(tidyr)
library(stringr)
library(SpiecEasi)
library(igraph)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)


setwd("/home/fiona_callahan/eDNA_sims_code")

#./get_ROC_stats.R multiSim_100sp_random_moreSamples 1 100 0

mode <- "ignore_direction" # ignore_sign ignore_direction cluster cluster_cov
ratio_of_avg <- FALSE #do we compute the average of the ratio or ratio of averages
numRuns <- 100
numSamples <- 100
logi <- TRUE
saveRes <- TRUE
covMode <- "noCount" # "all" "noCov" "cov" "covNoCount" "noCount"
dirName <- c("multiSim_10sp")
plot <- TRUE
seMode <- "negOnly" # negOnly both

mode <- args[2]
ratio_of_avg <- FALSE #do we compute the average of the ratio or ratio of averages
numRuns <- 100
numSamples <- args[3]
logi <- as.numeric(args[4]) == 1
saveRes <- TRUE
covMode <- "noCount" # "all" "noCov" "cov" "covNoCount" "noCount"
dirName <- c(args[1])
plot <- FALSE
seMode <- args[5] # "posOnly" # negOnly both

if (logi) {
  filtered <- FALSE
} else {
  filtered <- TRUE
}

#dirName <- c("multiSim_test2x10sp")
multiSimRes <- data.frame()
#resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
#resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
#              "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")
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


resNames <- c(
            paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
            paste0(seName2, "_infResGathered_", numRuns, "sims.csv")
            )


thisDir <-  paste0("/space/s1/fiona_callahan/", dirName, "/")
# load results into list
multiSimResL <- list()
for (i in seq_along(resNames)) {
  if(file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))) {
    #print(resNames[i])
    thisMultiSimRes <- fread(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]), header = TRUE)

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
source("./confusion_stats.R")

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
        if (str_detect(names(multiSimResL)[i], "bonferroni")) {
          modelSelect[i] <- TRUE
        } else {
          modelSelect[i] <- FALSE
        }
      } else {
        methods[i] <- "logistic"
        if (str_detect(names(multiSimResL)[i], "bonferroni")) {
          modelSelect[i] <- TRUE
        } else {
          modelSelect[i] <- FALSE
        }
      }
    } else if (str_detect(names(multiSimResL)[i], "linear")) {
      if (str_detect(names(multiSimResL)[i], "noCov")) {
        methods[i] <- "linear_noCov"
        if (str_detect(names(multiSimResL)[i], "bonferroni")) {
          modelSelect[i] <- TRUE
        } else {
          modelSelect[i] <- FALSE
        }
      } else {
        methods[i] <- "linear_cov"
        if (str_detect(names(multiSimResL)[i], "bonferroni")) {
          modelSelect[i] <- TRUE
        } else {
          modelSelect[i] <- FALSE
        }
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
      if(str_detect(names(multiSimResL)[i], "readAbd")) {
        methods[i] <- paste0(methods[i], "_readAbd")
      }
      if (str_detect(names(multiSimResL)[i], "cutoff")) {
        modelSelect[i] <- FALSE
      } else {
        modelSelect[i] <- TRUE
      }
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
      if (str_detect(names(multiSimResL)[i], "_cutoff0.3_")) {
        modelSelect[i] <- TRUE
      } else {
        modelSelect[i] <- FALSE
      }
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

    if(ratio_of_avg) {
      # TPR <- TP / (TP + FN)
      avg_TPR[i] <- mean(get_TPR(data = multiSimRes, mode = mode, return_components = TRUE)$TP, na.rm = TRUE) / 
                    (mean(get_TPR(data = multiSimRes, mode = mode, return_components = TRUE)$TP, na.rm = TRUE) + 
                          mean(get_TPR(data = multiSimRes, mode = mode, return_components = TRUE)$FN, na.rm = TRUE))
      # FPR <- FP / (FP + TN)
      avg_FPR[i] <- mean(get_FPR(data = multiSimRes, mode = mode, return_components = TRUE)$FP, na.rm = TRUE) / 
                    (mean(get_FPR(data = multiSimRes, mode = mode, return_components = TRUE)$FP, na.rm = TRUE) + 
                          mean(get_FPR(data = multiSimRes, mode = mode, return_components = TRUE)$TN, na.rm = TRUE))
      # prec <- TP / (TP + FP)
      avg_precision[i] <-  mean(get_precision(data = multiSimRes, mode = mode, return_components = TRUE)$TP, na.rm = TRUE) / 
                    (mean(get_precision(data = multiSimRes, mode = mode, return_components = TRUE)$TP, na.rm = TRUE) + 
                          mean(get_precision(data = multiSimRes, mode = mode, return_components = TRUE)$FP, na.rm = TRUE))
      # recall <- TP / (TP + FN)
      avg_recall[i] <-  mean(get_recall(data = multiSimRes, mode = mode, return_components = TRUE)$TP, na.rm = TRUE) / 
                    (mean(get_recall(data = multiSimRes, mode = mode, return_components = TRUE)$TP, na.rm = TRUE) + 
                          mean(get_recall(data = multiSimRes, mode = mode, return_components = TRUE)$FN, na.rm = TRUE)) 
    } else {
      TPR_temp <- get_TPR(data = multiSimRes, mode = mode)
      avg_TPR[i] <- mean(TPR_temp, na.rm = TRUE)
      TPR_sd[i] <- sqrt(var(TPR_temp, na.rm = TRUE) / length(TPR_temp))
      FPR_temp <- get_FPR(data = multiSimRes, mode = mode)
      avg_FPR[i] <- mean(FPR_temp, na.rm = TRUE)
      FPR_sd[i] <- sqrt(var(FPR_temp, na.rm = TRUE) / length(FPR_temp))
      avg_precision[i] <- mean(get_precision(data = multiSimRes, mode = mode), na.rm = TRUE)
      avg_recall[i] <- mean(get_recall(data = multiSimRes, mode = mode), na.rm = TRUE)
    }
}

se_include <- TRUE
if(se_include) {

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
          actualBeta <- params$filteredBeta
        } else {
          params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
          actualAlpha <- params$alpha
          actualBeta <- params$beta
        }
        # direction of alpha
        if (seMode == "posOnly") {
            alphaG <- graph_from_adjacency_matrix(actualAlpha > 0, mode = "max")
        } else if (seMode == "negOnly") {
            alphaG <- graph_from_adjacency_matrix(actualAlpha < 0, mode = "max")
        } else {
            #alphaG <- graph_from_adjacency_matrix(actualAlpha < 0, mode = "max")
            alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
        }

        numSpecies <- dim(actualAlpha)[1]

        if (mode == "cluster") {
          connected_alpha_actual <- (distances(alphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
          alphaG <- graph_from_adjacency_matrix(connected_alpha_actual != 0, mode = "max")

          se_inferred_alpha <- lapply(se$est$path, FUN = function(inferredAlpha) {
            inferredAlphaG <- graph_from_adjacency_matrix(as.matrix((inferredAlpha != 0)), mode = "max")
            connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(inferredAlpha)[1], ncol = dim(inferredAlpha)[1]) == 0)
            return(connected_alpha_inferred)
          })
        } else if (mode == "cluster_cov") {
          connectedCov_alpha_actual <- (distances(alphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
          
          for (cov in 1:(dim(actualBeta)[2] - 1)) { # nolint
              for (sp1 in 1:numSpecies) {
                  for (sp2 in 1:numSpecies) {
                      if (sp1 != sp2) {
                          if (actualBeta[sp1, cov] != 0 && actualBeta[sp2, cov] != 0) {
                              connectedCov_alpha_actual[sp1, sp2] <- 1
                              connectedCov_alpha_actual[sp2, sp1] <- 1
                          } 
                      }
                  }
              }
          }
          
          alphaG <- graph_from_adjacency_matrix(connectedCov_alpha_actual != 0, mode = "max")

          # inferred alpha is the same as above in cluster mode
          se_inferred_alpha <- lapply(se$est$path, FUN = function(inferredAlpha) {
            inferredAlphaG <- graph_from_adjacency_matrix(as.matrix((inferredAlpha != 0)), mode = "max")
            connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(inferredAlpha)[1], ncol = dim(inferredAlpha)[1]) == 0)
            return(connected_alpha_inferred)
          })
        } else {
          se_inferred_alpha <- se$est$path
        }

        # theta = true_interactions
        se.roc <- huge::huge.roc(se_inferred_alpha, theta = alphaG, verbose = FALSE)
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
          actualBeta <- params$filteredBeta
        } else {
          params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
          actualAlpha <- params$alpha
          actualBeta <- params$beta
        }
        # direction
        if (seMode == "posOnly") {
            alphaG <- graph_from_adjacency_matrix(actualAlpha > 0, mode = "max")
        } else if (seMode == "negOnly") {
            alphaG <- graph_from_adjacency_matrix(actualAlpha < 0, mode = "max")
        } else {
            #alphaG <- graph_from_adjacency_matrix(actualAlpha < 0, mode = "max")
            alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "max")
        }

        numSpecies <- dim(actualAlpha)[1]

        if (mode == "cluster") {
          connected_alpha_actual <- (distances(alphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
          alphaG <- graph_from_adjacency_matrix(connected_alpha_actual != 0, mode = "max")

          se_inferred_alpha <- lapply(se$est$path, FUN = function(inferredAlpha) {
            inferredAlphaG <- graph_from_adjacency_matrix(as.matrix(inferredAlpha != 0), mode = "max")
            connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(inferredAlpha)[1], ncol = dim(inferredAlpha)[1]) == 0)
            return(connected_alpha_inferred)
          })
        } else if (mode == "cluster_cov") {
          connectedCov_alpha_actual <- (distances(alphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
          
          for (cov in 1:(dim(actualBeta)[2] - 1)) { # nolint
              for (sp1 in 1:numSpecies) {
                  for (sp2 in 1:numSpecies) {
                      if (sp1 != sp2) {
                          if (actualBeta[sp1, cov] != 0 && actualBeta[sp2, cov] != 0) {
                              connectedCov_alpha_actual[sp1, sp2] <- 1
                              connectedCov_alpha_actual[sp2, sp1] <- 1
                          } 
                      }
                  }
              }
          }
          
          alphaG <- graph_from_adjacency_matrix(connectedCov_alpha_actual != 0, mode = "max")

          # inferred alpha is the same as above in cluster mode
          se_inferred_alpha <- lapply(se$est$path, FUN = function(inferredAlpha) {
            inferredAlphaG <- graph_from_adjacency_matrix(as.matrix((inferredAlpha != 0)), mode = "max")
            connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:numSpecies, to = 1:numSpecies) != Inf) * 
                                                (diag(nrow = dim(inferredAlpha)[1], ncol = dim(inferredAlpha)[1]) == 0)
            return(connected_alpha_inferred)
          })
        } else {
          se_inferred_alpha <- se$est$path
        }

        # theta = true_interactions
        se.roc <- huge::huge.roc(se_inferred_alpha, theta = alphaG, verbose = FALSE)
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
} # end of if se include

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
if (plot) {
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
                                  "spiecEasi_glasso" = "red", 
                                  "spiecEasi_mb" = "deeppink3",
                                  "sparcc" = "darkred",
                                  "INLA_noCov" = "blue", 
                                  "INLA_cov" = "deepskyblue", 
                                  "ecoCopula_cov" = "green", 
                                  "ecoCopula_noCov" = "#b1ff8d", 
                                  "ecoCopula_cov_readAbd" = "darkgreen", 
                                  "ecoCopula_noCov_readAbd" = "#00a200", 
                                  "logistic_noCov" = "yellow",
                                  "logistic_cov" = "brown1",
                                  "logistic" = "darkorange",
                                  "linear_cov" = "grey",
                                  "linear_noCov" = "black",
                                  "JAGS" = "darkorchid")) + # Manual color scale for method
    theme(text = element_text(size = 24))  # Set the base size for all text elements
    ROC_plot

}
#ROC_data[ROC_data$method == "SpiecEasi_glasso",]

#PR_plot <- ggplot(ROC_data, aes(x = avg_recall, y = avg_precision, color = method, shape = method)) +
#  scale_shape_manual(values = 1:12) +
#  geom_point(size = 3) +
#  stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
#  labs(title = paste("P-R: mode = ", mode), x = "Recall = TP/(TP+FN)", y = "Precision = TP/(TP+FP)") +
#  geom_abline(intercept = 0.5, slope = 0, linetype = "dashed", color = "gray") +   # Add y=0.5 line (no skill)
#  lims(x = c(0, 1), y = c(0, 1)) + # Set x and y-axis limits
#  theme(text = element_text(size = 24))  # Set the base size for all text elements

if (saveRes) {
  if (ratio_of_avg) {
    saveName <- paste0("_", mode, "_ratioOfAv_", numSamples, "samples_", numRuns, "runs")
  } else {
    saveName <- paste0("_", mode, "_", numSamples, "samples_", numRuns, "runs")
  }

  if (logi) {
    saveName <- paste0(saveName, "_logi")
  }

    if (seMode == "posOnly") {
        saveName <- paste0(saveName, "_sePosOnly")
    } else if (seMode == "negOnly") {
        saveName <- paste0(saveName, "_seNegOnly")
    } else {
        saveName <- paste0(saveName, "_seBoth")
    }

  write.csv(ROC_data, file = paste0(thisDir, "ROC_data_", saveName, ".csv"))

  if (plot) {
    ggsave(ROC_plot, filename = paste0(thisDir, "ROC_plot_", saveName, ".png"))
    #ggsave(PR_plot, filename = paste0(thisDir, "PR_plot_", saveName, ".png"))
  }
}

