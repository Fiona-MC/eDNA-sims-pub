library(ggplot2)
library(tidyr)
library(stringr)
library(igraph)

cluster <- FALSE
numRuns <- 100
numSamples <- 10000
logi <- TRUE
saveRes <- FALSE

dirNames <- c("multiSim_10sp", "multiSim_50sp", "multiSim_100sp", "multiSim_10sp_random_moreSamples", "multiSim_100sp_random_moreSamples")

# if you need to check these later, look at obsidian ROC curve note
# obsidian://open?vault=simulations&file=ROC%20curves
get_TPR <- function(data, mode = "ignore_sign", return_components = FALSE) {
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

get_FPR <- function(data, mode = "ignore_sign", return_components = FALSE) {
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

get_precision <- function(data, mode = "ignore_sign", return_components = FALSE) {
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

get_recall <- function(data, mode = "ignore_sign", return_components = FALSE) {
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

get_falseDiscovery <- function(data, mode = "ignore_sign", return_components = FALSE) {
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

#dirName <- c("multiSim_test2x10sp")
#resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
#resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
#              "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")

#ls /space/s1/fiona_callahan/multiSim_100
logistic_cutoffs <- c(0.05)
#log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _cov_2runs_cutoff", x, ".csv")})
#log_resnames_noCov <- 
#   sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _noCov_2runs_cutoff", x, ".csv")})


log_resnames_cov_logi <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_covNoCount_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
log_resnames_noCov_logi <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_noCov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})

log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_covNoCount_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_noCov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})

#if (numRuns == 100 && numSamples == 100 && dirName == c("multiSim_10sp")) {
#  log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
#                            paste0("logistic_mistakes_sampled", numSamples, "_cov_", numRuns, "runs_cutoff", x, ".csv")})
#  log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
#                            paste0("logistic_mistakes_sampled", numSamples, "_noCov_", numRuns, "runs_cutoff", x, ".csv")})
#}

#inla_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_infResGathered_cutoff", x, ".csv")})

#inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paperSep_infResGathered_cutoff", x, ".csv")})

#inla_cutoffs <- c(0.05)
#if (logi) {
#  inla_resnames_cov <- sapply(X = inla_cutoffs, FUN = function(x) {
#                        paste0("INLA_res_paperSep_sampled", numSamples, "_cov_logi_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
#  inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
#                        paste0("INLA_res_paperSep_sampled", numSamples, "_noCov_logi_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
#} else {
#  inla_resnames_cov <- sapply(X = inla_cutoffs, FUN = function(x) {
#                        paste0("INLA_res_paperSep_sampled", numSamples, "_cov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
#  inla_resnames_noCov <- sapply(X = inla_cutoffs, FUN = function(x) {
#                        paste0("INLA_res_paperSep_sampled", numSamples, "_noCov_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})
#}

#inla_cutoffs <- c(0, 0.0000001, .3, .5, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_sampled500_infResGathered_cutoff", x, ".csv")})

#inla_cutoffs <- c(0, 0.0000001, .3, .5, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
#inla_resnames_500 <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_sampled500_infResGathered_cutoff", x, ".csv")})

seName1_logi <- paste0("spiecEasi_res_sampled", numSamples, "_mb_logi")
seName2_logi <- paste0("spiecEasi_res_sampled", numSamples, "_glasso_logi")
seName3_logi <- paste0("spiecEasi_res_sampled", numSamples, "_sparcc_logi")

seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb_filtered")
seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso_filtered")
seName3 <- paste0("spiecEasi_res_sampled", numSamples, "_sparcc_filtered")



ecName1_logi <- paste0("ecoCopula_res_sampled", numSamples, "_noCov_logi")
ecName2_logi <- paste0("ecoCopula_res_sampled", numSamples, "_covNoCount_logi")

ecName1 <- paste0("ecoCopula_res_sampled", numSamples, "_noCov_filtered")
ecName2 <- paste0("ecoCopula_res_sampled", numSamples, "_covNoCount_filtered")



inla1_logi <- paste0("INLA_res_paperSep_sampled", numSamples, "_logi_infResGathered_", numRuns, "sims.csv")
jags1_logi <- paste0("JAGS_infResGathered_sampled", numSamples, "_logi_", numRuns, "sims.csv")

inla1 <- paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_infResGathered_", numRuns, "sims.csv")
jags1 <- paste0("JAGS_infResGathered_sampled", numSamples, "_filtered_", numRuns, "sims.csv")

if(logi) {
resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName3, "_infResGathered_", numRuns, "sims.csv"),
              inla1, 
              #inla_resnames_noCov,
              #inla_resnames_cov,
              log_resnames_cov,
              log_resnames_noCov,
              jags1,
              paste0(ecName1_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName3_logi, "_infResGathered_", numRuns, "sims.csv"),
              inla1_logi, 
              #inla_resnames_noCov,
              #inla_resnames_cov,
              log_resnames_cov_logi,
              log_resnames_noCov_logi,
              jags1_logi)
} else {
  resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName3, "_infResGathered_", numRuns, "sims.csv"),
              inla1, 
              #inla_resnames_noCov,
              #inla_resnames_cov,
              log_resnames_cov,
              log_resnames_noCov,
              jags1)
}



########################################
fdrL <- c()
fprL <- c()
tpL <- c()
fpL <- c()
tnL <- c()

fileNameL <- c()
simL <- c()

for (dirName in dirNames) { #iterate through sims
  # check if runs exist
  thisDir <-  paste0("/space/s1/fiona_callahan/", dirName, "/")
  for (i in seq_along(resNames)) {
    if(!file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))) {
      print("This file does not exist:")
      print(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))
    }
  }

  for (i in seq_along(resNames)) {
    if(file.exists(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))) {
      thisMultiSimRes <- read.csv(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]), header = TRUE)
      if (str_detect(resNames[i], "JAGS")) {
        cutoffIndex <- 9 # corresponds to 0.05
        thisMultiSimRes <- thisMultiSimRes[thisMultiSimRes$cutoff_val == cutoffIndex, ]
      } 
      thisFDR <- get_falseDiscovery(data = thisMultiSimRes, mode = "ignore_sign", return_components = TRUE)
      thisTPR <- get_FPR(data = thisMultiSimRes, mode = "ignore_sign", return_components = TRUE)
      tpL <- c(tpL, mean(thisFDR$TP, na.rm = TRUE))
      fpL <- c(fpL, mean(thisFDR$FP, na.rm = TRUE))
      tnL <- c(tnL, mean(thisTPR$TN, na.rm = TRUE))
      fdrL <- c(fdrL, mean(thisFDR$FDR, na.rm = TRUE))
      fprL <- c(fprL, mean(thisTPR$FPR, na.rm = TRUE))
      fileNameL <- c(fileNameL, resNames[i])
      
      if (str_detect(resNames[i], "_logi")) {
        thisSim <- paste0(dirName, "_logi")
      } else {
        thisSim <- dirName
      }
      simL <- c(simL, thisSim)
    }
  }
}

methodL <- c()
for (fileName in fileNameL) {
  baseName <- str_split(fileName, "_")[[1]][1]
  thisMethod <- baseName
  if (str_detect(fileName, "mb")) {
    thisMethod <- paste0(baseName, "_mb")
  } else if (str_detect(fileName, "glasso")) {
    thisMethod <- paste0(baseName, "_glasso")
  } else if (str_detect(fileName, "sparcc")) {
    thisMethod <- paste0(baseName, "_sparcc")
  }
  if (str_detect(fileName, "cov")) {
    thisMethod <- paste0(baseName, "_cov")
  } else if (str_detect(fileName, "noCov")) {
    thisMethod <- paste0(baseName, "_noCov")
  }
  methodL <- c(methodL, thisMethod)
}

simL[simL == "multiSim_10sp_random_moreSamples"] <- "multiSim_10sp_random"
simL[simL == "multiSim_100sp_random_moreSamples"] <- "multiSim_100sp_random"

FDRs <- data.frame(File = fileNameL, Method = methodL, Simulation = simL, FDR = fdrL)
FPRs <- data.frame(File = fileNameL, Method = methodL, Simulation = simL, FPR = fprL)

#FDRs <- data.frame(Method = fileNameL, Simulation = simL, FDR = fpL / (fpL + tpL))
#FPRs <- data.frame(Method = fileNameL, Simulation = simL, FPR = fpL / (fpL + tnL))

head(FPRs)

FDRs_filtered <- FDRs[!is.na(FDRs$FDR) & !is.nan(FDRs$FDR), ]
FPRs_filtered <- FPRs[!is.na(FPRs$FPR) & !is.nan(FPRs$FPR), ]

# Create the heatmap
ggplot(FDRs, aes(x = Simulation, y = Method, fill = FDR)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", na.value = "gray") + # You can adjust colors as needed
  geom_text(aes(label = sprintf("%.2f", FDR)), color = "black", size = 5) + # Display FDR values as text
  theme_minimal() + # You can change the theme as per your preference
  labs(x = "Simulation",
       y = "Method",
       fill = "FDR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Adjust the size value as needed
        axis.text.y = element_text(size = 14),  # Adjust the size value as needed
        axis.title = element_text(size = 14))   # Adjust the size value as needed



#ggplot(FPRs_filtered, aes(x = Simulation, y = Method, fill = FPR)) +
#  geom_tile() +
#  scale_fill_gradient(low = "gray", high = "blue") + # You can adjust colors as needed
#  geom_text(aes(label = sprintf("%.2f", FPR)), color = "black") + # Display FDR values as text
#  theme_minimal() + # You can change the theme as per your preference
#  labs(x = "Simulation",
#       y = "Method",
#       fill = "FPR")
