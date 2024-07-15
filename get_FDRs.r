library(ggplot2)
library(tidyr)
library(stringr)
library(igraph)

setwd("/home/fiona_callahan/eDNA_sims_code")

numRuns <- 100
numSamples <- 10000
logi <- TRUE # include logi
ratio_of_avg <- TRUE
bonferroni <- TRUE
mode <- "cluster_cov" # "ignore_sign" "cluster" "ignore_direction" "cluster_cov"

dirNames <- c("multiSim_10sp", "multiSim_100sp", "multiSim_10sp_random_moreSamples", "multiSim_100sp_random_moreSamples")

source("./confusion_stats.R")

#dirName <- c("multiSim_test2x10sp")
#resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
#resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
#              "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")

#ls /space/s1/fiona_callahan/multiSim_100
if (bonferroni) {
  logistic_cutoffs <- c("bonferroni")
} else {
  logistic_cutoffs <- c(0.05)
}
#log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _cov_2runs_cutoff", x, ".csv")})
#log_resnames_noCov <- 
#   sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _noCov_2runs_cutoff", x, ".csv")})


log_resnames_cov_logi <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_cov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
log_resnames_noCov_logi <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_noCov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})

log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_cov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
log_resnames_noCov <- sapply(X = logistic_cutoffs, FUN = function(x) {
                            paste0("logistic_mistakes_sampled", numSamples, "_noCov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})


if (bonferroni) {
  linear_cutoffs <- c("bonferroni")
} else {
  linear_cutoffs <- c(0.05)
}
#log_resnames_cov <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _cov_2runs_cutoff", x, ".csv")})
#log_resnames_noCov <- 
#   sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_sampled", numSamples, _noCov_2runs_cutoff", x, ".csv")})


lin_resnames_cov_logi <- sapply(X = linear_cutoffs, FUN = function(x) {
                            paste0("linearReg_mistakes_sampled", numSamples, "_cov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
lin_resnames_noCov_logi <- sapply(X = linear_cutoffs, FUN = function(x) {
                            paste0("linearReg_mistakes_sampled", numSamples, "_noCov_logi_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})

lin_resnames_cov <- sapply(X = linear_cutoffs, FUN = function(x) {
                            paste0("linearReg_mistakes_sampled", numSamples, "_cov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
lin_resnames_noCov <- sapply(X = linear_cutoffs, FUN = function(x) {
                            paste0("linearReg_mistakes_sampled", numSamples, "_noCov_filtered_", numRuns, "runs_cutoff", x, "_", numRuns, "sims.csv")})
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

seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb_filtered")
seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso_filtered")


sparcc_cutoffs <- c("pval_bootstrap")

sparcc_resnames_logi <- sapply(X = sparcc_cutoffs, FUN = function(x) {
                      paste0("sparcc_res_sampled", numSamples, "_logi_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})

sparcc_resnames <- sapply(X = sparcc_cutoffs, FUN = function(x) {
                    paste0("sparcc_res_sampled", numSamples, "_filtered_infResGathered_cutoff", x, "_", numRuns, "sims.csv")})





ecName1_logi <- paste0("ecoCopula_res_sampled", numSamples, "_noCov_logi")
ecName2_logi <- paste0("ecoCopula_res_sampled", numSamples, "_cov_logi")

ecName1 <- paste0("ecoCopula_res_sampled", numSamples, "_noCov_filtered")
ecName2 <- paste0("ecoCopula_res_sampled", numSamples, "_cov_filtered")

ecName3_logi <- paste0("ecoCopula_res_readAbd_sampled", numSamples, "_noCov_logi")
ecName4_logi <- paste0("ecoCopula_res_readAbd_sampled", numSamples, "_cov_logi")

ecName3 <- paste0("ecoCopula_res_readAbd_sampled", numSamples, "_noCov_filtered")
ecName4 <- paste0("ecoCopula_res_readAbd_sampled", numSamples, "_cov_filtered")



inla1_logi <- paste0("INLA_res_paperSep_sampled", numSamples, "_logi_infResGathered_", numRuns, "sims.csv")
jags1_logi <- paste0("JAGS_infResGathered_sampled", numSamples, "_logi_", numRuns, "sims.csv")

inla1 <- paste0("INLA_res_paperSep_sampled", numSamples, "_filtered_infResGathered_", numRuns, "sims.csv")
jags1 <- paste0("JAGS_infResGathered_sampled", numSamples, "_filtered_", numRuns, "sims.csv")

if(logi) {
resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName3, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName4, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              sparcc_resnames,
              inla1, 
              #inla_resnames_noCov,
              #inla_resnames_cov,
              log_resnames_cov,
              log_resnames_noCov,
              lin_resnames_cov,
              lin_resnames_noCov,
              jags1,
              paste0(ecName1_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName3_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName4_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1_logi, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2_logi, "_infResGathered_", numRuns, "sims.csv"), 
              sparcc_resnames_logi,
              inla1_logi, 
              #inla_resnames_noCov,
              #inla_resnames_cov,
              log_resnames_cov_logi,
              log_resnames_noCov_logi,
              lin_resnames_cov_logi,
              lin_resnames_noCov_logi,
              jags1_logi)
} else {
  resNames <- c(paste0(ecName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName2, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName3, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(ecName4, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName1, "_infResGathered_", numRuns, "sims.csv"), 
              paste0(seName2, "_infResGathered_", numRuns, "sims.csv"), 
              sparcc_resnames,
              inla1, 
              #inla_resnames_noCov,
              #inla_resnames_cov,
              log_resnames_cov,
              log_resnames_noCov,
              lin_resnames_cov,
              lin_resnames_noCov,
              jags1)
}



########################################
fdrL <- c()
fprL <- c()
tpL <- c()
fpL <- c()
tnL <- c()
fnL <- c()

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
      #print(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]))
      thisMultiSimRes <- read.csv(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]), header = TRUE)
      if (str_detect(resNames[i], "JAGS")) {
        cutoffIndex <- 9 # corresponds to 0.05
        thisMultiSimRes <- thisMultiSimRes[thisMultiSimRes$cutoff_val == cutoffIndex, ]
      } 


      thisFDR <- get_falseDiscovery(data = thisMultiSimRes, mode = mode, return_components = TRUE)
      thisTPR <- get_TPR(data = thisMultiSimRes, mode = mode, return_components = TRUE)
      thisFPR <- get_FPR(data = thisMultiSimRes, mode = mode, return_components = TRUE)
    

      tpL <- c(tpL, mean(thisFDR$TP, na.rm = TRUE))
      fpL <- c(fpL, mean(thisFDR$FP, na.rm = TRUE))
      tnL <- c(tnL, mean(thisFPR$TN, na.rm = TRUE))
      fnL <- c(fnL, mean(thisTPR$FN, na.rm = TRUE))
      fdrL <- c(fdrL, mean(thisFDR$FDR, na.rm = TRUE))
      fprL <- c(fprL, mean(thisFPR$FPR, na.rm = TRUE))
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
  } 
  if (str_detect(fileName, "cov")) {
    thisMethod <- paste0(baseName, "_cov")
  } else if (str_detect(fileName, "noCov")) {
    thisMethod <- paste0(baseName, "_noCov")
  }
  if (str_detect(fileName, "readAbd")) {
    thisMethod <- paste0(thisMethod, "_readAbd")
  }
  methodL <- c(methodL, thisMethod)
}

simL[simL == "multiSim_10sp_random_moreSamples"] <- "multiSim_10sp_random"
simL[simL == "multiSim_100sp_random_moreSamples"] <- "multiSim_100sp_random"


if (ratio_of_avg) {
  FDRs <- data.frame(File = fileNameL, Method = methodL, Simulation = simL, FDR = fpL / (fpL + tpL), FP = fpL, TP = tpL)
  FPRs <- data.frame(File = fileNameL, Method = methodL, Simulation = simL, FPR = fpL / (fpL + tnL))
} else {
  FDRs <- data.frame(File = fileNameL, Method = methodL, Simulation = simL, FDR = fdrL)
  FPRs <- data.frame(File = fileNameL, Method = methodL, Simulation = simL, FPR = fprL)
}

head(FPRs)

FDRs_filtered <- FDRs[!is.na(FDRs$FDR) & !is.nan(FDRs$FDR), ]
FPRs_filtered <- FPRs[!is.na(FPRs$FPR) & !is.nan(FPRs$FPR), ]

custom_labels_simulation <- c("multiSim_10sp" = "Set parms, 10sp",
                              "multiSim_10sp_logi" = "Covariance Mx, 10sp",
                              "multiSim_100sp" = "Set parms, 100sp",
                              "multiSim_100sp_logi" = "Covariance Mx, 100sp",
                              "multiSim_10sp_random" = "Random parms, 10sp",
                              "multiSim_100sp_random" = "Random parms, 100sp")

custom_labels_method <- c("ecoCopula_noCov" = "ecoCopula (noCov, pres-abs)",
                          "ecoCopula_cov" = "ecoCopula (cov, pres-abs)",
                          "ecoCopula_noCov_readAbd" = "ecoCopula (noCov, readAbd)",
                          "ecoCopula_cov_readAbd" = "ecoCopula (cov, readAbd)",
                          "spiecEasi_mb" = "spiecEasi (mb)",
                          "spiecEasi_glasso" = "spiecEasi (glasso)",      
                          "spiecEasi_sparcc" = "SPARCC",
                          "INLA" = "SDM-INLA",              
                          "logistic_cov" = "logistic (cov)",          
                          "logistic_noCov" = "logistic (noCov)",
                          "linearReg_cov" = "linear (cov)",
                          "linearReg_noCov" = "linear (noCov)",     
                          "JAGS" = "JSDM-MCMC")


if(logistic_cutoffs[1] == 0.05) {
  write.csv(FDRs, file = paste0("/space/s1/fiona_callahan/FDR_stats_", numSamples, "samples_", mode, ".csv"))
} else {
  write.csv(FDRs, file = paste0("/space/s1/fiona_callahan/FDR_stats_", numSamples, "samples_cutoff", logistic_cutoffs[1], "_", mode, ".csv"))
}

# Create the heatmap
FDRplot <- ggplot(FDRs, aes(x = Simulation, y = Method, fill = FDR)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", na.value = "gray", limits = c(0, 1)) + # You can adjust colors as needed
  geom_text(aes(label = sprintf("%.2f", FDR)), color = "black", size = 5) + # Display FDR values as text
  theme_minimal() + # You can change the theme as per your preference
  labs(x = "Simulation",
       y = "Method",
       fill = "FDR") +
  scale_x_discrete(labels = custom_labels_simulation) +  # Custom labels for the x-axis
  scale_y_discrete(labels = custom_labels_method) +      # Custom labels for the y-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Adjust the size value as needed
        axis.text.y = element_text(size = 14),  # Adjust the size value as needed
        axis.title = element_text(size = 14))   # Adjust the size value as needed

FDRplot

if(logistic_cutoffs[1] == 0.05) {
  ggsave(filename = paste0("/space/s1/fiona_callahan/FDR_plot_", numSamples, "samples_", mode, ".pdf"), plot = FDRplot)
} else {
  ggsave(filename = paste0("/space/s1/fiona_callahan/FDR_plot_", numSamples, "samples_cutoff", logistic_cutoffs[1], "_", mode, ".pdf"), plot = FDRplot)
}
#ggplot(FPRs_filtered, aes(x = Simulation, y = Method, fill = FPR)) +
#  geom_tile() +
#  scale_fill_gradient(low = "gray", high = "blue") + # You can adjust colors as needed
#  geom_text(aes(label = sprintf("%.2f", FPR)), color = "black") + # Display FDR values as text
#  theme_minimal() + # You can change the theme as per your preference
#  labs(x = "Simulation",
#       y = "Method",
#       fill = "FPR")

FDRs$totalDiscoveries <- FDRs$FP + FDRs$TP

DRplot <- ggplot(FDRs, aes(x = Simulation, y = Method, fill = totalDiscoveries)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", na.value = "gray") + # You can adjust colors as needed
  geom_text(aes(label = sprintf("%.2f", totalDiscoveries)), color = "black", size = 5) + 
  theme_minimal() + # You can change the theme as per your preference
  labs(x = "Simulation",
       y = "Method",
       fill = "Average Discoveries") +
  scale_x_discrete(labels = custom_labels_simulation) +  # Custom labels for the x-axis
  scale_y_discrete(labels = custom_labels_method) +      # Custom labels for the y-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Adjust the size value as needed
        axis.text.y = element_text(size = 14),  # Adjust the size value as needed
        axis.title = element_text(size = 14))   # Adjust the size value as needed


if(logistic_cutoffs[1] == 0.05) {
  ggsave(filename = paste0("/space/s1/fiona_callahan/DR_plot_", numSamples, "samples_", mode, ".pdf"), plot = FDRplot)
} else {
  ggsave(filename = paste0("/space/s1/fiona_callahan/DR_plot_", numSamples, "samples_cutoff", logistic_cutoffs[1], "_", mode, ".pdf"), plot = FDRplot)
}

if(logistic_cutoffs[1] == 0.05) {
  write.csv(FDRs, file = paste0("/space/s1/fiona_callahan/DR_stats_", numSamples, "samples_", mode, ".csv"))
} else {
  write.csv(FDRs, file = paste0("/space/s1/fiona_callahan/DR_stats_", numSamples, "samples_cutoff", logistic_cutoffs[1], "_", mode, ".csv"))
}
#DRplot