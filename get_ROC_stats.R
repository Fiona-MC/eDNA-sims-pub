library(ggplot2)
library(tidyr)
library(stringr)

dirName <- c("multiSim_100")
multiSimRes <- data.frame()
resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
              "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")

#ls /space/s1/fiona_callahan/multiSim_100
logistic_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
log_resnames <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_cutoff", x, ".csv")})
resNames <- c("ecoCopula_res_noCov_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
            "logistic_mistakes.csv", "INLA_infResGathered.csv", log_resnames)

#logistic_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
#log_resnames <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_dumb_cutoff", x, ".csv")})
resNames <- log_resnames

thisDir <-  paste0("/space/s1/fiona_callahan/", dirName, "/")
# load results into list
multiSimResL <- c()
for (i in seq_along(resNames)) {
  thisMultiSimRes <- read.csv(paste0("/space/s1/fiona_callahan/", dirName, "/", resNames[i]), header = TRUE)
  thisMultiSimRes$totalMistakes <- thisMultiSimRes$num_incorrectInferences + thisMultiSimRes$num_missedEffectsL
  thisMultiSimRes$fp_fp_tp <- thisMultiSimRes$num_incorrectInferences / 
                              (thisMultiSimRes$num_correctInferences + thisMultiSimRes$num_incorrectInferences)
  #thisMultiSimRes$fp_fp_tp_cluster <- thisMultiSimRes$num_incorrect_cluster / 
  #                            (thisMultiSimRes$num_correct_cluster + thisMultiSimRes$num_incorrect_cluster)
  actual_mean_fpr <- rep(NA, times = dim(thisMultiSimRes)[1])
  for (row in seq_len(dim(thisMultiSimRes)[1])) { #1:dim(multiSimRes)[1]
    if(is.na(thisMultiSimRes$fpr.mode[row])) {
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
  thisMultiSimRes$actual_mean_fpr <- actual_mean_fpr

  multiSimResL[[resNames[i]]] <- thisMultiSimRes
}

multiSimRes <- multiSimResL$logistic_mistakes.csv

# if you need to check these later, look at obsidian ROC curve note
# obsidian://open?vault=simulations&file=ROC%20curves
get_TPR <- function(data = multiSimRes, mode = "ignore_sign") {
    if (mode == "cluster") {
        TP <- data$num_correct_cluster 
        FN <- data$num_missed_cluster
    } else if (mode == "ignore_sign") {
        TP <- data$TP_ignoreSign
        FN <- data$FN_ignoreSign
    } else {
        # tp/(tp+fn)
        TP <- data$num_correctInferences
        FN <- data$num_missedEffectsL
    } 
    TPR <- TP / (TP + FN)
    return(TPR)
}

get_FPR <- function(data = multiSimRes, mode = "ignore_sign") {
    if (mode == "cluster") {
      FP <- data$num_incorrect_cluster
      actual_neg <- data$num_possibleEffectsL - (data$num_correct_cluster + data$num_missed_cluster)
      TN <- actual_neg - data$num_incorrect_cluster
    } else if (mode == "ignore_sign") {
      FP <- data$FP_ignoreSign
      TN <- data$TN_ignoreSign
    } else {
        FP <- data$num_incorrectInferences
        actual_neg <- data$num_possibleEffectsL - data$num_actualEffects
        TN <- actual_neg - data$num_incorrectInferences
    } 
    FPR <- FP / (FP + TN)
    return(FPR)
}

# method | threshold | avg_TPR | avg_FPR 
file <- rep(NA, times = length(multiSimResL))
methods <- rep(NA, times = length(multiSimResL))
thresholds <- rep(NA, times = length(multiSimResL))
avg_TPR <- rep(NA, times = length(multiSimResL))
avg_FPR <- rep(NA, times = length(multiSimResL))

for (i in seq_along(multiSimResL)) {
    multiSimRes <- multiSimResL[[i]]
    file[i] <- names(multiSimResL)[i]
    if (str_detect(names(multiSimResL)[i], "logistic")) {
      methods[i] <- "logistic"
    } else if (str_detect(names(multiSimResL)[i], "INLA")) {
      methods[i] <- "INLA"
    } else if (str_detect(names(multiSimResL)[i], "ecoCopula")) {
      methods[i] <- "ecoCopula"
    } else if (str_detect(names(multiSimResL)[i], "spiecEasi")) {
      methods[i] <- "spiecEasi"
    } 
    #if (str_detect(names(multiSimResL)[i], "cutoff")) {
    #  thresholds[i] <- str_split(c(names(multiSimResL)[i]), pattern = ".")
    #}
    avg_TPR[i] <- mean(get_TPR(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
    avg_FPR[i] <- mean(get_FPR(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
}

ROC_data <- data.frame(file = file, method = methods, threshold = thresholds, avg_FPR = as.numeric(avg_FPR), avg_TPR = as.numeric(avg_TPR))

ROC_data$avg_FPR <- as.numeric(ROC_data$avg_FPR)
ROC_data$avg_TPR <- as.numeric(ROC_data$avg_TPR)

ROC_plot <- ggplot(ROC_data, aes(x = avg_FPR, y = avg_TPR, color = method)) +
  geom_point(size = 3) +
  stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
  labs(title = "ROC", x = "FPR", y = "TPR") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +   # Add y=x line
  lims(x = c(0, 1), y = c(0, 1))  # Set x and y-axis limits


ggsave(ROC_plot, filename = paste0(thisDir, "ROC_plot.png"))
