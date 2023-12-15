library(ggplot2)
library(tidyr)
library(stringr)

cluster <- FALSE

dirName <- c("multiSim_100")
multiSimRes <- data.frame()
resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", 
              "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")

#ls /space/s1/fiona_callahan/multiSim_100
logistic_cutoffs <- c(0, 1e-128, 1e-64, 1e-32, 1e-16, 1e-8, 1e-4, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15,
             0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
log_resnames <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_cutoff", x, ".csv")})

logistic_cutoffs_noCov <- c(0, 1e-64, 1e-32, 1e-16, 1e-8, 1e-4, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15,
             0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
log_resnames_noCov <- sapply(X = logistic_cutoffs_noCov, FUN = function(x) {paste0("logistic_mistakes_noCov_cutoff", x, ".csv")})

#inla_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_infResGathered_cutoff", x, ".csv")})

inla_cutoffs <- c(0, 0.0000000000001, 0.0000001, 0.00001, .3, .5, .7, .9, 1, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)
inla_resnames <- sapply(X = inla_cutoffs, FUN = function(x) {paste0("INLA_res_paper_infResGathered_cutoff", x, ".csv")})

resNames <- c("ecoCopula_res_noCov_infResGathered.csv", 
              "spiecEasi_res_mb_infResGathered.csv", 
              "INLA_res_paper_infResGathered.csv", 
              inla_resnames,
              log_resnames,
              log_resnames_noCov)

#logistic_cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
#log_resnames <- sapply(X = logistic_cutoffs, FUN = function(x) {paste0("logistic_mistakes_dumb_cutoff", x, ".csv")})
#resNames <- log_resnames

thisDir <-  paste0("/space/s1/fiona_callahan/", dirName, "/")
# load results into list
multiSimResL <- list()
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
        TP <- data$TP_cluster 
        FN <- data$FN_cluster
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
      FP <- data$FP_cluster
      TN <- data$TN_cluster
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

get_precision <- function(data = multiSimRes, mode = "ignore_sign") {
    if (mode == "cluster") {
        TP <- data$TP_cluster 
        FP <- data$FP_cluster
    } else if (mode == "ignore_sign") {
        TP <- data$TP_ignoreSign
        FP <- data$FP_ignoreSign
    } 
    prec <- TP / (TP + FP)
    return(prec)
}

get_recall <- function(data = multiSimRes, mode = "ignore_sign") {
    if (mode == "cluster") {
      TP <- data$TP_cluster
      FN <- data$FN_cluster
    } else if (mode == "ignore_sign") {
      TP <- data$TP_ignoreSign
      FN <- data$FN_ignoreSign
    } 
    recall <- TP / (TP + FN)
    return(recall)
}

# method | threshold | avg_TPR | avg_FPR 
file <- rep(NA, times = length(multiSimResL))
methods <- rep(NA, times = length(multiSimResL))
thresholds <- rep(NA, times = length(multiSimResL))
avg_TPR <- rep(NA, times = length(multiSimResL))
avg_FPR <- rep(NA, times = length(multiSimResL))
avg_precision <- rep(NA, times = length(multiSimResL))
avg_recall <- rep(NA, times = length(multiSimResL))

for (i in seq_along(multiSimResL)) {
    multiSimRes <- multiSimResL[[i]]
    file[i] <- names(multiSimResL)[i]
    if (str_detect(names(multiSimResL)[i], "logistic")) {
      if (str_detect(names(multiSimResL)[i], "noCov")) {
        methods[i] <- "logistic_noCov"
      } else {
        methods[i] <- "logistic"
      }
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
    if (cluster) {
      avg_TPR[i] <- mean(get_TPR(data = multiSimRes, mode = "cluster"), na.rm = TRUE)
      avg_FPR[i] <- mean(get_FPR(data = multiSimRes, mode = "cluster"), na.rm = TRUE)
      avg_precision[i] <- mean(get_precision(data = multiSimRes, mode = "cluster"), na.rm = TRUE)
      avg_recall[i] <- mean(get_recall(data = multiSimRes, mode = "cluster"), na.rm = TRUE)
    } else {
      avg_TPR[i] <- mean(get_TPR(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
      avg_FPR[i] <- mean(get_FPR(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
      avg_precision[i] <- mean(get_precision(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
      avg_recall[i] <- mean(get_recall(data = multiSimRes, mode = "ignore_sign"), na.rm = TRUE)
    }
}

library(SpiecEasi)
library(igraph)
multiples <- FALSE
if(multiples) {
if (!cluster) {
  # for spiec-easi, just put separate lines for a couple of individual sims rather than average
  for (run in sample(1:100, size = 5)) {
    subdir <- paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/spiecEasi_res_mb/trial1/")
    se <- readRDS(paste0(subdir, "se_rawRes.Rdata"))
    params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
    actualAlpha <- params$alpha
    alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "undirected")
    # theta = true_interactions
    se.roc <- huge::huge.roc(se$est$path, theta = alphaG, verbose = FALSE)
    avg_TPR <- c(avg_TPR, se.roc$tp)
    avg_FPR <- c(avg_FPR, se.roc$fp)
    avg_precision <- c(avg_precision, rep(NA, times = length(se.roc$fp)))
    avg_recall <- c(avg_recall, rep(NA, times = length(se.roc$fp)))
    methods <- c(methods, rep(paste0("SpiecEasi", run), times = length(se.roc$fp)))
    file <- c(file, rep(NA, times = length(se.roc$fp)))
    thresholds <- c(thresholds, rep(NA, times = length(se.roc$fp)))
  } 
}
}

if (!cluster) {
  TPRs <- matrix(nrow = 100, ncol = 100) # rows are lambda values, columns are runs
  FPRs <- matrix(nrow = 100, ncol = 100)
  # for spiec-easi, just put separate lines for a couple of individual sims rather than average
  for (run in 1:100) {
    subdir <- paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/spiecEasi_res_mb/trial1/")
    se <- readRDS(paste0(subdir, "se_rawRes.Rdata"))
    params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))
    actualAlpha <- params$alpha
    alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "undirected")
    # theta = true_interactions
    se.roc <- huge::huge.roc(se$est$path, theta = alphaG, verbose = FALSE)
    TPRs[, run] <- se.roc$tp
    FPRs[, run] <- se.roc$fp
  } 
  this_avg_tpr <- apply(TPRs, MARGIN = 1, FUN = mean)
  this_avg_fpr <- apply(FPRs, MARGIN = 1, FUN = mean)
  avg_TPR <- c(avg_TPR, this_avg_tpr)
  avg_FPR <- c(avg_FPR, this_avg_fpr)
  avg_precision <- c(avg_precision, rep(NA, times = length(se.roc$fp)))
  avg_recall <- c(avg_recall, rep(NA, times = length(se.roc$fp)))
  methods <- c(methods, rep(paste0("SpiecEasi_avg"), times = length(se.roc$fp)))
  file <- c(file, rep(NA, times = length(se.roc$fp)))
  thresholds <- c(thresholds, rep(NA, times = length(se.roc$fp)))
}

ROC_data <- data.frame(file = file, method = methods, threshold = thresholds, 
                      avg_FPR = as.numeric(avg_FPR), avg_TPR = as.numeric(avg_TPR),
                      avg_precision = as.numeric(avg_precision), avg_recall = as.numeric(avg_recall))

ROC_data$avg_FPR <- as.numeric(ROC_data$avg_FPR)
ROC_data$avg_TPR <- as.numeric(ROC_data$avg_TPR)

ROC_plot <- ggplot(ROC_data, aes(x = avg_FPR, y = avg_TPR, color = method)) +
  geom_point(size = 3) +
  stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
  labs(title = paste("ROC: cluster = ", cluster), x = "FPR = FP/(FP+TN)", y = "TPR = TP/(TP+FN)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +   # Add y=x line (no skill)
  lims(x = c(0, 1), y = c(0, 1)) + # Set x and y-axis limits
  theme(text = element_text(size = 24))  # Set the base size for all text elements

PR_plot <- ggplot(ROC_data, aes(x = avg_recall, y = avg_precision, color = method)) +
  geom_point(size = 3) +
  stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
  labs(title = paste("P-R: cluster = ", cluster), x = "Recall = TP/(TP+FN)", y = "Precision = TP/(TP+FP)") +
  geom_abline(intercept = 0.5, slope = 0, linetype = "dashed", color = "gray") +   # Add y=0.5 line (no skill)
  lims(x = c(0, 1), y = c(0, 1)) + # Set x and y-axis limits
  theme(text = element_text(size = 24))  # Set the base size for all text elements


if (cluster) {
  ggsave(ROC_plot, filename = paste0(thisDir, "ROC_plot_cluster.png"))
} else {
  ggsave(ROC_plot, filename = paste0(thisDir, "ROC_plot.png"))
}

if (cluster) {
  ggsave(PR_plot, filename = paste0(thisDir, "PR_plot_cluster.png"))
} else {
  ggsave(PR_plot, filename = paste0(thisDir, "PR_plot.png"))
}
