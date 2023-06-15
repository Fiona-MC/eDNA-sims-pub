library(tidyr)
library(ggplot2)
library(gridExtra)


dirNums <- c("_ParmSet2")

multiSimLogistic <- data.frame()
dirNumL <- c()
for (dirNum in dirNums) {
  thisDir <- paste0("/space/s1/fiona_callahan/multiSim", dirNum, "/")
  thisMultiSimLogistic <- read.csv(paste0(thisDir, "logistic_mistakes.csv"), header = TRUE)
  thisMultiSimLogistic <- thisMultiSimLogistic[!is.na(thisMultiSimLogistic$sim_run), ]
  dirNumL <- c(dirNumL, rep(dirNum, times = length(thisMultiSimLogistic$sim_run)))
  multiSimLogistic <- plyr::rbind.fill(multiSimLogistic, thisMultiSimLogistic)
}
multiSimLogistic$dirNum <- dirNumL
multiSimLogistic$totalMistakes <- multiSimLogistic$num_incorrectInferences + multiSimLogistic$num_missedEffectsL
multiSimLogistic$fp_fp_tp <- multiSimLogistic$num_incorrectInferences / 
                              (multiSimLogistic$num_correctInferences + multiSimLogistic$num_incorrectInferences)

#write.csv(multiSimLogistic, file = "/home/fiona_callahan/eDNA_sims_code/multiSimLogistic_7_8_10_11_11a.csv")

actual_mean_fpr <- rep(NA, times = dim(multiSimLogistic)[1])
for (row in seq_len(dim(multiSimLogistic)[1])) { #1:dim(multiSimRes)[1]
    if(multiSimLogistic$fpr.mode[row] == "constant") {
      thisFPR <- multiSimLogistic$fpr.constant_fpr[row]
    } else if (multiSimLogistic$fpr.mode[row] == "none") {
      thisFPR <- 0
    } else if (multiSimLogistic$fpr.mode[row] == "dependent_sp") {
      thisFPR <- multiSimLogistic$fpr.mean_fpr[row]
    }
    actual_mean_fpr[row] <- thisFPR
}
multiSimLogistic$actual_mean_fpr <- actual_mean_fpr

hist(multiSimLogistic$totalMistakes)
mean(multiSimLogistic$totalMistakes)
hist(multiSimLogistic$num_incorrectInferences)
hist(multiSimLogistic$fp_fp_tp)
