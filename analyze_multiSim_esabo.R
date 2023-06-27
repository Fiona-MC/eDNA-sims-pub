library(ranger)
library(dplyr)
library(pdp)
library(ggplot2)
library(plyr)
library(gridExtra)
library(tidyr)

dirNums <- c(11)
multiSimRes <- data.frame()

dirNumL <- c()
for (dirNum in dirNums) {
  thisDir <- paste0("/space/s1/fiona_callahan/multiSim", dirNum, "/")
  thisMultiSimRes <- read.csv(paste0(thisDir, "ESABO_mistakes.csv"), header = TRUE)
  thisMultiSimRes <- thisMultiSimRes[!is.na(thisMultiSimRes$sim_run), ]
  dirNumL <- c(dirNumL, rep(dirNum, times = length(thisMultiSimRes$sim_run)))
  multiSimRes <- plyr::rbind.fill(multiSimRes, thisMultiSimRes)
}
multiSimRes$dirNum <- dirNumL

multiSimRes$totalMistakes <- multiSimRes$num_incorrectInferences + multiSimRes$num_missedEffectsL

actual_mean_fpr <- rep(NA, times = dim(multiSimRes)[1])
for (row in seq_len(dim(multiSimRes)[1])) { #1:dim(multiSimRes)[1]
    if(multiSimRes$fpr.mode[row] == "constant") {
      thisFPR <- multiSimRes$fpr.constant_fpr[row]
    } else if (multiSimRes$fpr.mode[row] == "none") {
      thisFPR <- 0
    } else if (multiSimRes$fpr.mode[row] == "dependent_sp") {
      thisFPR <- multiSimRes$fpr.mean_fpr[row]
    } else if (multiSimRes$fpr.mode[row] == "independent") {
      thisFPR <- multiSimRes$fpr.mean_fpr[row]
    }
    actual_mean_fpr[row] <- thisFPR
}
multiSimRes$actual_mean_fpr <- actual_mean_fpr

multiSimRes$fp_fp_tp <- multiSimRes$num_incorrectInferences / (multiSimRes$num_correctInferences + multiSimRes$num_incorrectInferences)

hist(multiSimRes$fp_fp_tp)
mean(multiSimRes$fp_fp_tp, na.rm = TRUE) # 0.596 for ~1000 sims in folder 11




# null expectations
# expectation of fp_tp_fp alpha
sum <- 0
num_possible_tp <- 2
num_possible_fp <- 4
level <- pnorm(-1) * 2
for(fp in 1:num_possible_fp) {
  for(tp in 1:num_possible_tp) {
    # binomial probability that fp=fp times prob(tp=tp|fp=fp)
    thisProb <- dbinom(size = num_possible_fp, x = fp, prob = level) * dbinom(size = num_possible_tp, x = tp, prob = 1 - level) 
    thisExp <- (fp / (fp + tp)) * thisProb
    sum <- sum + thisExp
  }
}
sum # 0.357 for esabo with reject level z = \pm 1

# expectation of fp_tp_fp alpha
sum <- 0
num_possible_tp <- 2
num_possible_fp <- 4
level <- pnorm(-2) * 2
for(fp in 1:num_possible_fp) {
  for(tp in 1:num_possible_tp) {
    # binomial probability that fp=fp times prob(tp=tp|fp=fp)
    thisProb <- dbinom(size = num_possible_fp, x = fp, prob = level) * dbinom(size = num_possible_tp, x = tp, prob = 1 - level) 
    thisExp <- (fp / (fp + tp)) * thisProb
    sum <- sum + thisExp
  }
}
sum # 0.061 for esabo with reject level z = \pm 2