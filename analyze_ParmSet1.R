library(tidyr)
library(ggplot2)
library(gridExtra)

dirName <- c("dumbSim2")
dirNums <- dirName

# INLA
multiSimRes <- data.frame()

dirNumL <- c()
for (dirNum in dirNums) {
  thisDir <- paste0("/space/s1/fiona_callahan/", dirName, "/")
  thisMultiSimRes <- read.csv(paste0(thisDir, "infResGathered.csv"), header = TRUE)
  dirNumL <- c(dirNumL, rep(dirNum, times = length(thisMultiSimRes$RunNum)))
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
    }
    actual_mean_fpr[row] <- thisFPR
}
multiSimRes$actual_mean_fpr <- actual_mean_fpr

multiSimRes$fp_fp_tp <- multiSimRes$num_incorrectInferences / (multiSimRes$num_correctInferences + multiSimRes$num_incorrectInferences)

hist(multiSimRes$totalMistakes)
hist(multiSimRes$num_incorrectInferences)
hist(multiSimRes$fp_fp_tp)
hist(multiSimRes$num_incorrect_alpha)
hist(multiSimRes$num_incorrect_beta)

mean(multiSimRes$fp_fp_tp, na.rm = TRUE) # with no violated assumptions should be 0.082



#Logistic
multiSimLogistic <- data.frame()
dirNumL <- c()
for (dirNum in dirNums) {
  thisDir <- paste0("/space/s1/fiona_callahan/", dirName, "/")
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
hist(multiSimLogistic$num_incorrectInferences)
hist(multiSimLogistic$fp_fp_tp)
hist(multiSimLogistic$num_incorrect_alpha)
hist(multiSimLogistic$num_incorrect_beta)

mean(multiSimLogistic$fp_fp_tp) # with no violated assumptions --should be 0.082

multiSimLogistic$fp_fp_tp_beta <- multiSimLogistic$num_incorrect_beta / 
                              ((3 - multiSimLogistic$num_missedEffects_beta) + multiSimLogistic$num_incorrect_beta)
mean(multiSimLogistic$fp_fp_tp_beta) # expectation is 0.075 if all assumptions are met
mean(multiSimLogistic$totalMistakes)



# expectation of fp_tp_fp
sum <- 0
for(fp in 1:10) {
  for(tp in 1:5) {
    # binomial probability that fp=fp times prob(tp=tp|fp=fp)
    thisProb <- dbinom(size = 10, x = fp, prob = 0.05) * dbinom(size = 5, x = tp, prob = 0.95) 
    thisExp <- (fp / (fp + tp)) * thisProb
    sum <- sum + thisExp
  }
}
sum # 0.082

# expectation of fp_tp_fp beta
sum <- 0
for(fp in 1:6) {
  for(tp in 1:3) {
    # binomial probability that fp=fp times prob(tp=tp|fp=fp)
    thisProb <- dbinom(size = 6, x = fp, prob = 0.05) * dbinom(size = 3, x = tp, prob = 0.95) 
    thisExp <- (fp / (fp + tp)) * thisProb
    sum <- sum + thisExp
  }
}
sum # 0.075
mean(multiSimLogistic$fp_fp_tp_beta) # expectation is 0.075 if all assumptions are met
