# sampling with or without replacement 
pool_of_sample <- c(rep("b", times = 2), rep("r", times = 118))
reps <- 10000
nReps <- 1
blueCountMeanL <- rep(NA, times = reps)
for(jj in 1:reps) {
  blue_count <- rep(NA, times = nReps)
  for(iii in 1:nReps) {
    thisSample <- sample(x = pool_of_sample, size = 50, replace = FALSE)
    blue_count[iii] <- sum(thisSample == "b")
  }
  blueCountMeanL[jj] <- mean(blue_count) / 50
}
hist(blueCountMeanL)
mean(blueCountMeanL)
2 / 120