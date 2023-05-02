library(ranger)
library(dplyr)
library(pdp)
library(ggplot2)

thisDir <- "/space/s1/fiona_callahan/multiSim7/"
multiSimRes <- read.csv(paste0(thisDir, "infResGathered.csv"), header = TRUE)

# put sum of two types of mistakes in a column
multiSimRes$totalMistakes <- multiSimRes$num_incorrectInferences + multiSimRes$num_missedEffectsL

lowMistakes_res_lt2 <- multiSimRes[multiSimRes$totalMistakes < 2 & !is.na(multiSimRes$totalMistakes), ]
lowMistakes_res_eq2 <- multiSimRes[multiSimRes$totalMistakes == 2 & !is.na(multiSimRes$totalMistakes), ]

#multiSimRes %>% 
#  summarise_if(is.numeric, c(mean, sd), na.rm = TRUE)

summary(multiSimRes)

parmNames <- readRDS(paste0(thisDir, "parmNames.Rdata"))
# this removes duplicates and makes it so that none of the parm names are not in the actual data set
parmNames <- names(multiSimRes[, parmNames])
# limit list to the parameters that I actually allowed to vary in the sim
#TODO some of these parms are functions of other parms (moment matching)
parmsVary <- as.vector(sapply(lapply(multiSimRes[, parmNames], unique), length) > 1)
parmNames <- parmNames[as.vector(parmsVary)]
#TODO these parm names need to be adjusted because some of them were inside other lists
#TODO c("mean_fpr", "N_50", "mean_mig_rate") were not in the list of parms because they are fcns of other parms -- calculate if they don't exist in list
randomParms <- c("num_samples_time", "num_samples_space", "radius", "fpr.mean_fpr", "fpr.mode", "covNoise_sd", "covMeasureNoise_sd", "r", "sigma", "N_50", "mean_mig_rate")
emergentParms <- c("Sp1PercentPresence", "Sp2PercentPresence", "Sp3PercentPresence")
# make formula from these parms
indep_var <- "totalMistakes"
#indep_var <- "finished_INLA"
fmla <- as.formula(paste(indep_var, " ~ ", paste(c(randomParms,emergentParms), collapse = "+")))

multiSimRes_na.rm <- multiSimRes[!is.na(multiSimRes[indep_var]), ]

#linear model just to look at
print(fmla)
lm_res <- lm(data = multiSimRes, formula = fmla)
summary(lm_res)

# random forest res
rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "impurity") 
summary(rf_res)
sort(importance(rf_res))

# pdp plotting
pdp_res <- pdp::partial(object = rf_res, pred.var = "num_samples_time", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "num_samples_space", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "N_50", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "r", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "Sp3PercentPresence", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "Sp2PercentPresence", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "Sp1PercentPresence", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "sigma", grid.resolution = 50, plot = TRUE)
pdp_res <- pdp::partial(object = rf_res, pred.var = "covMeasureNoise_sd", grid.resolution = 50, plot = TRUE)


#plot(multiSimRes$N_50, multiSimRes$totalMistakes)
multiSimRes_na.rm$totalMistakes.factor <- as.factor(multiSimRes_na.rm$totalMistakes)
ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = N_50)) +
  geom_violin(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_incorrectInferences)))+
  coord_flip() 

ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = N_50)) +
  geom_boxplot(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_incorrectInferences)))+
  coord_flip() +
  theme(legend.position = "none")

ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = sigma)) +
  geom_violin(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_incorrectInferences)))+
  coord_flip() 

ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = num_samples_space)) +
  geom_violin(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_incorrectInferences)))+
  coord_flip() 

hist(multiSimRes$num_samples_space)
hist(multiSimRes$num_samples_space[!multiSimRes$finished_INLA])

total_samples <- multiSimRes$num_samples_space * multiSimRes$num_samples_time
hist(total_samples)
total_samples_notFinished <- multiSimRes$num_samples_space[!multiSimRes$finished_INLA] * multiSimRes$num_samples_time[!multiSimRes$finished_INLA]
hist(total_samples_notFinished)

#killed_runNums <- c(447, 254, 487, 169, 426, 143, 102, 385, 373, 335)
#multiSimRes$killed <- multiSimRes$RunNum %in% killed_runNums
#multiSimRes[multiSimRes$killed, ]

multiSimResMeans <- multiSimRes %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
multiSimResMeans_killed <- multiSimRes[multiSimRes$RunNum %in% killed_runNums, ] %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
mean_df <- data.frame(varNames = names(multiSimResMeans), 
            allMeans = format(as.vector(unlist(multiSimResMeans[1,])), scientific = FALSE), 
            multiSimResMeans_killed = format(as.vector(unlist(multiSimResMeans_killed)), scientific = FALSE))
mean_df[mean_df$allMeans != mean_df$multiSimResMeans_killed, ]
# multiSimRes$X.1  #this is just 1:2000
hist(multiSimRes$totalMistakes)

mean_df$varNames[mean_df$allMeans != mean_df$multiSimResMeans_killed]
var <- "sigma"
var <- "r"
var <- "fpr.alpha_fpr"
var <- "covNoise_sd" #!!!!! this could honestly just be by chance though -- a mystery
var <- "num_samples_space"
var <- "radius"
means_random <- replicate(1000, mean(sample(multiSimRes[, var], 10, replace = FALSE)))
hist(means_random)
abline(v = mean(multiSimRes[multiSimRes$RunNum %in% killed_runNums, var]), col = 'red', lwd = 3, lty = 'dashed')


#TODO
#1 figure out how to look at variable importance and threshold values
#2 figure out if there are other parms I should set for Ranger rf fucntion

#3 figure out how to do basically the same thing with lasso regression
