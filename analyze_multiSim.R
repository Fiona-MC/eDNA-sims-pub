library(ranger)
library(dplyr)
library(pdp)
library(ggplot2)
library(plyr)
library(gridExtra)

dirNums<-c(7, 8, 10, 11)
multiSimRes <- data.frame()

dirNumL <- c()
for(dirNum in dirNums) {
  thisDir <- paste0("/space/s1/fiona_callahan/multiSim", dirNum, "/")
  thisMultiSimRes <- read.csv(paste0(thisDir, "infResGathered.csv"), header = TRUE)
  dirNumL <- c(dirNumL, rep(dirNum, times = length(thisMultiSimRes$RunNum)))
  multiSimRes <- plyr::rbind.fill(multiSimRes, thisMultiSimRes)
}

multiSimRes$dirNum <- dirNumL

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
randomParms <- c("num_samples_time", "num_samples_space", "radius", "fpr.mean_fpr", "fpr.mode", "covNoise_sd", "covMeasureNoise_sd", "r", "sigma", "N_50", "mean_mig_rate", "c2")
randomParmsVary <- as.vector(sapply(lapply(multiSimRes[, randomParms], unique), length) > 1)
randomParms <- randomParms[as.vector(randomParmsVary)] # take out any that I've changed to not be random

emergentParms <- c("Sp1PercentPresence", "Sp2PercentPresence", "Sp3PercentPresence")
# make formula from these parms
indep_var <- "totalMistakes"
#indep_var <- "finished_INLA"
fmlaParms<-c(randomParms,emergentParms)
fmla <- as.formula(paste(indep_var, " ~ ", paste(fmlaParms, collapse = "+")))

multiSimRes_na.rm <- multiSimRes[!is.na(multiSimRes[indep_var]), ]

#linear model just to look at
print(fmla)
lm_res <- lm(data = multiSimRes, formula = fmla)
summary(lm_res)

# random forest res
rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "impurity") 
summary(rf_res)
sort(importance(rf_res))
importanceDF<-as.data.frame(sort(importance(rf_res)))
importanceDF$Variable<-row.names(importanceDF)
names(importanceDF) <- c("RF_Importance", "Parameter")

better_parmNames <- data.frame(Parameter = fmlaParms, Parameter_name = c("n_time", "n_space", "neighborhood radius", "false detection rate", "false detect correlation", 
                    "covariate process noise", "covariate measurement noise", "r", "sigma", 
                    "false absence fcn parm", "migration rate", "inter-species affect scaling", 
                    "Sp1 Percent Presence", "Sp2 Percent Presence", "Sp3 Percent Presence"))

importanceDF <- merge(x = importanceDF, y = better_parmNames, by = "Parameter")

rf_importance_plot <- ggplot(importanceDF, aes(x = reorder(Parameter_name, RF_Importance), y=RF_Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(axis.text.y = element_text(hjust = 1, size = 24), axis.title.x = element_text(size = 24), axis.title.y = element_blank()) +
  coord_flip() +
  labs(y = "Random Forest Importance") 
ggsave(rf_importance_plot, file = "/space/s1/fiona_callahan/rf_importance.png", height = 10, width = 10)

# pdp plotting
pdp_resL <- list()

#for(parm in parms){
for (parm in fmlaParms) {
  pdp_res <- pdp::partial(object = rf_res, pred.var = parm, plot = FALSE)
  pdp_resL[[parm]] <- pdp_res
}

####### plot ##########
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  pdp_res <- pdp_resL[[parm]]
  p <- ggplot(pdp_res, aes(!!sym(parm), yhat)) +
    geom_point(size = 3) +
    geom_line() +
    labs(x = better_parmNames$Parameter_name[better_parmNames$Parameter == parm], y = paste0("Predicted ", indep_var)) +
    theme(axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 18))
  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
rf_marginal <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(rf_marginal, file = "/space/s1/fiona_callahan/rf_marginal.png", width = 40, height = 8)





#plot(multiSimRes$N_50, multiSimRes$totalMistakes)
multiSimRes_na.rm$totalMistakes.factor <- factor(multiSimRes_na.rm$totalMistakes, levels = 11:0)
multiSimRes_na.rm$total_samples <- multiSimRes_na.rm$num_samples_space * multiSimRes_na.rm$num_samples_time

fmlaParms
parm = "total_samples"
#violin
ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = !!sym(parm))) +
  geom_violin(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_missedEffects_alpha)))+
  coord_flip() 

# vertically aligned hists
#colored by total missed effects
ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), fill = as.factor(num_missedEffectsL))) + 
  geom_histogram() +
  facet_grid(rows = vars(totalMistakes.factor)) +
  scale_fill_viridis_d(name = "Missed effects") +
  labs(y = "Number of simulations") +
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = NULL, name = "Number of Mistakes in Inference (Type 1 and 2)")) 

# colored by alpha missed effects
ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), fill = as.factor(num_missedEffects_alpha))) + 
  geom_histogram() +
  facet_grid(rows = vars(totalMistakes.factor)) +
  scale_fill_viridis_d(name = "Missed effects Alpha") +
  labs(y = "Number of simulations") +
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = NULL, name = "Number of Mistakes in Inference (Type 1 and 2)")) 

# ONLY false inferences (type 1 errors)
ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), fill = as.factor(num_missedEffectsL))) + 
  geom_histogram() +
  facet_grid(rows = vars(factor(num_incorrectInferences, levels = 10:0))) +
  scale_fill_viridis_d(name = "Missed effects") +
  labs(y = "Number of simulations") +
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = NULL, name = "Number of Mistakes in Inference (Type 1)")) 


# violin with points
ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = !!sym(parm))) +
  geom_violin(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_incorrectInferences)))+
  coord_flip() 

# box plots with points
ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = !!sym(parm))) +
  geom_boxplot(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_incorrectInferences)))+
  coord_flip() 


# histograms of variables
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
