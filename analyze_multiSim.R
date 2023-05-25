library(ranger)
library(dplyr)
library(pdp)
library(ggplot2)
library(plyr)
library(gridExtra)

dirNums <- c(7, 8, 10, 11, "11a")
multiSimRes <- data.frame()

dirNumL <- c()
for (dirNum in dirNums) {
  thisDir <- paste0("/space/s1/fiona_callahan/multiSim", dirNum, "/")
  thisMultiSimRes <- read.csv(paste0(thisDir, "infResGathered.csv"), header = TRUE)
  dirNumL <- c(dirNumL, rep(dirNum, times = length(thisMultiSimRes$RunNum)))
  multiSimRes <- plyr::rbind.fill(multiSimRes, thisMultiSimRes)
}
multiSimRes$dirNum <- dirNumL

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

# put sum of two types of mistakes in a column
multiSimRes$totalMistakes <- multiSimRes$num_incorrectInferences + multiSimRes$num_missedEffectsL

lowMistakes_res_lt2 <- multiSimRes[multiSimRes$totalMistakes < 2 & !is.na(multiSimRes$totalMistakes), ]
lowMistakes_res_eq2 <- multiSimRes[multiSimRes$totalMistakes == 2 & !is.na(multiSimRes$totalMistakes), ]

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
#TODO c("mean_fpr", "N_50", "mean_mig_rate") were not in the list of parms because they are fcns of other parms 
# -- calculate if they don't exist in list
randomParms <- c("num_samples_time", "num_samples_space", "radius", 
                "actual_mean_fpr", "fpr.mode", "covNoise_sd", "covMeasureNoise_sd", 
                "r", "sigma", "N_50", "mean_mig_rate", "c2")
randomParmsVary <- as.vector(sapply(lapply(multiSimRes[, randomParms], unique), length) > 1)
randomParms <- randomParms[as.vector(randomParmsVary)] # take out any that I've changed to not be random

emergentParms <- c("Sp1PercentPresence", "Sp2PercentPresence", "Sp3PercentPresence")
# make formula from these parms
indep_var <- "totalMistakes"
#indep_var <- "fp_fp_tp"
#indep_var <- "finished_INLA"
fmlaParms <- c(randomParms)
fmla <- as.formula(paste(indep_var, " ~ ", paste(fmlaParms, collapse = "+")))

#remove NA's in response (run didnt finish or whatever)
multiSimRes_na.rm <- multiSimRes[!is.na(multiSimRes[[indep_var]]), ]

# inla average success
mean(multiSimRes_na.rm$fp_fp_tp, na.rm = TRUE)
quantile(multiSimRes_na.rm$fp_fp_tp, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
quantile(multiSimRes_na.rm$totalMistakes, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
mean(multiSimRes_na.rm$num_correctInferences + multiSimRes_na.rm$num_incorrectInferences, na.rm = TRUE)
mean(multiSimRes_na.rm$num_incorrectInferences, na.rm = TRUE)

# compare INLA to logistic
mean(multiSimLogistic$fp_fp_tp, na.rm = TRUE)
quantile(multiSimLogistic$fp_fp_tp, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
quantile(multiSimLogistic$totalMistakes, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
mean(multiSimLogistic$num_correctInferences + multiSimLogistic$num_incorrectInferences, na.rm = TRUE)
mean(multiSimLogistic$num_incorrectInferences, na.rm = TRUE)

multiSimResMerged <- merge(x = multiSimLogistic, y = multiSimRes_na.rm, 
                    by.x = c("sim_run", "dirNum", "trial"), by.y = c("RunNum", "dirNum", "trial"), suffixes = c("_logistic", "_INLA"))
mean(multiSimResMerged$num_incorrectInferences_logistic - multiSimResMerged$num_incorrectInferences_INLA, na.rm = TRUE)
mean(multiSimResMerged$totalMistakes_logistic - multiSimResMerged$totalMistakes_INLA, na.rm = TRUE)
mean(multiSimResMerged$fp_fp_tp_logistic - multiSimResMerged$fp_fp_tp_INLA, na.rm = TRUE)

t.test(x = multiSimResMerged$num_incorrectInferences_logistic, y = multiSimResMerged$num_incorrectInferences_INLA, paired = TRUE, na.rm = TRUE)
t.test(x = multiSimResMerged$totalMistakes_logistic, y = multiSimResMerged$totalMistakes_INLA, paired = TRUE, na.rm = TRUE)
t.test(x = multiSimResMerged$fp_fp_tp_logistic, y = multiSimResMerged$fp_fp_tp_INLA, paired = TRUE, na.rm = TRUE)

############ REMOVE fpr.mode = "constant" and "none" --correlation issues
### vvvvv this made no difference to the order of the "permutation" variable importances
#multiSimRes_na.rm <- multiSimRes_na.rm[multiSimRes_na.rm$fpr.mode != "none" & multiSimRes_na.rm$fpr.mode != "constant", ]

#linear model just to look at
print(fmla)
lm_res <- lm(data = multiSimRes_na.rm, formula = fmla)
lmSum <- summary(lm_res)
# order by pvalue
lmSum$coefficients[order(lmSum$coefficients[, 4]), ]

######## lm normalized ############
multiSimRes_na.rm.columnsNormalized <- multiSimRes_na.rm
for (colNum in seq_len(dim(multiSimRes_na.rm)[2])) { #1:dim(multiSimRes_na.rm)[2]
  if (names(multiSimRes_na.rm)[colNum] %in% fmlaParms && is.numeric(multiSimRes_na.rm[, colNum])) {
    var <- multiSimRes_na.rm.columnsNormalized[, colNum]
    var_normalized <- (var - mean(var)) / sd(var)
    multiSimRes_na.rm.columnsNormalized[, names(multiSimRes_na.rm)[colNum]] <- var_normalized
  }
}

lm_res_norm <- lm(data = multiSimRes_na.rm.columnsNormalized, formula = fmla)
lmSum <- summary(lm_res_norm)
# order by pvalue
lmSum$coefficients[order(lmSum$coefficients[, 4]), ]

# Extract coefficient estimates and standard errors
coef_data <- data.frame(
  coefficient = names(coef(lm_res_norm)),
  estimate = coef(lm_res_norm),
  std_error = summary(lm_res_norm)$coef[, "Std. Error"]
)
coef_data <- coef_data[coef_data$coefficient != "(Intercept)", ]

# Plot betas with error bars
lm_betaPlot <- ggplot(coef_data, aes(x = coefficient, y = estimate, ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  labs(x = "Coefficient", y = "Estimate") +
  coord_flip() +
  geom_hline(yintercept = 0, color = "red") +
  theme(text = element_text(size = 36))  

ggsave(lm_betaPlot, file = "/space/s1/fiona_callahan/lm_betaPlot.png", height = 12, width = 12)



# random forest res
#rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "impurity") 
rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "permutation") 
summary(rf_res)
sort(importance(rf_res))
importanceDF <- as.data.frame(sort(importance(rf_res)))
importanceDF$Variable <- row.names(importanceDF)
names(importanceDF) <- c("RF_Importance", "Parameter")


better_parmNames <- data.frame(Parameter = fmlaParms, Parameter_name = c("n_time", "n_space", 
                    "neighborhood radius", "false detection rate", "false detect correlation", 
                    "covariate process noise", "covariate measurement noise", "r", "sigma", 
                    "false absence parameter", "migration rate", "inter-species affect scaling"))#, 
                    #"Sp1 Percent Presence", "Sp2 Percent Presence", "Sp3 Percent Presence"))

importanceDF <- merge(x = importanceDF, y = better_parmNames, by = "Parameter")

rf_importance_plot <- ggplot(importanceDF, aes(x = reorder(Parameter_name, RF_Importance), y = RF_Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme(axis.text.y = element_text(hjust = 1, size = 24), axis.title.x = element_text(size = 24), axis.title.y = element_blank()) +
  coord_flip() +
  labs(y = "Random Forest Importance", title = paste("Random Forest Predictor Importance for", indep_var)) 
ggsave(rf_importance_plot, file = "/space/s1/fiona_callahan/rf_importance_noemergent_perm.png", height = 10, width = 10)

########### test set ###############
########### test set ###############
########### test set ###############
thisDir <- paste0("/space/s1/fiona_callahan/multiSim11_testSet/")
multiSimRes_testData <-  read.csv(paste0(thisDir, "infResGathered.csv"), header = TRUE)
# put sum of two types of mistakes in a column
multiSimRes_testData$totalMistakes <- multiSimRes_testData$num_incorrectInferences + multiSimRes_testData$num_missedEffectsL

actual_mean_fpr <- rep(NA, times = dim(multiSimRes_testData)[1])
for(row in seq_len(dim(multiSimRes_testData)[1])) { #1:dim(multiSimRes_testData)[1]
    if(multiSimRes_testData$fpr.mode[row] == "constant") {
      thisFPR <- multiSimRes_testData$fpr.constant_fpr[row]
    } else if(multiSimRes_testData$fpr.mode[row] == "none") {
      thisFPR <- 0
    } else if(multiSimRes_testData$fpr.mode[row] == "dependent_sp") {
      thisFPR <- multiSimRes_testData$fpr.mean_fpr[row]
    }
    actual_mean_fpr[row] <- thisFPR
}
multiSimRes_testData$actual_mean_fpr <- actual_mean_fpr

multiSimRes_testData$fp_fp_tp <- multiSimRes_testData$num_incorrectInferences / 
                                (multiSimRes_testData$num_correctInferences + multiSimRes_testData$num_incorrectInferences)
#remove NA's in response (run didnt finish or whatever)
multiSimRes_testData_na.rm <- multiSimRes_testData[!is.na(multiSimRes_testData[[indep_var]]), ]


naive_rmse <- sqrt(sum((4 - multiSimRes_testData_na.rm$totalMistakes)^2) / length(multiSimRes_testData_na.rm$totalMistakes))

lm_predictions <- predict(object = lm_res, newdata = multiSimRes_testData_na.rm)
lm_rmse <- sqrt(sum((lm_predictions - multiSimRes_testData_na.rm$totalMistakes)^2) / length(multiSimRes_testData_na.rm$totalMistakes))

rf_predictions <- predict(object = rf_res, data = multiSimRes_testData_na.rm)
rf_rmse <- sqrt(sum((rf_predictions$predictions - multiSimRes_testData_na.rm$totalMistakes)^2) / length(multiSimRes_testData_na.rm$totalMistakes))

(rf_rmse - naive_rmse) / naive_rmse
(rf_rmse - lm_rmse) / lm_rmse

########### end test set ###############
########### end test set ###############
########### end test set ###############

# pdp plotting
pdp_resL <- list()
ice_permL <- list()
cice_permL <- list()

#for(parm in parms){
for (parm in fmlaParms) {
  pdp_res <- pdp::partial(object = rf_res, pred.var = parm, plot = FALSE)
  pdp_resL[[parm]] <- pdp_res
  # ice plots
  ice_perm <- partial(object = rf_res, ice = TRUE, pred.var = parm, plot = FALSE)
  ice_permL[[parm]] <- ice_perm
  #cice plots
  cice_perm <- partial(object = rf_res, ice = TRUE, center = TRUE, pred.var = parm, plot = FALSE)
  cice_permL[[parm]] <- cice_perm
}

#ICE plots
ICEplotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  ice_perm <- ice_permL[[parm]]
  ice_perm_subset <- ice_perm[ice_perm$yhat.id %% 30 == 1, ]
  p <- ggplot(ice_perm_subset, aes(x = !!sym(parm), y = yhat, group = yhat.id)) +
      geom_line() +
      coord_cartesian(ylim = c(2, 8)) +
      #coord_cartesian(ylim = c(.2, .6)) +
      labs(x = better_parmNames$Parameter_name[better_parmNames$Parameter == parm], y = paste0("Predicted ", indep_var)) +
      theme(axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 18))

  ICEplotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
rf_marginal_ice <- grid.arrange(grobs = ICEplotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(rf_marginal_ice, file = "/space/s1/fiona_callahan/rf_marginal_ice.png", width = 40, height = 12)

#cice
ICEplotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  ice_perm <- cice_permL[[parm]]
  ice_perm_subset <- ice_perm[ice_perm$yhat.id %% 30 == 1, ]
  p <- ggplot(ice_perm_subset, aes(x = !!sym(parm), y = yhat, group = yhat.id)) +
      geom_line() +
      #coord_cartesian(ylim = c(-.2, .3)) +
      coord_cartesian(ylim = c(-2, 3)) +
      labs(x = better_parmNames$Parameter_name[better_parmNames$Parameter == parm], y = paste0("Predicted ", indep_var, "(centered)")) +
      theme(axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 18))

  ICEplotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
rf_marginal_cice <- grid.arrange(grobs = ICEplotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(rf_marginal_cice, file = "/space/s1/fiona_callahan/rf_marginal_cice.png", width = 40, height = 12)


####### plot marginal rf ##########
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  pdp_res <- pdp_resL[[parm]]
  p <- ggplot(pdp_res, aes(!!sym(parm), yhat)) +
    geom_point(size = 3) +
    geom_line() +
    coord_cartesian(ylim = c(3.5, 4.75)) +
    #coord_cartesian(ylim = c(.3, .5)) +
    labs(x = better_parmNames$Parameter_name[better_parmNames$Parameter == parm], y = paste0("Predicted ", indep_var)) +
    theme(axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 36), axis.text.x = element_text(size = 20))
  if (parm == "fpr.mode") { # fix angle
    p <- ggplot(pdp_res, aes(!!sym(parm), yhat)) +
      geom_point(size = 3) +
      #geom_line() +
      coord_cartesian(ylim = c(3.5, 4.75)) +
      #coord_cartesian(ylim = c(.3, .5)) +
      labs(x = better_parmNames$Parameter_name[better_parmNames$Parameter == parm], y = paste0("Predicted ", indep_var)) +
      theme(axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 36), axis.text.x = element_text(size = 20, angle = 10))
  }
  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
rf_marginal <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(rf_marginal, file = "/space/s1/fiona_callahan/rf_marginal_noemergent_perm.png", width = 40, height = 12)



############ plot marginal data ##############

#plot(multiSimRes$N_50, multiSimRes$totalMistakes)
multiSimRes_na.rm$totalMistakes.factor <- factor(multiSimRes_na.rm$totalMistakes, levels = 11:0)
multiSimRes_na.rm$total_samples <- multiSimRes_na.rm$num_samples_space * multiSimRes_na.rm$num_samples_time

fmlaParms
parm <- "fpr.mean_fpr"
#violin
ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = !!sym(parm))) +
  geom_violin(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_missedEffects_alpha))) +
  coord_flip() 

# vertically aligned hists
#colored by total missed effects
ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), fill = as.factor(num_missedEffectsL))) + 
  geom_histogram(bins = 10) +
  facet_grid(rows = vars(totalMistakes.factor)) +
  scale_fill_viridis_d(name = "Missed effects") +
  labs(y = "Number of simulations") +
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = NULL, name = "Number of Mistakes in Inference (Type 1 and 2)")) 

#as frequencies
ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), y = ..density.., fill = as.factor(num_missedEffectsL))) + 
  geom_histogram(bins = 10) +
  facet_grid(rows = vars(totalMistakes.factor)) +
  scale_fill_viridis_d(name = "Missed effects") +
  labs(y = "Number of simulations") +
  scale_y_continuous(sec.axis = sec_axis(~ ., breaks = NULL, name = "Number of Mistakes in Inference (Type 1 and 2)")) 


# colored by alpha missed effects
ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), fill = as.factor(num_missedEffects_alpha))) + 
  geom_histogram(bins = 15) +
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
  geom_point(aes(color = as.factor(num_incorrectInferences))) +
  coord_flip() 

# box plots with points
ggplot(multiSimRes_na.rm, aes(x = totalMistakes.factor, y = !!sym(parm))) +
  geom_boxplot(position = position_nudge()) +
  geom_point(aes(color = as.factor(num_incorrectInferences))) +
  coord_flip() 


plotL <- list()
orderImportances <- c()
for(parm in fmlaParms) {
  # hist of num_incorrectInferences
  if(class(multiSimRes_na.rm[[parm]]) == "character") {
    multiSimRes_na.rm$quantiles <- multiSimRes_na.rm[[parm]]
  } else {
    multiSimRes_na.rm$quantiles <- cut(multiSimRes_na.rm[[parm]], breaks = 2, labels = FALSE)
  }
  p <- ggplot(multiSimRes_na.rm, aes(x = num_incorrectInferences, fill = as.factor(quantiles))) + 
    geom_histogram(bins = 11, position = "dodge") +
    scale_fill_viridis_d(name = paste("Quantiles of", parm)) +
    labs(y = "Number of simulations", x = "Number of incorrect inferences (Type 1)") 

  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
hist_2q <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(hist_2q, file = "/space/s1/fiona_callahan/hist_2q.png", width = 40, height = 8)

#density
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms) {
  # hist of num_incorrectInferences
  if(class(multiSimRes_na.rm[[parm]]) == "character") {
    multiSimRes_na.rm$quantiles <- multiSimRes_na.rm[[parm]]
  }else{
    multiSimRes_na.rm$quantiles <- cut(multiSimRes_na.rm[[parm]], breaks = 3, labels = FALSE)
  }
  p <- ggplot(multiSimRes_na.rm, aes(x = num_incorrectInferences, y = ..density.., fill = as.factor(quantiles))) + 
    geom_histogram(bins = 11, position = "dodge") +
    scale_fill_viridis_d(name = paste("Quantiles of", parm)) +
    labs(x = "Number of incorrect inferences (Type 1)") 

  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
hist_3qd <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(hist_3qd, file = "/space/s1/fiona_callahan/hist_3qd.png", width = 40, height = 8)


#density for total mistakes
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms) {
  # hist of num_incorrectInferences
  if(class(multiSimRes_na.rm[[parm]]) == "character") {
    multiSimRes_na.rm$quantiles <- multiSimRes_na.rm[[parm]]
  } else {
    multiSimRes_na.rm$quantiles <- cut(multiSimRes_na.rm[[parm]], breaks = 3, labels = FALSE)
  }

  p <- ggplot(multiSimRes_na.rm, aes(x = totalMistakes, y = ..density.., fill = as.factor(quantiles))) + 
    geom_histogram(bins = 6, position = "dodge") +
    scale_fill_viridis_d(name = paste("Quantiles of", parm)) +
    labs(x = "Number of Total mistakes (Type 1 and 2)") 

  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
hist_3qd_totalmistakes <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(hist_3qd_totalmistakes, file = "/space/s1/fiona_callahan/hist_3qd_totalmistakes.png", width = 40, height = 8)


#density for total mistakes --with fill option
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms) {
  # hist of num_incorrectInferences
  if(class(multiSimRes_na.rm[[parm]]) == "character") {
    multiSimRes_na.rm$quantiles <- multiSimRes_na.rm[[parm]]
  } else {
    multiSimRes_na.rm$quantiles <- cut(multiSimRes_na.rm[[parm]], breaks = 3, labels = FALSE)
  }

  p <- ggplot(multiSimRes_na.rm, aes(x = totalMistakes, fill = as.factor(quantiles))) + 
    geom_histogram(bins = 6, position = "fill") +
    scale_fill_viridis_d(name = paste("Quantiles of", parm)) +
    labs(x = "Number of Total mistakes (Type 1 and 2)") 

  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
hist_3qd_totalmistakes_fill <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(hist_3qd_totalmistakes_fill, file = "/space/s1/fiona_callahan/hist_3qd_totalmistakes_fill.png", width = 40, height = 8)

#density for total mistakes --with fill option
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms) {
  # hist of num_incorrectInferences
  if(class(multiSimRes_na.rm[[parm]]) == "character") {
    multiSimRes_na.rm$quantiles <- multiSimRes_na.rm[[parm]]
  } else {
    multiSimRes_na.rm$quantiles <- cut(multiSimRes_na.rm[[parm]], breaks = 4, labels = FALSE)
  }

  p <- ggplot(multiSimRes_na.rm, aes(x = totalMistakes, fill = as.factor(quantiles))) + 
    geom_histogram(bins = 6, position = "fill") +
    scale_fill_viridis_d(name = paste("Quantiles of", parm)) +
    labs(x = "Number of Total mistakes (Type 1 and 2)") 

  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
hist_4qd_totalmistakes_fill <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(hist_4qd_totalmistakes_fill, file = "/space/s1/fiona_callahan/hist_4qd_totalmistakes_fill.png", width = 40, height = 8)


#stack
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms) {
  # hist of num_incorrectInferences
  if(class(multiSimRes_na.rm[[parm]]) == "character") {
    multiSimRes_na.rm$quantiles <- multiSimRes_na.rm[[parm]]
  } else {
    multiSimRes_na.rm$quantiles <- cut(multiSimRes_na.rm[[parm]], breaks = 3, labels = FALSE)
  }

  p <- ggplot(multiSimRes_na.rm, aes(x = totalMistakes, y = ..density.., fill = as.factor(quantiles))) + 
    geom_histogram(bins = 6, position = "stack") +
    scale_fill_viridis_d(name = paste("Quantiles of", parm)) +
    labs(x = "Number of Total mistakes (Type 1 and 2)") 

  plotL[parm] <- list(p)
  orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
}
hist_3qd_totalmistakes_stack <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(hist_3qd_totalmistakes_stack, file = "/space/s1/fiona_callahan/hist_3qd_totalmistakes_stack.png", width = 40, height = 8)






############## invert these plots #################
#density for total mistakes
plotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  if(class(multiSimRes_na.rm[[parm]]) != "character") {
    multiSimRes_na.rm$quantiles_response <- cut(multiSimRes_na.rm$totalMistakes, breaks = 3, labels = FALSE)

    p <- ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), y = ..density.., fill = as.factor(quantiles_response))) + 
      geom_histogram(bins = 6, position = "dodge") +
      scale_fill_brewer(palette = "Spectral", name = paste("Quantiles of Num mistakes")) 

    plotL[parm] <- list(p)
    orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
  }
}
histInverted_totalmistakes <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(histInverted_totalmistakes, file = "/space/s1/fiona_callahan/histInverted_totalmistakes.png", width = 40, height = 8)

plotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  if(class(multiSimRes_na.rm[[parm]]) != "character") {
    multiSimRes_na.rm$quantiles_response <- cut(multiSimRes_na.rm$totalMistakes, breaks = 3, labels = FALSE)

    p <- ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), y = ..density.., fill = as.factor(quantiles_response))) + 
      geom_histogram(bins = 6, position = "fill") +
      scale_fill_brewer(palette = "Spectral", name = paste("Quantiles of Num mistakes")) 

    plotL[parm] <- list(p)
    orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
  }
}
histInverted_totalmistakes <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(histInverted_totalmistakes, file = "/space/s1/fiona_callahan/histInverted_totalmistakes_stack.png", width = 40, height = 8)


plotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  if(class(multiSimRes_na.rm[[parm]]) != "character") {
    multiSimRes_na.rm$quantiles_response <- cut(multiSimRes_na.rm$totalMistakes, breaks = 6, labels = FALSE)

    p <- ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), y = ..density.., fill = as.factor(quantiles_response))) + 
      geom_histogram(bins = 6, position = "fill") +
      scale_fill_brewer(palette = "Spectral", name = paste("Quantiles of Num mistakes")) 

    plotL[parm] <- list(p)
    orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
  }
}
histInverted_totalmistakes <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(histInverted_totalmistakes, file = "/space/s1/fiona_callahan/histInverted_totalmistakes_stack.png", width = 40, height = 8)



plotL <- list()
orderImportances <- c()
for(parm in fmlaParms){
  if(class(multiSimRes_na.rm[[parm]]) != "character") {
    p <- ggplot(multiSimRes_na.rm, aes(x = !!sym(parm), y = ..density.., fill = totalMistakes.factor)) + 
      geom_histogram(bins = 6, position = "fill") +
      scale_fill_brewer(palette = "PRGn", name = paste("Num mistakes")) 

    plotL[parm] <- list(p)
    orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
  }
}
histInverted_totalmistakes <- grid.arrange(grobs = plotL[order(orderImportances, decreasing = TRUE)], nrow = 2)
ggsave(histInverted_totalmistakes, file = "/space/s1/fiona_callahan/histInverted_totalmistakes_stack2.png", width = 40, height = 8)











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
            allMeans = format(as.vector(unlist(multiSimResMeans[1, ])), scientific = FALSE), 
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
