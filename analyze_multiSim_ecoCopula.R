library(ranger)
library(dplyr)
library(pdp)
library(ggplot2)
library(plyr)
library(gridExtra)
library(tidyr)

dirName <- c("multiSim_100")
multiSimRes <- data.frame()
resNames <- c("ecoCopula_res_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "INLA_infResGathered.csv", "logistic_mistakes.csv")
resNames <- c("spiecEasi_res_glasso_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "spiecEasi_res_sparcc_infResGathered.csv", "ecoCopula_res_noCov_infResGathered.csv")

#ls /space/s1/fiona_callahan/multiSim_100
resNames <- c("ecoCopula_res_noCov_infResGathered.csv", "spiecEasi_res_mb_infResGathered.csv", "logistic_mistakes.csv", "INLA_infResGathered.csv")

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
    if(is.na(thisMultiSimRes$fpr.mode[row])){
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


multiSimRes <- multiSimResL$spiecEasi_res_mb_infResGathered.csv
multiSimRes <- multiSimResL$ecoCopula_res_infResGathered.csv
multiSimRes <- multiSimResL$INLA_infResGathered.csv
multiSimRes <- multiSimResL$logistic_mistakes.csv

multiSimRes <- multiSimResL$ecoCopula_res_noCov_infResGathered.csv
multiSimRes <- multiSimResL$spiecEasi_res_mb_infResGathered.csv

mean(multiSimRes$alpha_incorrect_undirected + multiSimRes$alpha_correct_undirected)
mean(multiSimRes$num_incorrectInferences + multiSimRes$num_correctInferences, na.rm = TRUE)
mean(multiSimRes$num_incorrect_cluster + multiSimRes$num_correct_cluster)

summary(multiSimRes$num_incorrect_alpha)
summary(multiSimRes$num_correctInferences)

mean(multiSimRes$fp_fp_tp, na.rm = TRUE)
mean(multiSimRes$alpha_incorrect_undirected / (multiSimRes$alpha_incorrect_undirected + multiSimRes$alpha_correct_undirected), na.rm = TRUE)
mean(multiSimRes$num_incorrect_cluster / (multiSimRes$num_incorrect_cluster + multiSimRes$num_correct_cluster), na.rm = TRUE)





###############################





#for (multiSimRes in multiSimResL) {
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

  emergentParms <- c("Sp1PercentPresence", "Sp2PercentPresence", "Sp3PercentPresence", "Sp4PercentPresence", "Sp5PercentPresence")
  # make formula from these parms
  indep_var <- "totalMistakes"
  #indep_var <- "fp_fp_tp"
  #indep_var <- "finished_INLA"
  fmlaParms <- c(randomParms)
  #fmlaParms <- c(randomParms[randomParms!="r"])

  fmla <- as.formula(paste(indep_var, " ~ ", paste(fmlaParms, collapse = "+")))

  #remove NA's in response (run didnt finish or whatever)
  multiSimRes_na.rm <- multiSimRes[!is.na(multiSimRes[[indep_var]]), ]

  # inla average success
  mean(multiSimRes_na.rm$fp_fp_tp, na.rm = TRUE)
  quantile(multiSimRes_na.rm$fp_fp_tp, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  quantile(multiSimRes_na.rm$totalMistakes, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  mean(multiSimRes_na.rm$num_correctInferences + multiSimRes_na.rm$num_incorrectInferences, na.rm = TRUE)
  mean(multiSimRes_na.rm$num_incorrectInferences, na.rm = TRUE)

  sum(multiSimRes$num_incorrectInferences, na.rm = TRUE) / 
      (sum(multiSimRes$num_correctInferences, na.rm = TRUE) + sum(multiSimRes$num_incorrectInferences, na.rm = TRUE))

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

  #ggsave(lm_betaPlot, file = "/space/s1/fiona_callahan/lm_betaPlot.png", height = 12, width = 12)

  # random forest res
  #rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "impurity") 
  rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "permutation") 
  summary(rf_res)
  sort(importance(rf_res))
  importanceDF <- as.data.frame(sort(importance(rf_res)))
  importanceDF$Variable <- row.names(importanceDF)
  names(importanceDF) <- c("RF_Importance", "Parameter")


  #better_parmNames <- data.frame(Parameter = fmlaParms, Parameter_name = c("n_time", "n_space", 
  #                    "neighborhood radius", "false detection rate", "false detect correlation", 
  #                    "covariate process noise", "covariate measurement noise", "r", "sigma", 
  #                    "false absence parameter", "migration rate", "inter-species affect scaling"))#, 
  #                    #"Sp1 Percent Presence", "Sp2 Percent Presence", "Sp3 Percent Presence"))

  importanceDF <- merge(x = importanceDF, y = better_parmNames, by = "Parameter")

  rf_importance_plot <- ggplot(importanceDF, aes(x = reorder(Parameter, RF_Importance), y = RF_Importance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme(axis.text.y = element_text(hjust = 1, size = 24), axis.title.x = element_text(size = 24), axis.title.y = element_blank()) +
    coord_flip() +
    labs(y = "Random Forest Importance", title = paste("Random Forest Predictor Importance for", indep_var)) 
  #ggsave(rf_importance_plot, file = "/space/s1/fiona_callahan/rf_importance.png", height = 10, width = 10)

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
  #ggsave(rf_marginal_ice, file = "/space/s1/fiona_callahan/rf_marginal_ice.png", width = 40, height = 12)

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
  #ggsave(rf_marginal_cice, file = "/space/s1/fiona_callahan/rf_marginal_cice.png", width = 40, height = 12)


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
  #ggsave(rf_marginal, file = "/space/s1/fiona_callahan/rf_marginal.png", width = 40, height = 12)

#}
