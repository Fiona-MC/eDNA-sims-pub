library(ranger)
library(dplyr)
library(pdp)
library(ggplot2)
library(plyr)
library(gridExtra)
library(tidyr)
library(stringr)
library(latex2exp)

setwd("~/eDNA_sims_code")

dirNames <- c("multiSim_10sp_random_1000")

#resNames <- c("logistic_mistakes_sampled10000_noCov_filtered_100runs_cutoff0.05_100sims.csv",
            #"logistic_mistakes_sampled10000_covNoCount_filtered_100runs_cutoff0.05_100sims.csv",
#            "spiecEasi_res_sampled10000_mb_filtered_infResGathered_100sims.csv",
#            "sparcc_res_sampled10000_filtered_infResGathered_cutoff0.3_100sims.csv",
#            "ecoCopula_res_sampled10000_noCov_filtered_infResGathered_100sims.csv",
            #"INLA_res_paperSep_sampled10000_filtered_noCov_infResGathered_cutoff0.05_100sims.csv",
#            "INLA_res_paperSep_sampled10000_filtered_covNoCount_infResGathered_cutoff0.05_100sims.csv")

numSamples <- 10000
bonferroni <- TRUE
mode <- "ignore_direction"
if (!bonferroni) {
  resNames <- c(paste0("logistic_mistakes_sampled", numSamples, "_noCov_filtered_1000runs_cutoff0.05_1000sims.csv"), 
              paste0("linearReg_mistakes_sampled", numSamples, "_noCov_filtered_1000runs_cutoff0.05_1000sims.csv"))
} else {
  resNames <- c(paste0("logistic_mistakes_sampled", numSamples, "_noCov_filtered_1000runs_cutoffbonferroni_1000sims.csv"), 
              paste0("linearReg_mistakes_sampled", numSamples, "_noCov_filtered_1000runs_cutoffbonferroni_1000sims.csv"))
}

source("./confusion_stats.R")

resName <- resNames[2]
plotL <- list()
for(resName in resNames) {
    method <- str_split(resName, pattern = "_")[[1]][1]

    if (str_detect(resName, "cov")){
        method <- paste0(method, "_cov")
    }

    multiSimRes <- data.frame()

    dirNameL <- c()
    for (dirName in dirNames) {
    thisDir <- paste0("/space/s1/fiona_callahan/", dirName, "/")
    thisMultiSimRes <- read.csv(paste0(thisDir, resName), header = TRUE)
    dirNameL <- c(dirNameL, rep(dirName, times = length(thisMultiSimRes$sim_run)))
    multiSimRes <- plyr::rbind.fill(multiSimRes, thisMultiSimRes)
    }
    multiSimRes$dirName <- dirNameL

    filter_criteria <- !str_detect(names(multiSimRes), "alpha") &
                        !str_detect(names(multiSimRes), "beta") &
                        !str_detect(names(multiSimRes), "names_species")
    multiSimRes <- multiSimRes[, filter_criteria]

    randomParms <- c("radius", "covMeasureNoise_sd", "r", "sigma", 
                    "c2", "mean_mig_rate", "indivSampleProb", "readSampleRate", 
                    "numCovs", "covar_scale_space", "readThreshold")
    randomParmsVary <- as.vector(sapply(lapply(multiSimRes[, randomParms], unique), length) > 1)
    randomParms <- randomParms[as.vector(randomParmsVary)] # take out any that I've changed to not be random

    multiSimRes$FDR <- get_falseDiscovery(data = multiSimRes, mode = mode, return_components = FALSE)
    multiSimRes$FPR <- get_FPR(data = multiSimRes, mode = mode, return_components = FALSE)
    multiSimRes$TPR <- get_TPR(data = multiSimRes, mode = mode, return_components = FALSE)

    indep_var <- "FDR"

    fmlaParms <- c(randomParms)
    #fmlaParms <- c(randomParms[randomParms!="r"])

    fmla <- as.formula(paste(indep_var, " ~ ", paste(fmlaParms, collapse = "+")))

    #remove NA's in response (run didnt finish or whatever)
    multiSimRes_na.rm <- multiSimRes[!is.na(multiSimRes[[indep_var]]), ]

    # random forest res
    #rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "impurity") 
    rf_res <- ranger(data = multiSimRes_na.rm, formula = fmla, importance = "permutation") 
    summary(rf_res)
    sort(importance(rf_res))
    importanceDF <- as.data.frame(sort(importance(rf_res)))
    importanceDF$Variable <- row.names(importanceDF)
    names(importanceDF) <- c("RF_Importance", "Parameter")

    # Create a named vector for LaTeX mappings
  latex_mappings <- c(
    "radius" = "$d_{mig}$",
    "covMeasureNoise_sd" = "$\\sigma_{cov}$",
    "r" = "$r$",
    "sigma" = "$\\sigma$",
    "c2" = "$c_2$",
    "mean_mig_rate" = "$mean\\_mig\\_rate$",
    "indivSampleProb" = "$pSample$",
    "readSampleRate" = "$lambdaReads$",
    "numCovs" = "$N_{cov}$",
    "covar_scale_space" = "$covar\\_scale\\_space$",
    "readThreshold" = "$R_{detect}$"
  )

  importanceDF$Parameter_TeX <- sapply(importanceDF$Parameter, function(x) latex_mappings[x])
    #better_parmNames <- data.frame(Parameter = fmlaParms, Parameter_name = c("n_time", "n_space", 
    #                    "neighborhood radius", "false detection rate", "false detect correlation", 
    #                    "covariate process noise", "covariate measurement noise", "r", "sigma", 
    #                    "false absence parameter", "migration rate", "inter-species affect scaling"))#, 
                        #"Sp1 Percent Presence", "Sp2 Percent Presence", "Sp3 Percent Presence"))

    #importanceDF <- merge(x = importanceDF, y = better_parmNames, by = "Parameter")


    rf_importance_plot <- ggplot(importanceDF, aes(x = reorder(Parameter_TeX, RF_Importance), y = RF_Importance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme(axis.text.y = element_text(hjust = 1, size = 24), axis.title.x = element_text(size = 24), 
              axis.title.y = element_blank(), plot.margin = unit(c(1, 1, 1, 2), "cm")) +
      coord_flip() +
      labs(y = "Random Forest Importance", title = paste(indep_var, "for", method)) +
      scale_x_discrete(labels = function(x) TeX(x))

    
    #rf_importance_plot
    plotL[[method]] <- rf_importance_plot

    # how well does the rf predict FDR (training set)
    data <- multiSimRes_na.rm
    naive_rmse <- sqrt(sum((median(data[[indep_var]]) - data[[indep_var]])^2, na.rm = TRUE) / length(data[[indep_var]]))
    rf_predictions <- predict(object = rf_res, data = data)
    rf_rmse <- sqrt(sum((rf_predictions$predictions - data[[indep_var]])^2) / length(data[[indep_var]]))
    (rf_rmse - naive_rmse) / naive_rmse # training


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
        coord_cartesian() +
        #coord_cartesian(ylim = c(.2, .6)) +
        labs(x = TeX(latex_mappings[parm]), y = paste0("Predicted ", indep_var)) +
        theme(axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 18))

    ICEplotL[parm] <- list(p)
    orderImportances <- c(orderImportances, importanceDF$RF_Importance[importanceDF$Parameter == parm])
  }
  rf_marginal_ice <- grid.arrange(grobs = ICEplotL[order(orderImportances, decreasing = TRUE)])

  if (bonferroni) {
    ggsave(filename = paste0("/space/s1/fiona_callahan/", dirName, "/RF_ice_", method, "_bonferroni_sampled", numSamples, ".png"), 
                plot = rf_marginal_ice, width = 16, height = 16)
  } else {
    ggsave(filename = paste0("/space/s1/fiona_callahan/", dirName, "/RF_ice_", method, "_cutoff0.05_sampled", numSamples, ".png"), 
                plot = rf_marginal_ice, width = 16, height = 16)
  }

}


if(bonferroni) {
ggsave(filename = paste0("/space/s1/fiona_callahan/", dirName, "/RF_importance_allmethods_bonferroni_sampled", numSamples, ".pdf"), 
                plot = grid.arrange(grobs = plotL), width = 16, height = 16)
} else {
ggsave(filename = paste0("/space/s1/fiona_callahan/", dirName, "/RF_importance_allmethods_cutoff0.05_sampled", numSamples, ".pdf"), 
                plot = grid.arrange(grobs = plotL), width = 16, height = 16)
}



###### test set ########
dirNames <- c("multiSim_10sp_random_testSet")
if (bonferroni) {
  resName <- paste0("logistic_mistakes_sampled", numSamples, "_noCov_filtered_100runs_cutoffbonferroni_100sims.csv")
  resName <- paste0("linearReg_mistakes_sampled", numSamples, "_noCov_filtered_100runs_cutoffbonferroni_100sims.csv")

} else {
  resName <- paste0("logistic_mistakes_sampled", numSamples, "_noCov_filtered_100runs_cutoff0.05_100sims.csv")
  resName <- paste0("linearReg_mistakes_sampled", numSamples, "_noCov_filtered_100runs_cutoff0.05_100sims.csv")
}

method <- str_split(resName, pattern = "_")[[1]][1]

if (str_detect(resName, "cov")) {
    method <- paste0(method, "_cov")
}

multiSimRes_test <- data.frame()

dirNameL <- c()
for (dirName in dirNames) {
thisDir <- paste0("/space/s1/fiona_callahan/", dirName, "/")
thisMultiSimRes <- read.csv(paste0(thisDir, resName), header = TRUE)
dirNameL <- c(dirNameL, rep(dirName, times = length(thisMultiSimRes$sim_run)))
multiSimRes_test <- plyr::rbind.fill(multiSimRes_test, thisMultiSimRes)
}
multiSimRes_test$dirName <- dirNameL

filter_criteria <- !str_detect(names(multiSimRes_test), "alpha") &
                    !str_detect(names(multiSimRes_test), "beta") &
                    !str_detect(names(multiSimRes_test), "names_species")
multiSimRes_test <- multiSimRes_test[, filter_criteria]


randomParms <- c("radius", "covNoise_sd", "covMeasureNoise_sd", "r", "sigma", 
                "c2", "mean_mig_rate", "indivSampleProb", "readSampleRate", 
                "numCovs", "covar_scale_space", "readThreshold")
randomParmsVary <- as.vector(sapply(lapply(multiSimRes_test[, randomParms], unique), length) > 1)
randomParms <- randomParms[as.vector(randomParmsVary)] # take out any that I've changed to not be random

multiSimRes_test$FDR <- get_falseDiscovery(data = multiSimRes_test, mode = "ignore_sign", return_components = FALSE)
multiSimRes_test$FPR <- get_FPR(data = multiSimRes_test, mode = "ignore_sign", return_components = FALSE)
multiSimRes_test$TPR <- get_TPR(data = multiSimRes_test, mode = "ignore_sign", return_components = FALSE)

indep_var <- "FDR"

fmlaParms <- c(randomParms)
#fmlaParms <- c(randomParms[randomParms!="r"])

fmla <- as.formula(paste(indep_var, " ~ ", paste(fmlaParms, collapse = "+")))

#remove NA's in response (run didnt finish or whatever)
multiSimRes_test_na.rm <- multiSimRes_test[!is.na(multiSimRes_test[[indep_var]]), ]

rf_predictions <- predict(object = rf_res, data = multiSimRes_test_na.rm)
rf_rmse <- sqrt(sum((rf_predictions$predictions - multiSimRes_test_na.rm$FDR)^2) / length(multiSimRes_test_na.rm$FDR))
naive_rmse <- sqrt(sum((median(multiSimRes_test_na.rm[[indep_var]]) - multiSimRes_test_na.rm[[indep_var]])^2, na.rm = TRUE) / 
              length(multiSimRes_test_na.rm[[indep_var]]))


#for test set
resName
rf_rmse
naive_rmse
(rf_rmse - naive_rmse) / naive_rmse



 

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




