library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)

nSp <- 100
dir_list <- c(paste0("/space/s1/fiona_callahan/sim_paper_stuff/multiSim_", nSp, "sp/"), 
            paste0("/space/s1/fiona_callahan/sim_paper_stuff/multiSim_", nSp, "sp_random_moreSamples/"),
            paste0("/space/s1/fiona_callahan/sim_paper_stuff/multiSim_", nSp, "sp_revision2/"),
            paste0("/space/s1/fiona_callahan/sim_paper_stuff/multiSim_", nSp, "sp_revision3/"))
logiL <- c(FALSE, FALSE, TRUE, TRUE)
samplesL <- c(100, 250, 10000)

noGlasso <- TRUE
no_PgeN <- TRUE
logistic_error_delete <- TRUE
correct_method <- "bh" # "bh" or "bonferroni"
numRuns <- 100
mode <- "ignore_direction" # "cluster" "ignore_sign" "ignore_direction" "cluster_cov"
ratio_of_avg <- FALSE

full_ROC <- data.frame()
simSetL <- c()
nSamplesL <- c()

for (i in seq_along(dir_list)) {
    thisDir <- dir_list[i]
    logi <- logiL[i]

    for (numSamples in samplesL) {
        if (ratio_of_avg) {
            saveName <- paste0("_", mode, "_ratioOfAv_", numSamples, "samples_", numRuns, "runs")
        } else {
            saveName <- paste0("_", mode, "_", numSamples, "samples_", numRuns, "runs")
        }


        if (str_detect(thisDir, "random")) {
            simSet <- "ecological (random params)"
        } else {
            simSet <- "ecological (set params)"
        }

        if (str_detect(thisDir, "revision2")) {
            simSet <- "covariance mx (no cov effect)"
        }
        if (str_detect(thisDir, "revision3")) {
            simSet <- "covariance mx (+cov effect)"
        }

        if (logi) {
            #simSet <- "covariance matrix"
            saveName <- paste0(saveName, "_logi")
        }

        if (file.exists(paste0(thisDir, "ROC_data_", saveName, ".csv"))) {
            ROC_data <- read.csv(file = paste0(thisDir, "ROC_data_", saveName, ".csv"))
            file_info <- file.info(paste0(thisDir, "ROC_data_", saveName, ".csv"))
            print(thisDir)
            print(logi)
            print(numSamples)
            print(paste0(thisDir, "ROC_data_", saveName, ".csv"))
            print(file_info$mtime)
            full_ROC <- rbind(full_ROC, ROC_data)
            simSetL <- c(simSetL, rep(simSet, times = dim(ROC_data)[1]))
            nSamplesL <- c(nSamplesL, rep(numSamples, times = dim(ROC_data)[1]))
        }
    }
}

full_ROC$sim_set <- simSetL
full_ROC$Num_Samples <- nSamplesL

if (noGlasso) {
    if (!sum(is.na(full_ROC$method)) == 0) {
        print("WARNING: may be filtering because some of the full_ROC$method are NA")
    }
    full_ROC <- full_ROC[full_ROC$method != "spiecEasi_glasso" & full_ROC$method != "SpiecEasi_glasso", ]
}

if (no_PgeN) { # delete regression runs for 100 species and 100 samples bc p>=n is not allowed
    if (!sum(is.na(full_ROC$method)) == 0) {
        print("WARNING: may be filtering because some of the full_ROC$method are NA")
    }
    if (!sum(is.na(full_ROC$Num_Samples)) == 0) {
        print("WARNING: may be filtering because some of the full_ROC$Num_Samples are NA")
    }
    regressionNames <- c("logistic_noCov", "logistic", "linear_noCov", "linear_cov")
    # remember for each run of this file, nSp is constant
    deleteRegression <- nSp >= as.numeric(full_ROC$Num_Samples) 
    delete <- deleteRegression & (full_ROC$method %in% regressionNames)
    full_ROC <- full_ROC[!delete, ]
    if (logistic_error_delete) {
        ## see 2025-01-20 -- glm.fit: fitted probabilities numerically 0 or 1 occurred
        regressionNames <- c("logistic_noCov", "logistic")
        # remember for each run of this file, nSp is constant
        deleteRegression <- (nSp == 100) & (full_ROC$Num_Samples == 250) 
        delete <- deleteRegression & (full_ROC$method %in% regressionNames)
        full_ROC <- full_ROC[!delete, ]
    }
}

custom_labeller <- function(variable, value) {
  if (variable == "Num_Samples") {
    return(paste(value, "samples"))
  }
  return(value)
}

full_ROC$method[full_ROC$method == "spiecEasi_glasso"] <- "SpiecEasi_glasso"
full_ROC$method[full_ROC$method == "spiecEasi_mb"] <- "SpiecEasi_mb"

full_ROC$file[is.na(full_ROC$file)] <- "not recorded"

if (correct_method == "bonferroni") {
    # delete bh corrected lines (keep bonferroni)
    full_ROC <- full_ROC[!str_detect(full_ROC$file, "bh"), ]
} else if (correct_method == "bh") {
    # delete bonferroni corrected lines (keep bh)
    full_ROC <- full_ROC[!str_detect(full_ROC$file, "bonferroni"), ]
}

full_ROC[full_ROC$method == "ecoCopula_noCov", ]
full_ROC <- full_ROC %>%
  arrange(sim_set, Num_Samples, method, avg_TPR)

ROC_plot <- ggplot(full_ROC, aes(x = avg_FPR, y = avg_TPR, color = method, group = method)) +
        scale_shape_manual(values = 1:12) +
        #geom_point(size = 5, aes(shape = modelSelect)) +
        geom_point(data = subset(full_ROC, modelSelect == TRUE), size = 5) +  # Show points only when modelSelect = TRUE
        stat_summary(data = subset(full_ROC, modelSelect == FALSE), aes(group = method), fun.y = mean, geom = "line", size = 1) +
        labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +   # Add y=x line (no skill)
        #geom_errorbar(aes(ymin = pmax(0, avg_TPR - TPR_sd), ymax = pmin(1, avg_TPR + TPR_sd)), width = 0.03) +  # Add TPR error bars
        #geom_errorbarh(aes(xmin = pmax(0, avg_FPR - FPR_sd), xmax = pmin(1, avg_FPR + FPR_sd)), height = 0.03) +  # Add FPR error bars
        lims(x = c(0, 1), y = c(0, 1)) + # Set x and y-axis limits
        scale_color_manual(values = c("SpiecEasi_glasso" = "red", 
                                "SpiecEasi_mb" = "deeppink3", 
                                "spiecEasi_glasso" = "red", 
                                "spiecEasi_mb" = "deeppink3",
                                "sparcc" = "darkred",
                                "INLA_noCov" = "blue", 
                                "INLA_cov" = "deepskyblue", 
                                "ecoCopula_cov" = "green", 
                                "ecoCopula_noCov" = "#b1ff8d", 
                                "ecoCopula_cov_readAbd" = "darkgreen", 
                                "ecoCopula_noCov_readAbd" = "#00a200", 
                                "logistic_noCov" = "yellow",
                                "logistic_cov" = "brown1",
                                "logistic" = "darkorange",
                                "linear_cov" = "grey",
                                "linear_noCov" = "black",
                                "JAGS" = "darkorchid"),
                                labels = c(
                                "SpiecEasi_glasso" = "SpiecEasi (glasso)", 
                                "SpiecEasi_mb" = "SpiecEasi (mb)", 
                                "spiecEasi_glasso" = "SpiecEasi (glasso)", 
                                "spiecEasi_mb" = "SpiecEasi (mb)",
                                "sparcc" = "SPARCC",
                                "INLA_noCov" = "SDM-INLA (no cov)", 
                                "INLA_cov" = "SDM-INLA (with cov)", 
                                "ecoCopula_cov" = "ecoCopula (with cov, pres-abs)", 
                                "ecoCopula_noCov" = "ecoCopula (no cov, pres-abs)", 
                                "ecoCopula_cov_readAbd" = "ecoCopula (with cov, read abd)", 
                                "ecoCopula_noCov_readAbd" = "ecoCopula (no cov, read abd)", 
                                "logistic_noCov" = "Logistic (no cov)",
                                "logistic_cov" = "Logistic (with cov)",
                                "logistic" = "Logistic (with cov)",
                                "linear_cov" = "Linear (with cov)",
                                "linear_noCov" = "Linear (no cov)",
                                "JAGS" = "JSDM-MCMC"
                                )) + # Manual color scale for method
        theme(text = element_text(size = 18))  + # Set the base size for all text elements
        facet_grid(sim_set ~ Num_Samples, labeller = custom_labeller)

ROC_plot

mean(full_ROC$FPR_sd[full_ROC$sim_set == "random parameters"], na.rm = TRUE)

mean(full_ROC$FPR_sd[full_ROC$sim_set == "set parameters"], na.rm = TRUE)


ggsave(paste0("/space/s1/fiona_callahan/sim_paper_stuff/facet_wrapped_plots_", nSp, "sp_", mode, ".pdf"), ROC_plot, width = 16, height = 18)

