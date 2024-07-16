library(ggplot2)
library(gridExtra)
library(stringr)

nSp <- 10
dir_list <- c(paste0("/space/s1/fiona_callahan/multiSim_", nSp, "sp/"), 
            paste0("/space/s1/fiona_callahan/multiSim_", nSp, "sp/"))
logiL <- c(TRUE, FALSE)
numSamples <- 10000
se_modeL <- seMode <- c("posOnly", "negOnly", "both")


numRuns <- 100
mode <- "ignore_direction" # "cluster" "ignore_sign" "ignore_direction" "cluster_cov"
ratio_of_avg <- FALSE

full_ROC <- data.frame()
simSetL <- c()
nSamplesL <- c()
seModeL <- c()

seSupp <- TRUE


for (i in seq_along(dir_list)) {
    thisDir <- dir_list[i]
    logi <- logiL[i]

    for (seMode in se_modeL) {
        if (ratio_of_avg) {
            saveName <- paste0("_", mode, "_ratioOfAv_", numSamples, "samples_", numRuns, "runs")
        } else {
            saveName <- paste0("_", mode, "_", numSamples, "samples_", numRuns, "runs")
        }


        if (str_detect(thisDir, "random")) {
            simSet <- "random parameters"
        } else {
            simSet <- "set parameters"
        }

        if (logi) {
            simSet <- "covariance matrix"
            saveName <- paste0(saveName, "_logi")
        }


        if (seMode == "posOnly") {
            saveName <- paste0(saveName, "_sePosOnly")
        } else if (seMode == "negOnly") {
            saveName <- paste0(saveName, "_seNegOnly")
        } else {
            saveName <- paste0(saveName, "_seBoth")
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
            seModeL <- c(seModeL, rep(seMode, times = dim(ROC_data)[1]))
        }
    }
}

full_ROC$sim_set <- simSetL
full_ROC$Num_Samples <- nSamplesL
full_ROC$se_mode <- seModeL

full_ROC$method[full_ROC$method == "spiecEasi_glasso"] <- "SpiecEasi_glasso"
full_ROC$method[full_ROC$method == "spiecEasi_mb"] <- "SpiecEasi_mb"

ROC_plot <- ggplot(full_ROC, aes(x = avg_FPR, y = avg_TPR, color = method, group = method)) +
        scale_shape_manual(values = 1:12) +
        #geom_point(size = 5, aes(shape = modelSelect)) +
        geom_point(data = subset(full_ROC, modelSelect == TRUE), size = 5) +  # Show points only when modelSelect = TRUE
        stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
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
                                "INLA_noCov" = "SDM-INLA (noCov)", 
                                "INLA_cov" = "SDM-INLA (cov)", 
                                "ecoCopula_cov" = "ecoCopula (cov)", 
                                "ecoCopula_noCov" = "ecoCopula (noCov)", 
                                "ecoCopula_cov_readAbd" = "ecoCopula (cov readAbd)", 
                                "ecoCopula_noCov_readAbd" = "ecoCopula (noCov readAbd)", 
                                "logistic_noCov" = "Logistic (noCov)",
                                "logistic_cov" = "Logistic (cov)",
                                "logistic" = "Logistic (cov)",
                                "linear_cov" = "Linear (cov)",
                                "linear_noCov" = "Linear (noCov)",
                                "JAGS" = "JSDM-MCMC"
                                )) + # Manual color scale for method
        theme(text = element_text(size = 18))  + # Set the base size for all text elements
        facet_grid(sim_set ~ se_mode)

ROC_plot

mean(full_ROC$FPR_sd[full_ROC$sim_set == "random parameters"], na.rm = TRUE)

mean(full_ROC$FPR_sd[full_ROC$sim_set == "set parameters"], na.rm = TRUE)


ggsave(paste0("/space/s1/fiona_callahan/facet_wrapped_plots_", nSp, "sp_", mode, "_seSupp.pdf"), ROC_plot)

