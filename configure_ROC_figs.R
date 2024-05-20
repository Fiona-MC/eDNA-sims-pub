library(ggplot2)
library(gridExtra)
library(stringr)

dir_list <- c("/space/s1/fiona_callahan/multiSim_10sp/", 
            "/space/s1/fiona_callahan/multiSim_10sp/", 
            "/space/s1/fiona_callahan/multiSim_10sp_random_moreSamples/")
logiL <- c(TRUE, FALSE, FALSE)
samplesL <- c(100, 1000, 10000)

numRuns <- 100

full_ROC <- data.frame()
simSetL <- c()
nSamplesL <- c()

for (i in seq_along(dir_list)) {
    thisDir <- dir_list[i]
    logi <- logiL[i]

    for (numSamples in samplesL) {
        saveName <- paste0("_", numSamples, "samples_", numRuns, "runs")

        if (str_detect(thisDir, "random")) {
            simSet <- "random parameters"
        } else {
            simSet <- "set parameters"
        }

        if (logi) {
            simSet <- paste0(simSet, " logi")
            saveName <- paste0(saveName, "_logi")
        }

        if (file.exists(paste0(thisDir, "ROC_data_", saveName, ".csv"))) {
            ROC_data <- read.csv(file = paste0(thisDir, "ROC_data_", saveName, ".csv"))
            full_ROC <- rbind(full_ROC, ROC_data)
            simSetL <- c(simSetL, rep(simSet, times = dim(ROC_data)[1]))
            nSamplesL <- c(nSamplesL, rep(numSamples, times = dim(ROC_data)[1]))
        }
    }
}

full_ROC$sim_set <- simSetL
full_ROC$Num_Samples <- nSamplesL

ROC_plot <- ggplot(full_ROC, aes(x = avg_FPR, y = avg_TPR, color = method, group = method)) +
        scale_shape_manual(values = 1:12) +
        #geom_point(size = 5, aes(shape = modelSelect)) +
        geom_point(data = subset(full_ROC, modelSelect == TRUE), size = 5) +  # Show points only when modelSelect = TRUE
        stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
        labs(x = "FPR", y = "TPR") +
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
                                        "ecoCopula_noCov" = "darkgreen", 
                                        "logistic_noCov" = "yellow",
                                        "logistic_cov" = "brown1",
                                        "logistic" = "darkorange",
                                        "JAGS" = "darkorchid")) + # Manual color scale for method
        theme(text = element_text(size = 18))  + # Set the base size for all text elements
        facet_grid(sim_set ~ Num_Samples)

ROC_plot

mean(full_ROC$FPR_sd[full_ROC$sim_set == "random parameters"], na.rm = TRUE)

mean(full_ROC$FPR_sd[full_ROC$sim_set == "set parameters"], na.rm = TRUE)


ggsave("/space/s1/fiona_callahan/facet_wrapped_plots_100sp.png", ROC_plot)
