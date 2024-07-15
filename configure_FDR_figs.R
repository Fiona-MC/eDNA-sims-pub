library(ggplot2)
library(gridExtra)
library(stringr)

mode_list <- c("ignore_sign", "ignore_direction", "cluster")
samples_list <- c(100, 10000)

numRuns <- 100
ratio_of_avg <- TRUE

full_FDR <- data.frame()
full_DR <- data.frame()
modeL <- c()
nSamplesL <- c()

for (i in seq_along(mode_list)) {
    mode <- mode_list[i]

    for (numSamples in samples_list) {

        FDRs <- read.csv(file = paste0("/space/s1/fiona_callahan/FDR_stats_", numSamples, "samples_cutoffbonferroni_", mode, ".csv"))
        DRs <- read.csv(file = paste0("/space/s1/fiona_callahan/DR_stats_", numSamples, "samples_cutoffbonferroni_", mode, ".csv"))

        full_FDR <- rbind(full_FDR, FDRs)
        full_DR <- rbind(full_DR, DRs)

        modeL <- c(modeL, rep(mode, times = dim(FDRs)[1]))
        nSamplesL <- c(nSamplesL, rep(numSamples, times = dim(FDRs)[1]))
    }
}

nSamplesL <- sapply(nSamplesL, FUN = function(x) {paste(x, "samples")})
full_FDR$mistake_mode <- modeL
full_FDR$Num_Samples <- as.factor(nSamplesL)

full_FDR$mistake_mode[full_FDR$mistake_mode == "cluster"] <- "indirect"
full_FDR$mistake_mode[full_FDR$mistake_mode == "ignore_sign"] <- "direct"

full_FDR$mistake_mode <- as.factor(full_FDR$mistake_mode)



full_DR$mistake_mode <- modeL
full_DR$Num_Samples <- as.factor(nSamplesL)

full_DR$mistake_mode[full_DR$mistake_mode == "cluster"] <- "indirect"
full_DR$mistake_mode[full_DR$mistake_mode == "ignore_sign"] <- "direct"

full_DR$mistake_mode <- as.factor(full_DR$mistake_mode)


custom_labels_simulation <- c("multiSim_10sp" = "Set parms, 10sp",
                              "multiSim_10sp_logi" = "Covariance Mx, 10sp",
                              "multiSim_100sp" = "Set parms, 100sp",
                              "multiSim_100sp_logi" = "Covariance Mx, 100sp",
                              "multiSim_10sp_random" = "Random parms, 10sp",
                              "multiSim_100sp_random" = "Random parms, 100sp")

custom_labels_method <- c("ecoCopula_noCov" = "ecoCopula (noCov, pres-abs)",
                          "ecoCopula_cov" = "ecoCopula (cov, pres-abs)",
                          "ecoCopula_noCov_readAbd" = "ecoCopula (noCov, readAbd)",
                          "ecoCopula_cov_readAbd" = "ecoCopula (cov, readAbd)",
                          "spiecEasi_mb" = "spiecEasi (mb)",
                          "spiecEasi_glasso" = "spiecEasi (glasso)",      
                          "spiecEasi_sparcc" = "SPARCC",
                          "INLA" = "SDM-INLA",              
                          "logistic_cov" = "logistic (cov)",          
                          "logistic_noCov" = "logistic (noCov)",
                          "linearReg_cov" = "linear (cov)",
                          "linearReg_noCov" = "linear (noCov)",     
                          "JAGS" = "JSDM-MCMC")



# Create the heatmap
# Create the heatmap
FDRplot <- ggplot(full_FDR, aes(x = Simulation, y = Method, fill = FDR)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", na.value = "gray", limits = c(0, 1)) + # Adjust colors as needed
  geom_text(aes(label = sprintf("%.2f", FDR)), color = "black", size = 5) + # Display FDR values as text
  theme_minimal() + # Change the theme as per your preference
  labs(x = "Simulation",
       y = "Method",
       fill = "FDR") +
  scale_x_discrete(labels = custom_labels_simulation) +  # Custom labels for the x-axis
  scale_y_discrete(labels = custom_labels_method) +      # Custom labels for the y-axis
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 14, vjust = 1),  # Adjust the size value as needed
        axis.text.y = element_text(size = 14),  # Adjust the size value as needed
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10)) +   # Adjust the size value as needed
  facet_grid(mistake_mode ~ Num_Samples)  # Facet by mode and Num_Samples

print(FDRplot)


# Create the heatmap
# Create the heatmap
DRplot <- ggplot(full_DR, aes(x = Simulation, y = Method, fill = totalDiscoveries)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", na.value = "gray") + # Adjust colors as needed
  geom_text(aes(label = sprintf("%.2f", totalDiscoveries)), color = "black") + # Display FDR values as text
  theme_minimal() + # Change the theme as per your preference
  labs(x = "Simulation",
       y = "Method",
       fill = "Average discoveries per run") +
  scale_x_discrete(labels = custom_labels_simulation) +  # Custom labels for the x-axis
  scale_y_discrete(labels = custom_labels_method) +      # Custom labels for the y-axis
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 14, vjust = 1),  # Adjust the size value as needed
        axis.text.y = element_text(size = 14),  # Adjust the size value as needed
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10)) +   # Adjust the size value as needed
  facet_grid(mistake_mode ~ Num_Samples)  # Facet by mode and Num_Samples

print(DRplot)
#ggsave("/space/s1/fiona_callahan/FDR_faceted_3.pdf")

ggsave("/space/s1/fiona_callahan/DR_faceted_3.pdf", width = 16, height = 12)
