library(ggplot2)
library(gridExtra)
library(stringr)

mode_list <- c("ignore_sign", "ignore_direction", "cluster")
samples_list <- c(100, 250, 10000)
correct_method <- "bh"

numRuns <- 100
ratio_of_avg <- FALSE

full_FDR <- data.frame()
full_DR <- data.frame()
modeL <- c()
nSamplesL <- c()

for (i in seq_along(mode_list)) {
    mode <- mode_list[i]

    for (numSamples in samples_list) {

        FDRs <- read.csv(file = paste0("/space/s1/fiona_callahan/sim_paper_stuff/FDR_stats_", numSamples, "samples_cutoff", correct_method, "_", mode, ".csv"))
        DRs <- read.csv(file = paste0("/space/s1/fiona_callahan/sim_paper_stuff/DR_stats_", numSamples, "samples_cutoff", correct_method, "_", mode, ".csv"))

        full_FDR <- rbind(full_FDR, FDRs)
        full_DR <- rbind(full_DR, DRs)

        modeL <- c(modeL, rep(mode, times = dim(FDRs)[1]))
        nSamplesL <- c(nSamplesL, rep(numSamples, times = dim(FDRs)[1]))
    }
}

nSamplesL <- sapply(nSamplesL, FUN = function(x) {
                                      paste(x, "samples")
                                      })
nSamplesL_factor <- factor(nSamplesL, levels = c("100 samples", "250 samples", "10000 samples"))

full_FDR$mistake_mode <- modeL
full_FDR$Num_Samples <- nSamplesL_factor

full_FDR$mistake_mode[full_FDR$mistake_mode == "cluster"] <- "indirect"
full_FDR$mistake_mode[full_FDR$mistake_mode == "ignore_direction"] <- "direct symmetric"
full_FDR$mistake_mode[full_FDR$mistake_mode == "ignore_sign"] <- "direct causal"

full_FDR$mistake_mode <- as.factor(full_FDR$mistake_mode)



full_DR$mistake_mode <- modeL
full_DR$Num_Samples <- nSamplesL_factor

full_DR$mistake_mode[full_DR$mistake_mode == "cluster"] <- "indirect"
full_DR$mistake_mode[full_DR$mistake_mode == "ignore_direction"] <- "direct symmetric"
full_DR$mistake_mode[full_DR$mistake_mode == "ignore_sign"] <- "direct causal"

full_DR$mistake_mode <- as.factor(full_DR$mistake_mode)


custom_labels_simulation <- c("multiSim_10sp" = "Ecological (Set, 10 species)",
                              "multiSim_10sp_logi" = "Covariance Mx alt, (10 species)",
                              "multiSim_100sp" = "Ecological (Set, 100 species)",
                              "multiSim_100sp_logi" = "Covariance Mx alt, 100 species",
                              "multiSim_10sp_random" = "Ecological (Random, 10 species)",
                              "multiSim_100sp_random" = "Ecological (Random, 100 species)",
                              "multiSim_10sp_revision2_logi" = "Covariance Mx (no cov effect, 10 species)",
                              "multiSim_100sp_revision2_logi" = "Covariance Mx (no cov effect, 100 species)",
                              "multiSim_10sp_revision3_logi" = "Covariance Mx (+cov effect, 10 species)",
                              "multiSim_100sp_revision3_logi" = "Covariance Mx (+cov effect, 100 species)")

if (correct_method == "bonferroni") {
  custom_labels_method <- c("ecoCopula_noCov" = "ecoCopula (no cov, pres-abs)",
                            "ecoCopula_cov" = "ecoCopula (with cov, pres-abs)",
                            "ecoCopula_noCov_readAbd" = "ecoCopula (no cov, read abd)",
                            "ecoCopula_cov_readAbd" = "ecoCopula (with cov, read abd)",
                            "spiecEasi_mb" = "spiecEasi (mb)",
                            "spiecEasi_glasso" = "spiecEasi (glasso)",      
                            "spiecEasi_sparcc" = "SPARCC",
                            "INLA" = "SDM-INLA",              
                            "logistic_cov" = "logistic (with cov, Bonferroni adj)",          
                            "logistic_noCov" = "logistic (no cov, Bonferroni adj)",
                            "linearReg_cov" = "linear (with cov, Bonferroni adj)",
                            "linearReg_noCov" = "linear (no cov, Bonferroni adj)",     
                            "JAGS" = "JSDM-MCMC")
} else if (correct_method == "bh") {
  custom_labels_method <- c("ecoCopula_noCov" = "ecoCopula (no cov, pres-abs)",
                            "ecoCopula_cov" = "ecoCopula (with cov, pres-abs)",
                            "ecoCopula_noCov_readAbd" = "ecoCopula (no cov, read abd)",
                            "ecoCopula_cov_readAbd" = "ecoCopula (with cov, read abd)",
                            "spiecEasi_mb" = "spiecEasi (mb)",
                            "spiecEasi_glasso" = "spiecEasi (glasso)",      
                            "spiecEasi_sparcc" = "SPARCC",
                            "INLA" = "SDM-INLA",              
                            "logistic_cov" = "logistic (with cov, BH adjusted)",          
                            "logistic_noCov" = "logistic (no cov, BH adjusted)",
                            "linearReg_cov" = "linear (with cov, BH adjusted)",
                            "linearReg_noCov" = "linear (no cov, BH adjusted)",     
                            "JAGS" = "JSDM-MCMC")
} else {
    custom_labels_method <- c("ecoCopula_noCov" = "ecoCopula (no cov, pres-abs)",
                            "ecoCopula_cov" = "ecoCopula (with cov, pres-abs)",
                            "ecoCopula_noCov_readAbd" = "ecoCopula (no cov, read abd)",
                            "ecoCopula_cov_readAbd" = "ecoCopula (with cov, read abd)",
                            "spiecEasi_mb" = "spiecEasi (mb)",
                            "spiecEasi_glasso" = "spiecEasi (glasso)",      
                            "spiecEasi_sparcc" = "SPARCC",
                            "INLA" = "SDM-INLA",              
                            "logistic_cov" = "logistic (with cov)",          
                            "logistic_noCov" = "logistic (no cov)",
                            "linearReg_cov" = "linear (with cov)",
                            "linearReg_noCov" = "linear (no cov)",     
                            "JAGS" = "JSDM-MCMC")
}

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

#print(DRplot)

ggsave(plot = FDRplot, file = paste0("/space/s1/fiona_callahan/sim_paper_stuff/FDR_faceted_3_", correct_method, ".pdf"), width = 20, height = 12)

ggsave(plot = DRplot, file = paste0("/space/s1/fiona_callahan/sim_paper_stuff/DR_faceted_3_", correct_method, ".pdf"), width = 24, height = 12)
