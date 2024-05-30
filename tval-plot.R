
library(ggplot2)

data_dir <- "/space/s1/fiona_callahan/multiSim_100sp/"
numRuns <- as.numeric(100)
outName <- "linearReg_mistakes_sampled100_noCov_filtered_100runs"

t_vals_df <- read.csv(paste0(data_dir, "tValsDistribution_", outName, "_", numRuns, "sims.csv"))
t_vals_df$actual_beta <- as.factor(t_vals_df$actual_beta)

plot1 <- ggplot(t_vals_df, aes(x = tValues, fill = actual_beta, color = actual_beta)) +
    geom_density(alpha = 0.5) +  # Add transparency for better visibility of overlapping densities
    #coord_cartesian(xlim = c(-1, 1)) +  # Set x-axis limits
    labs(title = paste0("Density Plot of tValues: ", outName),
        x = "tValues",
        y = "Density") +
    scale_fill_manual(values = rainbow(length(unique(t_vals_df$actual_beta)))) +  # Adjust colors
    scale_color_manual(values = rainbow(length(unique(t_vals_df$actual_beta)))) +  # Adjust colors
    theme_minimal()  # Optional: Change the theme
    
plot1
#ggsave(filename = paste0(data_dir, "tValsDistribution_", outName, "_", numRuns, "sims.png"), plot = plot1)

summary(t_vals_df)
