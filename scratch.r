


plotZvals <- TRUE
if (plotZvals) {
    library(ggplot2)
    z.values_0 <- z.values_0[!is.na(z.values_0)]
    z.values_pos <- z.values_pos[!is.na(z.values_pos)]
    z.values_neg <- z.values_neg[!is.na(z.values_neg)]
    z_vals_df <- data.frame(actual_beta = as.factor(c(rep(0, times = length(z.values_0)), 
                                        rep(1, times = length(z.values_pos)), 
                                        rep(-1, times = length(z.values_neg)))), 
                        zValues = c(z.values_0, z.values_pos, z.values_neg))

    normal_density <- function(x) dnorm(x, mean = 0, sd = 1)

    plot1 <- ggplot(z_vals_df, aes(x = zValues, fill = actual_beta, color = actual_beta)) +
        geom_density(alpha = 0.5) +  # Add transparency for better visibility of overlapping densities
        #coord_cartesian(xlim = c(-1, 1)) +  # Set x-axis limits
        labs(title = paste0("Density Plot of zValues: ", outName),
            x = "zValues",
            y = "Density") +
        #stat_function(fun = normal_density, aes(color = "Standard Normal"), size = 1) +  # Add standard normal density
        scale_fill_manual(values = rainbow(length(unique(z_vals_df$actual_beta)))) +  # Adjust colors
        scale_color_manual(values = rainbow(length(unique(z_vals_df$actual_beta)))) +  # Adjust colors
        theme_minimal()  # Optional: Change the theme
    
    ggsave(filename = paste0(data_dir, "zValsDistribution_", outName, "_", numRuns, "sims.png"), plot = plot1)

}
