#!/usr/bin/env Rscript

# this is not working due to gifski or somethign -- need to configure something for gif rendering

# run with:
#Rscript plots_of_sims.R /home/fiona_callahan/simData/FALSE_POS2/fpr_mean_0.1/dependent_FPR/run1/sitetab_longer.csv /home/fiona_callahan/simData/FALSE_POS2/fpr_mean_0.1/dependent_FPR/run1/

library(ggplot2)
library(reshape)
library(gridExtra)
#install.packages(c('gapminder','gganimate','gifski'))
library(gapminder)
library(gganimate)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("sitetab and output files need to be supplied", call. = FALSE)
} 

sitetab_loc <- args[1]
sim_sitetab <- read.csv(sitetab_loc)
subdir <- args[2]

# plot covs
cov_animate <- ggplot(sim_sitetab, aes(Lat, Long, size = abs(Cov1), color = sign(Cov1))) +
  geom_point(alpha = 0.7, show.legend = TRUE) +
  #scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  # Animating the plot
  labs(title = "Year: {frame_time}", x = "Latitude", y = "Longitude") +
  transition_time(Age) +
  ease_aes("linear")

animate(cov_animate)

anim_save(paste0(subdir, "animated_cov1.gif"), width = 800, height = 800)

cov_animate <- ggplot(sim_sitetab, aes(Lat, Long, size = abs(Cov2), color = sign(Cov2))) +
  geom_point(alpha = 0.7, show.legend = TRUE) +
  # scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  # Animating the plot
  labs(title = "Year: {frame_time}", x = "Latitude", y = "Longitude") +
  transition_time(Age) +
  ease_aes("linear")

animate(cov_animate)

anim_save(paste0(subdir, "animated_cov2.gif"), width = 800, height = 800)
