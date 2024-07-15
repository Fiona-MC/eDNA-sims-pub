# Load necessary packages
library(pwr)
library(ggplot2)

# Function to calculate sample size for linear regression
calculate_sample_size <- function(alpha, rho, power = 0.8) {
  # Calculate required sample size
  result <- pwr.r.test(r = rho, sig.level = alpha, power = power)
  return(result$n)
}

# Define parameters
num_species_values <- c(10, 100)
alpha_values <- 0.05 / num_species_values^2
power <- 0.8  # Desired power
rho_values <- seq(0.1, 0.9, by = 0.1)  # Range of correlation coefficients

# Calculate sample sizes for both num_species values
sample_sizes <- data.frame()
for (i in seq_along(num_species_values)) {
  num_species <- num_species_values[i]
  alpha <- alpha_values[i]
  for (rho in rho_values) {
    n <- calculate_sample_size(alpha, rho, power)
    sample_sizes <- rbind(sample_sizes, data.frame(num_species = num_species, rho = rho, sample_size = n))
  }
}

# Plotting with ggplot2
ggplot(sample_sizes, aes(x = rho, y = sample_size, color = as.factor(num_species))) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray") +
  labs(x = "Correlation",
       y = "Required Sample Size for power = 0.8",
       color = "Number of Species") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave("/space/s1/fiona_callahan/sample_size_versus_correlation.pdf")


















################
## EMPIRICAL CORRELATIONS FOR DATA SETS
###############

library(SpiecEasi)
data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
#> dim(amgut1.filt.n)
#[1] 289 127 (samples, taxa)
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

empirical_cor_amgut1 <- c()
for (sp1 in 1:dim(amgut1.filt)[2]) { # nolint: seq_linter.
  for (sp2 in 1:dim(amgut1.filt)[2]) { # nolint: seq_linter.
    empirical_cor_amgut1 <- c(empirical_cor_amgut1, cor(amgut1.filt[, sp1], amgut1.filt[, sp2]))
  }
}
hist(empirical_cor_amgut1)

amgut1.filt.pa <- (amgut1.filt > 5) * 1
empirical_cor_amgut1.pa <- c()
for (sp1 in 1:dim(amgut1.filt.pa)[2]) { # nolint: seq_linter.
  for (sp2 in 1:dim(amgut1.filt.pa)[2]) { # nolint: seq_linter.
    empirical_cor_amgut1.pa <- c(empirical_cor_amgut1.pa, cor(amgut1.filt.pa[, sp1], amgut1.filt.pa[, sp2]))
  }
}
hist(empirical_cor_amgut1.pa)

graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

empirical_cor_sim <- c()
for (sp1 in 1:dim(X)[2]) { # nolint: seq_linter.
  for (sp2 in 1:dim(X)[2]) { # nolint: seq_linter.
    empirical_cor_sim <- c(empirical_cor_sim, cor(X[, sp1], X[, sp2]))
  }
}
hist(empirical_cor_sim)


library(ecoCopula)
data(spider)
#> dim(spider$abund)
#[1] 28 12

empirical_cor_spider <- c()
for (sp1 in 1:dim(spider$abund)[2]) { # nolint: seq_linter.
  for (sp2 in 1:dim(spider$abund)[2]) { # nolint: seq_linter.
    empirical_cor_spider <- c(empirical_cor_spider, cor(spider$abund[, sp1], spider$abund[, sp2]))
  }
}
hist(empirical_cor_spider)


arctic_sitetab <- read.table("/space/s1/fiona_callahan/Arctic_original/sitetab.tsv", sep = "\t", header = TRUE)
names_animals <- names(arctic_sitetab)[12:25]
arctic_sitetab_animals <- arctic_sitetab[, names_animals]
arctic_sitetab_animals <- na.omit(arctic_sitetab_animals)

empirical_cor_arctic <- c()
for (sp1 in 1:dim(arctic_sitetab_animals)[2]) { # nolint: seq_linter.
  for (sp2 in 1:dim(arctic_sitetab_animals)[2]) { # nolint: seq_linter.
    empirical_cor_arctic <- c(empirical_cor_arctic, cor(arctic_sitetab_animals[, sp1], arctic_sitetab_animals[, sp2]))
  }
}
hist(empirical_cor_arctic)


sim_sitetab_sampled <- read.csv("/space/s1/fiona_callahan/multiSim_100sp/randomRun1/sim_sitetab_sampled10000_filtered.csv")
empirical_cor_sim <- c()
for (sp1 in 1:dim(sim_sitetab_sampled)[2]) { # nolint: seq_linter.
  for (sp2 in 1:dim(sim_sitetab_sampled)[2]) { # nolint: seq_linter.
    empirical_cor_sim <- c(empirical_cor_sim, cor(sim_sitetab_sampled[, sp1], sim_sitetab_sampled[, sp2]))
  }
}
hist(empirical_cor_sim)