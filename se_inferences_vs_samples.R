# this is for supplemental figure

library(igraph)
library(SpiecEasi)
library(ggplot2)

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

graph <- SpiecEasi::make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
######## ADDED ########
#PrecNoDiag <- Prec
#PrecNoDiag[diag(nrow = dim(Prec)[1]) == 1] <- 0
#graphTest <- igraph::graph_from_adjacency_matrix(PrecNoDiag != 0, mode = "undirected")
#########################
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

#se does it infer the same appx num of interactions regardless of num samples
sitetab <- read.csv("/space/s1/fiona_callahan/multiSim_100sp/randomRun1/sim_sitetab_readAbd_sampled10000_filtered.csv")
names_species <- names(sitetab)[grep("Sp", names(sitetab))]
sitetab <- sitetab[, names_species]

nInferredL <- c()
nSamplesL <- c()
dataSetL <- c()

for (n_samples in c(10, 25, 50, 75, 100, 125, 150, 175, 200, 250, 275)) {
  ### se data ###
  thisDataSet <- "se.synth"
  X_temp <- X[1:n_samples, ]
  se_temp <- spiec.easi(X_temp, method='mb', lambda.min.ratio=1e-2, nlambda=20)
  res.graph <- getRefit(se_temp)
  n_inferred <- sum(as.matrix(res.graph))

  nInferredL <- c(nInferredL, n_inferred)
  nSamplesL <- c(nSamplesL, n_samples)
  dataSetL <- c(dataSetL, thisDataSet)

  ### My data ###
  thisDataSet <- "multiSim_100sp"
  X_sitetab <- as.matrix(sitetab[1:n_samples, ])
  se_sitetab <- spiec.easi(X_sitetab, method='mb', lambda.min.ratio=1e-2, nlambda=20)
  res.graph_sitetab <- getRefit(se_sitetab)
  n_inferred_sitetab <- sum(as.matrix(res.graph_sitetab))

  nInferredL <- c(nInferredL, n_inferred_sitetab)
  nSamplesL <- c(nSamplesL, n_samples)
  dataSetL <- c(dataSetL, thisDataSet)

    ### real data ###
  thisDataSet <- "amgut1.filt"
  X_real <- amgut1.filt[1:n_samples, ]
  se_real <- spiec.easi(X_real, method='mb', lambda.min.ratio=1e-2, nlambda=20)
  res.graph_real <- getRefit(se_real)
  n_inferred_real <- sum(as.matrix(res.graph_real))

  nInferredL <- c(nInferredL, n_inferred_real)
  nSamplesL <- c(nSamplesL, n_samples)
  dataSetL <- c(dataSetL, thisDataSet)
}

data <- data.frame(NumSamples = nSamplesL, NumInferences = nInferredL, dataset = dataSetL)


# Create the ggplot
ggplot(data, aes(x = NumSamples, y = NumInferences, color = dataset)) +
  geom_point() +       # Add points
  geom_line() +        # Add lines
  labs(
    x = "Number of Samples",
    y = "Number of Inferences"
  ) +
  ylim(0, NA) +   
  theme_minimal()      # Use a minimal theme

ggsave("/space/s1/fiona_callahan/se_inferences_vs_samples.pdf")
