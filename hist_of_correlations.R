library(ggplot2)

# two uses -- 1 uses a list of num species and makes  wrapped plot just looking at the overall pattern in covariances
# second use is to look at the covariances for one folder (either 10 or 100 sp) and split by whether an interactuon was inferred

num_species <- c(10, 100)
#num_species <- c(10)
numRuns <- 100
regType <- "logistic" #"linearReg" "logistic"
numSamples <- 100

nonzero_covars_L <- c()
runs_L <- c()
nspL <- c()
covar_inferredTrueL <- c()
covar_inferredFalseL <- c()
for (nsp in num_species) {
    sim_dir <- paste0("/space/s1/fiona_callahan/sim_paper_stuff/multiSim_", nsp, "sp_revision2/")
    nonzero_covars <- c()
    runs <- c()
    for (run in 1:numRuns) {
        covarMx_L <- readRDS(paste0(sim_dir, "randomRun", run, "/rev_covar_mx.Rdata"))
        inferred_parms <- readRDS(paste0(sim_dir, "randomRun", run, "/", regType, "_mistakes_sampled", numSamples, "_noCov_logi_100runs/inferenceRes_cutoffbh.Rdata"))
        
        covarMx <- covarMx_L$covarMx
        diag(covarMx) <- 0
        these_nonzero <- covarMx[covarMx != 0]
        nonzero_covars <- c(nonzero_covars, these_nonzero)
        runs <- c(runs, rep(run, times = length(these_nonzero)))

        covar_inferredTrueL <- c(covar_inferredTrueL, covarMx[inferred_parms$alphaInferred != 0])
        covar_inferredFalseL <- c(covar_inferredFalseL, covarMx[inferred_parms$alphaInferred == 0])
    }
    nonzero_covars_L <- c(nonzero_covars_L, nonzero_covars)
    runs_L <- c(runs_L, runs)
    nspL <- c(nspL, rep(nsp, times = length(nonzero_covars)))
    print(sim_dir)
    min(nonzero_covars)
    max(nonzero_covars)
}

df <- data.frame(run = runs_L, nonzero_covars = nonzero_covars_L, nSpecies = nspL)

ggplot(df, aes(x = nonzero_covars, fill = as.factor(nSpecies))) + 
  geom_histogram(binwidth = .1, color = "black", alpha = 0.7) +
  facet_wrap(~ nSpecies, scales = "free_y") +
  labs(
    x = "Non-Zero Covariance Values",
    y = "Count",
    fill = "Number of Species"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("/space/s1/fiona_callahan/sim_paper_stuff/revision2_covariances.pdf")

df_inferred <- data.frame(covariance = c(covar_inferredTrueL, covar_inferredFalseL), 
                inferred = c(rep(TRUE, times = length(covar_inferredTrueL)), rep(FALSE, times = length(covar_inferredFalseL))))

ggplot(df_inferred, aes(x = covariance, fill = inferred)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.1, alpha = 0.7, position = "dodge") +
  labs(
    title = paste0(regType, ", ", numSamples, " samples", num_species[1], "species"),
    x = "Covariance",
    y = "Density",
    fill = "Inferred"
  ) +
  theme_minimal()


df_inferred_nozeros <- df_inferred[df_inferred$covariance != 0, ]

ggplot(df_inferred_nozeros, aes(x = covariance, fill = inferred)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.1, alpha = 0.7, position = "dodge") +
  labs(
    title = paste0(regType, ", ", numSamples, " samples", num_species[1], "species -- no 0's"),
    x = "Covariance",
    y = "Density",
    fill = "Inferred"
  ) +
  theme_minimal()

