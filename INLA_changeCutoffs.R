library(data.table)

args <- commandArgs(trailingOnly = TRUE)

sim_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun2/"
cutoff <- 1
save_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun2/INLA_res_faster/"
ROC_mode <- "noModelSelect"

sim_dir <- args[1]
cutoff <- as.numeric(args[2])
save_dir <- args[3]
ROC_mode <- args[4] #"noModelSelect" or "modelSelect"


trial <- 1
subdir <- paste0(save_dir, "trial", trial, "/")
params <- readRDS(paste0(sim_dir, "params.Rdata"))
sim_lists <- readRDS(paste0(subdir, "sim_lists-", trial, ".Rdata"))

modelcands <- names(sim_lists)

names_species <- params$names_species
names_cov <- params$names_cov

# waic
if (length(modelcands) > 1) {
    waictab <- sapply(names_species, function(response) {
        return(sapply(modelcands, function(modelname) {
                    sim_lists[[modelname]][[response]]$waic$waic}))
    })
    rownames(waictab) <- modelcands
}

# get inferred alpha and beta (just 1's -1's and 0's )
inferenceRes <- list()
inferenceRes$alphaInferred <- matrix(data = 0, nrow = length(names_species), ncol = length(names_species))
inferenceRes$betaInferred <- matrix(data = 0, nrow = length(names_species), ncol = params$numCovs)
for (sp_index in seq_along(names_species)) {
    if (ROC_mode == "noModelSelect") {
        if (length(modelcands) > 1) {
            this_best_model <- "Cov+Animal"
        } else {
            this_best_model <- modelcands[1]
        }
    } else {
        if(length(modelcands) > 1) {
            this_best_model <- modelcands[which.min(waictab[, sp_index])]
        } else {
            this_best_model <- modelcands[1]
        }
    }
    best_model_res <- sim_lists[[this_best_model]][[names_species[sp_index]]]
    #lower_bound <- best_model_res$summary.fixed$`0.025quant`
    #upper_bound <- best_model_res$summary.fixed$`0.975quant`
    z_score <- qnorm(1 - cutoff / 2)
    lower_bound <- best_model_res$summary.fixed$mean - z_score * best_model_res$summary.fixed$sd
    upper_bound <- best_model_res$summary.fixed$mean + z_score * best_model_res$summary.fixed$sd
    significant_pos_vars <- row.names(best_model_res$summary.fixed)[lower_bound > 0]
    significant_neg_vars <- row.names(best_model_res$summary.fixed)[upper_bound < 0]
    for (sp_index2 in seq_along(names_species)) { # put 1's and -1's into alpha and beta inferred based on 95% CI
    if (names_species[sp_index2] %in% significant_pos_vars) {inferenceRes$alphaInferred[sp_index, sp_index2] <- 1}
    if (names_species[sp_index2] %in% significant_neg_vars) {inferenceRes$alphaInferred[sp_index, sp_index2] <- -1}
    }
    for (cov_index in seq_along(names_cov)) {
    if (names_cov[cov_index] %in% significant_pos_vars) {inferenceRes$betaInferred[sp_index, cov_index] <- 1}
    if (names_cov[cov_index] %in% significant_neg_vars) {inferenceRes$betaInferred[sp_index, cov_index] <- -1}
    }
}

saveRDS(inferenceRes, paste0(subdir, "inferenceRes_cutoff", cutoff, ".Rdata"))