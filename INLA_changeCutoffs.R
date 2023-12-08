library(data.table)
numRuns <- 100
numTrials <- 2
sim_dirs <- sapply(X = 1:numRuns, FUN = function(run) {paste0(sim_dir, "randomRun", run, "/")})
cutoffs <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15)

#sim_dir = sim_dirs[2]
#trial = 1
#cutoff = cutoffs[1]

for (sim_dir in sim_dirs) {
    for (trial in 1:numTrials) {
        for (cutoff in cutoffs) {
            save_dir <- paste0(sim_dir, "INLA_res_paper/")
            subdir <- paste0(save_dir, "trial", trial, "/")
            params <- readRDS(paste0(sim_dir, "params.Rdata"))
            #save results for later
            sim_lists <- readRDS(paste0(subdir, "sim_lists-", trial, ".Rdata"))

            modelcands <- names(sim_lists)
       
            names_species <- params$names_species
            names_cov <- params$names_cov
            # waic
            waictab <- sapply(names_species, function(response) {
                return(sapply(modelcands, function(modelname) {
                            sim_lists[[modelname]][[response]]$waic$waic}))
            })
            rownames(waictab) <- modelcands
            
            # get inferred alpha and beta (just 1's -1's and 0's )
            inferenceRes <- list()
            inferenceRes$alphaInferred <- matrix(data = 0, nrow = length(names_species), ncol = length(names_species))
            inferenceRes$betaInferred <- matrix(data = 0, nrow = length(names_species), ncol = params$numCovs)
            for (sp_index in seq_along(names_species)) {
                this_best_model <- modelcands[which.min(waictab[, sp_index])]
                best_model_res <- sim_lists[[this_best_model]][[names_species[sp_index]]]
                significant_pos_vars <- row.names(best_model_res$summary.fixed)[best_model_res$summary.fixed$`0.025quant` > 0]
                significant_neg_vars <- row.names(best_model_res$summary.fixed)[best_model_res$summary.fixed$`0.975quant` < 0]
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
        }
    }
}