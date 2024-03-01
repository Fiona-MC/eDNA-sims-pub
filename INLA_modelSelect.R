source("./Arctic_eDNA_2021/script/INLA_ST_functions.R")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input and output files need to be supplied", call. = FALSE)
} 

#data_dir <- "/space/s1/fiona_callahan/multiSim_2x50sp/randomRun1/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_2x50sp/randomRun1/INLA_res_paperSep_sampled100/"

data_dir <- args[1]
save_dir <- args[2]

modelParmsL <- c("none", "cov", "sp", "spCov")
numTrials <- 1
plot <- FALSE

for (trial in 1:numTrials) {
    subdir <- paste0(save_dir, "trial", trial, "/")

    # load inla results
    sim_lists <- list()
    for (modelParms in modelParmsL) {
        sim_lists[[modelParms]] <- readRDS(paste0(subdir, "resList_", modelParms, ".Rdata"))
    }

    # CPO 
    modelcands <- names(sim_lists)

    params <- readRDS(paste0(data_dir, "params.Rdata"))
    names_species <- params$names_species
    names_cov <- params$names_cov
    
    cpotab <- sapply(names_species, function(response) {
        return(sapply(modelcands, function(modelname) { 
            -sum(log(sim_lists[[modelname]][[response]]$cpo$cpo), na.rm = TRUE) }))
    })
    
    rownames(cpotab) <- modelcands

    cpotabmelted <- melt(cpotab)
    colnames(cpotabmelted) <- c("Model", "Animal", "CPO")
    cpoplot <- ggplot(cpotabmelted) +
        geom_point(aes(x = Model, y = CPO, group = Animal, colour = Animal)) +
        geom_line(aes(x = Model, y = CPO, group = Animal, colour = Animal)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.title.x = element_blank(),
            text = element_text(size = 10)
        )
    
    # waic
    waictab <- sapply(names_species, function(response) {
        return(sapply(modelcands, function(modelname) {
                    sim_lists[[modelname]][[response]]$waic$waic}))
    })
    rownames(waictab) <- modelcands
    waictabmelted <- melt(waictab)
    colnames(waictabmelted) <- c("Model", "Animal", "WAIC")

    if (plot) { 
        waicplot <- ggplot(waictabmelted) +
        geom_point(aes(x = Model, y = WAIC, group = Animal, colour = Animal)) +
        geom_line(aes(x = Model, y = WAIC, group = Animal, colour = Animal)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title.x = element_blank(),
                text = element_text(size = 10)
        )
        
        
        cpoplotlet <- arrangeGrob(grobs = list(cpoplot), top = textGrob("A", x = unit(0, "npc"), 
                        y   = unit(1, "npc"), just = c("left", "top"), gp = gpar(col = "black", fontsize = 22)))
        waicplotlet <- arrangeGrob(grobs = list(waicplot), top = textGrob("B", x = unit(0, "npc"), 
                        y   = unit(1, "npc"), just = c("left", "top"), gp = gpar(col = "black", fontsize = 22)))
        model_plot <- grid.arrange(grobs = list(cpoplotlet, waicplotlet), 
                        ncol = 2, vp = viewport(width = 0.98, height = 0.98))
        ggsave(model_plot, file = paste0(subdir, "model_comparisons_", trial, ".png"))
    }
    
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

    saveRDS(inferenceRes, paste0(subdir, "inferenceRes.Rdata"))

    if (plot) {
        # Best model
        covarplots <- CreateCovarPlots(scoretab = waictab, names = names_species, modellist = sim_lists)
        best_WAIC <- arrangeGrob(grobs = covarplots) 
        ggsave(best_WAIC, file = paste0(subdir, "bestWAIC_INLA_", trial, ".png"))
        
        # Model without all cov NOT animal
        covarplots <- CreateCovarPlots(waictab, names_species, sim_lists, forcemodel = "cov")
        cov_only <- arrangeGrob(grobs = covarplots)
        ggsave(cov_only, file = paste0(subdir, "covOnly_INLA_", trial, ".png"))

        # Model without all cov NOT animal
        covarplots <- CreateCovarPlots(waictab, names_species, sim_lists, forcemodel = "sp")
        cov_only <- arrangeGrob(grobs = covarplots)
        ggsave(cov_only, file = paste0(subdir, "spOnly_INLA_", trial, ".png"))
        
        # Model with cov and animal
        covarplots <- CreateCovarPlots(waictab, names_species, sim_lists, forcemodel = "spCov")
        animal_cov <- arrangeGrob(grobs = covarplots)
        ggsave(animal_cov, file = paste0(subdir, "sp_cov_INLA_", trial, ".png"))
    }
}