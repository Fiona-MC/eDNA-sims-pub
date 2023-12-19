#!/usr/bin/env Rscript

# to run
# Rscript INLA_simAnalysis_paper.R /space/s1/fiona_callahan/multiSim_rw3/randomRun1/ /space/s1/fiona_callahan/multiSim_rw3/randomRun1/INLAres/

# NOTE this is the original unedited version!!!
source("/home/fiona_callahan/Arctic_eDNA_2021/script/INLA_ST_functions.R")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("input and output files need to be supplied", call. = FALSE)
} 

#data_dir <- "/home/fiona_callahan/simData/testing/run1/"
#save_dir <- "/home/fiona_callahan/simData/testing/INLAres/"
data_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun10/"
save_dir <- "/space/s1/fiona_callahan/multiSim_100/randomRun10/INLA_res_faster/"
sitetabName <- "sim_sitetab_sampled.csv"
#${folder}/ ${folder}/INLA_res_${INLA_type}/#
data_dir <- args[1]
save_dir <- args[2]
dir.create(save_dir)
sitetabName <- args[3]


plot <- FALSE
proj <- FALSE
numTrials <- 1

# load data
# sim_data_raw <- readRDS(paste0(data_dir, "sim_data.Rdata"))
sitetab <- read.csv(paste0(data_dir, sitetabName))
#locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))


names_cov <- params$names_cov
names_species <- params$names_species

# run analysis three times on the same data
for (trial in 1:numTrials){
  subdir <- paste0(save_dir, "trial", trial, "/")
  dir.create(subdir)
  # Select what you'll be working with
  k <- 11 # k is number of time knots
  #k <- 4
  # Prepare spatial and temporal meshes
  # alltidx <- seq(1,k)
  #reducedtidx <- seq(2,k,2)
  tknots <- seq(min(sitetab$Age), max(sitetab$Age), length = k)
  #edge case, this is probably the wrong way to handle this
  #(fixes bug if only one time pt)
  if(min(sitetab$Age) == max(sitetab$Age)) {
    tknots <- seq(min(sitetab$Age) - 1, max(sitetab$Age), length = k)
  }
  #reducedtknots <- tknots[seq(2,k,2)]
  mesh_t <- inla.mesh.1d(tknots)

  xytab <- cbind(sitetab$Long, sitetab$Lat)
  #xytab <- LonLatToAzi(xytab)
  n <- dim(sitetab)[1]
  
  bound <- inla.nonconvex.hull(as.matrix(xytab))
  #TODO what are these params??
  #mesh_s <- inla.mesh.2d(boundary = bound, max.edge = c(200, 400), max.n=200, cutoff = 500)
  mesh_s <- inla.mesh.2d(boundary = bound, max.edge = c(10, 25))#, max.edge = c(2, 4), max.n=200, cutoff = 500)
  plot(mesh_s)
  spde_s <- inla.spde2.matern(mesh_s)
  
  res_all <- list()
  for (response in names_species){ 
    # ADD MODELS
    # Run INLA 
    # Run with other animals and environment as covariates
    all_names_cov <- c(names_cov, names_species)
    all_names_cov <- all_names_cov[all_names_cov != response]
    res_all[[response]] <- RunInlaBin(response, sitetab, xytab, mesh_s, mesh_t, 1, 
                              namescov = all_names_cov, normcov = 2, project = proj)
  }
  
  sim_lists <- list()
  # ADD MODELS
  sim_lists[["Cov+Animal"]] <- res_all
  
  # save results for later
  saveRDS(sim_lists, paste0(subdir, "sim_lists-", trial, ".Rdata"))
  #readRDS(paste0(subdir, "sim_lists-", trial, ".Rdata"))
  run_data <- list()
  run_data$params <- params
  # run_data$sim_data_raw <- sim_data_raw
  run_data$sitetab <- sitetab
  saveRDS(run_data, paste0(subdir, "rundata", trial, ".Rdata"))
  #readRDS(paste0(subdir, "rundata", trial, ".Rdata"))

  # CPO 
  modelcands <- names(sim_lists)
  #modelcands <- c("None", "AllCov", "Cov+Animal")
  
  # waic
  waictab <- sapply(names_species, function(response) {
    return(sapply(modelcands, function(modelname) {
                sim_lists[[modelname]][[response]]$waic$waic}))
  })
  #rownames(waictab) <- modelcands
  #waictabmelted <- melt(waictab)
  #colnames(waictabmelted) <- c("Model", "Animal", "WAIC")

  if (plot) { #this doesnt work for only one model
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
    if (length(modelcands) > 1) {
      this_best_model <- modelcands[which.min(waictab[, sp_index])]
    } else {
      this_best_model <- modelcands[1]
    }
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
    covarplots <- CreateCovarPlots(waictab, names_species, sim_lists, forcemodel = "AllCov")
    cov_only <- arrangeGrob(grobs = covarplots)
    ggsave(cov_only, file = paste0(subdir, "covonly_INLA_", trial, ".png"))
    
    # Model with cov and animal
    covarplots <- CreateCovarPlots(waictab, names_species, sim_lists, forcemodel = "Cov+Animal")
    animal_cov <- arrangeGrob(grobs = covarplots)
    ggsave(animal_cov, file = paste0(subdir, "animal_cov_INLA_", trial, ".png"))
  }
}
