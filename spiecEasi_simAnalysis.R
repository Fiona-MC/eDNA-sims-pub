#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
library(igraph)

# to run
# Rscript spiecEasi_simAnalysis.R /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun2/ /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun2/spiecEasi_res/ glasso 1

#note -- so far not controlling for covs
args <- commandArgs(trailingOnly = TRUE)
data_dir <- "/space/s1/fiona_callahan/multiSim_manySp_testing2/randomRun3/"
save_dir <- "/space/s1/fiona_callahan/multiSim_manySp_testing2/randomRun3/spiecEasi_res_dumb_mb/"
data_dir <- args[1]
save_dir <- args[2]
se.method <- "mb"
se.method <- args[3]
numTrials <- args[4]
#se.method <- "mb"

# load data
locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))
sitetab_data <- read.csv(paste0(data_dir, "sitetab_abd_dumb.csv"))
sitetab_data <- read.csv(paste0(data_dir, "sim_sitetab_readAbd_sampled.csv"))

names_cov <- params$names_cov
names_species <- params$names_species


for (trial in 1:numTrials) {
    subdir <- paste0(save_dir, "trial", trial, "/")
    dir.create(subdir)
    # example code from: https://github.com/zdk123/SpiecEasi#basic-usage
    #data(amgut1.filt)
    #"Also, we should use a larger number of stars repetitions  for real data" (meaning >50)
    if (se.method == "mb") {
        res <- spiec.easi(as.matrix(sitetab_data[, params$names_species]), method = 'mb', lambda.min.ratio = 1e-3,
                            nlambda = 20, pulsar.params = list(rep.num = 500))
        res.graph <- getRefit(res)
    } else if (se.method == "glasso") {
        res <- spiec.easi(as.matrix(sitetab_data[, params$names_species]), method = 'glasso', lambda.min.ratio = 1e-3,
                            nlambda = 20, pulsar.params = list(rep.num = 500))
        res.graph <- getRefit(res)
    } else if (se.method == "sparcc") {
        res <- sparcc(as.matrix(sitetab_data[, params$names_species]))
        ## Define arbitrary threshold for SparCC correlation matrix for the graph
        res.graph <- abs(res$Cor) >= 0.3 #fiona edited
        diag(res.graph) <- 0
        res.graph <- Matrix(res.graph, sparse = TRUE)
    } 

    # get inferred alpha and beta (just 1's -1's and 0's )
    inferenceRes <- list()
    inferenceRes$alphaInferred <- as.matrix(res.graph)
    inferenceRes$betaInferred <- NA
    # put 0s on diagonal
    inferenceRes$alphaInferred <- inferenceRes$alphaInferred * as.integer(diag(nrow = params$numSpecies, ncol = params$numSpecies) == 0)

    saveRDS(inferenceRes, paste0(subdir, "inferenceRes.Rdata"))
}