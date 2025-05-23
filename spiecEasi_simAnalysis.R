#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
library(igraph)

# to run
# Rscript spiecEasi_simAnalysis.R /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun2/ /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun2/spiecEasi_res/ glasso 1
#note -- so far not controlling for covs
args <- commandArgs(trailingOnly = TRUE)
#data_dir <- "/space/s1/fiona_callahan/multiSim_manySp_testing2/randomRun4/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_manySp_testing2/randomRun4/spiecEasi_res_logiDir_mb/"
#data_dir <- "/space/s1/fiona_callahan/multiSim_2sp_test/randomRun1/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_2sp_test/randomRun1/spiecEasi_res_mb/"

#data_dir <- "/space/s1/fiona_callahan/multiSim_10sp_dep/randomRun1/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_10sp_dep/randomRun1/spiecEasi_res_test/"

data_dir <- args[1]
save_dir <- args[2]
se.method <- args[3]
numTrials <- args[4]
sitetab_name <- args[5]
#se.method <- "mb"
dir.create(save_dir)

# load data
#locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))
#sitetab_data <- read.csv(paste0(data_dir, "sitetab_abd_logi.csv"))
#sitetab_data <- read.csv(paste0(data_dir, "sitetab_abd_logi_dir.csv"))
sitetab_data <- read.csv(paste0(data_dir, sitetab_name))

names_cov <- params$names_cov
names_species <- names(sitetab_data)[grep("Sp", names(sitetab_data))]
num_species <- length(names_species)

for (trial in 1:numTrials) {
    subdir <- paste0(save_dir, "trial", trial, "/")
    dir.create(subdir)
    # example code from: https://github.com/zdk123/SpiecEasi#basic-usage
    #data(amgut1.filt)
    #"Also, we should use a larger number of stars repetitions  for real data" (meaning >50)
    if (se.method == "mb") {
        res <- spiec.easi(as.matrix(sitetab_data[, names_species]), method = 'mb', lambda.min.ratio = 1e-3,
                            nlambda = 100, pulsar.params = list(rep.num = 500))
        res.graph <- getRefit(res)
    } else if (se.method == "glasso") {
        res <- spiec.easi(as.matrix(sitetab_data[, names_species]), method = 'glasso', lambda.min.ratio = 1e-3,
                            nlambda = 100, pulsar.params = list(rep.num = 500))
        res.graph <- getRefit(res)
    } else if (se.method == "sparcc") {
        res <- sparcc(as.matrix(sitetab_data[, names_species]))
        ## Define arbitrary threshold for SparCC correlation matrix for the graph
        res.graph <- abs(res$Cor) >= 0.3 #fiona edited
        diag(res.graph) <- 0
        res.graph <- Matrix(res.graph, sparse = TRUE)
    } 

    saveRDS(res, file = paste0(subdir, "se_rawRes.Rdata"))

    # get inferred alpha and beta (just 1's -1's and 0's )
    inferenceRes <- list()
    inferenceRes$alphaInferred <- as.matrix(res.graph)
    inferenceRes$betaInferred <- NA
    # put 0s on diagonal
    inferenceRes$alphaInferred <- inferenceRes$alphaInferred * as.integer(diag(nrow = num_species, ncol = num_species) == 0)

    saveRDS(inferenceRes, paste0(subdir, "inferenceRes.Rdata"))
}
