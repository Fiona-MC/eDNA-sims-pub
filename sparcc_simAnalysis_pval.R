#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
library(igraph)

# to run
# Rscript sparcc_simAnalysis.r /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun2/ /space/s1/fiona_callahan/multiSim_5sp_testing/randomRun2/spiecEasi_res/ glasso 1
#note -- so far not controlling for covs
args <- commandArgs(trailingOnly = TRUE)

data_dir <- "/space/s1/fiona_callahan/multiSim_10sp/randomRun1/"
save_dir <- "/space/s1/fiona_callahan/multiSim_10sp/randomRun1/sparcc_res_sampled100_filtered/"
cutoff <- "pval_bootstrap"
numTrials <- "1"
sitetab_name <- "sim_sitetab_readAbd_sampled100_filtered.csv"
cpus <- "10"

data_dir <- args[1]
save_dir <- args[2]
cutoff <- args[3]
numTrials <- args[4]
sitetab_name <- args[5]
cpus <- as.numeric(args[6]) # only relevant for pval_bootstrap mode

dir.create(save_dir)

# load data
locList <- readRDS(paste0(data_dir, "locList.Rdata"))
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

    if (cutoff == "pval_bootstrap") {
        #time_before <- Sys.time()
        sparccbootres <- sparccboot(data = as.matrix(sitetab_data[, names_species]), R = 100, ncpus = cpus)
        #print(Sys.time() - time_before)
        pval_sparcc <- pval.sparccboot(x = sparccbootres)

        pval_mx <- triu2diag(pval_sparcc$pvals, diagval = 1)

        inferenceRes <- list()
        inferenceRes$alphaInferred <- (pval_mx < 0.05) * 1
        inferenceRes$betaInferred <- NA
        
        saveRDS(sparccbootres, file = paste0(subdir, "sparcc_rawRes", cutoff, ".Rdata"))

        #print(Sys.time() - time_before)
        # pval_sparcc$pvals res here is either upper or lower triangle of mx -- gotta fugure out order

    } else {
        # example code from: https://github.com/zdk123/SpiecEasi#basic-usage
        #data(amgut1.filt)
        #"Also, we should use a larger number of stars repetitions  for real data" (meaning >50)
        res <- sparcc(as.matrix(sitetab_data[, names_species]))
        ## Define arbitrary threshold for SparCC correlation matrix for the graph
        res.graph <- abs(res$Cor) >= as.numeric(cutoff) #fiona edited
        diag(res.graph) <- 0
        res.graph <- Matrix(res.graph, sparse = TRUE)

        saveRDS(res, file = paste0(subdir, "sparcc_rawRes", cutoff, ".Rdata"))

        # get inferred alpha and beta (just 1's -1's and 0's )
        inferenceRes <- list()
        inferenceRes$alphaInferred <- as.matrix(res.graph)
        inferenceRes$betaInferred <- NA
        # put 0s on diagonal
        inferenceRes$alphaInferred <- inferenceRes$alphaInferred * as.integer(diag(nrow = num_species, ncol = num_species) == 0)
    }

    saveRDS(inferenceRes, paste0(subdir, "inferenceRes_cutoff", cutoff, ".Rdata"))
}
