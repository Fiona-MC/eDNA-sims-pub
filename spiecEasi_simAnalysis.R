#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
#note -- so far not controlling for covs
args <- commandArgs(trailingOnly = TRUE)
data_dir <- "/space/s1/fiona_callahan/multiSim_5sp_testing/randomRun1/"
data_dir <- args[1]

# load data
locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))
sim_data_raw <- readRDS(paste0(data_dir, "sim_data_abd.Rdata"))
sitetab_data <- read.csv(paste0(data_dir, "sim_sitetab_readAbd_sampled.csv"))

# example code from: https://github.com/zdk123/SpiecEasi#basic-usage
#data(amgut1.filt)

depths <- rowSums(sitetab_data[, params$names_species]) #TODO filter out any rows where there are no reads
simData.n  <- t(apply(sitetab_data[depths != 0, params$names_species], 1, norm_to_total))
simData.cs <- round(simData.n * min(depths))

d <- ncol(simData.cs)
n <- nrow(simData.cs)
e <- d

#set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(simData.cs, mar = 2, distr = 'zinegbin', Sigma = Cor, n = n)

se <- spiec.easi(X, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 15)
huge::huge.roc(se$est$path, graph, verbose = FALSE)
stars.pr(getOptMerge(se), graph, verbose = FALSE)
# stars selected final network under: getRefit(se)