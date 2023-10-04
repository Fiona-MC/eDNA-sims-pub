#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)

data_dir <- "/space/s1/fiona_callahan/multiSim11a/randomRun30/"

# load data
locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))
sim_data_raw <- readRDS(paste0(data_dir, "sim_data_abd.Rdata"))


# example code from: https://github.com/zdk123/SpiecEasi#basic-usage
data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
huge::huge.roc(se$est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: getRefit(se)