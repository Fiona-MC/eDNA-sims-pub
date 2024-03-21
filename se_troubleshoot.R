library(igraph)
library(SpiecEasi)

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- SpiecEasi::make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
######## ADDED ########
PrecNoDiag <- Prec
PrecNoDiag[diag(nrow = dim(Prec)[1]) == 1] <- 0
graphTest <- igraph::graph_from_adjacency_matrix(PrecNoDiag != 0, mode = "undirected")
#########################
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)
# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done

huge::huge.roc(se$est$path, graph, verbose = FALSE)
############ ADDED #############
huge::huge.roc(se$est$path, graphTest, verbose = FALSE)
################################
#stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: getRefit(se)

sum(graph > 0) #254
sum(graph != 0) #254
sum(Cor < 0) #1142
sum(Cor > 0) #1317
sum(Cor < 0 & graph != 0) #128
sum(Cor > 0 & graph != 0) #126
sum(graphPos != 0) #126
sum(graphNeg != 0) #128

graphPos <- graph
graphPos[Cor < 0] <- 0
huge::huge.roc(se$est$path, graphPos, verbose=FALSE)

graphNeg <- graph
graphNeg[Cor > 0] <- 0
huge::huge.roc(se$est$path, graphNeg, verbose=FALSE)


graphPosPrec <- graph
graphPosPrec[Prec < 0] <- 0
huge::huge.roc(se$est$path, graphPosPrec, verbose=FALSE)

graphNegPrec <- graph
graphNegPrec[Prec > 0] <- 0
huge::huge.roc(se$est$path, graphNegPrec, verbose=FALSE)











# 10 species

data(amgut1.filt)
dim(amgut1.filt) # [1] 289 127
amgut1.filt <- amgut1.filt[, 1:10]
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- SpiecEasi::make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
PrecNoDiag <- Prec
PrecNoDiag[diag(nrow = dim(Prec)[1]) == 1] <- 0
graphTest <- graph_from_adjacency_matrix(PrecNoDiag != 0, mode = "undirected")
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-3, nlambda=100, pulsar.params = list(rep.num = 500))
# spiec.easi(as.matrix(sitetab_data[, params$names_species]), method = 'mb', lambda.min.ratio = 1e-3, nlambda = 100, pulsar.params = list(rep.num = 500))

# Applying data transformations...
# Selecting model with pulsar using stars...
# Fitting final estimate with mb...
# done

# these are the same as desired
huge::huge.roc(se$est$path, graph, verbose=FALSE)
huge::huge.roc(se$est$path, graphTest, verbose=FALSE)
#stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: getRefit(se)

sum(graph > 0) #20
sum(graph != 0) 
sum(Cor < 0) #24
sum(Cor > 0) #26
sum(Cor < 0 & graph != 0) #12
sum(Cor > 0 & graph != 0) #8

graphPos <- graph
graphPos[Cor < 0] <- 0
sum(graphPos != 0) #8
huge::huge.roc(se$est$path, graphPos, verbose=FALSE)

graphNeg <- graph
graphNeg[Cor > 0] <- 0
sum(graphNeg != 0) #12
huge::huge.roc(se$est$path, graphNeg, verbose=FALSE)





dirName <- c("multiSim_10sp")
run <- 1
logi <- TRUE
numSamples <- 1000
if(logi) {
  seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb_logi")
  seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso_logi")
  seName3 <- paste0("spiecEasi_res_sampled", numSamples, "_sparcc_logi")
} else {
  seName1 <- paste0("spiecEasi_res_sampled", numSamples, "_mb")
  seName2 <- paste0("spiecEasi_res_sampled", numSamples, "_glasso")
  seName3 <- paste0("spiecEasi_res_sampled", numSamples, "_sparcc")
}
subdir <- paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/", seName1, "/trial1/")
# spiec.easi(as.matrix(sitetab_data[, params$names_species]), method = 'mb', lambda.min.ratio = 1e-3, nlambda = 100, pulsar.params = list(rep.num = 500))
se <- readRDS(paste0(subdir, "se_rawRes.Rdata"))
params <- readRDS(paste0("/space/s1/fiona_callahan/", dirName, "/randomRun", run, "/params.Rdata"))

actualAlpha <- params$alpha
#alphaG <- graph_from_adjacency_matrix(actualAlpha < 0, mode = "undirected")
alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "undirected")
# theta = true_interactions
se.roc <- huge::huge.roc(se$est$path, theta = alphaG, verbose = FALSE)

