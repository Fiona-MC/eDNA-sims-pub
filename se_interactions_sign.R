# this is basically scratch coding space

runs <- 100
numPos <- 0
numNeg <- 0
for (run in 1:runs){
    inferredParms <- readRDS(paste0("/space/s1/fiona_callahan/multiSim_100/randomRun", run, "/spiecEasi_res_mb/trial1/inferenceRes.Rdata"))
    numPos <- numPos + sum(inferredParms$alphaInferred == 1)
    numNeg <- numNeg + sum(inferredParms$alphaInferred == -1)
}
numPos
numNeg
