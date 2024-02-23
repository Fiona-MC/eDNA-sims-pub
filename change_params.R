# the parameter num_samples_time was wrong

thisDir <- "/space/s1/fiona_callahan/multiSim_10sp/"
numRuns <- 100

for (run in 1:numRuns) {
    subdir <- paste0(thisDir, "randomRun", run, "/")
    oldParams <- readRDS(paste0(subdir, "params.Rdata"))
    newParams <- oldParams
    newParams$num_samples_time <- 500
    saveRDS(newParams, paste0(subdir, "paramsFixed.Rdata"))
}