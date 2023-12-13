# sample fewer data points randomly from sitetabs and save new
library(data.table)

sim_dir <- "/space/s1/fiona_callahan/multiSim_100/"
nSamples <- 500
numRuns <- 100

for (run in 1:numRuns) {
    save_dir <- paste0(sim_dir, "randomRun", run, "/")
    full_sitetab <- fread(file = paste0(save_dir, "sim_sitetab_sampled.csv"))
    full_sitetab_abd <- fread(file = paste0(save_dir, "sim_sitetab_readAbd_sampled.csv"))
    sample <- sample(x = 1:dim(full_sitetab)[1], size = nSamples, replace = FALSE) # nolint: seq_linter.
    sitetab_reSampled <- full_sitetab[sample, ]
    sitetab_reSampled_abd <- full_sitetab_abd[sample, ]
    fwrite(sitetab_reSampled, file = paste0(save_dir, "sim_sitetab_sampled", nSamples, ".csv"))
    fwrite(sitetab_reSampled_abd, file = paste0(save_dir, "sim_sitetab_readAbd_sampled", nSamples, ".csv"))
}
