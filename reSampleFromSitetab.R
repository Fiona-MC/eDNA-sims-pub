# sample fewer data points randomly from sitetabs and save new
library(data.table)

#sim_dir <- "/space/s1/fiona_callahan/multiSim_10sp/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_10x10sp/"
sim_dir <- "/space/s1/fiona_callahan/multiSim_100sp_random/"
save_dir <- "/space/s1/fiona_callahan/multiSim_100sp_random/"
numRuns <- 100

resample <- TRUE # do you want to just copy the stuff or resample it?
#nSamplesL <- c(50, 100, 500, 1000, 5000, 10000, 25000)
#nSamplesL <- c(50, 100, 500, 1000)
nSamplesL <- c(100)
dumb <- FALSE

if (dumb) {
    full_sitetabName <- "sitetab_dumb.csv"
    full_sitetabName_abd <- "sitetab_abd_dumb.csv"
} else {
    full_sitetabName <- "sim_sitetab_sampled.csv"
    full_sitetabName_abd <- "sim_sitetab_readAbd_sampled.csv"
}

dir.create(save_dir)

for (nSamples in nSamplesL) {
    for (run in 1:numRuns) {
        save_subdir <- paste0(save_dir, "randomRun", run, "/")
        sim_subdir <- paste0(sim_dir, "randomRun", run, "/")
        if(!dir.exists(save_subdir)) {dir.create(save_subdir)}
        file.copy(from = paste0(sim_subdir, "params.Rdata"), to = paste0(save_subdir, "params.Rdata"))
        file.copy(from = paste0(sim_subdir, "covList.Rdata"), to = paste0(save_subdir, "covList.Rdata"))
        file.copy(from = paste0(sim_subdir, "locList.Rdata"), to = paste0(save_subdir, "locList.Rdata"))
        full_sitetab <- fread(file = paste0(sim_subdir, full_sitetabName))
        full_sitetab_abd <- fread(file = paste0(sim_subdir, full_sitetabName_abd))
        if (dim(full_sitetab)[1] >= nSamples) {
            sample <- sample(x = 1:dim(full_sitetab)[1], size = nSamples, replace = FALSE) # nolint: seq_linter.
            sitetab_reSampled <- full_sitetab[sample, ]
            sitetab_reSampled_abd <- full_sitetab_abd[sample, ]
            if (dumb) {
                fwrite(sitetab_reSampled, file = paste0(save_subdir, "logiSim_sitetab_sampled", nSamples, ".csv"))
                fwrite(sitetab_reSampled_abd, file = paste0(save_subdir, "logiSim_sitetab_readAbd_sampled", nSamples, ".csv"))
            } else {
                fwrite(sitetab_reSampled, file = paste0(save_subdir, "sim_sitetab_sampled", nSamples, ".csv"))
                fwrite(sitetab_reSampled_abd, file = paste0(save_subdir, "sim_sitetab_readAbd_sampled", nSamples, ".csv"))
            }
        } else {
            print(paste("for run", run, ": dim(full_sitetab) =", dim(full_sitetab), "<", "nSamples=", nSamples))
        }
    }
}