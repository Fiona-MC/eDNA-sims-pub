# sample fewer data points randomly from sitetabs and save new
library(data.table)
library(stringr)

#sim_dir <- "/space/s1/fiona_callahan/multiSim_10sp/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_10x10sp/"
sim_dir <- "/space/s1/fiona_callahan/multiSim_10sp_random_moreSamples/"
numRuns <- 100
nSamplesL <- c(100, 1000, 10000)

# number of samples where species must be present
#nSamplesL <- c(50, 100, 500, 1000)
#nSamplesL <- c(1000, 10000)
logi <- FALSE
for (nSamples in nSamplesL) {
    presenceCutoff <- 0.1 * nSamples
    if (logi) {
        sitetabName <- paste0("logiSim_sitetab_sampled", nSamples, ".csv")
        sitetabName_abd <- paste0("logiSim_sitetab_readAbd_sampled", nSamples, ".csv")
    } else {
        sitetabName <- paste0("sim_sitetab_sampled", nSamples, ".csv")
        sitetabName_abd <- paste0("sim_sitetab_readAbd_sampled", nSamples, ".csv")
    }


    for (run in 1:numRuns) {
        sim_subdir <- paste0(sim_dir, "randomRun", run, "/")
        params <- readRDS(paste0(sim_subdir, "params.Rdata"))

        sitetab <- fread(file = paste0(sim_subdir, sitetabName))
        sitetab_abd <- fread(file = paste0(sim_subdir, sitetabName_abd))

        # remove species that have less than cutoff presences
        criteria1 <- colSums(sitetab) < presenceCutoff & str_detect(names(sitetab), "Sp")
        # absence cutoff
        criteria2 <- (nSamples - colSums(sitetab)) < presenceCutoff & str_detect(names(sitetab), "Sp")
        criteria <- criteria1

        #(data.table syntax here)
        sitetab_reSampled <- sitetab[, .SD, .SDcols = !criteria]
        sitetab_reSampled_abd <- sitetab_abd[, .SD, .SDcols = !criteria] 

        params$filteredSpNames <- names(sitetab_reSampled)[grep("Sp", names(sitetab_reSampled))]
        
        # Extract numeric parts using regular expressions
        numeric_parts <- gsub("\\D", "", params$filteredSpNames)
        # Convert to numeric
        numeric_species <- as.numeric(numeric_parts)

        actualBeta <- sign(params$beta)
        params$filteredBeta <- actualBeta[numeric_species, ]

        actualAlpha <- sign(params$alpha)
        params$filteredAlpha <- actualAlpha[numeric_species, numeric_species]

        saveRDS(params, file = paste0(sim_subdir, "paramsFiltered", nSamples, ".Rdata"))
        
        #print(sum(criteria1))
        #print(sum(criteria2))
        #print("")

        if (logi) {
            fwrite(sitetab_reSampled, file = paste0(sim_subdir, "logiSim_sitetab_sampled", nSamples, "_filtered.csv"))
            fwrite(sitetab_reSampled_abd, file = paste0(sim_subdir, "logiSim_sitetab_readAbd_sampled", nSamples, "_filtered.csv"))
        } else {
            fwrite(sitetab_reSampled, file = paste0(sim_subdir, "sim_sitetab_sampled", nSamples, "_filtered.csv"))
            fwrite(sitetab_reSampled_abd, file = paste0(sim_subdir, "sim_sitetab_readAbd_sampled", nSamples, "_filtered.csv"))
        }
    }
}