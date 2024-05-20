library(stats)
library(igraph)
library(stringr)

data_dir <- "/space/s1/fiona_callahan/multiSim_10sp_random_testSet/"
numRuns <- 100
covs <- (as.numeric("1") == 1)
logi <- (as.numeric("0") == 1)
sitetab_name <- "sim_sitetab_sampled100_filtered.csv"
outName <- "logistic_mistakes_sampled100_covNoCount_filtered_100runs"
countCovs <- FALSE

randomSim <- str_detect(data_dir, "random")

# to run
# Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim_moreSp.R /space/s1/fiona_callahan/multiSim_5sp_random/ 1000

runs <- 1:numRuns
numTrials <- 1
trials <- 1:1

for (run in runs) {
    for (trial in trials) { # basically ignore the trials thing -- I think this is deterministic so trials should be irrelevant
        if (file.exists(paste0(data_dir, "randomRun", run))) { # this is for the runs that were deleted
            print("run")
            print(run)
            
            sim_sitetab_sampled <- read.csv(file = paste0(data_dir, "randomRun", run, "/", sitetab_name), header = TRUE)

            speciesNames <- names(sim_sitetab_sampled)[grep("Sp", names(sim_sitetab_sampled))]
            numSpecies <- length(speciesNames)
            print(numSpecies)
        }
    }
}


