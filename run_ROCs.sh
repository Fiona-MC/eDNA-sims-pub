
#cluster <- as.numeric(args[2]) == 1
#ratio_of_avg <- FALSE #do we compute the average of the ratio or ratio of averages
#numRuns <- 100
#numSamples <- args[3]
#logi <- as.numeric(args[4]) == 1
#saveRes <- TRUE
#covMode <- "noCount" # "all" "noCov" "cov" "covNoCount" "noCount"

Rscript ./get_ROC_stats.R multiSim_10sp 1 100 0
Rscript ./get_ROC_stats.R multiSim_10sp 1 10000 0
Rscript ./get_ROC_stats.R multiSim_10sp 1 100 1
Rscript ./get_ROC_stats.R multiSim_10sp 1 10000 1
Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples 1 100 0
Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples 1 10000 0

Rscript ./get_ROC_stats.R multiSim_100sp 1 100 0
Rscript ./get_ROC_stats.R multiSim_100sp 1 10000 0
Rscript ./get_ROC_stats.R multiSim_100sp 1 100 1
Rscript ./get_ROC_stats.R multiSim_100sp 1 10000 1
Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples 1 100 0
Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples 1 10000 0