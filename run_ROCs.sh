
#cluster <- as.numeric(args[2]) == 1
#ratio_of_avg <- FALSE #do we compute the average of the ratio or ratio of averages
#numRuns <- 100
#numSamples <- args[3]
#logi <- as.numeric(args[4]) == 1
#saveRes <- TRUE
#covMode <- "noCount" # "all" "noCov" "cov" "covNoCount" "noCount"

#Rscript ./get_ROC_stats.R multiSim_10sp cluster_cov 100 0
#Rscript ./get_ROC_stats.R multiSim_10sp cluster_cov 10000 0
Rscript ./get_ROC_stats.R multiSim_10sp cluster_cov 100 1
Rscript ./get_ROC_stats.R multiSim_10sp cluster_cov 10000 1
#Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples cluster_cov 100 0
##Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples cluster_cov 10000 0

#Rscript ./get_ROC_stats.R multiSim_100sp cluster_cov 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp cluster_cov 10000 0
Rscript ./get_ROC_stats.R multiSim_100sp cluster_cov 100 1
Rscript ./get_ROC_stats.R multiSim_100sp cluster_cov 10000 1
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples cluster_cov 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples cluster_cov 10000 0

#Rscript ./get_ROC_stats.R multiSim_10sp cluster 100 0
#Rscript ./get_ROC_stats.R multiSim_10sp cluster 10000 0
Rscript ./get_ROC_stats.R multiSim_10sp cluster 100 1
Rscript ./get_ROC_stats.R multiSim_10sp cluster 10000 1
#Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples cluster 100 0
#Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples cluster 10000 0

#Rscript ./get_ROC_stats.R multiSim_100sp cluster 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp cluster 10000 0
Rscript ./get_ROC_stats.R multiSim_100sp cluster 100 1
Rscript ./get_ROC_stats.R multiSim_100sp cluster 10000 1
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples cluster 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples cluster 10000 0


#Rscript ./get_ROC_stats.R multiSim_10sp ignore_sign 100 0
#Rscript ./get_ROC_stats.R multiSim_10sp ignore_sign 10000 0
Rscript ./get_ROC_stats.R multiSim_10sp ignore_sign 100 1
Rscript ./get_ROC_stats.R multiSim_10sp ignore_sign 10000 1
#Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples ignore_sign 100 0
#Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples ignore_sign 10000 0

#Rscript ./get_ROC_stats.R multiSim_100sp ignore_sign 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp ignore_sign 10000 0
Rscript ./get_ROC_stats.R multiSim_100sp ignore_sign 100 1
Rscript ./get_ROC_stats.R multiSim_100sp ignore_sign 10000 1
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples ignore_sign 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples ignore_sign 10000 0


#Rscript ./get_ROC_stats.R multiSim_10sp ignore_direction 100 0
#Rscript ./get_ROC_stats.R multiSim_10sp ignore_direction 10000 0
Rscript ./get_ROC_stats.R multiSim_10sp ignore_direction 100 1
Rscript ./get_ROC_stats.R multiSim_10sp ignore_direction 10000 1
#Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples ignore_direction 100 0
#Rscript ./get_ROC_stats.R multiSim_10sp_random_moreSamples ignore_direction 10000 0

#Rscript ./get_ROC_stats.R multiSim_100sp ignore_direction 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp ignore_direction 10000 0
Rscript ./get_ROC_stats.R multiSim_100sp ignore_direction 100 1
Rscript ./get_ROC_stats.R multiSim_100sp ignore_direction 10000 1
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples ignore_direction 100 0
#Rscript ./get_ROC_stats.R multiSim_100sp_random_moreSamples ignore_direction 10000 0
