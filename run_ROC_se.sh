#multiSim_10sp, 10000 samp, neg only
Rscript ./ROC_se_troubleshoot.R multiSim_10sp ignore_direction 10000 0 negOnly
#multiSim_10sp, 10000 samp, pos only
Rscript ./ROC_se_troubleshoot.R multiSim_10sp ignore_direction 10000 0 posOnly
#multiSim_10sp, 10000 samp, both 
Rscript ./ROC_se_troubleshoot.R multiSim_10sp ignore_direction 10000 0 both

#multiSim_100sp, 10000 samp, neg only
Rscript ./ROC_se_troubleshoot.R multiSim_100sp ignore_direction 10000 0 negOnly
#multiSim_100sp, 10000 samp, pos only
Rscript ./ROC_se_troubleshoot.R multiSim_100sp ignore_direction 10000 0 posOnly
#multiSim_100sp, 10000 samp, both 
Rscript ./ROC_se_troubleshoot.R multiSim_100sp ignore_direction 10000 0 both

#multiSim_10sp, logi, 10000 samp, neg only
Rscript ./ROC_se_troubleshoot.R multiSim_10sp ignore_direction 10000 1 negOnly
#multiSim_10sp, logi, 10000 samp, pos only
Rscript ./ROC_se_troubleshoot.R multiSim_10sp ignore_direction 10000 1 posOnly
#multiSim_10sp, logi, 10000 samp, both 
Rscript ./ROC_se_troubleshoot.R multiSim_10sp ignore_direction 10000 1 both

#multiSim_100sp, logi, 10000 samp, neg only
Rscript ./ROC_se_troubleshoot.R multiSim_100sp ignore_direction 10000 1 negOnly
#multiSim_100sp, logi, 10000 samp, pos only
Rscript ./ROC_se_troubleshoot.R multiSim_100sp ignore_direction 10000 1 posOnly
#multiSim_100sp, logi, 10000 samp, both 
Rscript ./ROC_se_troubleshoot.R multiSim_100sp ignore_direction 10000 1 both
