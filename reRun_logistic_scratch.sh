#data_dir <- args[1]
#numRuns <- as.numeric(args[2])
#covs <- (as.numeric(args[3]) == 1)
#logi <- (as.numeric(args[4]) == 1)
#sitetab_name <- args[5]
#outName <- args[6]

#100sp
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 1 0 sim_sitetab_sampled100_filtered.csv logistic_mistakes_sampled100_covNoCount_filtered_100runs
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 0 0 sim_sitetab_sampled100_filtered.csv logistic_mistakes_sampled100_noCov_filtered_100runs 

Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 1 0 sim_sitetab_sampled1000_filtered.csv logistic_mistakes_sampled1000_covNoCount_filtered_100runs
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 0 0 sim_sitetab_sampled1000_filtered.csv logistic_mistakes_sampled1000_noCov_filtered_100runs

Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 1 0 sim_sitetab_sampled10000_filtered.csv logistic_mistakes_sampled10000_covNoCount_filtered_100runs
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 0 0 sim_sitetab_sampled10000_filtered.csv logistic_mistakes_sampled10000_noCov_filtered_100runs 

#100sp random
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples/ 100 1 0 sim_sitetab_sampled100_filtered.csv logistic_mistakes_sampled100_covNoCount_filtered_100runs
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples/ 100 0 0 sim_sitetab_sampled100_filtered.csv logistic_mistakes_sampled100_noCov_filtered_100runs 

Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples/ 100 1 0 sim_sitetab_sampled10000_filtered.csv logistic_mistakes_sampled10000_covNoCount_filtered_100runs
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples/ 100 0 0 sim_sitetab_sampled10000_filtered.csv logistic_mistakes_sampled10000_noCov_filtered_100runs 

#100sp logi
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 1 1 logiSim_sitetab_sampled100.csv logistic_mistakes_sampled100_covNoCount_logi_100runs
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 0 1 logiSim_sitetab_sampled100.csv logistic_mistakes_sampled100_noCov_logi_100runs 

Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 1 1 logiSim_sitetab_sampled10000.csv logistic_mistakes_sampled10000_covNoCount_logi_100runs
Rscript logistic_filtered.R /space/s1/fiona_callahan/multiSim_100sp/ 100 0 1 logiSim_sitetab_sampled10000.csv logistic_mistakes_sampled10000_noCov_logi_100runs 
