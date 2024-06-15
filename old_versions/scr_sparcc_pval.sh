#data_dir <- "/space/s1/fiona_callahan/multiSim_10sp/randomRun1/"
#save_dir <- "/space/s1/fiona_callahan/multiSim_10sp/randomRun1/sparcc_res_sampled100_filtered/"
#cutoff <- "pval_bootstrap"
#numTrials <- "1"
#sitetab_name <- "sim_sitetab_readAbd_sampled100_filtered.csv"
#cpus <- 1

# run all linear regression
#sim_dir=$1
#numRuns=$2
#numSamples=$3
#logi=$4
#filtered=$5 #100

# 10sp
#./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 100 0 1
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 10000 0 1

#logi
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 100 1 0
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 10000 1 0

# 100sp
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 100 0 1
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 10000 0 1

#logi
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 100 1 0
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 10000 1 0


# 10sp_random
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 100 0 1
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 10000 0 1

# 100sp_random
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 100 0 1
./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 10000 0 1

