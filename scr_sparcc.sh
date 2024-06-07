# run all linear regression
#sim_dir=$1
#numRuns=$2
#numSamples=$3
#logi=$4
#filtered=$5 #100

# 10sp
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 100 0 1
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 10000 0 1

#logi
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 100 1 0
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 10000 1 0

# 100sp
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 100 0 1
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 10000 0 1

#logi
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 100 1 0
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 10000 1 0


# 10sp_random
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 100 0 1
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples 100 10000 0 1

# 100sp_random
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 100 0 1
./run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples 100 10000 0 1

