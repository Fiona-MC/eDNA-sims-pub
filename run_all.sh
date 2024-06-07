#!/bin/bash
export OMP_NUM_THREADS=15

# ./run_all.sh /space/s1/fiona_callahan/multiSim_10sp 100 0 1

sim_dir=$1
numRuns=$2
logi=$3
filtered=$4

#cd /home/fiona_callahan/eDNA_sims_code

# where are they going when I send them to the background?? 

#for numSamples in 50 100 500 1000 10000 25000;
for numSamples in 100 10000;
do
#./runINLAsimAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} ${logi} ${filtered}

#./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 
#./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 

#./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 ${logi} ${filtered}
#./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 ${logi} ${filtered}

#./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 ${logi} ${filtered}
#./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 ${logi} ${filtered}

#./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} mb ${numSamples} ${logi} ${filtered}
#./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} glasso ${numSamples} ${logi} ${filtered}

./run_sparcc_simAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} ${logi} ${filtered}
done



