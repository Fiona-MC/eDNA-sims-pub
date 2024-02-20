#!/bin/bash
export OMP_NUM_THREADS=15

#for numSamples in 50 100 500 1000 5000 10000 25000; do
#./run_all_scr.sh /space/s1/fiona_callahan/multiSim_10sp 100 $numSamples
#done

sim_dir=$1
numRuns=$2
numSamples=$3

./runINLAsimAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} > ${sim_dir}/inlaCov${numSamples}Out.log 2> ${sim_dir}/inlaCov${numSamples}Err.log &

./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} > ${sim_dir}/logisticCov${numSamples}Out.log 2> ${sim_dir}/logisticCov${numSamples}Err.log &
./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} > ${sim_dir}/logistic${numSamples}Out.log 2> ${sim_dir}/logistic${numSamples}Err.log &

./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 > ${sim_dir}/ecCov${numSamples}Out.log 2> ${sim_dir}/ecCov${numSamples}Err.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 > ${sim_dir}/ec${numSamples}Out.log 2> ${sim_dir}/ec${numSamples}Err.log &

./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 > ${sim_dir}/ecCov_readAbd_${numSamples}Out.log 2> ${sim_dir}/ecCov_readAbd_${numSamples}Err.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 > ${sim_dir}/ec_readAbd_${numSamples}Out.log 2> ${sim_dir}/ec_readAbd_${numSamples}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} mb ${numSamples} > ${sim_dir}/seMb${numSamples}Out.log 2> ${sim_dir}/seMb${numSamples}Err.log &
./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} glasso ${numSamples} > ${sim_dir}/seGlasso${numSamples}Out.log 2> ${sim_dir}/seGlasso${numSamples}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} sparcc ${numSamples} > ${sim_dir}/seSparcc${numSamples}Out.log 2> ${sim_dir}/seSparcc${numSamples}Err.log &




