#!/bin/bash
export OMP_NUM_THREADS=15

#for numSamples in 50 100 500 1000 5000 10000 25000; do
#./run_all_scr.sh /space/s1/fiona_callahan/multiSim_10sp 100 $numSamples
#done

#./run_all_scr.sh /space/s1/fiona_callahan/multiSim_10sp 10 500

sim_dir=$1
numRuns=$2
numSamples=$3
logi=$4
logi=0

./runINLAsimAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} ${logi} > ${sim_dir}/inlaCov${numSamples}_${logi}Out.log 2> ${sim_dir}/inlaCov${numSamples}_${logi}Err.log &

./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} ${logi} > ${sim_dir}/logisticCov${numSamples}_${logi}Out.log 2> ${sim_dir}/logisticCov${numSamples}_${logi}Err.log &
./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} ${logi} > ${sim_dir}/logistic${numSamples}_${logi}Out.log 2> ${sim_dir}/logistic${numSamples}_${logi}Err.log &

./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 ${logi} > ${sim_dir}/ecCov${numSamples}_${logi}Out.log 2> ${sim_dir}/ecCov${numSamples}_${logi}Err.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 ${logi} > ${sim_dir}/ec${numSamples}_${logi}Out.log 2> ${sim_dir}/ec${numSamples}_${logi}Err.log &

./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 ${logi} > ${sim_dir}/ecCov_readAbd_${numSamples}_${logi}Out.log 2> ${sim_dir}/ecCov_readAbd_${numSamples}_${logi}Err.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 ${logi} > ${sim_dir}/ec_readAbd_${numSamples}_${logi}Out.log 2> ${sim_dir}/ec_readAbd_${numSamples}_${logi}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} mb ${numSamples} ${logi} > ${sim_dir}/seMb${numSamples}_${logi}Out.log 2> ${sim_dir}/seMb${numSamples}_${logi}Err.log &
./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} glasso ${numSamples} ${logi} > ${sim_dir}/seGlasso${numSamples}_${logi}Out.log 2> ${sim_dir}/seGlasso${numSamples}_${logi}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} sparcc ${numSamples} ${logi} > ${sim_dir}/seSparcc${numSamples}_${logi}Out.log 2> ${sim_dir}/seSparcc${numSamples}_${logi}Err.log &




