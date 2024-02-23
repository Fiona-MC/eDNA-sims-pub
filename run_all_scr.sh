#!/bin/bash
export OMP_NUM_THREADS=15

#for numSamples in 50 100 500 1000 5000 10000 25000; do
#./run_all_scr.sh /space/s1/fiona_callahan/multiSim_10sp 100 $numSamples
#done

#./run_all_scr.sh /space/s1/fiona_callahan/multiSim_10sp 10 500

sim_dir=$1
numRuns=$2
numSamples=$3

./runINLAsimAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} ${dumb} > ${sim_dir}/inlaCov${numSamples}_${dumb}Out.log 2> ${sim_dir}/inlaCov${numSamples}_${dumb}Err.log &

./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} ${dumb} > ${sim_dir}/logisticCov${numSamples}_${dumb}Out.log 2> ${sim_dir}/logisticCov${numSamples}_${dumb}Err.log &
./runLogisticSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} ${dumb} > ${sim_dir}/logistic${numSamples}_${dumb}Out.log 2> ${sim_dir}/logistic${numSamples}_${dumb}Err.log &

./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 ${dumb} > ${sim_dir}/ecCov${numSamples}_${dumb}Out.log 2> ${sim_dir}/ecCov${numSamples}_${dumb}Err.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 ${dumb} > ${sim_dir}/ec${numSamples}_${dumb}Out.log 2> ${sim_dir}/ec${numSamples}_${dumb}Err.log &

./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 ${dumb} > ${sim_dir}/ecCov_readAbd_${numSamples}_${dumb}Out.log 2> ${sim_dir}/ecCov_readAbd_${numSamples}vErr.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 ${dumb} > ${sim_dir}/ec_readAbd_${numSamples}_${dumb}Out.log 2> ${sim_dir}/ec_readAbd_${numSamples}_${dumb}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} mb ${numSamples} ${dumb} > ${sim_dir}/seMb${numSamples}_${dumb}Out.log 2> ${sim_dir}/seMb${numSamples}_${dumb}Err.log &
./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} glasso ${numSamples} ${dumb} > ${sim_dir}/seGlasso${numSamples}_${dumb}Out.log 2> ${sim_dir}/seGlasso${numSamples}_${dumb}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} sparcc ${numSamples} ${dumb} > ${sim_dir}/seSparcc${numSamples}_${dumb}Out.log 2> ${sim_dir}/seSparcc${numSamples}_${dumb}Err.log &




