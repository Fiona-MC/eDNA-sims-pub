#!/bin/bash
export OMP_NUM_THREADS=15

# ./run_all_noINLA.sh /space/s1/fiona_callahan/multiSim_10sp 100 0 1

sim_dir=$1
numRuns=$2
logi=$3
filtered=$4

#cd /home/fiona_callahan/eDNA_sims_code

# where are they going when I send them to the background?? 

#for numSamples in 10000 25000;
#for numSamples in 500 5000;
for numSamples in 100 1000 10000;
#for numSamples in 50 100 500 1000;
do
#./runINLAsimAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} > ${sim_dir}/inlaCov${numSamples}Out.log 2> ${sim_dir}/inlaCov${numSamples}Err.log &

./runLogisticSimAnalysis_filtered.sh ${sim_dir} ${numRuns} 1 ${numSamples} ${logi} ${filtered} > ${sim_dir}/logisticCov${numSamples}_${logi}Out.log 2> ${sim_dir}/logisticCov${numSamples}_${logi}Err.log &
./runLogisticSimAnalysis_filtered.sh ${sim_dir} ${numRuns} 0 ${numSamples} ${logi} ${filtered} > ${sim_dir}/logistic${numSamples}_${logi}Out.log 2> ${sim_dir}/logistic${numSamples}_${logi}Err.log &

#pres-abs
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 ${logi} ${filtered} > ${sim_dir}/ecCov${numSamples}_${logi}Out.log 2> ${sim_dir}/ecCov${numSamples}_${logi}Err.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 ${logi} ${filtered} > ${sim_dir}/ec${numSamples}_${logi}Out.log 2> ${sim_dir}/ec${numSamples}_${logi}Err.log &

#read abd
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 ${logi} ${filtered} > ${sim_dir}/ecCov_readAbd_${numSamples}_${logi}Out.log 2> ${sim_dir}/ecCov_readAbd_${numSamples}_${logi}Err.log &
./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 ${logi} ${filtered} > ${sim_dir}/ec_readAbd_${numSamples}_${logi}Out.log 2> ${sim_dir}/ec_readAbd_${numSamples}_${logi}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} mb ${numSamples} ${logi} ${filtered} > ${sim_dir}/seMb${numSamples}_${logi}Out.log 2> ${sim_dir}/seMb${numSamples}_${logi}Err.log &
./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} glasso ${numSamples} ${logi} ${filtered} > ${sim_dir}/seGlasso${numSamples}_${logi}Out.log 2> ${sim_dir}/seGlasso${numSamples}_${logi}Err.log &

./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} sparcc ${numSamples} ${logi} ${filtered} > ${sim_dir}/seSparcc${numSamples}_${logi}Out.log 2> ${sim_dir}/seSparcc${numSamples}_${logi}Err.log &
wait
done



