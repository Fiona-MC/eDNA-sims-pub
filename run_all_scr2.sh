#!/bin/bash
export OMP_NUM_THREADS=15

# ./run_all_scr2.sh

numRuns=100
logi=0
filtered=1

#cd /home/fiona_callahan/eDNA_sims_code

#for numSamples in 10000 25000;
#for numSamples in 500 5000;
for sim_dir in /space/s1/fiona_callahan/multiSim_10sp /space/s1/fiona_callahan/multiSim_100sp /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples;
do
    for numSamples in 100 10000;
    do
    #./runINLAsimAnalysis_scr.sh ${sim_dir} ${numRuns} ${numSamples} ${logi} ${filtered} &

    #./runLogisticSimAnalysis_filtered.sh ${sim_dir} ${numRuns} 1 ${numSamples} ${logi} ${filtered} > ${sim_dir}/logisticCov${numSamples}_${logi}Out.log 2> ${sim_dir}/logisticCov${numSamples}_${logi}Err.log &
    #./runLogisticSimAnalysis_filtered.sh ${sim_dir} ${numRuns} 0 ${numSamples} ${logi} ${filtered} > ${sim_dir}/logistic${numSamples}_${logi}Out.log 2> ${sim_dir}/logistic${numSamples}_${logi}Err.log &

    #pres-abs
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 ${logi} ${filtered} 
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 ${logi} ${filtered} 

    #read abd
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 ${logi} ${filtered} 
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 ${logi} ${filtered} 

    #./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} mb ${numSamples} ${logi} ${filtered} &
    #./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} glasso ${numSamples} ${logi} ${filtered} &

    #./run_sparcc_simAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} ${logi} ${filtered} &
    wait
    done
done




logi=1
filtered=0

#cd /home/fiona_callahan/eDNA_sims_code

#for numSamples in 10000 25000;
#for numSamples in 500 5000;
for sim_dir in /space/s1/fiona_callahan/multiSim_10sp /space/s1/fiona_callahan/multiSim_100sp;
do
    for numSamples in 100 10000;
    do
    #./runINLAsimAnalysis_scr.sh ${sim_dir} ${numRuns} ${numSamples} ${logi} ${filtered} &

    #./runLogisticSimAnalysis_filtered.sh ${sim_dir} ${numRuns} 1 ${numSamples} ${logi} ${filtered} > ${sim_dir}/logisticCov${numSamples}_${logi}Out.log 2> ${sim_dir}/logisticCov${numSamples}_${logi}Err.log &
    #./runLogisticSimAnalysis_filtered.sh ${sim_dir} ${numRuns} 0 ${numSamples} ${logi} ${filtered} > ${sim_dir}/logistic${numSamples}_${logi}Out.log 2> ${sim_dir}/logistic${numSamples}_${logi}Err.log &

    #pres-abs
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 ${logi} ${filtered} 
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 ${logi} ${filtered} 

    #read abd
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 ${logi} ${filtered} 
    ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 ${logi} ${filtered} 

    #./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} mb ${numSamples} ${logi} ${filtered} &
    #./run_spiecEasiSimAnalysis.sh ${sim_dir} ${numRuns} glasso ${numSamples} ${logi} ${filtered} &

    #./run_sparcc_simAnalysis.sh ${sim_dir} ${numRuns} ${numSamples} ${logi} ${filtered} &
    wait
    done
done
