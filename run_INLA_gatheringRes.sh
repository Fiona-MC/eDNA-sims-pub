#!/bin/bash
export OMP_NUM_THREADS=5

#./run_INLA_gatheringRes.sh /space/s1/fiona_callahan/sim_paper_stuff/multiSim_10sp_revision2 100 10000 1 0

sim_dir=$1
#sim_dir="/space/s1/fiona_callahan/savio/multiSim_10sp_random"
#numRuns=100
numRuns=$2
numSamples=$3
logi=$4
filtered=$5

echo $sim_dir
echo $numRuns
echo $numSamples


numTrials=1
INLA_type="paperSep"

if [ ${numSamples} == "None" ]
then
	resDirName=INLA_res_${INLA_type}
else
	resDirName=INLA_res_${INLA_type}_sampled${numSamples}
fi

scramble=0

#sitetab="sim_sitetab_sampled.csv"
if [ ${numSamples} == "None" ]
then
	sitetab="sim_sitetab_sampled.csv"
else
	sitetab=sim_sitetab_sampled${numSamples}.csv
fi

if [ ${filtered} == 1 ]
then
    sitetab=sim_sitetab_sampled${numSamples}_filtered.csv
    resDirName=${resDirName}_filtered
fi

if [ ${logi} == 1 ]
then
    sitetab=logiSim_sitetab_sampled${numSamples}.csv
    resDirName=${resDirName}_logi
fi

ROC_mode="noModelSelect" # this will mean there is no WAIC selection for the ones where the cutoff changes

#for cutoff in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15;
#for cutoff in 0.01;
for cutoff in 0 1 0.0000001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 .3 .5;
do
    saveDirName="${resDirName}_cov"
    Rscript ./gather_inferenceRes_general.R ${sim_dir}/ ${numRuns} ${numTrials} ${saveDirName} ${cutoff}

    saveDirName="${resDirName}_noCov"
    Rscript ./gather_inferenceRes_general.R ${sim_dir}/ ${numRuns} ${numTrials} ${saveDirName} ${cutoff}

    saveDirName="${resDirName}_covNoCount"
    Rscript ./gather_inferenceRes_general.R ${sim_dir}/ ${numRuns} ${numTrials} ${saveDirName} ${cutoff}
done

Rscript ./gather_inferenceRes_general.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName}

echo "all done"