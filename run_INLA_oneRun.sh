#!/bin/bash
export OMP_NUM_THREADS=5

#./run_INLA_oneRun.sh /space/s1/fiona_callahan/multiSim_10sp ${runNum} 100 0 1

sim_dir=$1
#sim_dir="/space/s1/fiona_callahan/multiSim_100"
runNum=$2
numSamples=$3
logi=$4
filtered=$5

echo "Starting INLA"
echo $sim_dir
echo $runNum
echo $numSamples

numTrials=1
INLA_type="paperSep"
timeout1=10
timeout2=24

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

folder=(${sim_dir}/randomRun${runNum})

echo "starting task $folder.."
mkdir "$folder/$resDirName/" 
for modelParms in none cov sp spCov; do
    # run INLA sim analysis
    echo $modelParms
    #timeout -k 10 ${timeout1}h Rscript ./INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ ${sitetab} ${modelParms} ${filtered}
done
Rscript ./INLA_modelSelect.R ${folder}/ ${folder}/${resDirName}/ ${filtered}
#./runINLA_checkAndReRun.sh ${sim_dir} ${resDirName} ${numRuns} 1 ${timeout2} ${INLA_type} ${sitetab}
Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ 0

for cutoff in 0 1 0.0000000000001 0.0000001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 .3 .5 .7 .9;
do
    saveDirName=${resDirName}_cov
    mkdir "$folder/$saveDirName/"
    Rscript ./INLA_changeCutoffs.R ${folder}/ ${cutoff} ${folder}/${saveDirName}/ ${folder}/${resDirName}/ ${ROC_mode} 1 ${filtered}
    Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${saveDirName}/ 1 ${cutoff}

    saveDirName=${resDirName}_noCov
    mkdir "$folder/$saveDirName/"
    Rscript ./INLA_changeCutoffs.R ${folder}/ ${cutoff} ${folder}/${saveDirName}/ ${folder}/${resDirName}/ ${ROC_mode} 0 ${filtered}
    Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${saveDirName}/ 0 ${cutoff}

    saveDirName=${resDirName}_covNoCount
    mkdir "$folder/$saveDirName/"
    Rscript ./INLA_changeCutoffs.R ${folder}/ ${cutoff} ${folder}/${saveDirName}/ ${folder}/${resDirName}/ ${ROC_mode} 1 ${filtered}
    Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${saveDirName}/ 0 ${cutoff}
done

echo task ${folder} done