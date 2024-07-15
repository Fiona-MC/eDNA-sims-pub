#!/bin/bash
export OMP_NUM_THREADS=5

#./runINLAsimAnalysis_scr.sh /space/s1/fiona_callahan/multiSim_100sp 100 100 0 1

sim_dir=$1
#sim_dir="/space/s1/fiona_callahan/multiSim_100"
#numRuns=100
numRuns=$2
numSamples=$3
logi=$4
filtered=$5

echo "Starting INLA"
echo $sim_dir
echo $numRuns
echo $numSamples
echo "filtered=" $filtered


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

ROC_mode="noModelSelect" # this will mean there is no WAIC selection for the ones where the cutoff changes

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

N=1 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)

folder=${sim_dir}/randomRun1


echo "starting task $folder.." >> ./INLA_time.txt
echo $(date) >> ./INLA_time.txt
for modelParms in none cov sp spCov; do
    # run INLA sim analysis
    echo $modelParms >> ./INLA_time.txt
    echo Rscript ./INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ ${sitetab} ${modelParms} ${filtered} >> INLA_time.txt

    Rscript ./INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ ${sitetab} ${modelParms} ${filtered}
    echo $(date) >> INLA_time.txt
done
echo Rscript ./INLA_modelSelect.R ${folder}/ ${folder}/${resDirName}/ ${filtered} >> INLA_time.txt
Rscript ./INLA_modelSelect.R ${folder}/ ${folder}/${resDirName}/ ${filtered}

echo $(date) >> INLA_time.txt
echo done >> INLA_time.txt
