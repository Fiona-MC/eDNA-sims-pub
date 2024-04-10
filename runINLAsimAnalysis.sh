#!/bin/bash
export OMP_NUM_THREADS=5

#./runINLAsimAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp_random 100 100 1
#./runINLAsimAnalysis.sh /global/scratch/users/fionacallahan/multiSim_10sp 100 ${numSamples}

sim_dir=$1
#sim_dir="/space/s1/fiona_callahan/multiSim_100"
#numRuns=100
numRuns=$2
numSamples=$3
filtered=$4

echo "Starting INLA"
echo $sim_dir
echo $numRuns
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

ROC_mode="noModelSelect" # this will mean there is no WAIC selection for the ones where the cutoff changes

if [ ${filtered} == 1 ]
then
    sitetab=sim_sitetab_sampled${numSamples}_filtered.csv
    resDirName=${resDirName}_filtered
fi
#INLA_type="faster"

#Rscript /home/fiona_callahan/eDNA_sims_code/filter_sims.R ${sim_dir}/ ${numRuns}
# Rscript /home/fiona_callahan/filter_sims.R /space/s1/fiona_callahan/multiSim3/

# this takes a minute (cp takes awhile -- if we just delete them it will be fast)
#mkdir ${sim_dir}/unrealisticRuns
#while IFS=',' read -r lineNum runNum reason; do
#    if [ -e "${sim_dir}/randomRun${runNum}" ]; then
    #    cp -pR "${sim_dir}/randomRun${runNum}" "${sim_dir}/unrealisticRuns/randomRun${runNum}"
    #    rm -r "${sim_dir}/randomRun${runNum}" #put this back when you edit the filtering criteria?
#    fi
#done < ${sim_dir}/unrealistic_runNums.csv


N=2 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)

folderNames=() # array of names of randomRun* so that diff numbers of runs can be done 
for ((i=1; i<=$numRuns; i++)); do
  folderNames+=(${sim_dir}/randomRun$i)
done

#folder=(${sim_dir}/randomRun1)

for folder in ${folderNames[@]}; do
    (
        #if test ! -d "${folder}/INLA_res_${INLA_type}/trial1" # if the folder is not already there 
        #then
            echo "starting task $folder.."
            mkdir "$folder/$resDirName/" 
            for modelParms in none cov sp spCov; do
                # run INLA sim analysis
                echo $modelParms
                timeout -k 10 ${timeout1}h Rscript ./INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ ${sitetab} ${modelParms} ${filtered}
            done
            Rscript ./INLA_modelSelect.R ${folder}/ ${folder}/${resDirName}/ ${filtered}
            #./runINLA_checkAndReRun.sh ${sim_dir} ${resDirName} ${numRuns} 1 ${timeout2} ${INLA_type} ${sitetab}
            Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ 1

            #for cutoff in 0.01;
            #for cutoff in 0 0.0000000000001 0.0000001 0.00001 .3 .5 .7 .9 1;
            for cutoff in 0 1 0.0000001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 .3 .5;
            do
                saveDirName=${resDirName}_cov
                mkdir "$folder/$saveDirName/"
                Rscript ./INLA_changeCutoffs.R ${folder}/ ${cutoff} ${folder}/${saveDirName}/ ${folder}/${resDirName}/ ${ROC_mode} 1 ${filtered}
                Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${saveDirName}/ 1 ${cutoff}

                saveDirName=${resDirName}_noCov
                mkdir "$folder/$saveDirName/"
                Rscript ./INLA_changeCutoffs.R ${folder}/ ${cutoff} ${folder}/${saveDirName}/ ${folder}/${resDirName}/ ${ROC_mode} 0 ${filtered}
                Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${saveDirName}/ 0 ${cutoff}
            done
            sleep $(( (RANDOM % 3) + 1)) # choose random number 1, 2, or 3 and sleep for that long -- no idea why
       #fi
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
done

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

#for cutoff in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15;
#for cutoff in 0.01;
for cutoff in 0 1 0.0000001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 .3 .5;
do
saveDirName=${resDirName}_cov
Rscript ./gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${saveDirName} ${cutoff}

saveDirName=${resDirName}_noCov
Rscript ./gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${saveDirName} ${cutoff}
done

Rscript ./gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName}

echo "all done"