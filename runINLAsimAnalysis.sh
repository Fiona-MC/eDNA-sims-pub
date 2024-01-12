#!/bin/bash
export OMP_NUM_THREADS=5

#./runINLAsimAnalysis.sh /space/s1/fiona_callahan/multiSim_100sp 100 500

sim_dir=$1
#sim_dir="/space/s1/fiona_callahan/multiSim_100sp"
#numRuns=100
numRuns=$2
numSamples=$3

numTrials=1
INLA_type="paper"

if [ ${numSamples} == "None" ]
then
	resDirName=INLA_res_${INLA_type}
else
	resDirName=INLA_res_${INLA_type}_sampled${numSamples}
fi

scramble=0
covs=1

#sitetab="sim_sitetab_sampled.csv"
if [ ${numSamples} == "None" ]
then
	sitetab="sim_sitetab_sampled.csv"
else
	sitetab=sim_sitetab_sampled${numSamples}.csv
fi

ROC_mode="noModelSelect" # this will mean there is no WAIC selection for the ones where the cutoff changes

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


N=5 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)

for folder in ${sim_dir}/randomRun*; do
    (
        #if test ! -d "${folder}/INLA_res_${INLA_type}/trial1" # if the folder is not already there 
        #then
            echo "starting task $folder.."
            mkdir "$folder/$resDirName/" 
            # run INLA sim analysis
            timeout -k 10 10h Rscript /home/fiona_callahan/eDNA_sims_code/INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ $sitetab
            #for cutoff in 0.01;
            #for cutoff in 0 0.0000000000001 0.0000001 0.00001 .3 .5 .7 .9 1;
            for cutoff in 0 1 0.0000001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 .3 .5;
            do
                Rscript /home/fiona_callahan/eDNA_sims_code/INLA_changeCutoffs.R ${folder}/ ${cutoff} ${folder}/${resDirName}/ ${ROC_mode}
                Rscript /home/fiona_callahan/eDNA_sims_code/count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ ${covs} ${cutoff}
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
Rscript /home/fiona_callahan/eDNA_sims_code/gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName} ${cutoff}
done

echo "all done"