#!/bin/bash
export OMP_NUM_THREADS=15

sim_dir="/space/s1/fiona_callahan/multiSim8"
numRuns=1000
numTrials=1

Rscript /home/fiona_callahan/eDNA_sims_code/filter_sims.R ${sim_dir}/
# Rscript /home/fiona_callahan/filter_sims.R /space/s1/fiona_callahan/multiSim3/

# this takes a minute (cp takes awhile -- if we just delete them it will be fast)
mkdir ${sim_dir}/unrealisticRuns
while IFS=',' read -r lineNum runNum; do
    if [ -e "${sim_dir}/randomRun${runNum}" ]; then
        cp -pR "${sim_dir}/randomRun${runNum}" "${sim_dir}/unrealisticRuns/randomRun${runNum}"
        rm -r "${sim_dir}/randomRun${runNum}"
    fi
done < ${sim_dir}/unrealistic_runNums.csv


N=8 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)

for folder in ${sim_dir}/randomRun*; do
    (
        echo "starting task $folder.."
        if test ! -d "${folder}/INLA_res_faster/trial1" # if the folder is not already there NOT WORKING
        then
            mkdir "$folder/INLA_res_faster/" 
            # run INLA sim analysis
            timeout -k 10 6h Rscript /home/fiona_callahan/eDNA_sims_code/INLA_simAnalysis_faster.R ${folder}/ ${folder}/INLA_res_faster/
            Rscript /home/fiona_callahan/eDNA_sims_code/count_INLAmistakes.R ${folder}/ ${folder}/INLA_res_faster/
        fi
        sleep $(( (RANDOM % 3) + 1)) # choose random number 1, 2, or 3 and sleep for that long -- no idea why
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

Rscript /home/fiona_callahan/eDNA_sims_code/gather_inferenceRes.R ${sim_dir}/ ${numRuns} ${numTrials}

echo "all done"
