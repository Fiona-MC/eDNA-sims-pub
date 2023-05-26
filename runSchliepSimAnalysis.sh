#!/bin/bash
export OMP_NUM_THREADS=5 # threads per job

# /home/fiona_callahan/eDNA_sims_code/runSchliepSimAnalysis.sh

sim_dir="/space/s1/fiona_callahan/multiSim11"
N=10 # number of jobs at a time

while IFS=',' read -r runNum; do
   (
        folder=${sim_dir}/randomRun${runNum}/
        echo "starting task $folder.."
        Rscript /home/fiona_callahan/eDNA_sims_code/schliep_sims.R ${folder} 10000 2000 
        sleep $(( (RANDOM % 3) + 1)) # choose random number 1, 2, or 3 and sleep for that long -- no idea why
    )&

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
done < ${sim_dir}/schliep_runNums.csv

wait
echo "all done"