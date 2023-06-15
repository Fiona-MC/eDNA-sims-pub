#!/bin/bash
export OMP_NUM_THREADS=15

sim_dir="/space/s1/fiona_callahan/multiSim_ParmSet1"
numRuns=100
numTrials=1

N=10 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)

for folder in ${sim_dir}/randomRun*; do
    (
    echo "starting task $folder.."  
    Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim.R ${sim_dir}/ ${numRuns}
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

echo "all done"
