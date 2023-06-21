#!/bin/bash
export OMP_NUM_THREADS=5 # threads per job

 #/home/fiona_callahan/eDNA_sims_code/runSchliepSimAnalysis.sh

# note: the file schliep_runNums.csv is produced in analyze_multiSim.R

sim_dir="/space/s1/fiona_callahan/multiSim11"
N=10 # number of jobs at a time
runNumFile="schliepRun20.csv"

while IFS=',' read -r runNum; do
   (
        folder=${sim_dir}/randomRun${runNum}/
        echo "starting task $folder.."
        #timeout -k 10 168h Rscript /home/fiona_callahan/eDNA_sims_code/schliep_sims.R ${folder} 50000 10000
        #echo ${folder}
        Rscript /home/fiona_callahan/eDNA_sims_code/schliep_sims.R ${folder} 100000 10000 
        sleep $(( (RANDOM % 3) + 1)) # choose random number 1, 2, or 3 and sleep for that long -- no idea why
    )&

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
done < ${sim_dir}/${runNumFile}

wait

# get file for analysis in analyze_multiSim.R and test convergence 
# .csv with convergence diagnostics will be produced for each folder
# the 2 here is the number of schliep trials -- which was hard coded into schliep_sims.R
Rscript /home/fiona_callahan/eDNA_sims_code/schliep_examineRes.R ${sim_dir}/ 2 ${runNumFile}

echo "all done"
