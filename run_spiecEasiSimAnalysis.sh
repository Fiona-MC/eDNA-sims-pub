#!/bin/bash
export OMP_NUM_THREADS=5

sim_dir="/space/s1/fiona_callahan/multiSim_10sp_indep"
numRuns=20
numTrials=1 # I think as this is implemented right now this needs to be 1
seMethod=mb
#seMethod=sparcc
#seMethod=glasso
random=1
plot=1

resDirName=spiecEasi_res_${seMethod}

#Rscript /home/fiona_callahan/eDNA_sims_code/filter_sims.R ${sim_dir}/ ${numRuns}
# Rscript /home/fiona_callahan/filter_sims.R /space/s1/fiona_callahan/multiSim3/

# this takes a minute (cp takes awhile -- if we just delete them it will be fast)
#mkdir ${sim_dir}/unrealisticRuns
#while IFS=',' read -r lineNum runNum reason; do
#    if [ -e "${sim_dir}/randomRun${runNum}" ]; then
#        cp -pR "${sim_dir}/randomRun${runNum}" "${sim_dir}/unrealisticRuns/randomRun${runNum}"
#        rm -r "${sim_dir}/randomRun${runNum}"
#    fi
#done < ${sim_dir}/unrealistic_runNums.csv



N=5 # 

for folder in ${sim_dir}/randomRun*; do 
    (
        #if test ! -d "${folder}/${resDirName}/trial1" # if the folder is not already there NOT WORKING
        #then
            echo "starting task $folder.."
            mkdir "${folder}/${resDirName}/" 
            # process as abundance
            # run INLA sim analysis
            Rscript spiecEasi_simAnalysis.R ${folder}/ ${folder}/${resDirName}/ ${seMethod} ${numTrials}
            # this ecoCopula one should work I think
            Rscript /home/fiona_callahan/eDNA_sims_code/countEcoCopulaMistakes.R ${folder}/ ${folder}/${resDirName}/ 
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

# this ecoCopula one should still work
Rscript /home/fiona_callahan/eDNA_sims_code/gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName}

echo "all done"
