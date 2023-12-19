#!/bin/bash
export OMP_NUM_THREADS=5

sim_dir="/space/s1/fiona_callahan/multiSim_100"
numRuns=100
numTrials=1 # I think as this is implemented right now this needs to be 1
scramble=0
covs=0
#sitetab_name="sim_sitetab_sampled.csv"
sitetab_name="sim_sitetab_sampled500.csv"

if [ ${covs} == 1 ]
then
	resDirName="ecoCopula_res_cov500"
else
	resDirName="ecoCopula_res_noCov500"
fi
#INLA_type="faster"

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



N=3 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)

for folder in ${sim_dir}/randomRun*; do
    (
        #if test ! -d "${folder}/${resDirName}/trial1" # if the folder is not already there NOT WORKING
        #then
            echo "starting task $folder.."
            mkdir "$folder/${resDirName}/" 
            # run sim analysis
            Rscript ecoCopula_simAnalysis.R ${folder}/ ${folder}/${resDirName}/ 0 ${sitetab_name}
            Rscript /home/fiona_callahan/eDNA_sims_code/count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ ${covs}
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

Rscript /home/fiona_callahan/eDNA_sims_code/gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName}

echo "all done"
