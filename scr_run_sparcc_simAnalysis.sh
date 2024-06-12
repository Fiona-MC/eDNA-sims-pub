#!/bin/bash
export OMP_NUM_THREADS=5

# ./scr_run_sparcc_simAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 100 0 1


sim_dir=$1
numRuns=$2
numSamples=$3
logi=$4
filtered=$5 #100

if [ ${numSamples} == "None" ]
then
	resDirName=sparcc_res
    sitetab_name=sim_sitetab_readAbd_sampled.csv
else
	resDirName=sparcc_res_sampled${numSamples}
    sitetab_name=sim_sitetab_readAbd_sampled${numSamples}.csv
fi


numTrials=1 # I think as this is implemented right now this needs to be 1
#seMethod=sparcc
#seMethod=glasso

plot=1
covs=0 # not a real option here

if [ ${logi} == 1 ]
then
    sitetab_name=logiSim_sitetab_readAbd_sampled${numSamples}.csv
    resDirName=${resDirName}_logi
fi

if [ ${filtered} == 1 ]
then
    sitetab_name=sim_sitetab_readAbd_sampled${numSamples}_filtered.csv
    resDirName=${resDirName}_filtered
fi


N=1 
folderNames=()

# Populate the array with the names
for ((i=1; i<=$numRuns; i++)); do
  folderNames+=(${sim_dir}/randomRun$i)
done

for folder in ${folderNames[@]}; do
    (
        #if test ! -d "${folder}/${resDirName}/trial1" # if the folder is not already there NOT WORKING
        #then
            echo "starting task $folder.."
            mkdir "${folder}/${resDirName}/" 
            # process as abundance
            # run INLA sim analysis
            #for cutoff in 0.0000000000000001 0.00000001 0.00001 0.0001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.9999 1 0.005 0.04 0.06 0.08 0.125 0.15 0.175 0.25;
            for cutoff in pval_bootstrap;
            do
            Rscript ./sparcc_simAnalysis_pval.R ${folder}/ ${folder}/${resDirName}/ ${cutoff} ${numTrials} ${sitetab_name} 10
            Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ ${covs} ${cutoff}
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

# this ecoCopula one should still work
#for cutoff in 0.0000000000000001 0.00000001 0.00001 0.0001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.9999 1 0.005 0.04 0.06 0.08 0.125 0.15 0.175 0.25;
for cutoff in pval_bootstrap;
do
Rscript ./gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName} ${cutoff}
done

echo "all done"

