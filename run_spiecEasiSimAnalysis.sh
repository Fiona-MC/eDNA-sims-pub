#!/bin/bash
export OMP_NUM_THREADS=5

# ./run_spiecEasiSimAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 10 glasso 100 0 1

#sim_dir="/space/s1/fiona_callahan/multiSim_100"
#numRuns=100
#seMethod=mb #glasso sparcc

sim_dir=$1
numRuns=$2
seMethod=$3
numSamples=$4
logi=$5
filtered=$6 #100

if [ ${numSamples} == "None" ]
then
	resDirName=spiecEasi_res
    sitetab_name=sim_sitetab_readAbd_sampled.csv
else
	resDirName=spiecEasi_res_sampled${numSamples}
    sitetab_name=sim_sitetab_readAbd_sampled${numSamples}.csv
fi


numTrials=1 # I think as this is implemented right now this needs to be 1
#seMethod=sparcc
#seMethod=glasso
random=1
plot=1
covs=0 # not a real option here

resDirName=${resDirName}_${seMethod}

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

N=1 # 

# Declare an array to store the names
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
            #mkdir "${folder}/${resDirName}/" 
            # process as abundance
            # run INLA sim analysis
            #Rscript spiecEasi_simAnalysis.R ${folder}/ ${folder}/${resDirName}/ ${seMethod} ${numTrials} ${sitetab_name}
            # this ecoCopula one should work I think
            Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ ${covs}
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
Rscript ./gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName}

echo "all done"

