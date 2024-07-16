#!/bin/bash
export OMP_NUM_THREADS=5

# ./run_ecoCopulaSimAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 0

#sim_dir="/space/s1/fiona_callahan/multiSim_100"
#numRuns=100
#covs=0

sim_dir=$1
numRuns=$2
covs=$3
numSamples=$4
readAbd=$5
logi=$6
filtered=$7

numTrials=1 # I think as this is implemented right now this needs to be 1
scramble=0

if [ ${numSamples} == "None" ]
then
    if [ ${readAbd} == 1 ]
    then
        resDirName=ecoCopula_readAbd_res
        sitetab_name="sim_sitetab_readAbd_sampled.csv"
    else 
        resDirName=ecoCopula_res
        sitetab_name="sim_sitetab_sampled.csv"
    fi
else
    if [ ${readAbd} == 1 ]
    then
        resDirName=ecoCopula_res_readAbd_sampled${numSamples}
        sitetab_name=sim_sitetab_readAbd_sampled${numSamples}.csv
    else 
        resDirName=ecoCopula_res_sampled${numSamples}
        sitetab_name=sim_sitetab_sampled${numSamples}.csv
    fi
fi

#cutoff=NA

if [ ${covs} == 1 ]
then
	resDirName=${resDirName}_cov
else
	resDirName=${resDirName}_noCov
fi

if [ ${logi} == 1 ]
then
    if [ ${readAbd} == 1 ]
    then
        sitetab_name=logiSim_sitetab_readAbd_sampled${numSamples}.csv
    else 
        sitetab_name=logiSim_sitetab_sampled${numSamples}.csv
    fi
    resDirName=${resDirName}_logi
fi

if [ ${filtered} == 1 ]
then
sitetab_name=sim_sitetab_sampled${numSamples}_filtered.csv
resDirName=${resDirName}_filtered
fi


N=3 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)


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
            mkdir "$folder/${resDirName}/" 
            # run sim analysis
            # 0 in here is
            Rscript ecoCopula_simAnalysis.R ${folder}/ ${folder}/${resDirName}/ ${scramble} ${sitetab_name} ${covs}
            
            for cutoff in 0 1 12 23 34 45 56 67 78 89 100;
            do
                Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ 0 ${cutoff} #recall that even if covs is true, beta is not inferred so this 0 is right
            done

            # covs here should be 0 because this is whether it infers the covariates, whereas above it just controls for them
            Rscript ./count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ 0 #recall that even if covs is true, beta is not inferred so this 0 is right

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

for cutoff in 0 1 12 23 34 45 56 67 78 89 100;
do
Rscript ./gather_inferenceRes_general.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName} ${cutoff}
done

Rscript ./gather_inferenceRes_general.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName}

echo "all done"
