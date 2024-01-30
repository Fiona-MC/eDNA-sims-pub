sim_dir=$1
resDirName=$2 # INLA_res_paperSep_sampled100
numRuns=$3 # 100
reRun=$4
timeout=$5

folderNames=() # array of names of randomRun* so that diff numbers of runs can be done 
for ((i=1; i<=$numRuns; i++)); do
  folderNames+=(${sim_dir}/randomRun$i)
done

for folder in ${folderNames[@]}; do
    for modelParms in none cov sp spCov; do
        thisRes=${folder}/${resDirName}/trial1/resList_${modelParms}.Rdata
        if ! test -f ${thisRes}; then
            echo ${thisRes} does not exist
            if [ ${reRun}==1 ]; then
            # run INLA sim analysis
            echo rerunning ${thisRes}...
            timeout -k 10 ${timeout}h Rscript /home/fiona_callahan/eDNA_sims_code/INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ ${sitetab} ${modelParms}
            fi
        fi
    done
    echo done with ${thisRes}.
    if [ ${reRun}==1 ]; then
    Rscript INLA_modelSelect.R ${folder}/ ${folder}/${resDirName}/
    fi
done

echo all done