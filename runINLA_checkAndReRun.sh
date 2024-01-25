sim_dir=$1
resDirName=INLA_res_paperSep_sampled100
for folder in ${sim_dir}/randomRun*; do
    for modelParms in none cov sp spCov; do
        thisRes=${folder}/${resDirName}/trial1/resList_${modelParms}.Rdata
        if ! test -f ${thisRes}; then
        echo ${thisRes} does not exist
        # run INLA sim analysis
        #timeout -k 10 10h Rscript /home/fiona_callahan/eDNA_sims_code/INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ ${sitetab} ${modelParms}
        #Rscript INLA_modelSelect.R ${folder}/ ${folder}/${resDirName}/
        fi
    done
done