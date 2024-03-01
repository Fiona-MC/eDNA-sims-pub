sim_dir=/space/s1/fiona_callahan/multiSim_10sp
numRuns=10
numTrials=1

for numSamples in 50 100 500 1000 5000 10000 25000; do
    #for resDirName in INLA_res_paperSep_sampled${numSamples} INLA_res_paperSep_sampled${numSamples}_cov INLA_res_paperSep_sampled${numSamples}_noCov ecoCopula_res_readAbd_sampled${numSamples}_cov ecoCopula_res_readAbd_sampled${numSamples}_noCov ecoCopula_res_sampled${numSamples}_cov ecoCopula_res_sampled${numSamples}_noCov spiecEasi_res_sampled${numSamples}_mb spiecEasi_res_sampled${numSamples}_glasso spiecEasi_res_sampled${numSamples}_sparcc; do
        echo $resDirName
        Rscript ./gather_inferenceRes_ecoCopula.R ${sim_dir}/ ${numRuns} ${numTrials} ${resDirName}
    done
done