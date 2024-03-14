dir=/space/s1/fiona_callahan

for folder in ${dir}/multiSim_10sp_random ${dir}/multiSim_50sp_random ${dir}/multiSim_100sp_random; do
    for run in {1..100}; do
        chmod -w $folder/randomRun${run}/params.Rdata
        chmod -w $folder/randomRun${run}/locList.Rdata
        chmod -w $folder/randomRun${run}/covList.Rdata
        chmod -w $folder/randomRun${run}/sim_sitetab_readAbd_sampled.csv
        chmod -w $folder/randomRun${run}/sim_sitetab_sampled.csv
        #for numSamples in 50 100 500 1000 10000 25000; do
        for numSamples in 100; do
            file1=$folder/randomRun${run}/sim_sitetab_readAbd_sampled${numSamples}.csv
            file2=$folder/randomRun${run}/sim_sitetab_sampled${numSamples}.csv
            chmod -w $file1
            chmod -w $file2
        done
    done
done