export OMP_NUM_THREADS=3


sim_dir=/space/s1/fiona_callahan/multiSim_10sp_random_moreSamples
numRuns=100 


numTrials=1

resDirName="/space/s1/KaiyuanLi/data/JAGS/results/filtered/sp10_100_rp"





sitetab='sim_sitetab_sampled100_filtered.csv'


N=10
for num in {1..100} 
do
    (
        #if test ! -d "${folder}/INLA_res_${INLA_type}/trial1" # if the folder is not already there 
        #then
            folder=${sim_dir}/randomRun${num} #no space after folder!Otherwise it's a function
            echo "starting task $folder.."
            mkdir "$resDirName/$num/" 
            # run INLA sim analysis
            timeout -k 10 10d Rscript /home/KaiyuanLi/current_code/JAGS/sim_data/sp10_ran.r ${folder}/ ${resDirName}/${num}/ ${sitetab}
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



echo "all done"