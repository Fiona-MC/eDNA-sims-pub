
numSpecies=10
numRuns=100
numSamples=10000

# /space/s1/KaiyuanLi/data/JAGS/results/filtered/sp10_100_filtered/mistakes.csv

cp /space/s1/KaiyuanLi/data/JAGS/results/filtered/sp${numSpecies}_${numSamples}_filtered/mistakes.csv /space/s1/fiona_callahan/multiSim_${numSpecies}sp/JAGS_infResGathered_sampled${numSamples}_filtered_${numRuns}sims.csv


# for random
numSpecies=10
numRuns=100
numSamples=10000

# /space/s1/KaiyuanLi/data/JAGS/results/filtered/sp10_100_filtered/mistakes.csv

cp /space/s1/KaiyuanLi/data/JAGS/results/filtered/sp${numSpecies}_${numSamples}_rp/mistakes.csv /space/s1/fiona_callahan/multiSim_${numSpecies}sp_random_moreSamples/JAGS_infResGathered_sampled${numSamples}_filtered_${numRuns}sims.csv




# for logi
numSpecies=10
numRuns=100
numSamples=100

# /space/s1/KaiyuanLi/data/JAGS/results/filtered/sp10_100_filtered/mistakes.csv

cp /space/s1/KaiyuanLi/data/JAGS/results/splogi${numSpecies}_${numSamples}/mistakes.csv /space/s1/fiona_callahan/multiSim_${numSpecies}sp/JAGS_infResGathered_sampled${numSamples}_logi_${numRuns}sims.csv