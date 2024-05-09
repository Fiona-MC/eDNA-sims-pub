# meant to be copy-pasted into savio or ponderosa command line -- dont run this as a script

# savio to ponderosa
folder=multiSim_10sp
subfolder1=INLA_res_paperSep_sampled100_logi_noCov
subfolder2=INLA_res_paperSep_sampled100_logi_cov
subfolder3=INLA_res_paperSep_sampled100_logi
subfolder4=INLA_res_paperSep_sampled100_logi_covNoCount
numRuns=100

for ((i = 1; i <= numRuns; i++)); do

scp -r -i ~/.ssh/id_rsa ./${folder}/randomRun${i}/${subfolder1} fiona_callahan@ponderosa.biol.berkeley.edu:/space/s1/fiona_callahan/${folder}/randomRun${i}/
scp -r -i ~/.ssh/id_rsa ./${folder}/randomRun${i}/${subfolder2} fiona_callahan@ponderosa.biol.berkeley.edu:/space/s1/fiona_callahan/${folder}/randomRun${i}/
scp -r -i ~/.ssh/id_rsa ./${folder}/randomRun${i}/${subfolder3} fiona_callahan@ponderosa.biol.berkeley.edu:/space/s1/fiona_callahan/${folder}/randomRun${i}/
scp -r -i ~/.ssh/id_rsa ./${folder}/randomRun${i}/${subfolder4} fiona_callahan@ponderosa.biol.berkeley.edu:/space/s1/fiona_callahan/${folder}/randomRun${i}/

done


# ponderosa to savio -- doesnt work to make shortcuts because of entering passwd 100 times

scp -r ./multiSim_10sp fionacallahan@dtn.brc.berkeley.edu:/global/scratch/users/fionacallahan