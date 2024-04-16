
folder=multiSim_10sp
subfolder1=INLA_res_paperSep_sampled10000_filtered_noCov
subfolder2=INLA_res_paperSep_sampled10000_filtered_cov
subfolder3=INLA_res_paperSep_sampled10000_filtered
numRuns=100

for ((i = 1; i <= numRuns; i++)); do
scp -r -i ~/.ssh/id_rsa ./${folder}/randomRun${i}/${subfolder1} fiona_callahan@ponderosa.biol.berkeley.edu:/space/s1/fiona_callahan/${folder}/randomRun${i}/
scp -r -i ~/.ssh/id_rsa ./${folder}/randomRun${i}/${subfolder2} fiona_callahan@ponderosa.biol.berkeley.edu:/space/s1/fiona_callahan/${folder}/randomRun${i}/
scp -r -i ~/.ssh/id_rsa ./${folder}/randomRun${i}/${subfolder3} fiona_callahan@ponderosa.biol.berkeley.edu:/space/s1/fiona_callahan/${folder}/randomRun${i}/
done

