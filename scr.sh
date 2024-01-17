#!/bin/bash
export OMP_NUM_THREADS=5

#./scr.sh /space/s1/fiona_callahan/multiSim_100sp

sim_dir=$1

numTrials=1
INLA_type="paper"

scramble=0

ROC_mode="noModelSelect" # this will mean there is no WAIC selection for the ones where the cutoff changes

folder=$sim_dir/randomRun2
numSamples=500
sitetab=sim_sitetab_sampled${numSamples}.csv
resDirName=INLA_test_2
echo "starting task $folder.."
mkdir "$folder/$resDirName/" 

Rscript /home/fiona_callahan/eDNA_sims_code/INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ $sitetab
      