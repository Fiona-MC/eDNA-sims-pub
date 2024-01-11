#!/bin/bash
export OMP_NUM_THREADS=5

#sim_dir=$1
sim_dir="/space/s1/fiona_callahan/multiSim_100sp"
numRuns=2
#numRuns=$2
numTrials=1
INLA_type="paper"
resDirName=INLA_res_${INLA_type}
scramble=0
covs=1
#sitetab="sim_sitetab_sampled.csv"
sitetab="sim_sitetab_sampled.csv"
ROC_mode="noModelSelect" # this will mean there is no WAIC selection for the ones where the cutoff changes


folder=$sim_dir/randomRun1

echo "starting task $folder..: $(date)"
mkdir "$folder/INLA_res_${INLA_type}/" 
# run INLA sim analysis
Rscript /home/fiona_callahan/eDNA_sims_code/INLA_simAnalysis_${INLA_type}.R ${folder}/ ${folder}/${resDirName}/ $sitetab
#for cutoff in 0.01;
#for cutoff in 0 0.0000000000001 0.0000001 0.00001 .3 .5 .7 .9 1;
for cutoff in 0 1 0.0000001 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 .3 .5;
do
    Rscript /home/fiona_callahan/eDNA_sims_code/INLA_changeCutoffs.R ${folder}/ ${cutoff} ${folder}/${resDirName}/ ${ROC_mode}
    Rscript /home/fiona_callahan/eDNA_sims_code/count_mistakes_general.R ${folder}/ ${folder}/${resDirName}/ ${covs} ${cutoff}
done

echo "finished task $folder..: $(date)"