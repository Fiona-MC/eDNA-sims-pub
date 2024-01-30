#!/bin/bash
export OMP_NUM_THREADS=15

# ./run_all.sh /space/s1/fiona_callahan/multiSim_10xtest 2

sim_dir=$1
numRuns=$2

#cd /home/fiona_callahan/eDNA_sims_code

# where are they going when I send them to the background?? 

#for numSamples in 50 100 500 1000 10000 25000;
for numSamples in 100;
do
./runINLAsimAnalysis.sh $1 $2 1 $numSamples > ${sim_dir}/inlaCov${numSamples}Out.log 2> ${sim_dir}/inlaCov${numSamples}Err.log &
#./runINLAsimAnalysis.sh $1 $2 0 $numSamples > ${sim_dir}/inla${numSamples}Out.log 2> ${sim_dir}/inla${numSamples}Err.log &

./runLogisticSimAnalysis.sh $1 $2 1 $numSamples > ${sim_dir}/logisticCov${numSamples}Out.log 2> ${sim_dir}/logisticCov${numSamples}Err.log &
./runLogisticSimAnalysis.sh $1 $2 0 $numSamples > ${sim_dir}/logistic${numSamples}Out.log 2> ${sim_dir}/logistic${numSamples}Err.log &

./run_ecoCopulaSimAnalysis.sh $1 $2 1 $numSamples > ${sim_dir}/ecCov${numSamples}Out.log 2> ${sim_dir}/ecCov${numSamples}Err.log &
./run_ecoCopulaSimAnalysis.sh $1 $2 0 $numSamples > ${sim_dir}/ec${numSamples}Out.log 2> ${sim_dir}/ec${numSamples}Err.log &

./run_spiecEasiSimAnalysis.sh $1 $2 mb $numSamples > ${sim_dir}/seMb${numSamples}Out.log 2> ${sim_dir}/seMb${numSamples}Err.log &
./run_spiecEasiSimAnalysis.sh $1 $2 glasso $numSamples > ${sim_dir}/seGlasso${numSamples}Out.log 2> ${sim_dir}/seGlasso${numSamples}Err.log &

./run_spiecEasiSimAnalysis.sh $1 $2 sparcc $numSamples > ${sim_dir}/seSparcc${numSamples}Out.log 2> ${sim_dir}/seSparcc${numSamples}Err.log &
done



