#!/bin/bash
export OMP_NUM_THREADS=15

# ./runLogisticSimAnalysis_filtered.sh /space/s1/fiona_callahan/sim_paper_stuff/multiSim_100sp_revision1 100 1 100 1 0

sim_dir="/space/s1/fiona_callahan/sim_paper_stuff/multiSim_100sp_revision1"
numRuns=100
covs=1
numSamples=10000
logi=1
filtered=0

sim_dir=$1
numRuns=$2
covs=$3
numSamples=$4
logi=$5
filtered=$6

echo "Starting logistic"
echo $sim_dir
echo $numRuns
echo $covs
echo $numSamples
echo logiSim? $logi

if [ ${numSamples} == "None" ]
then
	sitetab="sim_sitetab_sampled.csv"
	outname="logistic_mistakes"
else
	sitetab=sim_sitetab_sampled${numSamples}.csv
	outname=logistic_mistakes_sampled${numSamples}
fi

if [ ${covs} == 1 ]
then
	outname=${outname}_cov
else
	outname=${outname}_noCov
fi

if [ ${logi} == 1 ]
then
	sitetab=logiSim_sitetab_sampled${numSamples}.csv
	outname=${outname}_logi
fi

if [ ${filtered} == 1 ]
then
	sitetab=sim_sitetab_sampled${numSamples}_filtered.csv
	outname=${outname}_filtered
fi

outname=${outname}_${numRuns}runs

echo ${outname}

Rscript ./logistic_filtered.R ${sim_dir}/ ${numRuns} ${covs} ${logi} ${sitetab} ${outname}

echo "all done"
