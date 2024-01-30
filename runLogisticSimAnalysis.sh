#!/bin/bash
export OMP_NUM_THREADS=15

# ./runLogisticSimAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 0

#sim_dir="/space/s1/fiona_callahan/multiSim_100"
#numRuns=100
#covs=0

sim_dir=$1
numRuns=$2
covs=$3
numSamples=$4

echo "Starting logistic"
echo $sim_dir
echo $numRuns
echo $covs
echo $numSamples

if [ ${numSamples} == "None" ]
then
	sitetab="sim_sitetab_sampled.csv"
	outname="logistic_mistakes"
else
	sitetab=sim_sitetab_sampled${numSamples}.csv
	outname=logistic_mistakes_sampled${numSamples}
fi

dumb=0

if [ ${covs} == 1 ]
then
	outname=${outname}_cov
else
	outname=${outname}_noCov
fi

if [ ${dumb} == 1 ]
then
	outname=${outname}_dumb
fi

outname=${outname}_${numRuns}runs

Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim_moreSp.R ${sim_dir}/ ${numRuns} ${covs} ${dumb} ${sitetab} ${outname}

echo "all done"
