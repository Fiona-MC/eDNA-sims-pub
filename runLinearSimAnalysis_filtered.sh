#!/bin/bash
export OMP_NUM_THREADS=15

# ./runLinearSimAnalysis_filtered.sh /space/s1/fiona_callahan/sim_paper_stuff/multiSim_100sp_revision1 100 0 10000 1 0

sim_dir="/space/s1/fiona_callahan/sim_paper_stuff/multiSim_100sp_revision3"
numRuns=100
covs=0
numSamples=250
logi=1
filtered=0

sim_dir=$1
numRuns=$2
covs=$3
numSamples=$4
logi=$5
filtered=$6

echo "Starting linear regression"
echo $sim_dir
echo $numRuns
echo $covs
echo $numSamples
echo logiSim? $logi

if [ ${numSamples} == "None" ]
then
	sitetab="sim_sitetab_readAbd_sampled.csv"
	outname="linearReg_mistakes"
else
	sitetab=sim_sitetab_readAbd_sampled${numSamples}.csv
	outname=linearReg_mistakes_sampled${numSamples}
fi

if [ ${covs} == 1 ]
then
	outname=${outname}_cov
else
	outname=${outname}_noCov
fi

if [ ${logi} == 1 ]
then
	sitetab=logiSim_sitetab_readAbd_sampled${numSamples}.csv
	outname=${outname}_logi
fi

if [ ${filtered} == 1 ]
then
	sitetab=sim_sitetab_readAbd_sampled${numSamples}_filtered.csv
	outname=${outname}_filtered
fi

outname=${outname}_${numRuns}runs

echo ${outname}

Rscript ./linearRegression_simAnalysis.R ${sim_dir}/ ${numRuns} ${covs} ${logi} ${sitetab} ${outname}

echo "all done"
