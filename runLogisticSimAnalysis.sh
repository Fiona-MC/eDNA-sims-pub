#!/bin/bash
export OMP_NUM_THREADS=15

# ./runLogisticSimAnalysis.sh /space/s1/fiona_callahan/multiSim_10sp 100 0

#sim_dir="/space/s1/fiona_callahan/multiSim_100"
#numRuns=100
#covs=0

sim_dir=$1
numRuns=$2
covs=$3

sitetab_name="sim_sitetab_sampled.csv"
dumb=0

if [ ${covs} == 1 ]
then
	outname="logistic_mistakes_cov"
else
	outname="logistic_mistakes_noCov"
fi

if [ ${dumb} == 1 ]
then
	outname=${outname}_dumb
fi

if [ ${sitetab_name} == "sim_sitetab_sampled500.csv" ]
then
    outname=${outname}_500
fi

Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim_moreSp.R ${sim_dir}/ ${numRuns} ${covs} ${dumb} ${sitetab_name} ${outname}

echo "all done"
