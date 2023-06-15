#!/bin/bash
export OMP_NUM_THREADS=15
sim_dir="/space/s1/fiona_callahan/multiSim_ParmSet5"

numRuns=100
numTrials=1

#N=10 # N=10 resulted in average usage around 30 cores
# based on current rate with N=10 -- this should take ~6 days for 1000 runs (2 trials each)

Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim.R ${sim_dir}/ ${numRuns}

echo "all done"
