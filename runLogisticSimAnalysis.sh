#!/bin/bash
export OMP_NUM_THREADS=15
sim_dir="/space/s1/fiona_callahan/multiSim_100"
numRuns=100
dumb=0

Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim_moreSp.R ${sim_dir}/ ${numRuns} ${dumb}

echo "all done"
