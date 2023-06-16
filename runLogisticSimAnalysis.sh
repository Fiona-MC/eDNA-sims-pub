#!/bin/bash
export OMP_NUM_THREADS=15
sim_dir="/space/s1/fiona_callahan/dumbSim3"
numRuns=100

Rscript /home/fiona_callahan/eDNA_sims_code/logisticFromSim.R ${sim_dir}/ ${numRuns}

echo "all done"
