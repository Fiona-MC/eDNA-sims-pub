#!/bin/bash

# think about whether you wanna delete these
# I think I may have overwritten a bunch of these with the 500 samples but I'm not 100% sure
sim_dir="/space/s1/fiona_callahan/multiSim_100"

for folder in ${sim_dir}/randomRun*; do
   rm -r $folder/INLA_res_paper/
done

