#!/bin/bash
#SBATCH --job-name=INLA_test1
#SBATCH --partition=savio3_htc
#SBATCH --account=fc_popgen
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=2
#SBATCH --time=72:00:00
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=fiona_callahan@berkeley.edu
cd /global/scratch/users/fiona_callahan/eDNA-sims
module load r
module load rgdal
ht_helper.sh -t /global/scratch/users/fiona_callahan/eDNA-sims/INLA_test_tasks/INLA_test1.tasks