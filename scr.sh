export OMP_NUM_THREADS=5

sim_dir=$1
#sim_dir="/space/s1/fiona_callahan/multiSim_100sp"
#numRuns=100
numRuns=$2
numTrials=1
INLA_type="paper"
resDirName=INLA_res_${INLA_type}
scramble=0
covs=1
#sitetab="sim_sitetab_sampled.csv"
sitetab="sim_sitetab_sampled.csv"
ROC_mode="noModelSelect" 
echo $sim_dir
echo $numRuns