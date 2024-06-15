#sim_dir=$1
numRuns=100
#covs=$3
#numSamples=$4
#readAbd=$5
#logi=$6
#filtered=$7

for numSamples in 100 10000;
do
    for sim_dir in /space/s1/fiona_callahan/multiSim_10sp /space/s1/fiona_callahan/multiSim_100sp /space/s1/fiona_callahan/multiSim_10sp_random_moreSamples /space/s1/fiona_callahan/multiSim_100sp_random_moreSamples;
    do
        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 0 1
        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 0 1

        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 0 1
        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 0 1
    done
done


#logi
for numSamples in 100 10000;
do
    for sim_dir in /space/s1/fiona_callahan/multiSim_10sp /space/s1/fiona_callahan/multiSim_100sp;
    do
        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 0 1 0
        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 0 1 0

        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 1 ${numSamples} 1 1 0
        ./run_ecoCopulaSimAnalysis.sh ${sim_dir} ${numRuns} 0 ${numSamples} 1 1 0
    done
done