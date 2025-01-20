
nSpecies=100
simDir=/space/s1/fiona_callahan/sim_paper_stuff/multiSim_${nSpecies}sp_test_revision3
numRuns=100

mkdir ${simDir}
for run in $(seq 1 ${numRuns}); do
    mkdir ${simDir}/randomRun${run}/
    Rscript /home/fiona_callahan/eDNA_sims_code/revisions_code/logiSim_covMx_twostep.R ${simDir}/randomRun${run}/ ${nSpecies}
done

Rscript /home/fiona_callahan/eDNA_sims_code/reSampleFromSitetab.R ${simDir}/ ${simDir}/ ${numRuns}