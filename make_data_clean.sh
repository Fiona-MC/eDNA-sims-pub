nsp=10

# fixed parameters ecological sim
fullDir=/space/s1/fiona_callahan/sim_paper_stuff/multiSim_${nsp}sp
targetDir=/space/s1/fiona_callahan/sim_paper_stuff/clean_data/multiSim_${nsp}sp

for run in {1..100};
do
    for nsamp in 100 250 10000;
    do
        # Rscript make_parm_csv.R ${infile} ${outdir}
        Rscript /home/fiona_callahan/eDNA_sims_code/make_parm_csv.R ${fullDir}/randomRun${run}/paramsFiltered${nsamp}.Rdata ${targetDir}/randomRun${run}/

        cp ${fullDir}/randomRun${run}/sim_sitetab_readAbd_sampled${nsamp}_filtered.csv ${targetDir}/randomRun${run}/sim_sitetab_readAbd_sampled${nsamp}_filtered.csv
        cp ${fullDir}/randomRun${run}/sim_sitetab_sampled${nsamp}_filtered.csv ${targetDir}/randomRun${run}/sim_sitetab_sampled${nsamp}_filtered.csv
    done
    Rscript /home/fiona_callahan/eDNA_sims_code/make_parm_csv.R ${fullDir}/randomRun${run}/params.Rdata ${targetDir}/randomRun${run}/
done

# random parameters ecological sim
fullDir=/space/s1/fiona_callahan/sim_paper_stuff/multiSim_${nsp}sp_random_moreSamples
targetDir=/space/s1/fiona_callahan/sim_paper_stuff/clean_data/multiSim_${nsp}sp_random

for run in {1..100};
do
    for nsamp in 100 250 10000;
    do
        # Rscript make_parm_csv.R ${infile} ${outdir}
        Rscript /home/fiona_callahan/eDNA_sims_code/make_parm_csv.R ${fullDir}/randomRun${run}/paramsFiltered${nsamp}.Rdata ${targetDir}/randomRun${run}/

        cp ${fullDir}/randomRun${run}/sim_sitetab_readAbd_sampled${nsamp}_filtered.csv ${targetDir}/randomRun${run}/sim_sitetab_readAbd_sampled${nsamp}_filtered.csv
        cp ${fullDir}/randomRun${run}/sim_sitetab_sampled${nsamp}_filtered.csv ${targetDir}/randomRun${run}/sim_sitetab_sampled${nsamp}_filtered.csv
    done
    Rscript /home/fiona_callahan/eDNA_sims_code/make_parm_csv.R ${fullDir}/randomRun${run}/params.Rdata ${targetDir}/randomRun${run}/
done

# revision 2 is covariance mx without covariates
fullDir=/space/s1/fiona_callahan/sim_paper_stuff/multiSim_${nsp}sp_revision2
targetDir=/space/s1/fiona_callahan/sim_paper_stuff/clean_data/multiSim_${nsp}sp_revision2

for run in {1..100};
do
    for nsamp in 100 250 10000;
    do
        cp ${fullDir}/randomRun${run}/logiSim_sitetab_readAbd_sampled${nsamp}.csv ${targetDir}/randomRun${run}/logiSim_sitetab_readAbd_sampled${nsamp}.csv
        cp ${fullDir}/randomRun${run}/logiSim_sitetab_sampled${nsamp}.csv ${targetDir}/randomRun${run}/logiSim_sitetab_sampled${nsamp}.csv
    done
    # Rscript make_parm_csv.R ${infile} ${outdir}
    Rscript /home/fiona_callahan/eDNA_sims_code/make_parm_csv.R ${fullDir}/randomRun${run}/params.Rdata ${targetDir}/randomRun${run}/

    cp ${fullDir}/randomRun${run}/sitetab_abd_logi.csv ${targetDir}/randomRun${run}/sitetab_abd_logi.csv
    cp ${fullDir}/randomRun${run}/sitetab_logi.csv ${targetDir}/randomRun${run}/sitetab_logi.csv
done


# revision 3 is covariance mx with covariates
fullDir=/space/s1/fiona_callahan/sim_paper_stuff/multiSim_${nsp}sp_revision3
targetDir=/space/s1/fiona_callahan/sim_paper_stuff/clean_data/multiSim_${nsp}sp_revision3

for run in {1..100};
do
    for nsamp in 100 250 10000;
    do
        cp ${fullDir}/randomRun${run}/logiSim_sitetab_readAbd_sampled${nsamp}.csv ${targetDir}/randomRun${run}/logiSim_sitetab_readAbd_sampled${nsamp}.csv
        cp ${fullDir}/randomRun${run}/logiSim_sitetab_sampled${nsamp}.csv ${targetDir}/randomRun${run}/logiSim_sitetab_sampled${nsamp}.csv
    done
    # Rscript make_parm_csv.R ${infile} ${outdir}
    Rscript /home/fiona_callahan/eDNA_sims_code/make_parm_csv.R ${fullDir}/randomRun${run}/params.Rdata ${targetDir}/randomRun${run}/

    cp ${fullDir}/randomRun${run}/sitetab_abd_logi.csv ${targetDir}/randomRun${run}/sitetab_abd_logi.csv
    cp ${fullDir}/randomRun${run}/sitetab_logi.csv ${targetDir}/randomRun${run}/sitetab_logi.csv
done

