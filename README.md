# eDNA-sims

## Simulations
### Covariance Matrix simulation

Code for the covariance matrix simulations is in logiSim_covMx.R

It is set up to simulate interactions as they are set in the file params.Rdata, which can be generated by simulating data using the ecological simulation data.

The code to generate a list of parameters can be found in multiSimParms.R

(For a more parameter-independent version of these simulations, see logiSim_covMx_higherCorr.r)

Then, to create final sampled datasets: reSampleFromSitetab.R will re-sample randomly from the full simulated data.

### Ecological simulations

Main code for the ecological simulations can be found in multiSim_Mar2023.R

To run 100 simulations, with random parameters and 10 species: 
```
Rscript ./multiSim_Mar2023.R /path/to/output/directory/ 1 100 1 10
```

This depends on functions in the following files:

multiSimFunctions.R
multiSimParms.R

This also depends on R libraries:

library(ggplot2)
library(reshape)
library(gridExtra)
library(gapminder)
library(data.table)
library(fields)
library(MASS)

After running this, reSampleFromSitetab.R will re-sample randomly from the full simulated data.

Then, we filter for 10% presence of each species using filter_data_extinctions.R

This will output the final simulated data.

## Analyses

### Pipelines 
runLinearSimAnalysis_filtered.sh

runLogisticSimAnalysis_filtered.sh

runINLAsimAnalysis.sh (this analysis depends on scripts from https://github.com/wyc661217/Arctic_eDNA_2021)

Additonal note about INLA: For parallelized analysis using slurm scripts: 
run_INLA_oneRun.sh for each simulation run, followed by run_INLA_gatheringRes.sh

run_sparcc_simAnalysis.sh

run_spiecEasiSimAnalysis.sh

run_ecoCopulaSimAnalysis.sh


All of these files will output a csv file with all of the information to make ROC curves or get FDR.

### Curves and plots
get_ROC_stats.R will make a csv file with the statistics for all methods to make a ROC curve. 

get_FDRs.r will make tables and plots with statistics for false discovery rates and overall discovery rates

Both depend on confusion_stats.R

### Influence of specific parameters

random_forest_analysis.R

### Miscellaneous

check_nSpecies.R was used to delete simulations where too many species went extinct or were present at very low frequency.

compare_percent_presence.R was used to make figures showing percent presence of each species split by whether they had (inferred or actual) interactions.

configure_FDR_figs.R and configure_ROC_figs.R were used to make figures in the paper with multiple panels.

se_inferences_vs_samples.R was used to make supplemental figure about the number of inferences made by spiecEasi

