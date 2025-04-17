# Rscript make_parm_csv.R ${fullDir}/randomRun${run}/paramsFiltered${nsamp}.Rdata ${targetDir}/randomRun${run}/
# Rscript make_parm_csv.R ${fullDir}/randomRun${run}/params.Rdata ${targetDir}/randomRun${run}/

# Rscript make_parm_csv.R ${infile} ${outdir}

# fullDir="/space/s1/fiona_callahan/sim_paper_stuff/multiSim_10sp"
# targetDir="/space/s1/fiona_callahan/sim_paper_stuff/clean_data/multiSim_10sp"

# fullDir="/space/s1/fiona_callahan/sim_paper_stuff/multiSim_10sp_revision2"
# targetDir="/space/s1/fiona_callahan/sim_paper_stuff/clean_data/multiSim_10sp_revision2"

library(jsonlite)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outdir <- args[2]

#infile <- "/space/s1/fiona_callahan/sim_paper_stuff/multiSim_10sp_revision2/randomRun1/params.Rdata"
#outdir <- "/space/s1/fiona_callahan/sim_paper_stuff/clean_data/multiSim_10sp_revision2/randomRun1/"

if (!dir.exists(outdir)) {
    dir.create(outdir)
}

splitPath <- strsplit(infile, "/")
fileName <- splitPath[[1]][length(splitPath[[1]])]

params <- readRDS(infile)

if (str_detect(infile, "Filtered")) { # filtered in name
    filtered <- TRUE
} else {
    filtered <- FALSE
}

if (str_detect(infile, "revision")) { 
    base_path <- paste(splitPath[[1]][- length(splitPath[[1]])], collapse = "/")
    covar_parms <- readRDS(paste0(base_path, "/rev_covar_mx.Rdata"))
    params$betaCovar <- covar_parms$beta
    params$covarMx <- covar_parms$covarMx
    params$names_all_cov <- params$names_cov
}

write_json(params, path = paste0(outdir, strsplit(fileName, "\\.")[[1]][1], ".json"))

# params2 <- read_json(path = paste0(outdir, strsplit(fileName, "\\.")[[1]][1], ".json"), simplifyVector = TRUE)

if (filtered) {
    alpha_df <- as.data.frame(params$filteredAlpha)
    names(alpha_df) <- params$filteredSpNames
    rownames(alpha_df) <- params$filteredSpNames
    write.csv(x = alpha_df, file = paste0(outdir, "alpha_", strsplit(fileName, "\\.")[[1]][1], ".csv"))

    beta_df <- as.data.frame(params$filteredBeta)
    names(beta_df) <- params$names_all_cov
    rownames(beta_df) <- params$filteredSpNames
    write.csv(x = beta_df, file = paste0(outdir, "beta_", strsplit(fileName, "\\.")[[1]][1], ".csv"))
} else {
    alpha_df <- as.data.frame(params$alpha)
    names(alpha_df) <- params$names_species
    rownames(alpha_df) <- params$names_species
    write.csv(x = alpha_df, file = paste0(outdir, "alpha", ".csv"))

    beta_df <- as.data.frame(params$beta)
    names(beta_df) <- params$names_all_cov
    rownames(beta_df) <- params$names_species
    write.csv(x = beta_df, file = paste0(outdir, "beta", ".csv"))
}
