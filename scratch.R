thisMultiSimRes <- read.csv("/space/s1/fiona_callahan/multiSim_10sp/ecoCopula_res_readAbd_sampled100_noCov_filtered_infResGathered_100sims.csv")

source("~/eDNA_sims_code/confusion_stats.R")



names(multiSimResL)

for (i in 40:70) {
    name <- names(multiSimResL)[i]
    thisMultiSimRes <- multiSimResL[[name]]
    FPR <- get_FPR(thisMultiSimRes, mode = "cluster", return_components = TRUE)
    TPR <- get_TPR(thisMultiSimRes, mode = "cluster", return_components = TRUE)

    print(name)
    print(mean(FPR$FPR, na.rm = TRUE))
    print(mean(TPR$TPR, na.rm = TRUE))
}

inferenceRes0.1 <- readRDS("/space/s1/fiona_callahan/multiSim_100sp/randomRun13/sparcc_res_sampled10000_logi/trial1/inferenceRes_cutoff0.1.Rdata")

inferenceRes0.09 <- readRDS("/space/s1/fiona_callahan/multiSim_100sp/randomRun13/sparcc_res_sampled10000_logi/trial1/inferenceRes_cutoff0.09.Rdata")
