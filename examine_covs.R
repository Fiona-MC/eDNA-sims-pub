library(ggplot2)
library(gridExtra)
#TODO need to add this to the multiSimFunctions.R (in function makePlots)

runNum <- 30
thisDir <- paste0("/space/s1/fiona_callahan/multiSim11/randomRun", runNum,"/")
sitetab_longer <- read.csv(paste0(thisDir, "sitetab_longer.csv"))
head(sitetab_longer)

cov1plot <- ggplot(data = sitetab_longer[sitetab_longer$labID <= 8, ], aes(x = Age, y = Cov1, group = labID, color = labID)) +
    geom_point() +
    geom_line()
cov2plot <- ggplot(data = sitetab_longer[sitetab_longer$labID <= 8, ], aes(x = Age, y = Cov2, group = labID, color = labID)) +
    geom_point() +
    geom_line()
cov3plot <- ggplot(data = sitetab_longer[sitetab_longer$labID <= 8, ], aes(x = Age, y = Cov3, group = labID, color = labID)) +
    geom_point() +
    geom_line()
cov4plot <- ggplot(data = sitetab_longer[sitetab_longer$labID <= 8, ], aes(x = Age, y = Cov4, group = labID, color = labID)) +
    geom_point() +
    geom_line()

covplotsall <- arrangeGrob(cov1plot, cov2plot, cov3plot, cov4plot)
ggsave(file = "/space/s1/fiona_callahan/multiSim11/covPlots.png", covplotsall, width = 10, height = 8)



sitetab_sampled <- read.csv(paste0(thisDir, "sim_sitetab_sampled.csv"))

covPlotsL <- list()
for(cov in 1:4){
    covPlotsL[[cov]] <- ggplot(data = sitetab_sampled[sitetab_sampled$labID <= 8, ], aes(x = Age, y = !!sym(paste0("Cov", cov)), group = labID, color = labID)) + # nolint: line_length_linter.
        geom_point() +
        geom_line()
}
covplots <- arrangeGrob(grobs = covPlotsL)
ggsave(file = paste0("/space/s1/fiona_callahan/multiSim11/covPlots_sampled2_", runNum, ".png"), covplots, width = 10, height = 8)
