library(stringr)
library(ggplot2)

simdir <- "multiSim_10sp_random_moreSamples"
nSamples <- 100

nInteract <- c()
spNames <- c()
percentPresence <- c()
runsL <- c()

for (run in 1:100) {

sitetab <- read.csv(paste0("/space/s1/fiona_callahan/", simdir, "/randomRun", run, "/sim_sitetab_sampled", nSamples, "_filtered.csv"))

parms <- readRDS(paste0("/space/s1/fiona_callahan/", simdir, "/randomRun", run, "/paramsFiltered", nSamples, ".Rdata"))

actualAlpha <- parms$filteredAlpha

nInteract <- c(nInteract, apply(actualAlpha, MARGIN = 1, FUN = sum))
spNames <- c(spNames, parms$filteredSpNames)
percentPresence <- c(percentPresence, apply(sitetab[, str_detect(names(sitetab), "Sp")], MARGIN = 2, FUN = sum) / nSamples)
runsL <- c(runsL, rep(run, times = length(parms$filteredSpNames)))
}

percentPresenceDF <- data.frame(Run = runsL, Species = spNames, 
                                NumInteractions = nInteract, percentPresence = percentPresence, 
                                interact_bool = (nInteract != 0))


ggplot(percentPresenceDF, aes(x = percentPresence, fill = interact_bool)) +
  geom_histogram(position = "dodge") +
    labs(title = simdir) +
  #scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  theme_minimal()
