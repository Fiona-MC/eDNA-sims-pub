library(stringr)
library(ggplot2)

simdir <- "multiSim_100sp"
nSamples <- 100

nInteract <- c()
nInferredInteract <- c()
spNames <- c()
percentPresence <- c()
runsL <- c()

for (run in 1:100) {

sitetab <- read.csv(paste0("/space/s1/fiona_callahan/", simdir, "/randomRun", run, "/sim_sitetab_sampled", nSamples, "_filtered.csv"))

parms <- readRDS(paste0("/space/s1/fiona_callahan/", simdir, "/randomRun", run, "/paramsFiltered", nSamples, ".Rdata"))

actualAlpha <- parms$filteredAlpha

cutoff <- 0.05
inferredParms <- readRDS(paste0("/space/s1/fiona_callahan/", simdir, "/randomRun", run, "/logistic_mistakes_sampled100_noCov_filtered_100runs/inferenceRes_cutoff", cutoff, ".Rdata")) # nolint
alphaInferred <- inferredParms$alphaInferred
nInferredInteract <- c(nInferredInteract, apply(alphaInferred, MARGIN = 1, FUN = sum))

nInteract <- c(nInteract, apply(actualAlpha, MARGIN = 1, FUN = sum))
spNames <- c(spNames, parms$filteredSpNames)
percentPresence <- c(percentPresence, apply(sitetab[, str_detect(names(sitetab), "Sp")], MARGIN = 2, FUN = sum) / nSamples)
runsL <- c(runsL, rep(run, times = length(parms$filteredSpNames)))
if(length(nInferredInteract) != length(nInteract)) {
  print(run)
  print(length(nInferredInteract))
  print(length(nInteract))

}
}

percentPresenceDF <- data.frame(Run = runsL, Species = spNames, 
                                NumInteractions = nInteract, percentPresence = percentPresence, 
                                interact_bool = (nInteract != 0), 
                                nInferredInteract = nInferredInteract, inferredInteract_bool = (nInferredInteract != 0))


p <- ggplot(percentPresenceDF, aes(x = percentPresence, fill = interact_bool)) +
  geom_histogram(position = "dodge") +
    labs(title = simdir) +
    #scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    theme_minimal()
p 

p2 <- ggplot(percentPresenceDF, aes(x = percentPresence, fill = inferredInteract_bool)) +
  geom_histogram(position = "fill") +
    labs(title = simdir) +
    #scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
    theme_minimal()

p2
