# this is for after running runSchliepSimAnalysis.sh on a bunch of sims
# note: the making of the schliep_runNums.csv is in analyze_multiSim.R

dataDir <- "/space/s1/fiona_callahan/multiSim11/"

schliep_runNums <- read.csv(paste0(dataDir, "schliep_runNums.csv"), header = FALSE)

schliep_inferenceRes <- read.csv(paste0(dataDir, "schliep_inferenceRes_gathered.csv"), header = TRUE)
inla_multiSimRes <- read.csv(paste0(dataDir, "infResGathered.csv"), header = TRUE)

schliepResMerged <- merge(schliep_inferenceRes, inla_multiSimRes, by = "RunNum", all.x = TRUE, all.y = FALSE, suffixes = c("JSDM", "INLA"))

plot(schliepResMerged$totalMistakes_JSDM, schliepResMerged$totalMistakes_INLA, xlab = "totalMistakes_JSDM", ylab = "totalMistakes_INLA")


