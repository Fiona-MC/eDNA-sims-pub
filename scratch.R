args <- commandArgs(trailingOnly = TRUE)

print(as.numeric(args[1]))
print(typeof(as.numeric(args[1])))

is.na(as.numeric(args[1]))

is.na(as.numeric("NA"))

dmb <- read.csv("/space/s1/fiona_callahan/multiSim_100sp/randomRun1/sitetab_dumb.csv")
dim(dmb)