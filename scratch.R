args <- commandArgs(trailingOnly = TRUE)

print(as.numeric(args[1]))
print(typeof(as.numeric(args[1])))

is.na(as.numeric(args[1]))

is.na(as.numeric("NA"))
