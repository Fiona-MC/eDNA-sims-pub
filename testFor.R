args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(args[1])

stopifnot(i != 1)

Sys.sleep(3)
print("did not stop")
print(i)