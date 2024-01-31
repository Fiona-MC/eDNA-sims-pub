args <- commandArgs(trailingOnly = TRUE)
print("hello world")
i <- as.numeric(args[1])

stopifnot(i != 1)

Sys.sleep(3)
print("did not stop")
print(i)
