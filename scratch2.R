compare_append_c_preallocate <- function(n) {
  # Create an empty vector
  vec <- numeric(0)
  
  # Time for append
  append_time <- system.time(
    for (i in 1:n) {
      vec <- append(vec, i)
    }
  )[3]
  
  # Create an empty vector
  vec <- numeric(0)
  
  # Time for c()
  c_time <- system.time(
    for (i in 1:n) {
      vec <- c(vec, i)
    }
  )[3]
  
  # Pre-allocate space
  vec <- numeric(n)
  
  # Time for pre-allocating
  preallocate_time <- system.time(
    for (i in 1:n) {
      vec[i] <- i
    }
  )[3]
  
  results <- data.frame(Method = c("append", "c", "Pre-allocation"),
                        Time = c(append_time, c_time, preallocate_time))
  
  return(results)
}

# Test the function with a sample size
n <- 10000
results <- compare_append_c_preallocate(n)
print(results)
