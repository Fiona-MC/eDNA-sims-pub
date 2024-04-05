
# Create the matrix H
p <- 1000
diagonal <- 1-1/p
off_diagonal <- -1/p

G <- matrix(off_diagonal, nrow = p, ncol = p)
diag(G) <- diagonal
print(G)

det(G)

#G_inv <- solve(G)
# Print the inverse matrix H_inv
#print(G_inv)
G
testMx <- matrix(sample(-10:10, p * p, replace = TRUE), nrow = p, ncol = p)

testMx

G %*% testMx %*% G

diffs <- (G %*% testMx %*% G) - testMx

hist(diffs)
mean(diffs)
sd(diffs)

# for negative vals in omega, multiplying by G makes many of them positive and vice versa


# for omega drawn from -10 to 10

#p=100
#> mean(diffs)
#[1] -0.0039
#> sd(diffs)
#[1] 0.8180932

#p=10
#> mean(diffs)
#[1] -0.23
#> sd(diffs)
#[1] 2.428409

#p=1000
#> mean(diffs)
#[1] -0.013319
#> sd(diffs)
#[1] 0.2623388

