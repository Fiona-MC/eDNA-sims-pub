library(Matrix)
#library(highfrequency)
#library(magic)
# make small blocks fully connected with high correlations
make_blocks <- function(block_size, randomseed = block_size)
{
    set.seed(randomseed)
    block <- matrix(0, ncol = block_size, nrow = block_size)
    block <- block + diag(1, ncol = block_size, nrow = block_size)
    for(i in 1:(block_size - 1))
    {
        for(j in (i + 1):block_size)
        {
            cor_sign <- sample(c(-1, 1), 1)
            abs_element <- runif(1, 0.8, 1)
            block[i, j] <- cor_sign * abs_element
            block[j, i] <- block[i, j]
        }
    }
    new_block <- nearPD(block, corr = TRUE, keepDiag = TRUE)$mat
    return(new_block)
}
# example of make_blocks
bs1 <- 3
bs2 <- 2
result_matrix <- bdiag(make_blocks(bs1), make_blocks(bs2)) + matrix(0, nrow = bs1 + bs2, ncol = bs1 + bs2)
# adiag  useless function
CorrMx <- function(spec_num, block_num, block_size_lb, block_size_ub)
{
    block_size <- c(sample(block_size_lb:block_size_ub, 1))
    cor_Mx <- make_blocks(block_size[1])
    for(i in 1:block_num - 1)
    {
        temp_block_size <- sample(block_size_lb:block_size_ub, 1)
        cor_Mx <- bdiag(cor_Mx, make_blocks(temp_block_size, i))
        block_size <- c( block_size, temp_block_size)
    }
    cor_Mx <- bdiag(cor_Mx, matrix(0, nrow = spec_num - sum(block_size), ncol = spec_num - sum(block_size)))
    return(cor_Mx)
}
a <- CorrMx(spec_num = 100, block_num = 3, block_size_lb = 3, block_size_ub = 5)
a