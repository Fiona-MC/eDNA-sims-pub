
a1 <- rnorm(n = 100)
a2 <- a1 + rnorm(n = 100, sd = 0.00001)

y <- 2 * a1 + rnorm(n = 100, sd = 0.01)

df <- data.frame(y = y, a1 = a1, a2 = a2)

lm <- lm(data = df, formula = y ~ 1 + a1 + a2)

lm <- lm(data = df, formula = y ~ 1 + a2)

summary(lm)
