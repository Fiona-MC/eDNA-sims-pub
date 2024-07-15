
# Create synthetic data
n <- 12 # Number of samples
p <- 12   # Number of predictors

# Create predictor variables
predictors <- as.data.frame(matrix(rnorm(n * p, mean = 0, sd = 1), ncol = p))
names(predictors) <- paste0("x", 1:p)

y <- rnorm(n, 0, 1)  # Generate binary response variable

# Combine into a data frame
data <- data.frame(y, predictors)

# Fit the logistic regression model
model <- lm(y ~ ., data = data)

# View the summary of the model
model_summary <- summary(model)
print(model_summary)

model2 <- glm(y ~ ., data = data, family = gaussian)
model_summary2 <- summary(model2)
print(model_summary2)




# Create synthetic data
n <- 100 # Number of samples
p <- 100 # Number of predictors

# Create predictor variables
predictors <- as.data.frame(matrix(rnorm(n * p, mean = 0, sd = 1), ncol = p))
names(predictors) <- paste0("x", 1:p)

y <- rbinom(n, 1, 0.5)  # Generate binary response variable

# Combine into a data frame
data <- data.frame(y, predictors)

# Fit the logistic regression model
model <- glm(y ~ ., data = data, family = binomial)

# View the summary of the model
model_summary <- summary(model)
print(model_summary)














# Step 1: Create synthetic data
n <- 12 # Number of samples
p <- 10   # Number of predictors

# Create predictor variables
predictors <- as.data.frame(matrix(rnorm(n * p, mean = 0, sd = 1), ncol = p))
names(predictors) <- paste0("x", 1:p)

# Create continuous response variable
coefficients <- runif(p, min = 0.5, max = 1.5)  # Random coefficients for each predictor
y <- as.matrix(predictors) %*% coefficients + rnorm(n)  # Generate response variable with some noise

# Combine into a data frame
data <- data.frame(y, predictors)

# Step 2: Fit the linear regression model
model <- glm(y ~ ., data = data)

# Step 3: View the summary of the model
model_summary <- summary(model)
print(model_summary)







# Create synthetic data
n <- 10  # Number of samples
p <- 10   # Number of predictors

# Create predictor variables
predictors <- as.data.frame(matrix(rnorm(n * p, mean = 0, sd = 1), ncol = p))
names(predictors) <- paste0("x", 1:p)

# Create binary response variable
coefficients <- runif(p, min = -8, max = 8)  # Random coefficients for each predictor
linear_combination <- as.matrix(predictors) %*% coefficients
prob <- 1 / (1 + exp(-linear_combination))  # Convert to probabilities using the logistic function
y <- rbinom(n, 1, prob)  # Generate binary response variable

# Ensure there's variability in the response
table(y)

# Combine into a data frame
data <- data.frame(y, predictors)

# Fit the logistic regression model
model <- glm(y ~ ., data = data, family = binomial)

# View the summary of the model
model_summary <- summary(model)
print(model_summary)









# testing adding predictors
samp <- 100

wrongL <- c()
for (k in c(1, 10, 25, 50, 75, 95)) {
    wrong <- 0
    for (tests in 1:1000) {
        sp1 <- rnorm(n = samp, 0, 1)
        sp2 <- sp1 + rnorm(n = samp, 0, 2)
        sp3 <- rnorm(samp, 0, 1)

        cov <- matrix(data = rnorm(samp * k, 0, 1), nrow = samp, ncol = k)

        data <- data.frame(sp1 = sp1, sp2 = sp2, sp3 = sp3)
        data <- cbind(data, cov)

        cor(sp1, sp2)

        model <- lm(sp1 ~ ., data)
        model_summary <- data.frame(summary(model)$coefficients)

        if (model_summary[["Pr...t.."]][3] < model_summary[["Pr...t.."]][2]) {
            wrong <- wrong + 1
        }
    }
    wrongL <- c(wrongL, wrong)
}

data <- data.frame(k = c(1, 10, 25, 50, 75, 95), mistakes = wrongL / 1000)

ggplot(data = data, aes(x = k, y = mistakes)) +
    geom_point() + 
    geom_line()

ggsave("/space/s1/fiona_callahan/toy_lm_covVsMistakes.pdf")
