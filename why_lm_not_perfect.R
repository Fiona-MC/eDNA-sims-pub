# this is edited to work for filtered species names
# should actually work for both though

library(stats)
library(igraph)
library(stringr)

data_dir <- "/space/s1/fiona_callahan/multiSim_100sp_logi0.5/"
numRuns <- as.numeric(100)
covs <- (as.numeric("1") == 1)
logi <- (as.numeric("1") == 1)
sitetab_name <- "logiSim_sitetab_readAbd_sampled250.csv"
outName <- "linearReg_mistakes_sampled250_cov_logi_100runs"
countCovs <- FALSE

randomSim <- str_detect(data_dir, "random")
run <- 1
simParms <- readRDS(paste0(data_dir, "randomRun", run, "/params.Rdata"))

############### DO LOGISTIC REGRESSION ######################
# load sitetab
sim_sitetab_sampled <- read.csv(file = paste0(data_dir, "randomRun", run, "/", sitetab_name), header = TRUE)

sp_glm_L <- list()
speciesNames <- names(sim_sitetab_sampled)[grep("Sp", names(sim_sitetab_sampled))]
numSpecies <- length(speciesNames)
#print(numSpecies)
for (speciesName in speciesNames) {
    if (covs) {
        # glm
        fmla <- as.formula(paste(speciesName, "~ ", 
                            paste(c(speciesNames[speciesNames != speciesName]), collapse = "+"), "+",
                            paste(c(simParms$names_cov), collapse = "+")))
    } else {
        fmla <- as.formula(paste(speciesName, "~ ", 
                        paste(c(speciesNames[speciesNames != speciesName]), collapse = "+")))
    }
    model <- lm(fmla, data = sim_sitetab_sampled)
    model_summary <- data.frame(summary(model)$coefficients)
    sp_glm_L[[speciesName]] <- model_summary
}

minPval <- unlist(lapply(sp_glm_L, FUN = function(model_summary) {
    sort(model_summary[, "Pr...t.."])[1]
}))

second_pval <- unlist(lapply(sp_glm_L, FUN = function(model_summary) {
    sort(model_summary[, "Pr...t.."])[2]
}))

data.frame(y = names(minPval), min_pval = minPval, second_pval = second_pval)

# check if lowest pvalue is the correct one
miss <- 0
for (speciesName in speciesNames) {
    spNum <- as.numeric(gsub("[^0-9]", "", speciesName))
    if (spNum <= 50) {
        actualInteractSpNum <- spNum + 50
    } else {
        actualInteractSpNum <- spNum - 50
    }
    actualInteractSpName <- paste0("Sp", actualInteractSpNum)
    # get interaction with min pval
    model_summary <- sp_glm_L[[speciesName]]
    low_pval_interaction <- row.names(model_summary)[model_summary[["Pr...t.."]] == min(model_summary[["Pr...t.."]])]
    #print(spNum)
    #print(low_pval_interaction)
    #print(actualInteractSpName)
    if (low_pval_interaction != actualInteractSpName) {
        miss <- miss + 1
        #print(simParms$alpha[spNum, actualInteractSpNum])
    }
    pval_of_correct <- model_summary[["Pr...t.."]][row.names(model_summary) == actualInteractSpName]
    rank_of_correct <- which(sort(model_summary[["Pr...t.."]]) == pval_of_correct)
    min_pval <- min(model_summary[["Pr...t.."]])
    print(speciesName)
    print(rank_of_correct)
}



# hist of t-value by interaction
tvals_interact <- c()
tvals_noInteract <- c()
for (speciesName in speciesNames) {
    spNum <- as.numeric(gsub("[^0-9]", "", speciesName))
    if (spNum <= 50) {
        actualInteractSpNum <- spNum + 50
    } else {
        actualInteractSpNum <- spNum - 50
    }
    actualInteractSpName <- paste0("Sp", actualInteractSpNum)
    # get interaction with min pval
    model_summary <- sp_glm_L[[speciesName]]
    model_summary <- model_summary[grepl("^Sp", row.names(model_summary)), ]
    tvals_interact <- c(tvals_interact, model_summary[["t.value"]][row.names(model_summary) == actualInteractSpName])
    tvals_noInteract <- c(tvals_noInteract, model_summary[["t.value"]][row.names(model_summary) != actualInteractSpName])
}

tvals <- data.frame(
  value = c(tvals_interact, tvals_noInteract),
  group = c(rep("interact", length(tvals_interact)), 
            rep("no interact", length(tvals_noInteract)))
)

library(ggplot2)
# Plot the histograms
ggplot(tvals, aes(x = value, fill = group)) +
  geom_histogram(position = "dodge", alpha = 0.5, bins = 30) +
  labs(x = "t value", y = "Frequency") +
  theme_minimal()








# this part is to look at expected empirical correlations between samples by chance
run <- 1
sim_sitetab <- read.csv(file = paste0(data_dir, "randomRun", run, "/", "sitetab_abd_logi.csv"), header = TRUE)

n <- 250
corr_null <- c()
corr_alt <- c()
corr_alt_pos <- c()
reps <- 10000

corr_alt_119 <- c()

for (rep in 1:reps) {
    samp <- sample(1:length(sim_sitetab$Sp1), size = n) #nolint
    nullSamp1 <- sim_sitetab$Sp1[samp]
    nullSamp2 <- sim_sitetab$Sp2[samp]
    corr_null <- c(corr_null, cor(nullSamp1, nullSamp2))

    altSamp1 <- sim_sitetab$Sp3[samp]
    altSamp2 <- sim_sitetab$Sp53[samp]
    corr_alt <- c(corr_alt, cor(altSamp1, altSamp2))

    altSamp3 <- sim_sitetab$Sp50[samp]
    altSamp4 <- sim_sitetab$Sp100[samp]
    corr_alt_pos <- c(corr_alt_pos, cor(altSamp3, altSamp4))

    altSamp5 <- sim_sitetab$Sp27[samp]
    altSamp6 <- sim_sitetab$Sp77[samp]
    corr_alt_119 <- c(corr_alt_119, cor(altSamp5, altSamp6))
}

hist(corr_alt_119)

data <- data.frame(
  value = c(corr_null, corr_alt, corr_alt_pos, corr_alt_119),
  group = c(rep("corr_null", length(corr_null)), 
            rep("corr_alt", length(corr_alt)), 
            rep("corr_alt_pos", length(corr_alt_pos)), 
            rep("corr_alt_119", length(corr_alt_119)))
)

library(ggplot2)
# Plot the histograms
ggplot(data, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(x = "Values", y = "Frequency") +
  theme_minimal()

