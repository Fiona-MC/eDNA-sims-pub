# the idea here was to run a regression on the data simulated in the spiecEasi package 
# to see if that does better than SE
library(SpiecEasi)
library(stats)
library(igraph)
library(ggplot2)

mode <- "ignore_sign"
cluster <- FALSE

# get spiecEasi test data

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

graph <- SpiecEasi::make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

X_df <- data.frame(X)
names(X_df) <- sapply(seq_len(dim(X_df)[2]), FUN = function(num) {paste0("Sp", num)})
cutoffs <- c(0, 1e-128, 1e-64, 1e-32, 1e-16, 1e-8, 1e-4, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15,
             0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)


# run glm
simParms <- list()
simParms$names_species <- names(X_df)
simParms$numSpecies <- length(names(X_df))
covs <- FALSE
simParms$numCovs <- 1

sp_glm_L <- list()
for (speciesName in simParms$names_species) {
    if (covs) {
        # glm
        fmla <- as.formula(paste(speciesName, "~ ", 
                            paste(c(simParms$names_species[simParms$names_species != speciesName]), collapse = "+"), "+",
                            paste(c(simParms$names_cov), collapse = "+")))
    } else {
        fmla <- as.formula(paste(speciesName, "~ ", 
                        paste(c(simParms$names_species[simParms$names_species != speciesName]), collapse = "+")))
    }
    model <- glm(fmla, family = gaussian, data = X_df)
    model_summary <- data.frame(summary(model)$coefficients)
    sp_glm_L[[speciesName]] <- model_summary
}

i <- 1
TP_cluster <- rep(NA, times = length(cutoffs))
FP_cluster <- rep(NA, times = length(cutoffs))
FN_cluster <- rep(NA, times = length(cutoffs))
TN_cluster <- rep(NA, times = length(cutoffs))
TP_ignoreSign <- rep(NA, times = length(cutoffs))
FP_ignoreSign <- rep(NA, times = length(cutoffs))
FN_ignoreSign <- rep(NA, times = length(cutoffs)) 
TN_ignoreSign <- rep(NA, times = length(cutoffs))
num_incorrect_alphaL <- rep(NA, times = length(cutoffs))
num_correctInferencesL <- rep(NA, times = length(cutoffs))
num_incorrectInferencesL <- rep(NA, times = length(cutoffs))
num_missedEffects_alphaL <- rep(NA, times = length(cutoffs))
num_missedEffectsL <- rep(NA, times = length(cutoffs))
alpha_direction_mistakesL <- rep(NA, times = length(cutoffs))
num_correct_clusterL <- rep(NA, times = length(cutoffs))
num_incorrect_clusterL <- rep(NA, times = length(cutoffs))
num_missed_clusterL <- rep(NA, times = length(cutoffs))
alpha_incorrect_undirectedL <- rep(NA, times = length(cutoffs))
alpha_correct_undirectedL <- rep(NA, times = length(cutoffs))
num_actualEffects <- rep(NA, times = length(cutoffs))
num_actualEffectsL <- rep(NA, times = length(cutoffs))
num_possibleEffectsL <- rep(NA, times = length(cutoffs))

for (cutoff in cutoffs) {

    avg_alphInferred <- matrix(0, nrow = simParms$numSpecies, ncol = simParms$numSpecies)
    avg_betInferred <- matrix(0, nrow = simParms$numSpecies, ncol = simParms$numCovs - 1)
    nCompleteB <- 0
    nCompleteA <- 0

    # TODO get df of actual and predicted values to test if I am making the confusion mx wrong
    actual_predicted_vals <- data.frame()
    alpha_actualL <- c()
    alpha_predictedL <- c()
    beta_actualL <- c()
    beta_predictedL <- c()

    ############### LOAD ACTUAL PARAMS ######################

    actualBeta <- NA

    actualAlpha <- sign(Prec)
    for(iii in seq_len(dim(actualAlpha)[1])) {actualAlpha[iii, iii] <- 0}
    
    ############### GET INFERRED ALPHA AND BETA ######################
    betaInferred <- matrix(NA, nrow = simParms$numSpecies, ncol = (simParms$numCovs - 1))
    if (covs) {
        for (spNum in 1:simParms$numSpecies) {
            for (covNum in 1:(simParms$numCovs - 1)) {
                speciesName <- paste0("Sp", spNum)
                covName <- paste0("Cov", covNum)
                model_summary <- sp_glm_L[[speciesName]]
                betaInferred[spNum, covNum] <- (model_summary[covName, "Pr...t.."] < cutoff) * sign(model_summary[covName, "Estimate"])
            }
        }
        if(!is.na(sum(betaInferred))) {
            avg_betInferred <- avg_betInferred + betaInferred
            nCompleteB <- nCompleteB + 1
        }   
    }

    alphaInferred <- matrix(NA, nrow = simParms$numSpecies, ncol = simParms$numSpecies)
    for (spNum1 in 1:simParms$numSpecies) {
        for (spNum2 in 1:simParms$numSpecies) {
            speciesName1 <- paste0("Sp", spNum1)
            speciesName2 <- paste0("Sp", spNum2)
            model_summary <- sp_glm_L[[speciesName1]]
            if (spNum1 == spNum2) {
                alphaInferred[spNum1, spNum2] <- 0
            } else {
                alphaInferred[spNum1, spNum2] <- (model_summary[speciesName2, "Pr...t.."] < cutoff) * 
                                                sign(model_summary[speciesName2, "Estimate"])
            }
        }
    }


    if(!is.na(sum(alphaInferred))) {
        avg_alphInferred <- avg_alphInferred + alphaInferred
        nCompleteA <- nCompleteA + 1
    }

    ############### COMPARE INFERRECE RESULTS TO ACTUAL ALPHA AND BETA ######################
    inferredParms <- list()
    inferredParms$alphaInferred <- alphaInferred
    inferredParms$betaInferred <- betaInferred
    inferredParms$cutoff <- cutoff

    alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "undirected")
    inferredAlphaG <- graph_from_adjacency_matrix(inferredParms$alphaInferred != 0, mode = "undirected")
    connected_alpha_actual <- (distances(alphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) != Inf) * 
                                (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
    connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) != Inf) * 
                                (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
    # Note: inferredParms$betaInferred * actualBeta == 1 if and only if both are 1 or both are -1
    num_correct_alpha <- sum(inferredParms$alphaInferred * actualAlpha == 1)
    if (covs) {
        num_correct_beta <- sum(inferredParms$betaInferred * actualBeta == 1)
        num_correct <- num_correct_alpha + num_correct_beta
    } else {
        num_correct <- num_correct_alpha
    }

    num_correct_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == 1)

    # Note: here we want to add together inferrence in the wrong direction with saying there is an effect when it's actually 0
    # type 1 error -- inferring an effect where there is none or wrong direction of effect
    count_incorrectT1_alpha <- sum(inferredParms$alphaInferred * actualAlpha == -1) + 
                                    sum(actualAlpha == 0 & inferredParms$alphaInferred != 0)
    alpha_direction_mistakes <- sum(inferredParms$alphaInferred * actualAlpha == -1)

    if (covs) {
        count_incorrectT1_beta <- sum(inferredParms$betaInferred * actualBeta == -1) + 
                                    sum(actualBeta == 0 & inferredParms$betaInferred != 0)
        count_incorrectT1 <- count_incorrectT1_beta + count_incorrectT1_alpha
    } else {
        count_incorrectT1 <- count_incorrectT1_alpha
    }

    count_incorrect_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == -1) + 
                                    sum(connected_alpha_actual == 0 & connected_alpha_inferred != 0)

    if (covs) {
    TP_cluster[i] <- sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1) + 
                        sum(abs(betaInferred) * abs(actualBeta) == 1)
    FP_cluster[i] <- sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0)) + 
                        sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0))
    FN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1)) + 
                        sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1))
    TN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0)) + 
                        sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0))
    } else {
        TP_cluster[i] <- sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1) 
        FP_cluster[i] <- sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0)) 
        FN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1)) 
        TN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0)) 
    }
    
    # type 2 error
    # count missed effects (number of times that actual effect is 1 or -1 and inferred effect is 0)
    num_missedEffects_alpha <- sum(actualAlpha != 0 & inferredParms$alphaInferred == 0)

    if (covs) {
        num_missedEffects_beta <- sum(actualBeta != 0 & inferredParms$betaInferred == 0)
    }

    if (covs) {
        num_missedEffects <- num_missedEffects_alpha + num_missedEffects_beta
    } else {
        num_missedEffects <- num_missedEffects_alpha 
    }

    num_missedEffects_cluster <- sum(connected_alpha_actual != 0 & connected_alpha_inferred == 0)

    # undirected meaning that a-->b iff b-->a in "actual alpha". This also ignores the sign.
    undirected_alpha_actual <- distances(alphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) == 1
    undirected_alpha_inferred <- distances(inferredAlphaG, v = 1:simParms$numSpecies, to = 1:simParms$numSpecies) == 1

    alpha_incorrect_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == -1) + 
                                    sum(undirected_alpha_actual == 0 & undirected_alpha_inferred != 0)
    alpha_correct_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == 1)

    # add to running lists
    if (covs) {
        num_incorrect_betaL[i] <- count_incorrectT1_beta
        num_missedEffects_betaL[i] <- num_missedEffects_beta
    } 

    if (covs) {
        TP_ignoreSign[i] <- sum(abs(alphaInferred) * abs(actualAlpha) == 1) + 
                            sum(abs(betaInferred) * abs(actualBeta) == 1)
        FP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0)) + 
                            sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0))
        FN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1)) + 
                            sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1))
        TN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0)) + 
                            sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0))
    } else {
        TP_ignoreSign[i] <- sum(abs(alphaInferred) * abs(actualAlpha) == 1) 
        FP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0)) 
        FN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1)) 
        TN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0)) 
    }

    num_incorrect_alphaL[i] <- count_incorrectT1_alpha
    num_correctInferencesL[i] <- num_correct
    num_incorrectInferencesL[i] <- count_incorrectT1
    num_missedEffects_alphaL[i] <- num_missedEffects_alpha
    num_missedEffectsL[i] <- num_missedEffects
    alpha_direction_mistakesL[i] <- alpha_direction_mistakes
    num_correct_clusterL[i] <- num_correct_cluster
    num_incorrect_clusterL[i] <- count_incorrect_cluster
    num_missed_clusterL[i] <- num_missedEffects_cluster
    alpha_incorrect_undirectedL[i] <- alpha_incorrect_undirected
    alpha_correct_undirectedL[i] <- alpha_correct_undirected
    # count total number of actual effects in the model
    num_actualEffects <- sum(abs(actualAlpha))

    num_actualEffectsL[i] <- num_actualEffects
    num_possibleEffectsL[i] <- dim(actualAlpha)[1]^2 - dim(actualAlpha)[1]
    i <- i + 1
}

res_stats <- data.frame(
        TP_cluster = TP_cluster,
        FP_cluster = FP_cluster,
        FN_cluster = FN_cluster,
        TN_cluster = TN_cluster,
        TP_ignoreSign = TP_ignoreSign,
        FP_ignoreSign = FP_ignoreSign, 
        FN_ignoreSign = FN_ignoreSign,
        TN_ignoreSign = TN_ignoreSign,
        num_incorrect_alpha = num_incorrect_alphaL,
        num_correctInferences = num_correctInferencesL,
        num_incorrectInferences = num_incorrectInferencesL,
        num_missedEffects_alpha = num_missedEffects_alphaL,
        num_missedEffects = num_missedEffectsL,
        alpha_direction_mistakes = alpha_direction_mistakesL,
        num_correct_cluster = num_correct_clusterL,
        num_incorrect_cluster = num_incorrect_clusterL,
        num_missed_cluster = num_missed_clusterL,
        alpha_incorrect_undirected = alpha_incorrect_undirectedL,
        alpha_correct_undirected = alpha_correct_undirectedL,
        num_actualEffects = num_actualEffects,
        num_actualEffects = num_actualEffectsL,
        num_possibleEffects = num_possibleEffectsL
        )


# run spiecEasi on the same data
TPRs <- rep(NA, times = 100) 
FPRs <- rep(NA, times = 100)

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=100)



alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "undirected")
# theta = true_interactions
se.roc <- huge::huge.roc(se$est$path, theta = alphaG, verbose = FALSE)
se.TPRs <- se.roc$tp
se.FPRs <- se.roc$fp
se.methods <- rep("SpiecEasi", times = length(se.roc$fp))
se.thresholds <- se$lambda

length(se.TPRs) == length(se.FPRs) & length(se.TPRs) == length(se.thresholds)


# method | threshold | TPR | FPR 
lm.methods <- rep("linear regression", times = length(cutoffs))
lm.thresholds <- cutoffs
lm.TPR <- res_stats$TP_ignoreSign / (res_stats$TP_ignoreSign + res_stats$FN_ignoreSign) # TPR <- TP / (TP + FN)
lm.FPR <- res_stats$FP_ignoreSign / (res_stats$FP_ignoreSign + res_stats$TN_ignoreSign) # FPR <- FP / (FP + TN)

length(lm.TPR) == length(lm.FPR) & length(lm.TPR) == length(lm.thresholds)

ROC_data <- data.frame(
    method = c(se.methods, lm.methods),
    thresholds = c(se.thresholds, lm.thresholds),
    TPR = c(se.TPRs, lm.TPR),
    FPR = c(se.FPRs, lm.FPR)
)   


ROC_plot <- ggplot(ROC_data, aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 3) +
  stat_summary(aes(group = method), fun.y = mean, geom = "line", size = 1) +
  labs(title = paste("ROC: cluster = ", cluster), x = "FPR = FP/(FP+TN)", y = "TPR = TP/(TP+FN)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +   # Add y=x line (no skill)
  lims(x = c(0, 1), y = c(0, 1)) + # Set x and y-axis limits
  theme(text = element_text(size = 24))  # Set the base size for all text elements


# huge roc example
#generate data
#L = huge::huge.generator(d = 200, graph = "cluster", prob = 0.3)
#out1 = huge::huge(L$data)

#draw ROC curve
#Z1 = huge::huge.roc(out1$path,L$theta)

#Maximum F1 score
#max(Z1$F1)
