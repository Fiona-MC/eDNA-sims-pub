#readRDS('/space/s1/fiona_callahan/multiSim_100/randomRun1/spiecEasi_res_mb/trial1/inferenceRes.Rdata')
library(dplyr)
library(ggplot2)
.libPaths(c("/home/KaiyuanLi/R/x86_64-pc-linux-gnu-library/4.3/", .libPaths()))
library(igraph)
data_dir='/space/s1/KaiyuanLi/data/JAGS/results/filtered/sp'
args <- commandArgs(trailingOnly = TRUE)
sp_num=10
sample_num=1000
true_dir='/space/s1/fiona_callahan/multiSim_10sp/randomRun'
save_dir='/space/s1/KaiyuanLi/data/JAGS/results/filtered/sp'
covs <- FALSE
cutoff <- c(0, 1, 0.0000001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08 ,0.09, 0.1, 0.15,0.2, .3,.4, .5,.6,.7,.8,.9)
#load(paste0(data_dir,sp_num,'_',10000,'/',5,'/sim_rho.RData'))
#dim(Rho[[1]])
numRuns <- 100 # how many Randomruns
numcutoff <- length(cutoff) # the number of cutoffs
cluster <- FALSE #just controls cluster plotting
crediable<-function(rho,.n.species,cutoff)
{
    .n.species=length(rho[[1]][1,1,])
    #corMx<-SUMMARY(rho, median)
    inferred_corMx<-matrix(0,nrow=.n.species,ncol=.n.species)
    for(i in 1:(.n.species-1))
    {
        for(j in (i+1):.n.species)
        {
            element_temp<-c()
            n.chains<-length(rho)
            for(k in 1:n.chains)
            {
                element_temp<-c(element_temp,rho[[k]][,i,j])
            }
            if(quantile(element_temp,cutoff/2)>0 )
            {
                inferred_corMx[i,j]=1
                inferred_corMx[j,i]=1
            }
            if(quantile(element_temp,1-cutoff/2)<0)
            {
                inferred_corMx[i,j]=-1
                inferred_corMx[j,i]=-1
            }

        }
    }
    return(inferred_corMx)
}

#!/usr/bin/env Rscript


runL <- c()
trialL <- c()
num_correctInferencesL <- c() # true pos
num_incorrect_alphaL <- c() # species interaction incorrectInference
num_incorrectInferencesL <- c() # false pos (or wrong direction)
num_actualEffectsL <- c() # total actual effects
num_possibleEffectsL <- c() # n^2 + n*p -n (minus n because of the diagonal of alpha)
num_missedEffects_alphaL <- c()
num_missedEffects_betaL <- c()
num_missedEffectsL <- c()# false negatives
timeL <- c()
finished_trL <- TRUE # T or F whether this trial finished the INLA part
num_correct_clusterL <- c()
num_incorrect_clusterL <- c()
num_missed_clusterL <- c()
alpha_direction_mistakesL <- c()
alpha_incorrect_undirectedL <- c()
alpha_correct_undirectedL <- c()
num_incorrect_betaL <- c()
TP_ignoreSign <- c()
FP_ignoreSign <- c()
TN_ignoreSign <- c()
FN_ignoreSign <- c()
TP_cluster <- c()
FP_cluster <- c()
TN_cluster <- c()
FN_cluster <- c()
TP_sign <- c()
FP_sign <- c()
TN_sign <- c()
FN_sign <- c()
TP_ignoreDirection <- c()
FP_ignoreDirection <-c()
TN_ignoreDirection <-c()
FN_ignoreDirection <- c()

i <- 1
for (run in 1:numRuns) {
    alpha_path=paste0(data_dir,sp_num,'_',sample_num,'_filtered/',run,'/sim_rho.RData')
    beta_path=paste0(data_dir,sp_num,'_',sample_num,'_filtered/',run,'/sim_beta.RData')
    if(!file.exists(alpha_path))
    {
        next
    }  
    true_path=paste0(true_dir,run,'/paramsFiltered',sample_num,'.Rdata')
    simParms <- readRDS(true_path)
    actualAlpha <- sign(simParms$filteredAlpha)
    load(alpha_path)
    Beta=NA
    if (covs) { #haven't finished cov yet
            actualBeta <- sign(simParms$beta)
            covTypes <- unlist(lapply(simParms$covVars, FUN = function(x) {return(x[["type"]])}))
            if (dim(actualBeta)[2] == length(covTypes)) {
                actualBeta <- actualBeta[, covTypes != "constant"]
            }
            load(beta_path) 
        }
    for (trial in 1:numcutoff) {
        alphaInferred=crediable(Rho,sp_num,cutoff[trial])
        betaInferred=Beta
        inferredParms=list(alphaInferred=alphaInferred,betaInferred=betaInferred)
         # store run and trial info
        runL<- c(runL,run)
        trialL <- c(trialL,trial)
        if (finished_trL) {
            if (covs) { # take out constant covariate
                covTypes <- unlist(lapply(simParms$covVars, FUN = function(x) {return(x[["type"]])}))
                if (dim(inferredParms$betaInferred)[2] == length(covTypes)) {
                    inferredParms$betaInferred <- inferredParms$betaInferred[, covTypes != "constant"]
                }
            }
            alphaG <- graph_from_adjacency_matrix(actualAlpha != 0, mode = "undirected")
            inferredAlphaG <- graph_from_adjacency_matrix(inferredParms$alphaInferred != 0, mode = "undirected")
            connected_alpha_actual <- (distances(alphaG, v = 1:length(Rho[[1]][1,1,]), to = 1:length(Rho[[1]][1,1,])) != Inf) * 
                                        (diag(nrow = dim(actualAlpha)[1], ncol = dim(actualAlpha)[1]) == 0)
            connected_alpha_inferred <- (distances(inferredAlphaG, v = 1:length(Rho[[1]][1,1,]), to = 1:length(Rho[[1]][1,1,])) != Inf) * 
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
            count_incorrectT1_alpha <- sum(inferredParms$alphaInferred * actualAlpha == -1) + sum(actualAlpha == 0 & inferredParms$alphaInferred != 0)
            alpha_direction_mistakes <- sum(inferredParms$alphaInferred * actualAlpha == -1)

            if (covs) {
                count_incorrectT1_beta <- sum(inferredParms$betaInferred * actualBeta == -1) + sum(actualBeta == 0 & inferredParms$betaInferred != 0)
                count_incorrectT1 <- count_incorrectT1_beta + count_incorrectT1_alpha
            } else {
                count_incorrectT1 <- count_incorrectT1_alpha
            }

            count_incorrect_cluster <- sum(connected_alpha_inferred * connected_alpha_actual == -1) + 
                                            sum(connected_alpha_actual == 0 & connected_alpha_inferred != 0)

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

            alphaInferred <- inferredParms$alphaInferred
            betaInferred <- inferredParms$betaInferred
            if (covs) {
                TP_cluster[i] <- sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1) + 
                                    sum(abs(betaInferred) * abs(actualBeta) == 1)
                FP_cluster[i] <- sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0)) + 
                                    sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0))
                FN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1))
                TN_cluster[i] <- sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0)) -
                                    length(Rho[[1]][1,1,]) # subtract diagonal
                try(if (TP_cluster[i] + FP_cluster[i] + FN_cluster[i] + TN_cluster[i] != 
                        (length(Rho[[1]][1,1,])^2 - length(Rho[[1]][1,1,]) + length(Rho[[1]][1,1,]) * (simParms$numCovs - 1))) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx in cluster")})

            } else {
                TP_cluster<-c(TP_cluster, sum(abs(connected_alpha_inferred) * abs(connected_alpha_actual) == 1) )
                FP_cluster <-c(FP_cluster, sum((abs(connected_alpha_inferred) == 1) & (abs(connected_alpha_actual) == 0)) )
                FN_cluster<- c(FN_cluster,sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 1)) )
                TN_cluster <-c(TN_cluster, sum((abs(connected_alpha_inferred) == 0) & (abs(connected_alpha_actual) == 0)) -
                                    length(Rho[[1]][1,1,]) )# subtract diagonal
                #try(if (TP_cluster[i] + FP_cluster[i] + FN_cluster[i] + TN_cluster[i] != length(Rho[[1]][1,1,])^2 - length(Rho[[1]][1,1,])) 
                         #   {stop("count_mistakes_general.R something is wrong with the confusion mx in cluster")})
            }

            # undirected meaning that a-->b iff b-->a in "actual alpha". This also ignores the sign.
            undirected_alpha_actual <- distances(alphaG, v = 1:length(Rho[[1]][1,1,]), to = 1:length(Rho[[1]][1,1,])) == 1
            undirected_alpha_inferred <- distances(inferredAlphaG, v = 1:length(Rho[[1]][1,1,]), to = 1:length(Rho[[1]][1,1,])) == 1

            alpha_incorrect_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == -1) + 
                                            sum(undirected_alpha_actual == 0 & undirected_alpha_inferred != 0)
            alpha_correct_undirected <- sum(undirected_alpha_inferred * undirected_alpha_actual == 1)

            # add to running lists
            if (covs) {
                num_incorrect_betaL[i] <- count_incorrectT1_beta
                num_missedEffects_betaL[i] <- num_missedEffects_beta
            } 
            
            if (covs) {
                TP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 1)) + 
                                    sum((abs(betaInferred) == 1) & (abs(actualBeta) == 1))
                FP_ignoreSign[i] <- sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0)) + 
                                    sum((abs(betaInferred) == 1) & (abs(actualBeta) == 0))
                FN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 1))
                TN_ignoreSign[i] <- sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0)) + 
                                    sum((abs(betaInferred) == 0) & (abs(actualBeta) == 0)) -
                                    length(Rho[[1]][1,1,]) # subtract diagonal
                try(if (TN_ignoreSign[i] + FN_ignoreSign[i] + FP_ignoreSign[i] + TP_ignoreSign[i] != 
                        (length(Rho[[1]][1,1,])^2 - length(Rho[[1]][1,1,]) + length(Rho[[1]][1,1,]) * (simParms$numCovs - 1))) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx")})

            } else {
                TP_ignoreSign<- c(TP_ignoreSign,sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 1)))
                FP_ignoreSign<- c(FP_ignoreSign, sum((abs(alphaInferred) == 1) & (abs(actualAlpha) == 0)) )
                FN_ignoreSign<-  c(FN_ignoreSign,sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 1)) )
                TN_ignoreSign<-  c(TN_ignoreSign,sum((abs(alphaInferred) == 0) & (abs(actualAlpha) == 0)) -
                                        length(Rho[[1]][1,1,]) )# subtract diagonal
                # subtract diagonal
                TP_ignoreDirection <- c(TP_ignoreDirection,sum(abs(undirected_alpha_inferred) * abs(undirected_alpha_actual) == 1, na.rm = TRUE) )
                FP_ignoreDirection <- c(FP_ignoreDirection,sum((abs(undirected_alpha_inferred) == 1) & (abs(undirected_alpha_actual) == 0), na.rm = TRUE) )
                FN_ignoreDirection<- c(FN_ignoreDirection,sum((abs(undirected_alpha_inferred) == 0) & (abs(undirected_alpha_actual) == 1), na.rm = TRUE) )
                TN_ignoreDirection <- c( TN_ignoreDirection,sum((abs(undirected_alpha_inferred) == 0) & (abs(undirected_alpha_actual) == 0), na.rm = TRUE) - 
                                    sp_num) # su
                        #         {stop("count_mistakes_general.R something is wrong with the confusion mx")})
            }

            if (covs) {
                TP_sign[i] <- sum(alphaInferred == 1 & actualAlpha == 1) + 
                                    sum(betaInferred == 1 & actualBeta == 1) +
                                    sum(alphaInferred == -1 & actualAlpha == -1) + 
                                    sum(betaInferred == -1 & actualBeta == -1)
                FP_sign[i] <- sum(alphaInferred != 0 & actualAlpha == 0) + 
                                    sum(betaInferred != 0 & actualBeta == 0) +
                                    sum(alphaInferred == 1 & actualAlpha == -1) + 
                                    sum(betaInferred == 1 & actualBeta == -1) +
                                    sum(alphaInferred == -1 & actualAlpha == 1) + 
                                    sum(betaInferred == -1 & actualBeta == 1) 
                FN_sign[i] <- sum(alphaInferred == 0 & actualAlpha != 0) + 
                                    sum(betaInferred == 0 & actualBeta != 0)
                TN_sign[i] <- sum(alphaInferred == 0 & actualAlpha == 0) + 
                                    sum(betaInferred == 0 & actualBeta == 0) -
                                    length(Rho[[1]][1,1,]) # subtract diagonal
                try(if (TP_sign[i] + FP_sign[i] + FN_sign[i] + TN_sign[i] != 
                        (length(Rho[[1]][1,1,])^2 - length(Rho[[1]][1,1,]) + length(Rho[[1]][1,1,]) * (simParms$numCovs - 1))) 
                            {stop("count_mistakes_general.R something is wrong with the confusion mx")})

            } else {
                TP_sign <- c(TP_sign,sum(alphaInferred != 0 & actualAlpha == alphaInferred))
                FP_sign <- c(FP_sign,sum(alphaInferred != 0 & actualAlpha == 0) +
                                sum(alphaInferred == 1 & actualAlpha == -1) + 
                                sum(alphaInferred == -1 & actualAlpha == 1) )
                FN_sign <- c(FN_sign,sum(alphaInferred == 0 & actualAlpha != 0))
                TN_sign <- c(TN_sign,sum(alphaInferred == 0 & actualAlpha == 0) -
                                    length(Rho[[1]][1,1,])) # subtract diagonal
                #try(if (TP_sign[i] + FP_sign[i] + FN_sign[i] + TN_sign[i] != length(Rho[[1]][1,1,])^2 - length(Rho[[1]][1,1,])) 
                 #           {stop("count_mistakes_general.R something is wrong with the confusion mx")})
            }

            # add to running lists
            num_incorrect_alphaL<-c(num_incorrect_alphaL, count_incorrectT1_alpha)
            num_correctInferencesL<-c(num_correctInferencesL, num_correct)
            num_incorrectInferencesL<-c(num_incorrectInferencesL, count_incorrectT1)
            num_missedEffects_alphaL<-c(num_missedEffects_alphaL, num_missedEffects_alpha)
            num_missedEffectsL<-c(num_missedEffectsL, num_missedEffects)
            alpha_direction_mistakesL<-c(alpha_direction_mistakesL, alpha_direction_mistakes)
            num_correct_clusterL<-c(num_correct_clusterL, num_correct_cluster)
            num_incorrect_clusterL<-c(num_incorrect_clusterL, count_incorrect_cluster)
            num_missed_clusterL<-c(num_missed_clusterL, num_missedEffects_cluster)
            alpha_incorrect_undirectedL<-c(alpha_incorrect_undirectedL, alpha_incorrect_undirected)
            alpha_correct_undirectedL<-c(alpha_correct_undirectedL, alpha_correct_undirected)

            
        } 
        
        # count total number of actual effects in the model
        if (covs) {
            num_actualEffects <- sum(abs(actualAlpha)) + sum(abs(actualBeta))
        } else {
            num_actualEffects <- sum(abs(actualAlpha))
        }
        num_actualEffectsL<-c(num_actualEffectsL, num_actualEffects)

        if (covs) {
            num_possibleEffects <- dim(actualAlpha)[1]^2 - dim(actualAlpha)[1] + (dim(actualBeta)[1] * dim(actualBeta)[2])
        } else {
            num_possibleEffects <- dim(actualAlpha)[1]^2 - dim(actualAlpha)[1]
        }
        num_possibleEffectsL<-c(num_possibleEffectsL, num_possibleEffects)

        if (cluster) {
          layout <- layout_with_mds(graph = alphaG)
          plot(alphaG, main = "Actual Network", vertex.size = 0)
          png(filename = paste0(data_dir,sp_num,'_',sample_num,'/', "actualNetworkPlot.png"), height = 800, width = 800, units = "px")
          dev.off()
          layout <- layout_with_mds(graph = inferredAlphaG)
          plot(inferredAlphaG, main = paste0("Inferred Network: ", data_dir,sp_num,'_',sample_num,'/'), layout = layout, vertex.size = 0)
          png(filename = paste0(data_dir,sp_num,'_',sample_num,'/', "inferredNetworkPlot.png"), height = 800, width = 800, units = "px")
          dev.off()
        }
        i <- i + 1
    }
}

df <- data.frame(sim_run = runL, 
            cutoff_val = trialL, 
            num_incorrect_alpha = num_incorrect_alphaL,
            num_correctInferences = num_correctInferencesL, 
            num_incorrectInferences = num_incorrectInferencesL, 
            num_actualEffects = num_actualEffectsL,
            num_missedEffects_alpha = num_missedEffects_alphaL,
            num_missedEffectsL = num_missedEffectsL,
            num_possibleEffectsL = num_possibleEffectsL,
            num_correct_cluster = num_correct_clusterL,
            num_incorrect_cluster = num_incorrect_clusterL,
            num_missed_cluster = num_missed_clusterL,
            alpha_direction_mistakes = alpha_direction_mistakesL,
            alpha_incorrect_undirected = alpha_incorrect_undirectedL,
            alpha_correct_undirected = alpha_correct_undirectedL,
            TP_ignoreSign = TP_ignoreSign,
            FP_ignoreSign = FP_ignoreSign,
            TN_ignoreSign = TN_ignoreSign,
            FN_ignoreSign = FN_ignoreSign,
            TP_cluster = TP_cluster,
            FP_cluster = FP_cluster,
            TN_cluster = TN_cluster,
            FN_cluster = FN_cluster,
            TP_sign = TP_sign,
            FP_sign = FP_sign,
            TN_sign = TN_sign,
            FN_sign = FN_sign,
            TP_ignoreDirection = TP_ignoreDirection,
            FP_ignoreDirection = FP_ignoreDirection,
            TN_ignoreDirection = TN_ignoreDirection,
            FN_ignoreDirection = FN_ignoreDirection)

df
head(df)
# print(df)
#df <-df %>% arrange(trial)
#df<-cbind(df,cutoff)
roc_nosign<-df %>% group_by(cutoff_val) %>% summarise(TP=mean(TP_ignoreSign),FP=mean(FP_ignoreSign)
,FN=mean(FN_ignoreSign),TN=mean(TN_ignoreSign)) %>% dplyr::mutate(TPR=TP/(TP+FN),FPR=FP/(FP+TN))

roc_nosigned <- rbind(data.frame(FPR = 0, TPR = 0), data.frame(FPR =roc_nosign$FPR,TPR=roc_nosign$TPR))
ROC_plot <- ggplot(roc_nosigned, aes(x = FPR, y = TPR))+geom_line(col='red')+geom_point(size = 1)+geom_abline(,slope = 1, intercept = 0, color = "black", linetype = "dashed") 
print(ROC_plot)

ggsave(paste0(data_dir,sp_num,'_',sample_num,'_filtered/nosignfilted',sample_num,'.png'),ROC_plot)


write.csv(df, paste0(data_dir,sp_num,'_',sample_num,'_filtered/', "mistakes.csv"))



# for all the elements in one matrix, median<0, negative, t;median>0, positive,t
# for same elements in matrixs, t statistics.
# for (run in 1:numRuns) 
# {
#     alpha_path=paste0(data_dir,sp_num,'_',sample_num,'/',run,'/sim_rho.RData')
#     for(i in 1:(.n.species-1))
#     {
#         for(j in (i+1):.n.species)
#         {
#             positive_t<-c()
#             negative_t<-c()
#             element_temp<-c()
#             for(k in 1:n.chains)
#             {
#                 element_temp<-c(element_temp,Rho[[k]][,i,j])
#             }
#             element_p=t.test(element_temp,mu=0,alternative='two.sided')
#             temp_med=median(element_p)
#             if(temp_med<0)
#             {
#                 negative_t<-c(negative,element_p$p.value)
#             }
#             if(temp_med>0)
#             {
#                 positive_t<-c(negative,element_p$p.value)
#             }
#         }
#    }
# }
# t_zero<-function(rho,n.chains,.n.species)
# {
#     for(i in 1:(.n.species-1))
#     {
#         for(j in (i+1):.n.species)
#         {
#             positive_t<-c()
#             negative_t<-c()
#             element_temp<-c()
#             for(k in 1:n.chains)
#             {
#                 element_temp<-c(element_temp,rho[[k]][,i,j])
#             }
#             element_p=t.test(element_temp,mu=0,alternative='two.sided')
#             temp_med=median(element_p)
#             if(temp_med<0)
#             {
#                 negative_t<-c(negative,element_p$p.value)
#             }
#             if(temp_med>0)
#             {
#                 positive_t<-c(negative,element_p$p.value)
#             }
#         }
#    }
#     return(list(ipositive_t,negative_t))
# }


# crediable_CI<-function(rho,.n.species,cutoff)
# {
#     #corMx<-SUMMARY(rho, median)
#     positive_t<-c()
#     negative_t<-c()
#     zero_t<-c()
#     for(i in 1:(.n.species-1))
#     {
#         for(j in (i+1):.n.species)
#         {
#             element_temp<-c()
#             n.chains<-length(rho)
#             for(k in 1:n.chains)
#             {
#                 element_temp<-c(element_temp,rho[[k]][,i,j])
#             }
#             if(quantile(element_temp,cutoff/2)>0 )
#             {
#                 inferred_corMx[i,j]=median
#                 inferred_corMx[j,i]=1
#             }
#             if(quantile(element_temp,1-cutoff/2)<0)
#             {
#                 inferred_corMx[i,j]=-1
#                 inferred_corMx[j,i]=-1
#             }

#         }
#     }
#     return(inferred_corMx)
# }

# ###very confused about this one.
# positive_t<-c()
# negative_t<-c()
# .n.species=10
# zero_t<-c()
# for (run in 9:9) 
# {
#     alpha_path=paste0(data_dir,sp_num,'_',sample_num,'_filtered/',run,'/sim_rho.RData')
#     load(alpha_path)
#     .n.species=.n.species=length(Rho[[1]][1,1,])
#     for(i in 1:(.n.species-1))
#     {
#         for(j in (i+1):.n.species)
#         {
#             element_temp<-c()
#             n.chains<-length(Rho)
#             for(k in 1:n.chains)
#             {
#                 element_temp<-c(element_temp,Rho[[k]][,i,j])
#             }
#             temp_med=median(element_temp)
#             if(temp_med<0)
#             {
#                 element_p=t.test(element_temp,mu=0,alternative='less')
#                 negative_t<-c(negative_t,temp_med)
#             }
#             if(temp_med>0)
#             {
#                 element_p=t.test(element_temp,mu=0,alternative='greater')
#                 positive_t<-c(positive_t,temp_med)
#             }
#             if(quantile(element_temp,0.05/2)*quantile(element_temp,1-0.05/2) )
#             {
#                 zero_t<-c(zero_t,temp_med)
#             }
            
#         }
#    }
# }
# library(ggplot2)
# library(reshape2)

# # Example vectors
# positive <-positive_t  # Example data 1
# negative <- negative_t  # Example data 2
# zero <-zero_t

# # Create a combined long-format data frame
# data_long <- data.frame(
#   value = c(vector1, vector2, vector3),
#   variable = factor(rep(c("postive", "negative", "zero"), times = c(length(vector1), length(vector2), length(vector3))))
# )

# # Melt is not needed since data is already in long format
# # Plot using ggplot2
# library(ggplot2)
# ggplot(data_long, aes(x = value, fill = variable)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "Comparison of Distributions",
#        x = "Value",
#        y = "Density") +
#   theme_minimal()






# median<-matrix(0,ncol=45,nrow=100)
# .n.species=10
# .n.species*(i-1)+(j-1)
# for (run in 1:numRuns) 
# {
#     alpha_path=paste0(data_dir,sp_num,'_',sample_num,'/',run,'/sim_rho.RData')
#     load(alpha_path)
#     for(i in 1:(.n.species-1))
#     {
#         for(j in (i+1):.n.species)
#         {
#             element_temp<-c()
#             n.chains<-length(Rho)
#             for(k in 1:n.chains)
#             {
#                 element_temp<-c(element_temp,Rho[[k]][,i,j])
#             }
#             median(element_temp)
            
#         }
#    }
# }


# The distribution of the medians for single random run

# The distribution of one interaction over 100 simulations, just need to make sure that all the
## they are not the same

# store the crediable interval inferred matrix for every cutoff

# add more cutoff into the threshold vector

# Beta is env to species, try to figure out the filteredBeta!!!