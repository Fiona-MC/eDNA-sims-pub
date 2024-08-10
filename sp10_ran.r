setwd( "/home/KaiyuanLi/current_code/JAGS")
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
data_dir=args[1]
save_dir=args[2]
#dir.create(save_dir,recursive =TRUE)
sitetabName=args[3]
site_tab=read.csv(paste0(data_dir,sitetabName))



n_env_vars <-length(grep("^Cov", colnames(site_tab), value = TRUE))
n_sites  <-length(unique(site_tab$labID))
n_species <-length(grep("^Sp", colnames(site_tab), value = TRUE))


#Consider each observation as a single site

Occur <-as.matrix(site_tab %>%select(starts_with("Sp")))
colnames(Occur)<-NULL
X<-as.matrix(site_tab %>%select(starts_with("Cov")))
colnames(X)<-NULL




n.chains <- 3
n.iter <- 5000
n.burn <- 4000
n.thin <- 1


df <- 1

model_name <- 'sim_model'
source('fit_JSDM.r')

#Diagnose(Beta, 'rhat')
#save(file ='/home/KaiyuanLi/current code/JAGS/sim_ave.RData')
#load('/home/KaiyuanLi/current code/JAGS/sim_ave.RData')
save(time.taken,file=paste0(save_dir, 'time_sim.RData'))
save(Beta,file=paste0(save_dir, 'sim_beta.RData'))
save(Rho,file=paste0(save_dir, 'sim_rho.RData'))
save(Mu,file=paste0(save_dir, 'sim_Mu.RData'))
save(EnvRho,file=paste0(save_dir, 'sim_EnvRho.RData'))








