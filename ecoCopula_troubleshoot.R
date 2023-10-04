library(ecoCopula)

data_dir <- "/space/s1/fiona_callahan/multiSim11a/randomRun30/"

# load data
sim_data_raw <- readRDS(paste0(data_dir, "sim_data.Rdata"))
sitetab <- read.csv(paste0(data_dir, "sim_sitetab_sampled.csv"))
#locList <- readRDS(paste0(data_dir, "locList.Rdata"))
params <- readRDS(paste0(data_dir, "params.Rdata"))

names_cov <- params$names_cov
names_species <- params$names_species

pa_real <- as.matrix(sitetab[1:28, names_species])
X_real <- sitetab[1:28, names_cov]
fit0 <- stackedsdm(pa_real,~., data = X_real, family="binomial", ncores = 1) #eqiv. manyglm()
# ERROR!!! but sometimes not. 4 or fewer species doesnt work. 
# could maybe be fixed somehow with lambda -- not sure if it's worth doing this or just run it with bigger sim
cgr_sim <- cgr(fit0, method="AIC")
plot(cgr_sim, pad = 1)


#example that works
# spider data is stored in ecoCopula
data(spider)
X <- as.data.frame(spider$x) # environmental covariates
X3 <- X[,1:5] # environmental covariates
abund <- spider$abund
pa = (abund>0)*1 # presence-absence of spiders
pa3<-pa[,1:4]

# fit marginal model
spider_pa <- stackedsdm(pa3,~1, data = X, family="binomial",ncores = 2) #eqiv. manyglm()
# fit copula ordination 
#spid_lv=cord(spider_pa)
# biplot
#plot(spid_lv,biplot = TRUE)

# fit copula ordination 
spid_gr=cgr(spider_pa, seed=3)
# biplot
plot(spid_gr, pad=1)
