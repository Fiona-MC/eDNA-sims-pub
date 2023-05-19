# the purpose of this is to use the siumulation data in the code from Schliep paper  (/Users/fionacallahan/Documents/Nielsen_lab/schliep_code/RunMultivariateOrdinalSpatialAlgorithm.R)
# (DOI: 10.1111/geb.12666)

#install.packages(c("stats","ggplot2","coda","Rcpp","RcppArmadillo","inline","MCMCpack","fields","Matrix","spam","mvtnorm","bayesSurv","msm","fdrtool","tmvtnorm"))
library(stats)
library(ggplot2)
library(coda)

#MCMC Algorithms for model with extensions 

source("/home/fiona_callahan/schliep_code/MultivariateSpatialOrdinalAlgorithm.R")
library(Rcpp)
library(RcppArmadillo)
library(inline)
source("/home/fiona_callahan/schliep_code/MultivariateSpatialOrdinalAlgorithmExtensions.R")
source("/home/fiona_callahan/schliep_code/MultivariateExtensionsRcppCode.R")

numIters=5000
burn=500

configureDataSchliep<-function(subdir, sim_data_raw, locList, params, mode = "N", num_samples=num_samples){
  #timePts<-1:params$numGens #time points to "collect"
  timePts<-seq(from=1, to=params$numGens, by = 1000/num_samples)
  
  # Wrangle data
  locDF<-as.data.frame(matrix(unlist(locList), nrow=length(locList), ncol=2, byrow=T))
  distMx<-dist(locDF, diag=T, upper=T)
  
  Y<-list()
  X<-list()
  for(t_index in 1:length(timePts)){
    t<-timePts[t_index]
    if(mode == "X"){ # configure data differently for extended mode
      # Y: list containing t items where each item is also a list of length n, 
      # each of which is a matrix with S columns and J_(i,j) rows. The values are the ordinal abundance scores for each replicate at location i and time j for species i (column i). 
      Y[[t_index]]=lapply(X=sim_data_raw[[t]]$y, FUN=function(y_s){return(matrix(y_s, ncol=params$numSpecies, nrow = 1))})
    }else{
      #dataSchliep$Y=Y #list containing t items where each item is an n x S matrix of ordinal responses 
      Y[[t_index]]<-matrix(as.double(unlist(sim_data_raw[[t]]$y)), nrow=length(locList), ncol=params$numSpecies, byrow=T)
    }
    cov1<-sim_data_raw[[t]]$covs[[1]]
    cov2<-sim_data_raw[[t]]$covs[[2]]
    #cov3<-sim_data_raw[[t]]$covs[[3]]
    X[[t_index]]<-matrix(c(rep(1,times=length(cov1)),cov1,cov2), nrow=length(locList), ncol=params$numCovs+1, byrow=F)
  }
  
  dataSchliep=list()
  dataSchliep$t=length(timePts) #number of survey times
  dataSchliep$L=2 #number of ordinal categories (2= binary)
  dataSchliep$S=params$numSpecies #number of species (2 or more)
  dataSchliep$p=params$numCovs+1 #number of covariates
  dataSchliep$n=length(locList) #number of locations
  dataSchliep$dist=as.matrix(distMx) #distance matrix between locations
  dataSchliep$Y=Y #list containing t items where each item is an n x S matrix of ordinal responses 
  dataSchliep$X=X #list containing t items where each item is an n x p matrix of covariates
  if(mode == "X"){
    #J is matrix with n rows and t columns. Each J_(i,j) is the number of replicate observations for plot i and year j. Even if no replicates, this matrix needs to be provided with 1 (when the location is observed that year) or 0 (when the location is not observed that year)
    J = matrix(rep(x = as.integer(1), times = length(locList)*length(timePts)), nrow = length(locList), ncol = length(timePts))
    dataSchliep$J=J
  }
  
  # dataSchliep is list of the format of HWAdataM.Rdata (from paper supplement)
  #load("../schliep_code/HWAdataMBinary.Rdata")
  # below is from schliep supplement
  return(dataSchliep)
}

#priors=get.priors(HWAdataMBinary,temporalRE=T) 
#pars=get.startingValues(HWAdataMBinary,temporalRE=T) 
#run=Multivariate.Ordinal.Spatial.Model(HWAdataMBinary,pars,priors,temporalRE=T,iters=1000,print.out=100)

testData<-function(data, mode="N"){
  
  # Y: list containing t items where each item is also a list of length n, 
  # each of which is a matrix with S columns and J_(i,j) rows. The values are the ordinal abundance scores for each replicate at location i and time j for species i (column i). 
  
  if(length(data[["Y"]])!=data[["t"]]){
    print("t does not match length of Y")
  }
  if(length(data[["X"]])!=data[["t"]]){
    print("t does not match length of X")
  }
  for(time in 1:data[["t"]]){
    if(mode == "X"){
      if(length(unique(unlist(as.vector(data[["Y"]][[time]]))))!=data[["L"]]){
        print("number of ordinal categories (L) doesnt match values in Y")
      }
      if(dim(data[["Y"]][[time]][[1]])[1]!=data[["J"]][1,time]){
        print("dimension of Y not matched to J_ij")
      }
      if(dim(data[["Y"]][[time]][[1]])[2]!=data[["S"]]){
        print("dimension of Y not matched to S")
      }

      if(length(data[["Y"]][[time]])!=data[["n"]]){
        print("length of Y not matched to n")
      }
      if(typeof(data[["Y"]][[time]][[1]])!="integer"){
        print("type of Y[[t]] wrong")
      }
      if(dim(data[["J"]])[1] != data[["n"]]){
        print("J dimension wrong")
      }
      if(dim(data[["J"]])[2] != data[["t"]]){
        print("J dimension wrong")
      }
      if(typeof(data[["J"]])!="integer"){
        print("check typeof J" )
      }
      if(class(data[["J"]])[1]!="matrix"){
        print("check typeof J" )
      }
    } else {
      if(length(unique(as.vector(data[["Y"]][[time]])))!=data[["L"]]){
        print("number of ordinal categories (L) doesnt match values in Y")
      }
      if(dim(data[["Y"]][[time]])[1]!=data[["n"]]){
        print("dimension of Y not matched to n")
      }
      if(dim(data[["Y"]][[time]])[2]!=data[["S"]]){
        print("dimension of Y not matched to S")
      }
      if(typeof(data[["Y"]][[time]])!="double"){
        print("type of Y[[t]] wrong")
      }
    }
    if(dim(data[["X"]][[time]])[1]!=data[["n"]]){
      print("dimension of X not matched to n")
    } 
    if(dim(data[["X"]][[time]])[2]!=data[["p"]]){
      print("dimension of X not matched to p")
    } 
    if(typeof(data[["X"]][[time]])!="double"){
      print("type of X[[t]] wrong")
    }
  }
  if(typeof(data[["X"]])!="list"){
    print("X type wrong")
  }
  if(typeof(data[["Y"]])!="list"){
    print("Y type wrong")
  }
  if(typeof(data$dist)!="double"){
    print("dist type wrong")
  }
  if(!(all.equal(data[["dist"]],t(data[["dist"]])))){
    print("distance mx not symmetric")
  }
  if(dim(data[["dist"]])[1] != data[["n"]]){
    print("dist mx not the right dimensions")
  }
  if(dim(data[["dist"]])[2] != data[["n"]]){
    print("dist mx not the right dimensions")
  }
  if(class(data$dist)[1]!="matrix"){
    print("class of dist wrong")
  }
}

#############
# plot #
#############
plotPosteriorTable<-function(schliepOutput, ci_percent=0.95, mode = "N", burn_in=500){
  p=schliepOutput$data$p #number of covariates
  S=schliepOutput$data$S #number of species
  t=schliepOutput$data$t #number of time points
  n=schliepOutput$data$n #number of sites
  iter<-dim(schliepOutput$A)[3]
  
  credibleIntervalsL<-list()
  
  A=schliepOutput$A[,,burn_in:iter]
  for(species1 in 1:S){
    for(species2 in 1:S){
      varName = paste0("A",species1,species2)
      mean = mean(A[species1, species2, ])
      quants = quantile(x = A[species1, species2, ], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
      credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
    }
  }
  
  beta<-schliepOutput$beta
  if(mode == "X"){
    beta<-beta[burn_in:iter,]
    i=1
    for(species in 1:S){
      for(cov in 0:(p-1)){
        varName<-paste0("beta_cov",cov,"_sp",species)
        #### BUG FIXED
        mean = mean(beta[,i])
        quants = quantile(x = beta[,i], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
        credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
        i=i+1
      }
    }
  }else{
    beta<-beta[,,burn_in:iter]
    for(species in 1:S){
      for(cov in 0:(p-1)){
        varName<-paste0("beta_cov",cov,"_sp",species)
        mean = mean(beta[cov+1, species, ])
        quants = quantile(x = beta[cov+1, species, ], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
        credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
      }
    }
  }
  
  if(mode == "X"){
    rho <- schliepOutput$rho[,,burn_in:iter]
    for(species1 in 1:S){
      for(species2 in 1:S){
        varName = paste0("rho",species1,species2)
        mean = mean(rho[species1, species2, ])
        quants = quantile(x = rho[species1, species2, ], probs = c((1-ci_percent)/2,(1+ci_percent)/2))
        credibleIntervalsL[[varName]]<-c(as.numeric(mean), as.numeric(quants[1]), as.numeric(quants[2]))
      }
    }
  }
  
  credibleIntervals<-data.frame(t(data.frame(credibleIntervalsL)))
  names(credibleIntervals)<-c("mean", "lowBound", "highBound")
  credibleIntervals$varName = row.names(credibleIntervals)
  return(credibleIntervals)
}

############################### ############################### ############################### 
############################### RUN ############################### 
############################### ############################### ############################### 
mode = "X"
for(num_samples in c(100)){
for(trial in 1:2){
  for(run in 1:1){
    print(Sys.time())
    dataDir = paste0("./simData/2sp_2cov_competition/run1/")
    subdir = paste0("./2sp_2cov_competition/SchliepRes5000-500/trial",trial,"X/")
    dir.create(subdir)
    
    #inla_runData<-readRDS(paste0("./2sp_2cov_noInteractions/trial",trial,"/rundata-iid",trial,".Rdata"))
    
    sim_data_raw<-readRDS(paste0(dataDir,"sim_data.Rdata"))
    locList<-readRDS(paste0(dataDir,"locList.Rdata"))
    params<-readRDS(paste0(dataDir,"params.Rdata"))
    
    #put random 0s and 1s into a matrix to test Y
    #testY<-list()
    #for(time in 1:dataSchliepSim[["t"]]){
    #    set.seed(time)
    #    testY[[time]]<-matrix(data=rbern(n=dataSchliepSim[["n"]]*dataSchliepSim[["S"]], prob = .3), nrow = dataSchliepSim[["n"]], ncol = dataSchliepSim[["S"]])
    #}
    
    dataSchliepSim<-configureDataSchliep(subdir, sim_data_raw, locList, params, mode = mode, num_samples=num_samples)
    testData(dataSchliepSim, mode=mode)
    
    saveRDS(dataSchliepSim, paste0(subdir,"dataSchliep",trial,"Burn500_t",num_samples,".Rdata"))
    
    #################################################################
    ### General Bivariates spatio-temporal model for ordinal data ###
    #################################################################
    #required libraries: MCMCpack, fields, Matrix, spam, mvtnorm, bayesSurv, msm, fdrtool, tmvtnorm
    #MCMC Algorithmsn for general model:
    if(mode == "N"){
    
      #Get hyperparameters for priors; whether or not temporal random effects are included (True or False)
      priors=get.priors(dataSchliepSim,temporalRE=T) 
      
      #Get starting values for the algorithm; whether or not temporal random effects are included (True or False)
      pars=get.startingValues(dataSchliepSim,temporalRE=T) 
      
      #Run the model: Functions take the data, starting values for parameters, prior info, number of iterations to run, how often to print out the iteration number so you can monitor progress, and whether or not temporal random effects are included (True or False)
      schliepRes=Multivariate.Ordinal.Spatial.Model(dataSchliepSim,pars,priors,temporalRE=T,iters=numIters,print.out=500)
      
      #Output (in list form): 
      #data: returns the data that was input 
      #priors: returns the hyperparameters of the priors 
      #alpha: array of dimension (t x S x (iters+1)) of posterior samples (if temporalRE=T)
      #beta: array of dimension (p x S x (iters+1)) of posterior samples
      #lambda: if L>2, array of dimension (S x (L-1) x (iters+1)) of posterior samples
      #A: array of dimension (S x S x (iters+1)) of posterior samples
      #phi: matrix of dimension ((iters+1) x S) of poseterior samples
      #K: list of S arrays, each of dimension (n  x t x (iters+1)) of posterior samples
      
    }else if(mode == "X"){
      
      #Get hyperparameters for priors; whether or not temporal random effects are included (True or False), temporal dynamics (True or False), and replicate observations (True or False)
      priorsX=get.priorsX(dataSchliepSim,temporalRE=T,temporalDyn=T,replicates=F)
      
      #Get starting values for the algorithm; whether or not temporal random effects are included (True or False), temporal dynamics (True or False), and replicate observations (True or False)
      parsX=get.startingvaluesX(dataSchliepSim,temporalRE=T,temporalDyn=T,replicates=F)
      
      #Run the model: Functions take the data, starting values for parameters, prior info, number of iterations to run, how often to print out the iteration number so you can monitor progress, whether or not to include cross species spatial dependence, temporal random effects, temporal dynamics, cross species temporal dynamics, and replicate observations
      schliepRes=Multivariate.Ordinal.Spatial.ModelX(dataSchliepSim,parsX,priorsX,iters=numIters,print.out=500,spDep=T,temporalRE=T,temporalDyn=T,crossDyn=T,replicates=F)
      
      
      ########################################
      ### Compute Rank Probability Scores: ###
      ########################################
      
      #RPS takes the following arguments
      # 1) output from Bivariate.Ordinal.Spatial.ModelX
      # 2) beginning value of the chain to take posterior draws from, e.g., if b=100, first 99 draws will be burn-in
      # 3) ending value of the chain to take posterior darws from
      # 4) number of posterior samples to use to evaluate RPS
      
      #rps=RPS(run,b=5,e=10,iters=3)
      #schliepRes$RPS<-rps
      #output is an array of dimension n x t x S
  
      #Output (in list form): 
      #data: returns the data that was input 
      #priors: returns the hyperparameters of the priors 
      #alpha: array of dimension (t x S x (iters+1)) of posterior samples (if temporalRE=T)
      #beta: matrix of dimension (iters+1 x (p*S)) of posterior samples
      #lambda: if L>2, array of dimension (S x (L-1) x (iters+1)) of posterior samples
      #A: array of dimension (S x S x (iters+1)) of posterior samples
      #phi: matrix of dimension ((iters+1) x S) of posterior samples
      #K: matrix of dimension (n*S  x t x (iters+1)) of posterior samples
      #delta: matrix of dimension ((iters+1) x S) of posterior samples
      #rho: array of dimension (S x S x (iters+1)) of posterior samples
      #Omega: array of dimension (n x S^2 x (iters+1)) of posterior samples; element Omega[i,,j] is the vectorized covariance matrix for location i and iteration j.
      #RPS: array of dimension n x t x S
      
    }else{
      stop("mode not valid in Schliep_sim")
    }
    
    saveRDS(schliepRes, paste0(subdir,"SchliepRes", trial,"Burn500_t",num_samples,".Rdata"))
    
    #schliepRes<-readRDS("./fiona_simulation/inlaRun_SchliepcodeComparison/ScliepRes.Rdata")
    
    posteriorTable<-plotPosteriorTable(schliepRes, mode = mode, ci_percent = 0.95, burn_in =burn)
    posteriorTable$color <- sign(posteriorTable$lowBound) * sign(posteriorTable$highBound)
    
    p1 <- ggplot(data=posteriorTable) + 
      geom_point(aes(y=varName,x=mean, col=as.factor(color)), show.legend=F) +
      geom_errorbarh(aes(y=varName,xmin =lowBound,xmax=highBound, col=as.factor(color)), show.legend=F) +
      geom_vline(xintercept = 0) +
      scale_color_manual(values=c("#999999", "#999999", "#56B4E9")) +
      ggtitle(paste0("Credible Intervals for Parameter Posteriors"))

    ggsave(filename = paste0(subdir, "CI_plot_schliep",trial,"Burn500_t",num_samples,".png"), plot = p1)
    
  print(Sys.time()) 
  }
}
}


