library(ggplot2) # plotting

configureDataSchliep <- function(subdir, sim_data_raw, locList, params, mode = "N", num_samples=num_samples) {
  #timePts<-1:params$numGens #time points to "collect"
  #timePts <- seq(from = 1, to = params$numGens, by = 1000/num_samples)

  
  # Wrangle data
  locDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = T))
  distMx <- dist(locDF, diag = T, upper = T)
  
  Y <- list()
  X <- list()
  for (t_index in 1:length(timePts)) {
    t <- timePts[t_index]
    if (mode == "X") { # configure data differently for extended mode
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