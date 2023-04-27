source("./inla_code/INLA_ST_func.R")
num_samples = 100
#for(num_samples in c(10,20,50,100)){

for(trial in 1:3){
  subdir = paste0("./sim_INLAres_fpr/independent/trial", trial, "/")
  dataDir = "./simData/FALSE_POS/independent_FPR/run1/"
  dir.create(subdir)
  
  sim_data_raw<-readRDS(paste0(dataDir, "sim_data.Rdata"))
  locList<-readRDS(paste0(dataDir,"locList.Rdata"))
  params<-readRDS(paste0(dataDir,"params.Rdata"))
  
  #timePts<-1:params$numGens #time points to "collect"
  timePts<-seq(from=1, to=params$numGens, by = 1000/num_samples)
  
  # Wrangle data
  sitetab_data<-numeric(length(timePts)*length(locList)*4+params$numCovs+params$numSpecies)
  i=1
  for(t in timePts){
    for(s_index in 1:length(locList)){
      sitetab_data[i] <- s_index
      sitetab_data[i+1] <- params$numGens-t
      sitetab_data[i+2] <- locList[s_index][[1]][1]
      sitetab_data[i+3] <- locList[s_index][[1]][2]
      sitetab_data[i+4] <- sim_data_raw[[t]]$covs[[1]][s_index]
      sitetab_data[i+5] <- sim_data_raw[[t]]$covs[[2]][s_index]
      #sitetab_data[i+6] <- sim_data_raw[[t]]$covs[[3]][s_index]
      for(k in 1:params$numSpecies){
        sitetab_data[i+3+params$numCovs+k] <- sim_data_raw[[t]]$y[[s_index]][k]
      }
      i=i+4+params$numCovs+params$numSpecies
    }
  }
  
  sim_sitetab<-as.data.frame(matrix(data = sitetab_data, nrow=length(timePts)*length(locList), ncol=4+params$numCovs+params$numSpecies, byrow=T))
  names(sim_sitetab)<-c("labID", "Age", "Lat", "Long", "Cov1", "Cov2", "Sp1", "Sp2")# "Sp4", "Sp5")
  
  
  # Select what you'll be working with
  k <- 11 # k is number of time knots
  #k <- 4
  sitetab<-sim_sitetab
  write.csv(sitetab, paste0(subdir,"sitetab.csv"))
  # Prepare spatial and temporal meshes
  # alltidx <- seq(1,k)
  #reducedtidx <- seq(2,k,2)
  tknots <- seq(min(sitetab$Age), max(sitetab$Age), length = k)
  #reducedtknots <- tknots[seq(2,k,2)]
  mesh.t <- inla.mesh.1d(tknots)
  
  xytab <- cbind(sitetab$Long,sitetab$Lat)
  #xytab <- LonLatToAzi(xytab)
  n <- dim(sitetab)[1]
  
  bound <- inla.nonconvex.hull(as.matrix(xytab))
  #TODO what are these params??
  #mesh.s <- inla.mesh.2d(boundary = bound, max.edge = c(200, 400), max.n=200, cutoff = 500)
  mesh.s <- inla.mesh.2d(boundary = bound, max.edge = c(10, 25))#, max.edge = c(2, 4), max.n=200, cutoff = 500)
  plot(mesh.s)
  spde.s <- inla.spde2.matern(mesh.s)
  
  
  
  goodanimalnames<-c("Sp1", "Sp2")#, "Sp3")#, "Sp4", "Sp5")
  # this takes a long time even with only 10 time points
  # Model with no covariates
  reslist <- list()
  for(response in goodanimalnames){ # NOTE -- this takes a few minutes
    
    # Run INLA - Binomial model with number of trials = 1 (no covariates)
    reslist[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=c(),normcov=0,project = F)
  }
  
  #reslist1<-reslist
  # Define lists to store results
  # ADD MODELS
  resAll <- list()
  
  resCovOnly <- list()
  
  for(response in goodanimalnames){ 
    # ADD MODELS
    # Run INLA 
    #response="Sp1"
    namesCov<-c("Cov1","Cov2")#, "Cov3")
    resCovOnly[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesCov,normcov=0, project = F)
    
    # Run with other animals as covariates
    namesCov<-c(namesCov,goodanimalnames)
    namesCov<-namesCov[namesCov!=response]
    resAll[[response]] <- RunInlaBin(response,sitetab,xytab,mesh.s,mesh.t,1,namescov=namesCov,normcov=0, project = F)
  }
  
  
  
  simLists <- list()
  simLists[["None"]] <- reslist
  #simLists[["None2"]] <- reslist1
  # ADD MODELS
  simLists[["AllCov"]]<- resCovOnly
  simLists[["Cov+Animal"]]<- resAll
  
  #save results for later
  saveRDS(simLists, paste0(subdir, "simLists-iid-",num_samples,"-",trial,".Rdata"))
  runData<-list()
  runData$params<-params
  runData$sim_data_raw<-sim_data_raw
  runData$sitetab<-sitetab
  saveRDS(runData,paste0(subdir,"rundata-iid-",num_samples,"-",trial,".Rdata"))
  # CPO 
  modelcands <- names(simLists)
  #modelcands <- c("None", "AllCov", "Cov+Animal")
  
  # modelcands<- modelcands[c(F,T,T,T,T,T,T,T,T,T)]
  # next line not working -- I think because the no covariates model didnt work on all of the animals -- Not sure why though
  # only hare, wolf, reindeer, and horse
  # nevermind that didnt work
  # in goodanimalnames, there is caribou but the name in animallists is reindeer -- these are the same animal
  # animal table uses reindeer
  # fixed-- not sure how this happened, can't find anywhere we could have gotten caribou...
  cpotab <- sapply(goodanimalnames, function(response){
    return(sapply( modelcands, function(modelname){ -sum(log(simLists[[modelname]][[response]]$cpo$cpo),na.rm=TRUE) }))
  })
  
  ### debugging
  #for(response in goodanimalnames){
  ##  for(modelname in modelcands){
  #    print(modelname)
  #   print(response)
  #  -sum(log(animallists[[modelname]][[response]]$cpo$cpo),na.rm=TRUE)
  #}
  #}
  ###
  
  
  rownames(cpotab) <- modelcands
  ############ ??????????????
  #colnames(cpotab)[which(colnames(cpotab) == "Reindeer")] <- "Caribou"
  cpotabmelted <- melt(cpotab)
  colnames(cpotabmelted) <- c("Model","Animal","CPO")
  cpoplot <- ggplot(cpotabmelted) +
    geom_point(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
    geom_line(aes(x=Model,y=CPO,group=Animal,colour=Animal)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.title.x=element_blank(),
          text = element_text(size = 10)
    )
  
  # waic
  #modelcands <- names(animallists)
  waictab <- sapply(goodanimalnames, function(response){
    return(sapply( modelcands, function(modelname){simLists[[modelname]][[response]]$waic$waic}))
  })
  rownames(waictab) <- modelcands
  waictabmelted <- melt(waictab)
  colnames(waictabmelted) <- c("Model","Animal","WAIC")
  waicplot <- ggplot(waictabmelted) +
    geom_point(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
    geom_line(aes(x=Model,y=WAIC,group=Animal,colour=Animal)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x=element_blank(),
          text = element_text(size = 10)
    )
  
  
  cpoplotlet <- arrangeGrob(cpoplot, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                                    gp=gpar(col="black", fontsize=22)))
  waicplotlet <- arrangeGrob(waicplot, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),
                                                      gp=gpar(col="black", fontsize=22)))
  model_plot<-grid.arrange(cpoplotlet,waicplotlet, ncol=2, vp=viewport(width=0.98, height=0.98))
  ggsave(model_plot, file=paste0(subdir,"model_comparisons-iid-",num_samples,"-",trial,".png"))
  
  # Best model
  covarplots <- CreateCovarPlots(scoretab = waictab,names = goodanimalnames,modellist = simLists)
  bestWAIC<-grid.arrange(covarplots[[1]],covarplots[[2]], vp=viewport(width=0.98, height=0.98))
  ggsave(bestWAIC, file=paste0(subdir,"bestWAIC_INLA-iid-",num_samples,"-",trial,".png"))
  
  # Model without all cov NOT animal
  covarplots <- CreateCovarPlots(waictab,goodanimalnames,simLists,forcemodel="AllCov")
  covOnly<-grid.arrange(covarplots[[1]],covarplots[[2]], vp=viewport(width=0.98, height=0.98))
  ggsave(covOnly, file=paste0(subdir,"covonly_INLA-iid-",num_samples,"-",trial,".png"))
  
  # Model with cov and animal
  covarplots <- CreateCovarPlots(waictab,goodanimalnames,simLists,forcemodel="Cov+Animal")
  animalCov<-grid.arrange(covarplots[[1]],covarplots[[2]],vp=viewport(width=0.98, height=0.98))
  ggsave(animalCov, file=paste0(subdir,"animalCov_INLA-iid-",num_samples,"-",trial,".png"))
}
#}