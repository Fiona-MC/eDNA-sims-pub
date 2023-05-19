# the purpose of this is to use the siumulation data in the code from Schliep paper  (/Users/fionacallahan/Documents/Nielsen_lab/schliep_code/RunMultivariateOrdinalSpatialAlgorithm.R)
# (DOI: 10.1111/geb.12666)
library(stats)
library(ggplot2) # plotting
library(coda) # mcmc diagnostics

#MCMC Algorithms for model with extensions 
args <- commandArgs(trailingOnly = TRUE)

source("/home/fiona_callahan/schliep_code/MultivariateSpatialOrdinalAlgorithm.R")
library(Rcpp) #all three of these packages are for c++ to work
library(RcppArmadillo)
library(inline) 
source("/home/fiona_callahan/schliep_code/MultivariateSpatialOrdinalAlgorithmExtensions.R")
source("/home/fiona_callahan/schliep_code/MultivariateExtensionsRcppCode.R")
source("/home/fiona_callahan/eDNA_sims_code/Schliep_sims_fcns.R")


numIters <- 100
burn <- 20
dataDir <- "/space/s1/fiona_callahan/multiSim11/randomRun10/"

dataDir <- args[1]
numIters <- args[2]
burn <- args[3]

numTrials <- 2
mode <- "X"
numRums <- 1
plot <- TRUE

#Rscript /home/fiona_callahan/eDNA_sims_code/schliep_sims.R /space/s1/fiona_callahan/multiSim11/randomRun10/ 10000 2000

######### RUN #########
dir.create(paste0(dataDir, "SchliepRes/"))
sim_data_raw <- readRDS(paste0(dataDir, "sim_data.Rdata"))
INLA_sitetab <- read.csv(paste0(dataDir, "sim_sitetab_sampled.csv"), header = TRUE)
locList <- readRDS(paste0(dataDir, "locList.Rdata"))
params <- readRDS(paste0(dataDir, "params.Rdata"))

timePts <- sort(params$num_gens - unique(INLA_sitetab$Age))[-1] # take out first time point because N_0 is 10 and that is too low for detection in most cases

locDF <- unique(INLA_sitetab[, c("Lat", "Long")])
locList.asDF <- as.data.frame(matrix(unlist(locList), nrow = length(locList), ncol = 2, byrow = TRUE))
locIndices <- unname(apply(locDF, MARGIN = 1, FUN = function(latLon) {
                                              return(which(apply(locList.asDF, 1, function(row) {all(round(row, 2) == round(latLon, 2))} )))
                                              }))
# testing that above works
#locations <- locList.asDF[locIndices,]
#round(locations, 4)==round(locDF, 4)

#put random 0s and 1s into a matrix to test Y
#testY<-list()
#for(time in 1:dataSchliepSim[["t"]]){
#    set.seed(time)
#    testY[[time]]<-matrix(data=rbern(n=dataSchliepSim[["n"]]*dataSchliepSim[["S"]], prob = .3), nrow = dataSchliepSim[["n"]], ncol = dataSchliepSim[["S"]])
#}


for (trial in 1:numTrials) {
  outdir = paste0(dataDir, "SchliepRes/trial",trial,"/")
  dir.create(outdir)

  dataSchliepSim <- configureDataSchliep(outdir, sim_data_raw, locList, params, mode = mode, timePts = timePts, locIndices = locIndices)
  if (testData(dataSchliepSim, mode = mode) == FALSE) {
    print("BAD DATA")
  }

  saveRDS(dataSchliepSim, paste0(outdir,"dataSchliep.Rdata"))

  # save info about this run
  schliepParms <- list(dataDir = dataDir, numIters = numIters, burn = burn, mode = mode)
  saveRDS(schliepParms, paste0(outdir, "schliepParms.Rdata"))

  #################################################################
  ### General Bivariates spatio-temporal model for ordinal data ###
  #################################################################
  #required libraries: MCMCpack, fields, Matrix, spam, mvtnorm, bayesSurv, msm, fdrtool, tmvtnorm
  #MCMC Algorithmsn for general model:
  if (mode == "N") {
    #Get hyperparameters for priors; whether or not temporal random effects are included (True or False)
    priors = get.priors(dataSchliepSim,temporalRE=T) 
    
    #Get starting values for the algorithm; whether or not temporal random effects are included (True or False)
    pars = get.startingValues(dataSchliepSim,temporalRE=T) 
    
    #Run the model: Functions take the data, starting values for parameters, prior info, number of iterations to run, how often to print out the iteration number so you can monitor progress, and whether or not temporal random effects are included (True or False)
    schliepRes = Multivariate.Ordinal.Spatial.Model(dataSchliepSim,pars,priors,temporalRE=T,iters=numIters,print.out=100)
    
    #Output (in list form): 
    #data: returns the data that was input 
    #priors: returns the hyperparameters of the priors 
    #alpha: array of dimension (t x S x (iters+1)) of posterior samples (if temporalRE=T)
    #beta: array of dimension (p x S x (iters+1)) of posterior samples
    #lambda: if L>2, array of dimension (S x (L-1) x (iters+1)) of posterior samples
    #A: array of dimension (S x S x (iters+1)) of posterior samples
    #phi: matrix of dimension ((iters+1) x S) of poseterior samples
    #K: list of S arrays, each of dimension (n  x t x (iters+1)) of posterior samples
    
  } else if (mode == "X") {
    
    #Get hyperparameters for priors; whether or not temporal random effects are included (True or False), temporal dynamics (True or False), and replicate observations (True or False)
    priorsX = get.priorsX(dataSchliepSim,temporalRE=T,temporalDyn=T,replicates=F)
    
    #Get starting values for the algorithm; whether or not temporal random effects are included (True or False), temporal dynamics (True or False), and replicate observations (True or False)
    parsX = get.startingvaluesX(dataSchliepSim,temporalRE=T,temporalDyn=T,replicates=F)
    
    #Run the model: Functions take the data, starting values for parameters, prior info, number of iterations to run, how often to print out the iteration number so you can monitor progress, whether or not to include cross species spatial dependence, temporal random effects, temporal dynamics, cross species temporal dynamics, and replicate observations
    schliepRes = Multivariate.Ordinal.Spatial.ModelX(dataSchliepSim,parsX,priorsX,iters=numIters,print.out=100,spDep=T,temporalRE=T,temporalDyn=T,crossDyn=T,replicates=F)
    
    
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
    
  } else {
    stop("mode not valid in Schliep_sim")
  }
  
  saveRDS(schliepRes, paste0(outdir,"SchliepRes.Rdata"))
  
  if (plot) {
    posteriors <- plotPosteriorTable(schliepRes, mode = mode, ci_percent = 0.95, burn_in = burn)
    posteriorTable <- posteriors$credibleIntervals
    posteriorTable$color <- sign(posteriorTable$lowBound) * sign(posteriorTable$highBound)
    
    p1 <- ggplot(data = posteriorTable) + 
      geom_point(aes(y = varName,x = mean, col = as.factor(color)), show.legend = FALSE) +
      geom_errorbarh(aes(y = varName, xmin = lowBound, xmax = highBound, col = as.factor(color)), show.legend = FALSE) +
      geom_vline(xintercept = 0) +
      scale_color_manual(values = c("#999999", "#999999", "#56B4E9")) +
      ggtitle(paste0("Credible Intervals for Parameter Posteriors"))

    ggsave(filename = paste0(outdir, "CI_plot_schliep.png"), plot = p1)
  }
}

