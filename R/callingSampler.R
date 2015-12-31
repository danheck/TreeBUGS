# #' Calling the Sampler (JAGS or OpenBUGS or WinBUGS)
# #'
# #' @param model either "betaMPT", "traitMPT"
# #' @param mergedTree Ouptut of mergeBranches()
# #' @param data data frame with frequencies
# #' @param modelfile Full path to the model description file. (Made with makeModelDescription())
# #' @param S How many parameters does the model have
# #' not necessary:  @param parameters List with parameters that should be returned from the sampler.
# #' @param transformedPar vector with names of transformed parameter to be sampled (e.g., parameter differences)
# #' @param n.iter Number of iterations.
# #' @param n.burnin Burnin period.
# #' @param n.update Update paramater for JAGS
# #' @param n.thin Thinning rate.
# #' @param n.chains number of MCMC chains
# #' @param sampler The Sampler to be used. Options are JAGS, OpenBUGS or WinBUGS. Note that you need to install the sampler on # your computer.
# #' @param autojags whether to run JAGS until the MCMC chains converge (see \link{autojags}). Can take a lot of time for large # models.
# #' not necessary:  @param savetable name of the file to which the results should be saved on the hard drive.
# #' @param ... Arguments to be passed to other methods.
# #' @author Nina R. Arnold, Denis Arnold, Daniel W. Heck
# #' @export
callingSampler <- function(model,
                           mergedTree,
                           data,
                           modelfile,
                           S,
                           transformedPar=NULL,
                           covPars=NULL,
                           covData=NULL,
                           hyperpriors=NULL,
                           n.iter=100000,
                           n.burnin=NULL,
                           n.update= 10,
                           n.thin=2,
                           n.chains=3,
                           sampler="JAGS",
                           autojags=TRUE,
                           # savetable = NULL,
                           ...){

  if(is.na(n.burnin)){n.burnin=n.iter/2}
  parameters <- list("theta", "alph")
  if(model == "betaMPT"){
    parameters <- c(parameters, list("alph","bet", "mean", "sd"))
  }else{
    parameters <- c(parameters, list("mu", "mean", "rho", "sigma"))
  }

  responses <- data
  subjs <- nrow(responses)


  ### prepare data for WinBUGS

  treeNames=unique(mergedTree$Tree)
  NresponsesTree=vector("numeric",length=length(treeNames))

  for(i in 1:length(treeNames)){
    NresponsesTree[i]=length(which(mergedTree$Tree==treeNames[i]))
  }

  index=0
  for(i in 1:length(treeNames)){
    assign(paste("items",treeNames[i],sep="."),
           rowSums(responses[(index+1):(index+NresponsesTree[i])]))
    assign(paste("response",treeNames[i],sep="."),
           matrix(as.vector(t(responses[(index+1):(index+NresponsesTree[i])])),
                  ncol=NresponsesTree[i],nrow=subjs,byrow=TRUE))
    index=index+NresponsesTree[i]
  }

  # make data
  data <- c(paste("items",treeNames, sep="."),
            paste("response",treeNames, sep="."),
            "subjs", "S")
  if(model == "traitMPT"){
    data <- c(data, "V", "df")
    df <- hyperpriors$df
    V <- hyperpriors$V
  }
  if(!is.null(covData)){
    if(is.character(covData)){
      covData <- read.csv(covData, header=T, sep= ",",  strip.white = TRUE)
    }
    covSD <- apply(covData, 2, sd)
    data <- c(data, "covData", "covSD")
  }
#   for(i in 1:length(treeNames)){
#     data <- c(data,
#               paste("items",treeNames[i],sep="."),
#               paste("response",treeNames[i],sep="."))
#   }
  # data <- c(data, "subjs", "S")


  # call Sampler

  if(sampler%in%c("jags","JAGS")){
    # W=diag(S)
    # data <- c(data, "W")
    parametervector=c(unlist(parameters),
                      "T1.obs","T1.pred","p.T1","p.T1ind",
                      transformedPar, covPars)
    # random initial values
    inits <- function() list(theta=matrix(runif(subjs*S), S, subjs),
                             alpha=rgamma(S, 2, 1),
                             beta=rgamma(S, 2, 1))
    samples = jags.parallel(data,
                            inits=inits,
                            parameters.to.save=parametervector,
                            model.file = modelfile,
                            n.iter=n.iter,
                            n.burnin=n.burnin,
                            n.chains=n.chains,
                            DIC=T,
                            envir=environment(),
                            ...)
    if(autojags){
      recompile(samples)
      samples.upd <- autojags(samples, n.update = n.update)
      samples=samples.upd
    }

  }else{
    if(sampler%in%c("openbugs","OpenBUGS","winbugs","WinBUGS")){

      samples = bugs(data,
                     inits=NULL,
                     parameters,
                     model.file = modelfile,
                     n.chains=n.chains,
                     n.iter=n.iter,
                     n.burnin=n.burnin,
                     n.thin=n.thin,
                     DIC=T,
                     codaPkg=F,
                     debug=F,
                     program = sampler,
                     ...)
      # get coda samples (better to process afterwards)
      # tmp <- capture.output({samples <- read.bugs(out)})
    }else{
      print(paste("Unknown sampler:",sampler))
    }
  }

  return(samples)
}
