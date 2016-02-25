# #' Calling the Sampler (JAGS)
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
                           X_list=list(),   # list with design matrices for fixed effects
                           hyperpriors=NULL,
                           n.iter=100000,
                           n.burnin=NULL,
                           n.update= 3,
                           n.thin=2,
                           n.chains=3,
                           autojags=FALSE,
                           # savetable = NULL,
                           ...){

  if(is.na(n.burnin)){n.burnin=n.iter/2}

  if(model == "betaMPT"){
    parameters <- list("theta", "alph","bet", "mean", "sd")
  }else{
    parameters <- list("theta", "mu", "mean", "rho", "sigma")
  }

  responses <- data
  subjs <- nrow(responses)


  ### prepare data for WinBUGS

  treeNames=unique(mergedTree$Tree)
  NresponsesTree=vector("numeric",length=length(treeNames))
  parameters <- c(parameters, paste("response",treeNames, "pred.mean", sep="."))

  for(i in 1:length(treeNames)){
    NresponsesTree[i]=length(which(mergedTree$Tree==treeNames[i]))
  }

  # make character vector with data object names
  data <- c("subjs", "S")
  index=0
  for(i in 1:length(treeNames)){
    ### variable names in JAGS
    name.mean <-     paste("response",treeNames[i],"mean",sep=".")
    name.response <- paste("response",treeNames[i],sep=".")
    name.items <- paste("items",treeNames[i],sep=".")

    assign(name.items, rowSums(responses[(index+1):(index+NresponsesTree[i])]))

    # check whether any N=0
    if (any(rowSums(responses[(index+1):(index+NresponsesTree[i])]) == 0))
      warning("One or more responses do not have responses for tree",
              treeNames[i], ". As a solution, the critical participants might be excluded.")

    assign(name.response,  matrix(as.vector(t(responses[(index+1):(index+NresponsesTree[i])])),
                                  ncol=NresponsesTree[i],nrow=subjs,byrow=TRUE))

    # mean frequencies for T1 statistic
    assign(name.mean, colMeans(get(name.response)))

    index=index+NresponsesTree[i]

    data <- c(data, name.mean, name.response, name.items)
  }


  if(model == "traitMPT"){
    data <- c(data, "V", "df")
    df <- hyperpriors$df
    V <- hyperpriors$V

    if(length(X_list) != 0){
      for(pp in 1:length(X_list))
          assign(names(X_list)[pp], X_list[[pp]])
      data <- c(data, names(X_list))
    }
  }
  if(!is.null(covData)){
    covSD <- apply(covData, 2, sd)
    data <- c(data, "covData", "covSD")
  }

  # call Sampler

  # if(sampler %in% c("jags","JAGS")){
    parametervector=c(unlist(parameters),
                      "T1.obs","T1.pred","p.T1","p.T1ind",
                      transformedPar, covPars)
    # random initial values: required (otherwise T1 statistic results in errors: no variance!)
    if(model == "betaMPT"){
      inits <- function() list("theta"=matrix(runif(subjs*S), S, subjs)
                               )
    }else{
      inits <- function() list("delta.part.raw" = matrix(rnorm(subjs*S, 0,1), S, subjs),
                               "xi"=runif(S,.9,1.1),
                               "T.prec.part" = rWishart(1,df,V)[,,1]
                               )
    }
    samples <- jags.parallel(data,
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
      cat("#####################################\n
#### Autojags started. Might require some time
#### (use n.update to adjust maximum number of updates). See ?autojags\n
#####################################\n")
      recompile(samples, n.iter=n.iter)
      samples.upd <- autojags(samples, n.update = n.update)
      samples=samples.upd
    }

#   }else{
#       print(paste("Unknown sampler:",sampler))
#     }


  return(samples)
}
