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
                           fixedPar=NULL,
                           transformedPar=NULL,
                           covPars=NULL,
                           covData=NULL,
                           X_list=list(),   # list with design matrices for fixed effects
                           # groupT1=NULL,    # list with groupMatT1 und NgroupT1  for splitted T1 statistic
                           hyperpriors=NULL,
                           n.iter=20000,
                           n.adapt=2000,
                           n.burnin=2000,
                           n.update= 2000,
                           n.thin=5,
                           n.chains=3,
                           autojags=NULL,
                           # savetable = NULL,
                           ...){

  if(is.na(n.burnin)){n.burnin=n.iter/2}

  if(model == "betaMPT"){
    parameters <- list("theta", "alph","bet", "mean", "sd")
  }else{
    parameters <- list("theta", "mu", "mean", "rho", "sigma")
  }
  if(!is.null(fixedPar)){
    parameters <- c(parameters, "thetaFE")
  }

  responses <- data
  subjs <- nrow(responses)


  ### prepare data for JAGS
  treeNames=unique(mergedTree$Tree)
  NresponsesTree=vector("numeric",length=length(treeNames))

  for(i in 1:length(treeNames)){
    NresponsesTree[i]=length(which(mergedTree$Tree==treeNames[i]))
  }

  # make character vector with data object names
  data <- c("subjs", "S")
  index=0
  for(i in 1:length(treeNames)){

    ### variable names in JAGS
    name.response <- paste("response",treeNames[i],sep=".")
    name.items <- paste("items",treeNames[i],sep=".")

    assign(name.items, rowSums(responses[(index+1):(index+NresponsesTree[i])]))

    # check whether any N=0
    rowsums <- rowSums(responses[(index+1):(index+NresponsesTree[i])])
    if (any(is.na(rowsums)) || any(rowsums == 0))
      warning("One or more participants do not have responses for tree",
              treeNames[i], ". As a solution, the critical participants might be excluded.")

    assign(name.response,  matrix(as.vector(t(responses[(index+1):(index+NresponsesTree[i])])),
                                  ncol=NresponsesTree[i],nrow=subjs,byrow=TRUE))

    index=index+NresponsesTree[i]

    data <- c(data,
              name.response, name.items)
  }


  if(model == "traitMPT"){
    df <- hyperpriors$df
    V <- hyperpriors$V
    data <- c(data, "V", "df")

    if(length(X_list) != 0){
      for(pp in 1:length(X_list))
        assign(names(X_list)[pp], X_list[[pp]])
      data <- c(data, names(X_list))
    }
  }

  if(!is.null(covData) & !any(dim(covData) == 0)){

    covData <- as.matrix(covData)
    if(any(is.na(covData))){
      warning("Data frame with covariates contains missing values (NA).",
              "\n  This is likely to cause problems for JAGS.")
    }

    # covData only required for predictors/discrete factors:
    if(!is.null(covPars))
      data <- c(data, "covData")
  }


  # call Sampler
  parametervector=c(unlist(parameters),
                    transformedPar, covPars)

  # random initial values: required (otherwise T1 statistic results in errors: no variance!)
  if(model == "betaMPT"){
    inits <- function() list("theta"=matrix(runif(subjs*S), S, subjs)
    )
  }else{
    # draw appropriate random starting values:
    mu <- xi <- rep(NA, S)
    for(s in 1:S){
      tmp <- ifelse(length(hyperpriors$mu)<=1,hyperpriors$mu[1],hyperpriors$mu[s])
      mu[s] <- eval(parse(text=sub("d","r", sub("(","(1,", tmp,  fixed=TRUE))))
      tmp <- ifelse(length(hyperpriors$xi)<=1,hyperpriors$xi[1],hyperpriors$xi[s])
      xi[s] <- eval(parse(text=sub("d","r", sub("(","(1,", tmp,  fixed=TRUE))))
    }

    inits <- function() list("delta.part.raw" = matrix(rnorm(subjs*S, -1,1), S, subjs),
                             "xi"=xi,  "mu" = mu,
                             "T.prec.part" = as.matrix(rWishart(1,df,V)[,,1])
    )
  }
  inits.list <- replicate(n.chains, inits(), simplify=FALSE)
  for(i in 1:length(inits.list))
    inits.list[[i]]$.RNG.name <- c("base::Wichmann-Hill",
                                   "base::Marsaglia-Multicarry",
                                   "base::Super-Duper",
                                   "base::Mersenne-Twister")[1+ (i-1)%% 4]
  #         samples2 <- jags.parallel(data,
  #                             inits=inits,
  #                             parameters.to.save=parametervector,
  #                             model.file = modelfile,
  #                             n.iter=n.iter,
  #                             n.burnin=n.burnin,
  #                             n.chains=n.chains,
  #                             DIC=T,
  #                             envir=environment(),
  #                             ...)
  n.samples <- ceiling((n.iter-n.burnin)/n.thin)
  choice <- ""
  if(n.samples > 10000){
    choice <- readline(prompt =
                         paste0("\n############ Warning ################\n",
                                "Your present MCMC settings for n.burnin/n.iter/n.thin\n",
                                "imply that more than 10,000 samples are stored per parameter per chain.\n",
                                "This might result in problems due to an overload of your computers memory (RAM).\n",
                                "If you are sure you want to continue, press <RETURN>."))
  }
  if(choice != "")
    stop("Model fitting terminated by user.")


  data.list <-  lapply(data, get, envir=environment())
  names(data.list) <- data
  samples <- run.jags(model = modelfile,
                      monitor=c(parametervector, "deviance"),
                      data=data.list,
                      n.chains=n.chains,
                      inits=inits.list,
                      burnin=n.burnin,
                      adapt=n.adapt,
                      sample=n.samples,
                      thin=n.thin,
                      modules=c("dic","glm"),
                      summarise=FALSE,
                      method="parallel",
                      ...)


  if(!is.null(autojags)){
    cat("#####################################\n
        #### Autojags started. Might require some time
        #### (use max.time to adjust maximum time for updating). See ?autoextend.jags\n
        #####################################\n")
    #       recompile(samples, n.iter=n.iter)
    #       samples <- autojags(samples, n.update = n.update)
    samples <- do.call(autoextend.jags,
                       c(list(runjags.object = samples,
                              summarise=FALSE),
                         autojags))  # additional user arguments
  }

  return(samples)   # own summary
}
