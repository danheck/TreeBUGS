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

    data <- c(data, name.response, name.items)
  }


  if (model == "traitMPT"){
    df <- hyperpriors$df
    V <- hyperpriors$V
    data <- c(data, "V", "df")

    if (length(X_list) != 0){
      for (pp in 1:length(X_list))
        assign(names(X_list)[pp], X_list[[pp]])
      data <- c(data, names(X_list))
    }
  }

  if (!is.null(covData) & !any(dim(covData) == 0)){

    covData <- as.matrix(covData)
    if (anyNA(covData)){
      warning("Data frame with covariates contains missing values (NA).",
              "\n  This is likely to cause problems for JAGS.")
    }

    # covData only required for predictors/discrete factors:
    if (!is.null(covPars)){
      # standardization:
      covVar <- diag(cov(covData))
      # covData <- as.matrix(scale(covData))
      data <- c(data, "covData", "covVar")
    }
  }

  parametervector=c(unlist(parameters), transformedPar, covPars)

  ############################## starting values
  inits <- function() {
    if (model == "betaMPT"){
      ini <- list("theta"=matrix(runif(subjs*S), S, subjs))
    } else {
      # draw appropriate random starting values:
      mu <- xi <- rep(NA, S)
      for(s in 1:S){
        tmp <- ifelse(length(hyperpriors$mu)<=1,hyperpriors$mu[1],hyperpriors$mu[s])
        mu[s] <- eval(parse(text=sub("d","r", sub("(","(1,", tmp,  fixed=TRUE))))

        if(hyperpriors$xi[1] == "dunif(0,10)"){
          xi[s] <- runif(1,.2,1)             # less extreme starting values
        }else{
          tmp <- ifelse(length(hyperpriors$xi)<=1,hyperpriors$xi[1],hyperpriors$xi[s])
          xi[s] <- eval(parse(text=sub("d","r", sub("(","(1,", tmp,  fixed=TRUE))))
        }
      }
      ini <- list("delta.part.raw" = matrix(rnorm(subjs*S, -1,1), S, subjs),
                  "xi"=xi,"mu" = mu)

      # starts with small correlations and scaling parameters close to 1
      if (!anyNA(hyperpriors$V))
        ini$T.prec.part <- as.matrix(rWishart(1,df+30,V)[,,1])

      # check starting values:
      # hist(replicate(1000,cov2cor(solve(rWishart(1,4+1+30,diag(2))[,,1]))[1,2]))
      # hist(replicate(5000,runif(1,.2,1)*sqrt(solve(rWishart(1,4+1+30,diag(2))[,,1])[1,1])))
    }
    ini
  }
  inits.list <- replicate(n.chains, inits(), simplify=FALSE)
  for(i in 1:length(inits.list)){
    inits.list[[i]]$.RNG.name <- c("base::Wichmann-Hill",
                                   "base::Marsaglia-Multicarry",
                                   "base::Super-Duper",
                                   "base::Mersenne-Twister")[1+ (i-1)%% 4]
    inits.list[[i]]$.RNG.seed <-  sample.int(1e4, 1)
  }

  n.samples <- ceiling((n.iter-n.burnin)/n.thin)
  if(n.samples > 30000)
    warning("Note: Your present MCMC settings for n.burnin/n.iter/n.thin\n",
            "      imply that more than 30,000 samples are stored per parameter per chain.\n",
            "      This might result in problems due to an overload of your computers memory (RAM).")

  data.list <-  lapply(data, get, envir=environment())
  names(data.list) <- data
  if (any(c("initlist", "inits", "init") %in% names(list(...)))){
    samples <- run.jags(model = modelfile, monitor=c(parametervector, "deviance"),
                        data=data.list, n.chains=n.chains,
                        burnin=n.burnin, adapt=n.adapt,  sample=n.samples, thin=n.thin,
                        modules=c("dic","glm"), summarise=FALSE, method="parallel", ...)
  } else {
    samples <- run.jags(model = modelfile, monitor=c(parametervector, "deviance"),
                        data=data.list, inits=inits.list, n.chains=n.chains,
                        burnin=n.burnin, adapt=n.adapt,  sample=n.samples, thin=n.thin,
                        modules=c("dic","glm"), summarise=FALSE, method="parallel", ...)
  }

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

# data generation for categorical point-mass distribution
rcat <- function(n, const){
  rep(const, n)
}
