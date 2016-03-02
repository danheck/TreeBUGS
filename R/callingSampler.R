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
                           groupT1=NULL,    # list with groupMatT1 und NgroupT1  for splitted T1 statistic
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
    rowsums <- rowSums(responses[(index+1):(index+NresponsesTree[i])])
    if (any(is.na(rowsums)) || any(rowsums == 0))
      warning("One or more participants do not have responses for tree",
              treeNames[i], ". As a solution, the critical participants might be excluded.")

    assign(name.response,  matrix(as.vector(t(responses[(index+1):(index+NresponsesTree[i])])),
                                  ncol=NresponsesTree[i],nrow=subjs,byrow=TRUE))

    # mean frequencies for T1 statistic
    assign(name.mean, colMeans(get(name.response)))

    index=index+NresponsesTree[i]

    data <- c(data, name.mean, name.response, name.items)

    # add G x numCat matrix with mean frequencies (one line per group)
    if(!is.null(groupT1)){
      name.group.mean <- paste0("group.resp.",treeNames[i],".mean")

      mean.per.group <- c()
      for(g in 1:nrow(groupT1$groupMatT1)){
        idx <- groupT1$groupMatT1[1:groupT1$NgroupT1[g]]
        mean.per.group <- rbind(mean.per.group,
                                colMeans(get(name.response)[idx,,drop=FALSE]))

      }
      assign(name.group.mean, mean.per.group)
      data <- c(data, name.group.mean)
    }
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
  if(!is.null(covData) & !any(dim(covData) == 0)){
    covData <- as.matrix(covData)
    data <- c(data, "covData")
    if(any(grepl("cor_", covPars))){
      covSD <- apply(covData, 2, sd, na.rm=TRUE)
      data <- c(data, "covSD")
    }
  }

  # call Sampler
    parametervector=c(unlist(parameters),
                      "T1.obs","T1.pred","p.T1","p.T1ind",
                      transformedPar, covPars)

    if(!is.null(groupT1)){
      groupMatT1 <- groupT1$groupMatT1
      NgroupT1 <- groupT1$NgroupT1
      data <- c(data, "groupMatT1", "NgroupT1")
      parametervector <- c(parametervector, "T1.group.obs", "T1.group.pred", "p.T1.group",
                           paste0("group.resp.",treeNames, ".pred.mean"))
    }

    # random initial values: required (otherwise T1 statistic results in errors: no variance!)
    if(model == "betaMPT"){
      inits <- function() list("theta"=matrix(runif(subjs*S), S, subjs)
                               )
    }else{
      inits <- function() list("delta.part.raw" = matrix(rnorm(subjs*S, 0,1), S, subjs),
                               "xi"=runif(S,2/3,3/2), "mu" = rnorm(S, 0, .5),
                               "T.prec.part" = rWishart(1,2*df,V)[,,1]
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

    data.list <-  lapply(data, get, envir=environment())
    names(data.list) <- data
    samples <- run.jags(model = modelfile,
                        monitor=parametervector,
                        data=data.list,
                        n.chains=n.chains,
                        inits=inits.list,
                        burnin=n.burnin,
                        adapt=n.adapt,
                        sample=ceiling((n.iter-n.burnin)/n.thin),
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
      samples <- do.call(autoextend.jags, c(list(runjags.object = samples,
                                                 summarise=FALSE),
                                            autojags))  # additional user arguments
    }

  return(samples)   # own summary
}
