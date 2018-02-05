
# general model fitting function
# type: "traitMPT" or "betaMPT"
fitModel <- function(type,
                     eqnfile,  # statistical model stuff
                     data,
                     restrictions,
                     covData,
                     predStructure,  # predictor structure
                     predType,    # c("c," f", "r")
                     transformedParameters,
                     corProbit=TRUE,

                     hyperprior,

                     # MCMC stuff:
                     n.iter=20000,
                     n.adapt = 2000,
                     n.burnin=2000,
                     n.thin=5,
                     n.chains=3,
                     dic =FALSE,
                     ppp = 0,

                     # File Handling stuff:
                     modelfilename,
                     parEstFile,
                     posteriorFile,
                     autojags=NULL,
                     call = NULL,
                     ...){

  if(missing(restrictions)) restrictions <- NULL
  if(missing(covData)) covData <- NULL
  if(missing(transformedParameters)) transformedParameters <- NULL
  if(missing(predStructure)) predStructure <- NULL
  if(missing(predType)) predType <- NULL

  checkParEstFile(parEstFile)
  modelfilename <- checkModelfilename(modelfilename)
  data <- readData(data)

  if(n.iter <= n.burnin)
    stop("n.iter must be larger than n.burnin")

  # MPT structure for JAGS
  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  data <- readSubjectData(data, unique(mergedTree$Category))

  tHoutput <- thetaHandling(mergedTree,restrictions)
  SubPar <- tHoutput$SubPar
  mergedTree <- tHoutput$mergedTree
  fixedPar <- tHoutput$fixedPar

  thetaNames <- SubPar[,1:2]
  thetaUnique <- thetaNames[rownames(unique(thetaNames[2])),]$Parameter
  S <- max(SubPar$theta)
  isIdentifiable(S, mergedTree)

  # transformed parameters
  transformedPar <- getTransformed(thetaNames = thetaNames,
                                   transformedParameters = transformedParameters)
  N <- nrow(data)

  # (A) traitMPT: adjust predictor structure
  # (B) both: remove discrete variables (+ predictors) and correlate covData<->theta

  # covariate: reading + checking
  covData <- covDataRead(covData, N)
  predType <- predTypeDefault(covData, predType=predType)

  if(type == "traitMPT"){

    # inverse Wishart prior
    if (is.null(hyperprior$V)){
      hyperprior$V <- diag(S)
      if(is.null(hyperprior$df))
        hyperprior$df <- S+1

    # independent chi-square
    } else if (is.na(hyperprior$V)){
      if(is.null(hyperprior$df))
        hyperprior$df <- 1
      # hyperprior$V <- NULL
    }
    hyperprior$mu <- check.hyperprior(hyperprior$mu, thetaUnique, label="mu")

    ##################### TRAIT MPT
    covData <- covDataCenter(covData, predType=predType)

    # PREDICTORS: assign covariates to parameters and handle factor levels
    predTmp1 <- covHandling(covData, predStructure, N, thetaNames, predType=predType,
                            defaultExclude="ALL_COVARIATES")
    predFactorLevels <- predTmp1$predFactorLevels
    predTable <- predTmp1$covTable
    covDataNumeric <- predTmp1$covData
    # get model string: phi(....)
    predTmp2 <- covStringPredictor(predTable, S=S,
                                   predFactorLevels=predFactorLevels,
                                   IVprec=hyperprior$IVprec)
    predString <- predTmp2$modelString
    covPars <- predTmp2$covPars
    X_list <- predTmp2$X_list



  }else{
    ##################### BETA MPT
    hyperprior$alpha <- check.hyperprior(hyperprior$alpha, thetaUnique, label="alpha")
    hyperprior$beta <- check.hyperprior(hyperprior$beta, thetaUnique, label="beta")

    predString <- NULL
    X_list <- covDataNumeric <- covPars <- NULL
    predTable <- predFactorLevels <- NULL
  }

  makeModelFile(model = type, filename = modelfilename, mergedTree = mergedTree ,
                S = S, hyperprior = hyperprior, predString = predString,
                parString=transformedPar$modelstring, fixedPar=fixedPar)

  time0 <- Sys.time()
  cat("MCMC sampling started at ", format(time0), "\n")
  runjags <- callingSampler(model = type,
                            mergedTree = mergedTree,
                            data = data,
                            modelfile = modelfilename,
                            S = max(SubPar$theta),
                            fixedPar=fixedPar,
                            transformedPar = transformedPar$transformedParameters,
                            covPars=covPars,
                            covData=covDataNumeric,   # must be a numeric matrix!
                            X_list=X_list,
                            hyperpriors = hyperprior,
                            n.iter = n.iter,
                            n.adapt = n.adapt,
                            n.burnin = n.burnin,
                            n.thin = n.thin,
                            n.chains = n.chains,
                            autojags = autojags,
                            ...)
  time1 <- Sys.time()
  cat("MCMC sampling finished at", format(time1), "\n  ")
  print(time1-time0)

  # store details about model:
  mptInfo <- list(model=type,
                  thetaNames = thetaNames,
                  thetaUnique = thetaUnique,
                  thetaFixed = fixedPar$Parameter[!duplicated(fixedPar$theta)],
                  MPT=mergedTree,
                  eqnfile=eqnfile,
                  data=data,
                  restrictions=restrictions,
                  covData=covData,
                  corProbit=corProbit,
                  predTable=predTable,
                  predFactorLevels=predFactorLevels,
                  predType=predType,
                  transformedParameters=transformedPar$transformedParameters,
                  hyperprior=hyperprior)

  # correlation of posterior samples:
  if(!is.null(covData) | type == "betaMPT"){
    if(!is.null(predTable)){
      isPred <- (1:ncol(covData)) %in% predTable$covIdx
    }else{
      isPred <- rep(FALSE, length(predType))
    }

    sel <- predType == "c" & !isPred
    if(any(sel) | type == "betaMPT"){
      cdat <- covData[,sel,drop = FALSE]
      runjags$mcmc <- as.mcmc.list(
        lapply(runjags$mcmc, corSamples,
               covData=cdat, thetaUnique=thetaUnique,
               rho=ifelse(type == "betaMPT", TRUE, FALSE),
               corProbit = corProbit))
    }
  }


  # own summary (more stable than runjags)
  mcmc.summ <- summarizeMCMC(runjags$mcmc)
  summary <- summarizeMPT(mcmc = runjags$mcmc, mptInfo = mptInfo,
                          summ=mcmc.summ)
  summary$call <- call

  if(dic){
    summary$dic <-   extract(runjags, "dic", ...)
  }else{
    summary$dic <- NULL
  }


  # class structure for TreeBUGS
  fittedModel <- list(summary=summary,
                      mptInfo=mptInfo,
                      mcmc.summ = mcmc.summ,
                      runjags = runjags,
                      postpred=NULL,
                      call=call,
                      time=time1-time0)
  class(fittedModel) <- type


  fittedModel <- addPPP(fittedModel, M=ppp)

  # write results to file
  writeSummary(fittedModel, parEstFile)
  if(!missing(posteriorFile) && !is.null(posteriorFile))
    try(save(fittedModel, file=posteriorFile))
  gc(verbose=FALSE)

  fittedModel
}
