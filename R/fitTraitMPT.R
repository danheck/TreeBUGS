#' Fit a Hierarchical latent-trait MPT Model
#'
#' Fits a latent-trait MPT model (Klauer, 2010) based on a standard MPT model file (.eqn) and individual data table (.csv).
#'
#' @inheritParams betaMPT
#' @param predStructure  Similar to \code{covStructure}, but defines the mapping which variables in \code{covData} are predictors for which MPT parameters (whereas \code{covStructure} determines only sampled correlations). Default: No predictors.
#' @param predType a character vector specifying the type of continuous or discrete predictors in each column of \code{covData}: \code{"c"} = continuous covariate (which are centered to have a mean of zero); \code{"f"} = discrete predictor, fixed effect (default for character/factor variables); \code{"r"} = discrete predictor, random effect.
#' @param corProbit whether to use probit-transformed MPT parameters to compute correlations (the default for trait-MPT)
#' @param mu hyperprior for group means of probit-transformed parameters. Default is a standard normal distribution that implied a uniform distribution on the MPT probability parameters
#' @param xi hyperprior for scaling parameters of the group-level parameter variances. Default is a uniform distribution on the interval [0,100]
#' @param V  S x S matrix used as a hyperprior for the inverse-Wishart hyperprior parameters with as many rows and columns as there are core MPT parameters. Default is a diagonal matrix.
#' @param df degrees of freedom for the inverse-Wishart hyperprior for the individual parameters. Minimum is S+1, where S gives the number of core MPT parameters.
#'
#' @param IVprec hyperprior on the precision (i.e., the inverse of the variance) of the slope parameter for continuous independent variables. Default implies a Cauchy prior (Rouder et. al, 2012). To use a more-informative standard-normal prior, use \code{IVprec = 'dcat(1)'}.
#' @return a list of the class \code{traitMPT} with the objects:
#' \itemize{
#'  \item \code{summary}: MPT tailored summary. Use \code{summary(fittedModel)}
#'  \item \code{mptInfo}: info about MPT model (eqn and data file etc.)
#'  \item \code{mcmc}: the object returned from the MCMC sampler. Note that the object \code{fittedModel$mcmc} is an \link[runjags]{runjags} object, whereas \code{fittedModel$mcmc$mcmc} is a mcmc.list as used by the coda package (\link[coda]{mcmc})
#' }
#' @author Daniel Heck, Denis Arnold, Nina R. Arnold
#' @references
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70-98.
#'
#' Matzke, D., Dolan, C. V., Batchelder, W. H., & Wagenmakers, E.-J. (2015). Bayesian estimation of multinomial processing tree models with heterogeneity in participants and items. Psychometrika, 80, 205-235.
#'
#' Rouder, J. N., Morey, R. D., Speckman, P. L., & Province, J. M. (2012). Default Bayes factors for ANOVA designs. Journal of Mathematical Psychology, 56, 356-374.


#' @export

traitMPT <- function(eqnfile,  # statistical model stuff
                    data,
                    restrictions,
                    covData,
                    covStructure,   # correlation
                    predStructure,  # predictor structure
                    predType,    # c("c," f", "r")
                    transformedParameters,
                    T1group,
                    corProbit=TRUE,

                    # hyperpriors:
                    mu = "dnorm(0,1)",
                    xi = "dunif(0,100)",
                    V,
                    df,
                    IVprec = "dchisq(1)",  # change to "dcat(1)" to set beta ~ dnorm(0,1)

                    # MCMC stuff:
                    n.iter=20000,
                    n.adapt = 2000,
                    n.burnin=2000,
                    n.thin=5,
                    n.chains=3,
                    dic =FALSE,

                    # File Handling stuff:
                    modelfilename,
                    parEstFile,
                    autojags=NULL,
                    ...){
  if(missing(restrictions)) restrictions <- NULL
  if(missing(covData)) covData <- NULL
  if(missing(covStructure)) covStructure <- NULL
  if(missing(predStructure)) predStructure <- NULL
  if(missing(predType)) predType <- NULL
  if(missing(transformedParameters)) transformedParameters <- NULL
  if(missing(T1group)) T1group <- NULL

  checkParEstFile(parEstFile)
  modelfilename <- checkModelfilename(modelfilename)
  data <- readData(data)

  if(n.iter <= n.burnin)
    stop("n.iter must be larger than n.burnin")

  # MPT structure
  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  data <- readSubjectData(data, unique(mergedTree$Category))

  tHoutput <- thetaHandling(mergedTree,restrictions)
  SubPar <- tHoutput$SubPar
  mergedTree <- tHoutput$mergedTree

  thetaNames <- SubPar[,1:2]
  thetaUnique <- thetaNames[rownames(unique(thetaNames[2])),]$Parameter
  S <- max(SubPar$theta)
  isIdentifiable(S, mergedTree)

  # transformed parameters
  transformedPar <- getTransformed(model = "traitMPT",
                                   thetaNames = thetaNames,
                                   transformedParameters = transformedParameters)
  N <- nrow(data)


  # covariate: reading + checking
  covData <- covDataRead(covData, N)
  predType <- predTypeDefault(covData, predType=predType)
  groupT1 <- getGroupT1(covData, predType, T1group=T1group)
  covData <- covDataCenter(covData, predType=predType)


  # PREDICTORS: assign covariates to parameters and handle factor levels
  predTmp1 <- covHandling(covData, predStructure, N, thetaNames, predType=predType,
                          defaultExclude="ALL_COVARIATES", T1group=T1group)
  predFactorLevels <- predTmp1$predFactorLevels
  predTable <- predTmp1$covTable
  covData <- predTmp1$covData
  # get model string: phi(....)
  predTmp2 <- covStringPredictor(predTable, S=S,
                                 predFactorLevels=predFactorLevels, IVprec=IVprec)
  predString <- predTmp2$modelString
  predPars <- predTmp2$covPars

  # CORRELATIONS
  covTmp1 <- covHandling(covData, covStructure, N, thetaNames, predType=predType,
                         defaultExclude=unique(predTable$Covariate))
  covTable <- covTmp1$covTable
  covTmp2 <- covStringCorrelation(covTable, corProbit=corProbit)
  corString <- covTmp2$modelString
  corPars <- covTmp2$covPar

  # T1 per group split
  covData <- covData[,sapply(covData, class) %in% c("numeric", "integer"), drop=FALSE]

  if( any(covTable$Covariate %in% predTable$Covariate))
    warning("Some variables in covData appear both as predictors and as covariates:\n    ",
            paste(covTable$Covariate[covTable$Covariate %in% predTable$Covariate], collapse=", "))


  # hyperpriors
  if(missing(V) || is.null(V))
    V <- diag(S)
  if(missing(df) || is.null(df))
    df <- S+1

  makeModelFile(model = "traitMPT",
                filename = modelfilename,
                mergedTree = mergedTree ,
                S = S,
                hyperprior = list(mu = mu, xi = xi),
                predString = predString,
                corString = corString,
                parString=transformedPar$modelstring,
                groupMatT1=groupT1$groupMatT1)

  time0 <- Sys.time()
  cat("MCMC sampling started at ", format(time0), "\n")
  runjags <- callingSampler(model = "traitMPT",
                         mergedTree = mergedTree,
                         data = data,
                         modelfile = modelfilename,
                         S = max(SubPar$theta),
                         transformedPar = transformedPar$transformedParameters,
                         covPars=c(corPars, predPars),
                         covData=covData,
                         X_list=predTmp2$X_list,
                         groupT1=groupT1,
                         hyperpriors = list(V=V, df=df),
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

  # Beta MPT: rename parameters and get specific summaries
#   summary <- summarizeMPT(model = "traitMPT",
#                           mcmc = mcmc,
#                           thetaNames = thetaNames,
#                          # covIncluded = !is.null(covData),
#                           predFactorLevels=predFactorLevels,
#                           transformedParameters = transformedPar$transformedParameters,
#                           NgroupT1 = groupT1$NgroupT1)

  mptInfo <- list(model="traitMPT",
                  thetaNames = thetaNames,
                  thetaUnique = thetaUnique,
                  MPT=mergedTree,
                  eqnfile=eqnfile,
                  data=data,
                  restrictions=restrictions,
                  covData=covData,
                  covTable=covTable,
                  corProbit=corProbit,
                  predTable=predTable,
                  predFactorLevels=predFactorLevels,
                  transformedParameters=transformedPar$transformedParameters,
                  T1group=groupT1)

  # own summary (more stable than runjags)
  mcmc.summ <- summarizeMCMC(runjags$mcmc)
  summary <- summarizeMPT(mcmc = runjags$mcmc, mptInfo = mptInfo, summ=mcmc.summ)
  if(dic){
    summary$dic <-   extract(runjags, "dic", ...)
  }else{summary$dic <- NULL}


  # class structure for TreeBUGS
  # mcmc$BUGSoutput <- renameBUGSoutput(mcmc$BUGSoutput, thetaUnique, "traitMPT")
  fittedModel <- list(summary=summary,
                      mptInfo=mptInfo,
                      mcmc.summ = mcmc.summ,
                      runjags = runjags,
                      call=match.call(),
                      time=time1-time0)
  class(fittedModel) <- "traitMPT"

  # write results to file
  writeSummary(fittedModel, parEstFile)

  gc(verbose=FALSE)

  return(fittedModel)
}
