#' Fit a Hierarchical latent-trait MPT Model
#'
#' Fits a latent-trait MPT model (Klauer, 2010) based on a standard MPT model file (.eqn) and individual data table (.csv).
#'
#' @inheritParams betaMPT
#' @param covType a character vector specifying the type of covariates in each column of \code{covData}: \code{"c"} = continuous covariate (which is the default for continuous variables); \code{"f"} = fixed effect (default for character/factor variables); \code{"r"} = random effect.
#' @param mu hyperprior for group means of probit-transformed parameters. Default is a standard normal distribution.
#' @param xi hyperprior for scaling parameters of the group-level parameter variances. Default is a uniform distribution on the interval [0,100]
#' @param V  S x S matrix used as a hyperprior for the inverse-Wishart hyperprior parameters with as many rows and columns as there are core MPT parameters. Default is a diagonal matrix.
#' @param df degrees of freedom for the inverse-Wishart hyperprior for the individual parameters. Minimum is S+1, where S gives the number of core MPT parameters.
#'
#' @return a list of the class \code{traitMPT} with the objects:
#' \itemize{
#'  \item \code{summary}: MPT tailored summary. Use \code{summary(fittedModel)}
#'  \item \code{mptInfo}: info about MPT model (eqn and data file etc.)
#'  \item \code{mcmc}: the object returned from the MCMC sampler. Standard: An \code{\link{jags.parallel}} object. Note that the sample can be transformed into an \code{mcmc.list} for analysis with the \code{coda} package by \code{as.mcmc.list(fittedModel$mcmc$BUGSoutput)}
#'  \item \code{sampler}: the type of sampler used (standard: \code{"JAGS"})
#' }
#' @author Daniel Heck, Denis Arnold, Nina R. Arnold
#' @references
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70-98.
#'
#' Matzke, D., Dolan, C. V., Batchelder, W. H., & Wagenmakers, E.-J. (2015). Bayesian estimation of multinomial processing tree models with heterogeneity in participants and items. Psychometrika, 80, 205-235.

#' @export

traitMPT <- function(eqnfile,  # statistical model stuff
                    data,
                    restrictions,
                    covData,
                    covStructure,
                    covType,  # c("c," f", "r")
                    transformedParameters,

                    # hyperpriors:
                    mu = "dnorm(0,1)",
                    xi = "dunif(0,100)",
                    V,
                    df,

                    # MCMC stuff:
                    n.iter=50000,
                    n.burnin=5000,
                    n.thin=5,
                    n.chains=3,

                    # File Handling stuff:
                    modelfilename,
                    parEstFile,
                    sampler="JAGS",
                    autojags=FALSE,
                    ...){
  if(missing(restrictions)) restrictions <- NULL
  if(missing(covData)) covData <- NULL
  if(missing(covStructure)) covStructure <- NULL
  if(missing(transformedParameters)) transformedParameters <- NULL

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
  S <- max(SubPar$theta)
  isIdentifiable(S, mergedTree)

  # transformed parameters
  transformedPar <- getTransformed(model = "traitMPT",
                                   thetaNames = thetaNames,
                                   transformedParameters = transformedParameters)
  N <- nrow(data)


  # assign covariates to parameters and handle factor levels
  covTmp1 <- covHandling(covData, covStructure, N, thetaNames, covType=covType)
  covData <- covTmp1$covData
  covFactorLevels <- covTmp1$covFactorLevels
  covTable <- covTmp1$covTable
  covType <- covTmp1$covType

  if(!is.null(covData)){
    for(i in 1:ncol(covData)){
      if(covType[i] == "c" &  colnames(covData)[i] %in% covTable$Covariate){
        scaled <- scale(covData[,i], scale = FALSE)
        if(any(scaled != covData[,i])){
          warning("Predictor covariate ", colnames(covData)[i],
                  "is automatically centererd to a MEAN OF ZERO \n  (recommended for latent-trait MPT)!")
          covData[,i] <- scaled
        }
      }
    }
  }

  # get model string: phi(....)
  covTmp2 <- covStringTrait(covTable, S=S, covFactorLevels=covFactorLevels)
  covString <- covTmp2$modelString
  covPars <- covTmp2$covPars

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
                covString = covString,
                sampler=sampler,
                parString=transformedPar$modelstring)

  time0 <- Sys.time()
  cat("MCMC sampling started at ", format(time0), "\n")
  mcmc <- callingSampler(model = "traitMPT",
                         mergedTree = mergedTree,
                         data = data,
                         modelfile = modelfilename,
                         S = max(SubPar$theta),
                         transformedPar = transformedPar$transformedParameters,
                         covPars=covPars,
                         covData=covData,
                         X_list=covTmp2$X_list,
                         hyperpriors = list(V=V, df=df),
                         n.iter = n.iter,
                         n.burnin = n.burnin,
                         n.thin = n.thin,
                         n.chains = n.chains,
                         sampler = sampler,
                         autojags = autojags,
                         ...)
  time1 <- Sys.time()
  cat("MCMC sampling finished at", format(time1), "\n  ")
  print(time1-time0)

  # Beta MPT: rename parameters and get specific summaries
  summary <- summarizeMPT(model = "traitMPT",
                          mcmc = mcmc,
                          thetaNames = thetaNames,
                          sampler = sampler,
                          covIncluded = !is.null(covData),
                          covFactorLevels=covFactorLevels,
                          transformedParameters = transformedPar$transformedParameters)

  mptInfo <- list(thetaNames = thetaNames,
                  MPT=mergedTree,
                  eqnfile=eqnfile,
                  data=data,
                  restrictions=restrictions,
                  covData=covData,
                  covTable=covTable,
                  covFactorLevels=covFactorLevels,
                  transformedParameters=transformedPar$transformedParameters)

  # class structure for TreeBUGS
  fittedModel <- list(summary=summary,
                      mptInfo=mptInfo,
                      mcmc=mcmc,
                      sampler=sampler,
                      call=match.call(),
                      time=time1-time0)
  class(fittedModel) <- "traitMPT"

  # write results to file
  writeSummary(fittedModel, parEstFile)


  return(fittedModel)
}
