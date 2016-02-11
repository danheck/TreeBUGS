#' Fit a Hierarchical Beta-MPT Model
#'
#' Fits a Beta-MPT model (Smith & Batchelder, 2010) based on a standard MPT model file (.eqn) and individual data table (.csv).
#'
#' @param eqnfile The (full path to the) file that specifies the MPT model (standard .eqn syntax)
#' @param data The (full path to the) csv file with the data (semicolon separated; category labels in first row). Alternatively: a data frame or matrix (rows=individuals, columns = individual category frequencies, category labels as column names)
#' @param restrictions  Optional: Either the (full path to the) file that specifies which parameters should be constants and which should be equal; alternatively: a list of restrictions, e.g., \code{list("D1=D2","g=0.5")}
#' @param covData The path to the csv data file for the individual values on the covariates (semicolon-separated; rows=individuals in the same order as \code{data}, covariate labels in first row). Alternatively: a data frame or matrix (rows=individuals, columns = individual values on covariates, covariate labels as column names)
#' @param covStructure  Optional: Either the (full path to the) file that specifies the assigment of MPT parameters to covariates (that is, each row assigns one or more MPT parameters to one or more covariates, separated by a semicolon, e.g., \code{Do g; age extraversion}). Can also be provided as a list, e.g., \code{list("Do Dn ; age", "g ; extraversion"}). Default: All combinations included (not recommended, could be unstable).
#' @param transformedParameters list with parameter transformations that should be computed based on the posterior samples (e.g., for testing parameter differences: \code{list("diffD=Do-Dn")})
#' @param modelfilename Name that the modelfile that is made by the function to work with JAGS should get.
#'        Default is to write this information to the tempdir as required by CRAN standards.
#' @param alpha Hyperprior of for the alpha and beta parameters (default: uniform prior on the interval [1,5000]).
#' @param beta Second hyperprior, see \code{alpha}
#' @param parEstFile Name of the file to with the estimates should be stored (e.g., "parEstFile.txt")
#' @param n.iter Number of iterations (including burnin samples).
#' @param n.burnin Burnin period.
#' @param n.thin Thinning rate.
#' @param n.chains number of MCMC chains
#' @param autojags whether to run JAGS until the MCMC chains converge (see \link{autojags}).  Use \code{n.update=3} as an additional argument to control how often JAGS is rerun. Can take a lot of time for large models.
#' @param ... Arguments to be passed to the sampling function (default: \code{\link{jags.parallel}}.
#'
#' @details Note that, in the Beta-MPT model, correlations of individual MPT parameters with covariates are sampled. Hence, the covariates do not affect the estimation of the actual Beta-MPT parameters. Therefore, the correlation of covariates with the individual MPT parameters can equivalently be performed after fitting the model using the sampled posterior parameter values stored in \code{betaMPT$mcmc}
#'
#' @return a list of the class \code{betaMPT} with the objects:
#' \itemize{
#'  \item \code{summary}: MPT tailored summary. Use \code{summary(fittedModel)}
#'  \item \code{mptInfo}: info about MPT model (eqn and data file etc.)
#'  \item \code{mcmc}: the object returned from the MCMC sampler. Standard: An \code{\link{jags.parallel}} object. Note that the sample can be transformed into an \code{mcmc.list} for analysis with the \code{coda} package by \code{as.mcmc.list(fittedModel$mcmc$BUGSoutput)}
#' }
#' @author Nina R. Arnold, Denis Arnold, Daniel Heck
#' @references Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial processing tree models for addressing individual differences. Journal of Mathematical Psychology, 54, 167-183.
#' @export

betaMPT <- function(eqnfile,  # statistical model stuff
                    data,
                    restrictions,
                    covData,
                    covStructure,
                    transformedParameters,
                    alpha = "dunif(.01,5000)",
                    beta = "dunif(.01,5000)",

                    # MCMC stuff:
                    n.iter=50000,
                    n.burnin=5000,
                    n.thin=5,
                    n.chains=3,

                    # File Handling stuff:
                    modelfilename,
                    parEstFile,
                    autojags = FALSE,
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


  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  data <- readSubjectData(data,unique(mergedTree$Category))

  tHoutput <- thetaHandling(mergedTree,restrictions)
  SubPar <- tHoutput$SubPar
  mergedTree <- tHoutput$mergedTree

  thetaNames <- tHoutput[[1]][,1:2]
  S <- max(SubPar$theta)
  isIdentifiable(S, mergedTree)

  # transformed parameters
  transformedPar <- getTransformed(model = "betaMPT",
                                   thetaNames = thetaNames,
                                   transformedParameters = transformedParameters)

  N <- nrow(data)

  # covariate: get neat table and appropriate JAGS string
  covTmp1 <- covHandling(covData, covStructure, N, thetaNames)
  covData <- covTmp1$covData
  covTable <- covTmp1$covTable
  if( any(covTable$covType != "c"))
    stop("To compute correlations in betaMPT, only continuous covariates are allowed.")

  covTmp2 <- covStringBeta(covTable)
  covString <- covTmp2$modelString
  covPars <- covTmp2$covPar


  makeModelFile(model = "betaMPT",
                filename = modelfilename,
                mergedTree = mergedTree ,
                S = max(SubPar$theta),
                hyperprior = list(alpha=alpha, beta = beta),
                covString = covString,
                parString = transformedPar$modelstring)

  time0 <- Sys.time()
  cat("MCMC sampling started at ", format(time0), "\n")
  mcmc <- callingSampler(model = "betaMPT",
                         mergedTree = mergedTree,
                         data = data,
                         modelfile = modelfilename,
                         S = max(SubPar$theta),
                         transformedPar = transformedPar$transformedParameters,
                         covPars = covPars,
                         covData = covData,
                         # parameters,
                         n.iter = n.iter,
                         n.burnin = n.burnin,
                         n.thin = n.thin,
                         n.chains = n.chains,
                         autojags = autojags,
                         ...)
  time1 <- Sys.time()
  cat("MCMC sampling finished at", format(time1), "\n  ")
  print(time1-time0)

  # Beta MPT: rename parameters and get specific summaries
  summary <- summarizeMPT(model = "betaMPT",
                          mcmc = mcmc,
                          thetaNames = thetaNames,
                          covIncluded = !is.null(covData),
                          transformedParameters = transformedPar$transformedParameters)

  mptInfo <- list(thetaNames = thetaNames,
                  MPT=mergedTree,
                  eqnfile = eqnfile,
                  data = data,
                  restrictions = restrictions,
                  covData=covData,
                  covTable=covTable,
                  transformedParameters=transformedPar$transformedParameters)

  # class structure for TreeBUGS
  fittedModel <- list(summary = summary,
                      mptInfo = mptInfo,
                      mcmc = mcmc,
                      call = match.call(),
                      time = time1-time0)

  class(fittedModel) <- "betaMPT"

  # write results to file
  writeSummary(fittedModel, parEstFile)


  return(fittedModel)
}


