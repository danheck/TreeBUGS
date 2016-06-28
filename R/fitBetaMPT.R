#' Fit a Hierarchical Beta-MPT Model
#'
#' Fits a Beta-MPT model (Smith & Batchelder, 2010) based on a standard MPT model file (.eqn) and individual data table (.csv).
#'
#' @param eqnfile The (full path to the) file that specifies the MPT model (standard .eqn syntax). Note that category labels must start with a letter (different to multiTree) and match the column names of \code{data}
#' @param data The (full path to the) csv file with the data (comma separated; category labels in first row). Alternatively: a data frame or matrix (rows=individuals, columns = individual category frequencies, category labels as column names)
#' @param restrictions  Specifies which parameters should be (a) constant (e.g., \code{"a=b=.5"}) or (b) constrained to be identical (e.g., \code{"Do=Dn"}) or (c) treated as fixed effects (i.e., identical for all participants; \code{"a=b=FE"}). Either given as the path to a text file with restrictions per row or as a list of restrictions, e.g., \code{list("D1=D2","g=0.5")}
#' @param covData The path to the csv data file for the individual values on the covariates (comma-separated; rows=individuals in the same order as \code{data}, covariate labels in first row). Alternatively: a data frame or matrix (rows=individuals, columns = individual values on covariates, covariate labels as column names)
# @param covStructure  Optional: Which correlations to compute? Either the (full path to the) file that specifies the assigment of MPT parameters to covariates (that is, each row assigns one or more MPT parameters to one or more covariates, separated by a semicolon, e.g., \code{Do g; age extraversion}). Can also be provided as a list, e.g., \code{list("Do Dn ; age", "g ; extraversion"}). Default: All continuous variables in covData included.
#' @param transformedParameters list with parameter transformations that should be computed based on the posterior samples (e.g., for testing parameter differences: \code{list("diffD=Do-Dn")})
#' @param modelfilename Name that the modelfile that is made by the function to work with JAGS should get.
#'        Default is to write this information to the tempdir as required by CRAN standards.
#' @param corProbit whether to use probit-transformed MPT parameters to compute correlations (probit-values of \code{+Inf} are truncated to \code{max(5,max(probit))}; similarly for \code{-Inf}). Default for beta-MPT: MPT parameters are used on the probability scale [0,1].
#' @param alpha Hyperprior for the alpha and beta parameters in JAGS syntax (default: uniform prior on the interval [1,5000] for all parameters). A vector can be used to specify separate hyperpriors for each MPT parameter (to check the order of parameters, use \code{\link{readEQN}} with \code{paramOrder = TRUE}).
#' @param beta Second hyperprior, see \code{alpha}
#' @param parEstFile Name of the file to with the estimates should be stored (e.g., "parEstFile.txt")
#' @param n.iter Number of iterations per chain (including burnin samples). See \code{\link[runjags]{run.jags}} for details.
#' @param n.adapt number of adaption samples to adjust MCMC sampler in JAGS. The sampler will be more efficient if it is tuned well.
#' @param n.burnin Number of samples for burnin (samples will not be stored and removed from n.iter)
#' @param n.thin Thinning rate.
#' @param n.chains number of MCMC chains (sampled in parallel).
#' @param dic whether to compute DIC using \code{\link[runjags]{extract}}, which requires additional sampling. Can also be computed and added after fitting the model by \code{fittedModel$dic <- extract(fittedModel$runjags, "dic")}
#' @param ppp number of samples to compute  posterior predictive p-value (see \code{\link{posteriorPredictive}})
#' @param autojags if provided (as an empty list or with arguments passed to \link[runjags]{autoextend.jags}), JAGS runs repeatedly until the MCMC chains converges . E.g., use \code{list(max.time="30m")} to restrict sampling to 30 minutes (similarly for hours, days, and weeks)
#' @param ... Arguments to be passed to the JAGS sampling function (i.e., to \code{\link[runjags]{run.jags}}.
#'
#' @details Note that, in the Beta-MPT model, correlations of individual MPT parameters with covariates are sampled. Hence, the covariates do not affect the estimation of the actual Beta-MPT parameters. Therefore, the correlation of covariates with the individual MPT parameters can equivalently be performed after fitting the model using the sampled posterior parameter values stored in \code{betaMPT$mcmc}
#'
#' @return a list of the class \code{betaMPT} with the objects:
#' \itemize{
#'  \item \code{summary}: MPT tailored summary. Use \code{summary(fittedModel)}
#'  \item \code{mptInfo}: info about MPT model (eqn and data file etc.)
#'  \item \code{runjags}: the object returned from the MCMC sampler. Note that the object \code{fittedModel$runjags} is an \link[runjags]{runjags} object, whereas \code{fittedModel$runjags$mcmc} is a \code{mcmc.list} as used by the coda package (\link[coda]{mcmc})
#' }
#' @author Daniel Heck, Nina R. Arnold, Denis Arnold,
#' @references Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial processing tree models for addressing individual differences. Journal of Mathematical Psychology, 54, 167-183.
#' @export

betaMPT <- function(eqnfile, data, restrictions,
                    covData, #covStructure,
                    transformedParameters,
                    corProbit=FALSE,
                    alpha = "dunif(.01,5000)",
                    beta = "dunif(.01,5000)",

                    # MCMC stuff:
                    n.iter=20000, n.adapt=2000,
                    n.burnin=2000, n.thin=5,
                    n.chains=3, dic =FALSE,
                    ppp = 0,

                    # File Handling stuff:
                    modelfilename, parEstFile,
                    autojags = NULL,   ...){

  hyperprior <- list(alpha=alpha, beta=beta)

  fittedModel <- fitModel(type="betaMPT", eqnfile=eqnfile,
                          data=data,restrictions=restrictions,
                          covData=covData,#covStructure=covStructure,
                          transformedParameters=transformedParameters,
                          corProbit=corProbit, hyperprior=hyperprior,
                          n.iter=n.iter, n.adapt = n.adapt,
                          n.burnin=n.burnin, n.thin=n.thin,
                          n.chains=n.chains, dic =dic,  ppp = ppp,
                          modelfilename=modelfilename,
                          parEstFile=parEstFile,
                          autojags=autojags,
                          ...)
  fittedModel$call <- match.call()

  return(fittedModel)
}


######################## OLD DUPLICATE CODE (=> see fitModel)



# if(missing(restrictions)) restrictions <- NULL
# if(missing(covData)) covData <- NULL
# if(missing(covStructure)) covStructure <- NULL
# if(missing(transformedParameters)) transformedParameters <- NULL
#
# checkParEstFile(parEstFile)
# modelfilename <- checkModelfilename(modelfilename)
# data <- readData(data)
#
# if(n.iter <= n.burnin)
#   stop("n.iter must be larger than n.burnin")
#
#
# Tree <- readEQN(eqnfile)
# mergedTree <- mergeBranches(Tree)
# data <- readSubjectData(data,unique(mergedTree$Category))
#
# tHoutput <- thetaHandling(mergedTree,restrictions)
# SubPar <- tHoutput$SubPar
# mergedTree <- tHoutput$mergedTree
# fixedPar <- tHoutput$fixedPar
#
# thetaNames <- tHoutput[[1]][,1:2]
# thetaUnique <- thetaNames[rownames(unique(thetaNames[2])),]$Parameter
# S <- max(SubPar$theta)
# isIdentifiable(S, mergedTree)
#
# # transformed parameters
# transformedPar <- getTransformed(model = "betaMPT",
#                                  thetaNames = thetaNames,
#                                  transformedParameters = transformedParameters)
#
# N <- nrow(data)
#
# # covariate: reading + checking
# covData <- covDataRead(covData, N, binaryToNumeric=TRUE)
# predType <- predTypeDefault(covData, predType=NULL)
# # get neat table and appropriate JAGS string
# covTmp1 <- covHandling(covData, covStructure, N, thetaNames, predType=predType, defaultExclude=NULL)
#                        # onlyContinuous=TRUE)
# covTable <- covTmp1$covTable
# if( any(covTable$predType != "c"))
#   stop("To compute correlations in betaMPT, only continuous covariates are allowed.")
#
# # only allow correlations
# covTmp2 <- covStringCorrelation(covTable, corProbit=corProbit)
# corString <- covTmp2$modelString
# covPars <- covTmp2$covPar
# covData <- covData[,sapply(covData, class) %in% c("numeric", "integer"), drop=FALSE]
#
#
#
# makeModelFile(model = "betaMPT",
#               filename = modelfilename,
#               mergedTree = mergedTree ,
#               S = max(SubPar$theta),
#               hyperprior = list(alpha=alpha, beta = beta),
#               corString = corString,
#               parString = transformedPar$modelstring,
#               fixedPar=fixedPar)
#
# time0 <- Sys.time()
# cat("MCMC sampling started at ", format(time0), "\n")
# runjags <- callingSampler(model = "betaMPT",
#                        mergedTree = mergedTree,
#                        data = data,
#                        modelfile = modelfilename,
#                        S = max(SubPar$theta),
#                        fixedPar=fixedPar,
#                        transformedPar = transformedPar$transformedParameters,
#                        covPars = covPars,
#                        covData = covData,
#
#                        # parameters,
#                        n.iter = n.iter,
#                        n.adapt = n.adapt,
#                        n.burnin = n.burnin,
#                        n.thin = n.thin,
#                        n.chains = n.chains,
#                        autojags = autojags,
#                        ...)
# time1 <- Sys.time()
# cat("MCMC sampling finished at", format(time1), "\n  ")
# print(time1-time0)
#
# mptInfo <- list(model="betaMPT",
#                 thetaNames = thetaNames,
#                 thetaUnique = thetaUnique,
#                 thetaFixed = unique(fixedPar$Parameter),
#                 MPT=mergedTree,
#                 eqnfile = eqnfile,
#                 data = data,
#                 restrictions = restrictions,
#                 covData=covData,
#                 covTable=covTable,
#                 corProbit=corProbit,
#                 transformedParameters=transformedPar$transformedParameters,
#                 hyperprior=list(alpha=alpha, beta=beta))
#
#
# # own summary (more stable than runjags)
# mcmc.summ <- summarizeMCMC(runjags$mcmc)
# # Beta MPT: rename parameters and get specific summaries
# summary <- summarizeMPT(mcmc = runjags$mcmc, mptInfo = mptInfo, M=M.T1, summ=mcmc.summ)
# if(dic){
#   summary$dic <-   extract(runjags, "dic", ...)
# }else{summary$dic <- NULL}
#
#
# # class structure for TreeBUGS
# # mcmc$BUGSoutput <- renameBUGSoutput(mcmc$BUGSoutput, thetaUnique, "betaMPT")
# fittedModel <- list(summary = summary,
#                     mptInfo = mptInfo,
#                     mcmc.summ = mcmc.summ,
#                     runjags = runjags,
#                     call = match.call(),
#                     time = time1-time0)
#
# class(fittedModel) <- "betaMPT"
#
# # write results to file
# writeSummary(fittedModel, parEstFile)
# gc(verbose=FALSE)
