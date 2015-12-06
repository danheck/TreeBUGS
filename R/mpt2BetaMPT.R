#' Fit a Hierarchical Beta-MPT Model
#'
#' Fits a Beta-MPT model based on a standard MPT model file (.eqn) and individual data table (.csv).
#'
#' @param eqnfile The (full path to the) file that specifies the MPT model (standard .eqn syntax)
#' @param data The (full path to the) cvs file with the data (comma separated; category labels in first row)
#' @param restrictions  Optional: The (full path to the) file that specifies which parameters should be constants and which should be equal. Alternatively: a list of restrictions, e.g., \code{list("D1=D2","g=0.5")}
#' @param transformedParameters list with parameter transformations that should be computed based on the posterior samples (e.g., for testing parameter differences: \code{list("diffD=Do-Dn")})
#' @param modelfilename Name that the modelfile that is made by the function to work with WinBUGS should get.
#'        Default is to write this information to the tempdir as required by CRAN standards.
#' @param alpha Hyperprior of for the alpha and beta parameters (default: uniform prior on the interval [1,5000]).
#' @param beta Second hyperprior, see \code{alpha}
#' @param parEstFile Name of the file to with the estimates should be stored (e.g., "parEstFile.txt")
#' @param n.iter Number of iterations.
#' @param n.burnin Burnin period.
#' @param n.thin Thinning rate.
#' @param n.chains number of MCMC chains
#' @param sampler Which sampler should be used? Default is "JAGS". Further options are "OpenBUGS" and "WinBugs" (without MPT specific summary and plotting functions)
#' @param autojags whether to run JAGS until the MCMC chains converge (see \link{autojags}). Can take a lot of time for large models.
#' @param ... Arguments to be passed to the sampling function (default: \code{\link{jags.parallel}}.
#'
#' @return a list of the class \code{betaMPT} with the objects:
#' \itemize{
#'  \item \code{summary}: MPT tailored summary. Use \code{summary(fittedModel)}
#'  \item \code{mptInfo}: info about MPT model (eqn and data file etc.)
#'  \item \code{mcmc}: the object returned from the MCMC sampler. Standard: An \code{\link{jags.parallel}} object. Note that the sample can be transformed into an \code{mcmc.list} for analysis with the \code{coda} package by \code{as.mcmc.list(fittedModel$mcmc$BUGSoutput)}
#'  \item \code{sampler}: the type of sampler used (standard: \code{"JAGS"})
#' }
#' @author Nina R. Arnold, Denis Arnold, Daniel Heck
#' @export

mpt2BetaMPT<-function(eqnfile,  # statistical model stuff
                   data,
                   restrictions = NULL,
                   transformedParameters=NULL,
                   alpha="dunif(1,5000)",
                   beta="dunif(1,5000)",

                   # MCMC stuff:
                   n.iter=100000,
                   n.burnin=50000,
                   n.thin=2,
                   n.chains=3,

                   # File Handling stuff:
                   modelfilename=NULL,
                   parEstFile = NULL,
                   sampler="JAGS",
                   autojags=TRUE,
                   ...){

  if(is.null(modelfilename)){
    modelfilename=tempfile(pattern="MODELFILE",fileext=".txt")
  }

  ### read data
  if(is.matrix(data) | is.data.frame(data)){
    data <- as.data.frame(data)
  }else{
    data = read.csv(data, header=TRUE)
  }
  if(any(is.na(data))){
    warning("Check for missings in the data.")
  }

  Tree=readMultiTree(eqnfile)
  data=readSubjectData(data,unique(Tree$Answers))
  # reorder data:
  # data <- data[,sort(unique(Tree$Answers))] OLD Daniel fix

  Tree=mergeBranches(Tree) # OLD mergeBranches ,names(data))

  tHoutput=thetaHandling(Tree,restrictions)
  SubPar=tHoutput[[1]]
  Tree=tHoutput[[2]]

  data=data[,Tree$Answers] #ordering data according to Tree
  thetaNames <- tHoutput[[1]][,1:2]

  # transformed parameters
  transformedPar <- getTransformed(thetaNames, transformedParameters)

  makeModelDescription(modelfilename,Tree ,max(SubPar$theta),
                       alpha=alpha,beta = beta,
                       sampler=sampler,
                       parString=transformedPar$modelstring)
  mcmc <- callingBetaMPT(Tree,
                            data,
                            modelfile=modelfilename,
                            numberOfParameters=max(SubPar$theta),
                            transformedPar=transformedPar$transformedParameters,
                            # parameters,
                            n.iter=n.iter,
                            n.burnin=n.burnin,
                            n.thin=n.thin,
                            n.chains=n.chains,
                            sampler=sampler,
                            autojags=autojags,
                            ...)

  # Beta MPT: rename parameters and get specific summaries
  summary <- summarizeBetaMPT(mcmc, thetaNames, sampler=sampler,
                              transformedParameters=transformedPar$transformedParameters)

  mptInfo <- list(thetaNames = thetaNames,
                  eqnfile=eqnfile,
                  data=data,
                  restrictions=restrictions)

  # class structure for TreeBUGS
  fittedModel <- list(summary=summary,
                      mptInfo=mptInfo,
                      mcmc=mcmc,
                      sampler=sampler,
                      call=match.call()  )

  # write results
  if(!(missing(parEstFile) | is.null(parEstFile))){
    write.table(summary,  file=parEstFile, sep ="\t",
                na="NA",dec=".",row.names=T,col.names=T,quote=F)
  }

  class(fittedModel) <- "betaMPT"
  return(fittedModel)
}


