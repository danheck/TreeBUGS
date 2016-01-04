#' Fit a Hierarchical latent-trait MPT Model
#'
#' Fits a latent-trait MPT model (Klauer, 2010) based on a standard MPT model file (.eqn) and individual data table (.csv).
#'
#' @inheritParams betaMPT
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
#' @references Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70-98.
#' @export

traitMPT <- function(eqnfile,  # statistical model stuff
                    data,
                    restrictions = NULL,
                    transformedParameters=NULL,

                    # hyperpriors:
                    mu = "dnorm(0,1)",
                    xi = "dunif(0,100)",
                    V=NULL,
                    df=NULL,

                    # MCMC stuff:
                    n.iter=50000,
                    n.burnin=5000,
                    n.thin=5,
                    n.chains=3,

                    # File Handling stuff:
                    modelfilename=NULL,
                    parEstFile = NULL,
                    sampler="JAGS",
                    autojags=FALSE,
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

  Tree <- readEQN(eqnfile)
  # reorder data:
  # data <- data[,sort(unique(Tree$Answers))] OLD Daniel fix

  mergedTree <- mergeBranches(Tree) # OLD mergeBranches ,names(data))
  data <- readSubjectData(data, unique(mergedTree$Category))

  tHoutput <- thetaHandling(mergedTree,restrictions)
  SubPar <- tHoutput$SubPar
  mergedTree <- tHoutput$mergedTree

  data <- data[,mergedTree$Category] #ordering data according to Tree
  thetaNames <- SubPar[,1:2]
  S <- max(SubPar$theta)
  isIdentifiable(S, mergedTree)

  # transformed parameters
  transformedPar <- getTransformed(model = "traitMPT",
                                   thetaNames = thetaNames,
                                   transformedParameters = transformedParameters)


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
                sampler=sampler,
                parString=transformedPar$modelstring)

  mcmc <- callingSampler(model = "traitMPT",
                         mergedTree = mergedTree,
                         data = data,
                         modelfile = modelfilename,
                         S = max(SubPar$theta),
                         transformedPar = transformedPar$transformedParameters,
                         hyperpriors = list(V=V, df=df),
                         n.iter = n.iter,
                         n.burnin = n.burnin,
                         n.thin = n.thin,
                         n.chains = n.chains,
                         sampler = sampler,
                         autojags = autojags,
                         ...)


  # Beta MPT: rename parameters and get specific summaries
  summary <- summarizeMPT(model = "traitMPT",
                          mcmc = mcmc,
                          thetaNames = thetaNames,
                          sampler = sampler,
                          transformedParameters = transformedPar$transformedParameters)

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

  class(fittedModel) <- "traitMPT"
  return(fittedModel)
}
