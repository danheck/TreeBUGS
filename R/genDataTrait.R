
#' Generate Data for Trait MPT Models
#'
#' Generating a data file with known parameter structure using the Trait-MPT. Useful for simulations and robustness checks.
#' @inheritParams betaMPT
#' @inheritParams genMPT
#' @param N number of participants
#' @param mean Named vector of true group means of individual MPT parameters (probabilities in the interval [0,1]). If the vector is not named, the internal order of parameters is used (can be obtained using \code{\link{readEQN}}).
#' @param sigma (named) vector of group standard deviations of latent (!) individual MPT parameters. Default is zero (no person heterogeneity).
#' @param rho (named) correlation matrix for latent (!) individual MPT parameters. Must be symmetric and positive definite (e.g., no correlations of 1 or -1 allowed). Default: a diagonal matrix (i.e., zero correlations)
#'
#' Data are generated independently from the JAGS model files used for fitting the Trait-MPT model. If data for an equality-constrained version of the MPT model are required, the restrictions need to be hard-coded into the EQN-model file. Note that equal means still result in nonidentical MPT parameters on the individual level!
#'
#' @return a list including the generated frequencies (\code{data}) and the true, underlying parameters (\code{parameters})
#'
#' @examples
#' # Example: Standard Two-High-Threshold Model (2HTM)
#' EQNfile <- system.file("MPTmodels/2htm.eqn", package="TreeBUGS")
#' rho <- matrix(c(1,.8,.2,
#'                 .8,1,.1,
#'                 .2,.1,1), nrow=3)
#' colnames(rho) <- rownames(rho) <- c("Do","Dn","g")
#' genDat <- genTraitMPT(N = 100,
#'                       numItems = c(Target=250, Lure=250),
#'                       eqnfile = EQNfile,
#'                       mean = c(Do=.7, Dn=.7, g=.5),
#'                       sigma =   c(Do=.3, Dn=.3, g=.15),
#'                       rho = rho)
#' head(genDat$data, 3)
#' plotFreq(genDat$data, eqn=EQNfile)
#' @references Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70-98.
#' @seealso \code{\link{genMPT}}
#' @export
genTraitMPT <- function(N, numItems, eqnfile, restrictions,
                        mean=NULL, sigma=NULL, rho=NULL,
                        warning=TRUE){
  if(missing(restrictions))
    restrictions <- NULL
  Tree <- readEQN(eqnfile)
  # mergedTree <- mergeBranches(Tree)
  # thetaNames <- getParameter(mergedTree)
  mergedTree <- mergeBranches(Tree)
  Tree.restr <- thetaHandling(mergedTree,restrictions)
  thetaNames <-  Tree.restr$SubPar[,1:2]
  thetaNames <- thetaNames[rownames(unique(thetaNames[2])),]$Parameter
  treeLabels <- unique(mergedTree$Tree)
  S <- length(thetaNames)

  # default values:
  if(is.null(rho)){
    rho <- diag(S)
    dimnames(rho) <- list(thetaNames, thetaNames)
  }
  if(missing(sigma) || is.null(sigma)){
    sigma <- rep(0, S)
    names(sigma) <- thetaNames
  }
  sigma <- checkNaming(S, thetaNames, sigma, "sigma", warning=warning)
  rho <- checkNamingMatrix(S, thetaNames, rho, "rho", warning=warning)


  if(!is.null(mean)){
    mean <- checkNaming(S, thetaNames, mean, "mean", warning=warning)
  }else{
    stop("Either 'mean' must be provided.")
  }


  # generate multivariate normal using cholesky decomposition
  if(S>1){
    # iidNormal <- matrix(rnorm(N*S), nrow = N,
    #                     dimnames=list(NULL, thetaNames)) # i.i.d. standard normal
    covMatrix <- diag(sigma) %*% rho %*% diag(sigma)
    # cholesky <- chol(covMatrix)
    # thetaLatent <- matrix(qnorm(mean), N, S, byrow = TRUE) + iidNormal %*% cholesky
    thetaLatent <- mvrnorm(N, qnorm(mean), covMatrix)
  }else{
    thetaLatent <- matrix(rnorm(N, qnorm(mean), sigma), N, S)
  }
  colnames(thetaLatent) <- thetaNames
  theta <- pnorm(thetaLatent)

  # response frequencies:
  freq <- genMPT(theta, numItems, eqnfile, restrictions, warning=warning)

  res <- list(data = freq, parameters = list(theta=theta,
                                             thetaLatent = thetaLatent,
                                             mean=mean, sigma=sigma, rho=rho))
  return(res)
}


