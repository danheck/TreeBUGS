#' Generate Data for Latent-Trait MPT Models
#'
#' Generating a data set with known parameter structure using the Trait-MPT.
#' Useful for simulations and robustness checks.
#'
#' @inheritParams betaMPT
#' @inheritParams genMPT
#' @param N number of participants
#' @param mean named vector of data-generating group means of the individual MPT parameters
#'     on the probability scale.
#'     If the vector is not named, the internal order of parameters is used
#'     (can be obtained using \code{\link{readEQN}}).
#' @param mu an alternative way to define the group-level means on the latent-probit scale
#'     (i.e., \code{mu = qnorm(mean)} or equivalently, \code{mean = pnorm(mu)}).
#' @param sigma (named) vector of group standard deviations of
#'     individual MPT parameters on the latent probit scale.
#'     Default is zero (no person heterogeneity).
#' @param rho (named) correlation matrix for individual MPT parameters on the
#'     latent probit scale. Must be symmetric and positive definite (e.g., no
#'     correlations of 1 or -1 allowed).
#'     Default: a diagonal matrix (i.e., zero correlations).
#'
#' @details
#' This functions implements a two-step sampling procedure. First, the person
#' parameters on the latent probit-scale are sampled from the multivariate normal
#' distribution (based on the mean \code{mu = qnorm(mean)}, the standard deviations
#' \code{sigma}, and the correlation matrix \code{rho}).
#' These person parameters are then transformed to the probability scale using
#' the probit-link.
#' In a last step, observed frequencies are sampled for each person using the MPT equations.
#'
#' Note that the user can generate more complex structures for the latent person parameters,
#' and then supply these person parameters to the function \code{\link{genMPT}}.
#'
#' @return a list including the generated frequencies per person (\code{data})
#'     and the sampled individual parameters (\code{parameters}) on the probit
#'     and probability scale (\code{thetaLatent} and \code{theta}, respectively).
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
                        mean, mu, sigma, rho,
                        warning = TRUE){

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
  if (missing(rho) || is.null(rho)){
    rho <- diag(S)
    dimnames(rho) <- list(thetaNames, thetaNames)
  } else if(S==1){
    rho <- as.matrix(rho)
  }
  if (missing(sigma) || is.null(sigma)){
    sigma <- rep(0, S)
    names(sigma) <- thetaNames
  }
  sigma <- checkNaming(S, thetaNames, sigma, "sigma",
                       interval = c(0, Inf), warning=warning)
  rho <- checkNamingMatrix(S, thetaNames, rho, "rho", warning=warning)

  if (!missing(mu) && !is.null(mu)){
    if (!missing(mean) && !is.null(mean))
      stop ("Only one of 'mean' and 'mu' can be defined.")
    mu <- checkNaming(S, thetaNames, mu, "mu",
                      interval = c(-Inf, Inf), warning=warning)
    mean <- pnorm(mu)
  } else if (!missing(mean) && !is.null(mean)){
    mean <- checkNaming(S, thetaNames, mean, "mean",
                        interval = c(0, 1), warning=warning)
    mu <- qnorm(mean)
  } else {
    stop("Either 'mean' or 'mu' must be provided.")
  }


  # generate multivariate normal using cholesky decomposition
  if(S>1){
    # iidNormal <- matrix(rnorm(N*S), nrow = N,
    #                     dimnames=list(NULL, thetaNames)) # i.i.d. standard normal
    covMatrix <- diag(sigma) %*% rho %*% diag(sigma)
    # cholesky <- chol(covMatrix)
    # thetaLatent <- matrix(qnorm(mean), N, S, byrow = TRUE) + iidNormal %*% cholesky
    thetaLatent <- mvrnorm(N, mu, covMatrix)
  }else{
    thetaLatent <- rnorm(N, mu, sigma)
  }
  dim(thetaLatent) <- c(N, S)
  colnames(thetaLatent) <- thetaNames
  theta <- pnorm(thetaLatent)

  # response frequencies:
  freq <- genMPT(theta, numItems, eqnfile, restrictions, warning=warning)
  list(data = freq,
       parameters = list(theta=theta, thetaLatent = thetaLatent,
                         mean=mean, mu=mu, sigma=sigma, rho=rho))
}


