
#' Generate Data for Beta MPT Models
#'
#' Generating a data file with known parameter structure using the Beta-MPT. Useful for robustness checks.
#'
#' @inheritParams betaMPT
#' @param N number of participants
#' @param numItems number of responses per tree (a named vector with tree labels)
#' @param mean Named vector of true group means of individual MPT parameters. If the vector is not named, the internal order of parameters is used (can be obtained using \code{\link{readEQN}}).
#' @param sd named vector of group standard deviations of individual MPT parameters.
#' @param alpha Alternative specification of the group-level distribution using the shape parameters of the beta distribution (see \link{dbeta}).
#' @param beta see \code{alpha}
#'
#' Data are generated independently from the JAGS model files used for fitting the Beta-MPT model. If data for an equality-constrained version of the MPT model are required, the restrictions need to be hard-coded into the EQN-model file. Note that equal means still result in nonidentical MPT parameters on the individual level!
#'
#' @return a list including the generated frequencies (\code{data}) and the true, underlying parameters (\code{parameters})
#'
#' @examples
#' # Example: Standard Two-High-Threshold Model (2HTM)
#' EQNfile <- paste0(.libPaths()[1], "/TreeBUGS/MPTmodels/2htm.eqn")
#' genDat <- genBetaMPT(N = 100,
#'                      numItems = c(Target=250, Lure=250),
#'                      eqnfile = EQNfile,
#'                      mean = c(Do=.7, Dn=.7, g=.5),
#'                      sd =   c(Do=.1, Dn=.1, g=.05))
#' @importFrom stats  rbeta
#' @references Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial processing tree models for addressing individual differences. Journal of Mathematical Psychology, 54, 167-183.
#' @export
genBetaMPT <- function(N, numItems, eqnfile, mean=NULL, sd=NULL, alpha=NULL, beta=NULL){

  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  thetaNames <- getParameter(mergedTree)
  S <- length(thetaNames)


  if(!is.null(mean) & !is.null(sd)){
    mean <- checkNaming(S, thetaNames, mean, "mean")
    sd <- checkNaming(S, thetaNames, sd, "sd")

    alpha <- -(mean*(mean^2-mean+sd^2))/sd^2
    beta <- ((mean-1)*(mean^2-mean+sd^2))/sd^2
    if(any(alpha<=0) | any(beta<=0))
      stop("Check numerical values for mean and sd, result in negative alpha/beta parameters of beta-hyperprior distribution.")
  }else if(!is.null(alpha) & !is.null(beta)){
    alpha <- checkNaming(S, thetaNames, alpha, "alpha")
    beta <- checkNaming(S, thetaNames, beta, "beta")

  }else{
    stop("Either 'mean' and 'sd' or 'alpha' and 'beta' must be provided.")
  }

  individPar <- c()
  freq <- matrix(NA, N, nrow(mergedTree), dimnames=list(NULL, mergedTree$Category))
  for(s in 1:S){
    individPar <- cbind(individPar, rbeta(N, shape1 = alpha[s], shape2 = beta[s]))
  }
  colnames(individPar) <- thetaNames
  for(n in 1:N){
    mergedTree$prob <- sapply(mergedTree$Equation,
                              function(ff) eval(parse(text=ff),
                                                as.list(individPar[n,])))

    numTrees <- length(unique(mergedTree$Tree))
    for(k in 1:numTrees){
      sel <- mergedTree$Tree %in% names(numItems)[k]
      cat <- findInterval(runif(numItems[k]), cumsum(mergedTree$prob[sel]))+1
      catLabel <- mergedTree$Category[sel][cat]
      freq[n,mergedTree$Category[sel]] <- table(factor(catLabel,
                                                       levels=mergedTree$Category[sel]),
                                                exclude=NA)
    }
  }



  res <- list(data = freq, parameters = list(individPar=individPar, mean=mean,
                                            sd=sd, alpha=alpha, beta=beta))
  return(res)
}

#' Generate Data for Trait MPT Models
#'
#' Generating a data file with known parameter structure using the Trait-MPT. Useful for robustness checks.
#'
#' @inheritParams betaMPT
#' @param N number of participants
#' @param numItems number of responses per tree (a vector, possibly named with tree labels)
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
#' EQNfile <- paste0(.libPaths()[1], "/TreeBUGS/MPTmodels/2htm.eqn")
#' genDat <- genTraitMPT(N = 100,
#'                       numItems = c(Target=250, Lure=250),
#'                       eqnfile = EQNfile,
#'                       mean = c(Do=.7, Dn=.7, g=.5),
#'                       sigma =   c(Do=.3, Dn=.3, g=.15),
#'                       rho = matrix(c(1 , .8, .2,
#'                                      .8, 1 , .1,
#'                                      .2, .1, 1 ),
#'                                    nrow=3))
#'
#' @references Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70-98.
#' @export
genTraitMPT <- function(N, numItems, eqnfile, mean=NULL, sigma=NULL, rho=NULL){

  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  thetaNames <- getParameter(mergedTree)
  S <- length(thetaNames)

  if(is.null(rho)){
    rho <- diag(S)
    dimnames(rho) <- list(thetaNames, thetaNames)
  }
  if(is.null(sigma)){
    sigma <- rep(0, S)
    names(sigma) <- thetaNames
  }
  sigma <- checkNaming(S, thetaNames, sigma, "sigma")
  rho <- checkNamingMatrix(S, thetaNames, rho, "rho")


  if(!is.null(mean)){
    mean <- checkNaming(S, thetaNames, mean, "mean")
  }else{
    stop("Either 'mean' must be provided.")
  }

  iidNormal <- matrix(rnorm(N*S), nrow = N, dimnames=list(NULL, thetaNames))

  # generate multivariate normal using cholesky decomposition
  covMatrix <- diag(sigma) %*% rho %*% diag(sigma)
  cholesky <- chol(covMatrix)
  individParLatent <- mean + iidNormal %*% cholesky
  colnames(individParLatent) <- thetaNames
  individPar <- pnorm(individParLatent)

  freq <- matrix(NA, N, nrow(mergedTree), dimnames=list(NULL, mergedTree$Category))
  for(n in 1:N){
    mergedTree$prob <- sapply(mergedTree$Equation,
                              function(ff) eval(parse(text=ff),
                                                as.list(individPar[n,])))

    numTrees <- length(unique(mergedTree$Tree))
    for(k in 1:numTrees){
      sel <- mergedTree$Tree %in% names(numItems)[k]
      cat <- findInterval(runif(numItems[k]), cumsum(c(0, mergedTree$prob[sel], 1)))#+1
      catLabel <- mergedTree$Category[sel][cat]
      freq[n,mergedTree$Category[sel]] <- table(factor(catLabel,
                                                       levels=mergedTree$Category[sel]),
                                                exclude=NA)
    }
  }



  res <- list(data = freq, parameters = list(individPar=individPar, individParLatent = individParLatent,
                                             mean=mean, sigma=sigma, rho=rho))
  return(res)
}






################### HELPER FUNCTIONS ####################################


checkNamingMatrix <- function(S, thetaNames, matrix, matrixName = "rho"){
  if(any(S != dim(matrix)))
    stop("Dimensions of matrix '",matrixName,"' not correct, should be ", S, "x", S)
  if(is.null(dimnames(matrix))){
    warning("Matrix '",matrixName,"' not named. Internal order of parameters is used, see ?readMultiTree and check parameters by generatedData$parameters")
    dimnames(matrix) <- list(thetaNames, thetaNames)
  }else if(any(sort(thetaNames) != sort(rownames(matrix)))){
    stop("Row names of matrix '", matrixName,"' do not match parameter labels in eqn file.")
  }else if(any(sort(thetaNames) != sort(colnames(matrix)))){
    stop("Column names of matrix '", matrixName,"' do not match parameter labels in eqn file.")
  }else{
    matrix <- matrix[thetaNames,]
    matrix <- matrix[,thetaNames]
  }

  if(any(diag(matrix) != 1))
    stop("Diagonal must contain ones!")

  if(any(matrix != t(matrix)))
    stop("Matrix ", matrixName, " must be symmetric!")

  if(any(matrix < -1 | matrix >1))
    stop("'",matrixName,"' cannot be negative!")

  return(matrix)
}

checkNaming <- function(S, thetaNames, vector, vectorName){
  if(S != length(vector))
    stop("Length of '",vectorName,"' not correct, should be ", S)
  if(is.null(names(vector))){
    warning("Vector '",vectorName,"' not named. Internal order of parameters is used, see ?readMultiTree and check parameters by generatedData$parameters")
    names(vector) <- thetaNames
#     cat(vectorName, ":\n")
#     print(vector)
  }else if(any(sort(thetaNames) != sort(names(vector)))){
    stop("Parameter names of vector '",vectorName,"' do not match parameter labels in eqn file.")
  }else{
    vector <- vector[thetaNames]
  }

  if(any(vector<0))
    stop("'",vectorName,"' cannot be negative!")

  return(vector)
}
