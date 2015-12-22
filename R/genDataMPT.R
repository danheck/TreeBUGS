
#' Simulate Data for Beta MPT Models
#'
#' Generating a data file with known parameter structure using the Beta-MPT. Useful for robustness checks.
#'
#' @inheritParams betaMPT
#' @param N number of participants
#' @param numItems number of responses per tree (possibly a named vector with tree labels)
#' @param mean Named vector of true group means of individual MPT parameters. If the vector is not named, the internal order of parameters is used (can be obtained using \code{\link{readEQN}}).
#' @param sd named vector of group standard deviations of individual MPT parameters.
#' @param alpha Alternative specification of the group-level distribution using the shape parameters of the beta distribution (see \link{dbeta}).
#' @param beta see \code{alpha}
#'
#' Data are generated independently from the JAGS model files used for fitting the Beta-MPT model.
#' @return a list including the generated frequencies (\code{data}) and the true, underlying parameters (\code{parameters})
#'
#' @examples
#' \dontrun{
#' genDat <- genBetaMPT(N = 100,
#'                      numItems = c(Target=250, Lure=250),
#'                      eqnfile = "2htm.eqn",
#'                      mean = c(Do=.7, Dn=.7, g=.5),
#'                      sd =   c(Do=.1, Dn=.1, g=.05))
#'                      }
#' @importFrom stats  rbeta
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


checkNaming <- function(S, thetaNames, vector, vectorName){
  if(S != length(vector))
    stop("Length of '",vectorName,"' not correct, should be ", S)
  if(is.null(names(vector))){
    warning("Vector '",vectorName,"' not named. Internal order of parameters is used, see ?readMultiTree")
  }else if(any(sort(thetaNames) != sort(names(vector)))){
    stop("Parameter names of vector '",vectorName,"' do not match parameter labels in eqn file.")
  }else{
    vector <- vector[thetaNames]
  }

  if(any(vector<0))
    stop("'",vectorName,"' cannot be negative!")

  return(vector)
}
