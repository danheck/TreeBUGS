
################# posterior predicted frequencies ###########################


#' Get Posterior Predictive Samples
#'
#' Draw frequencies based on posterior distribution of individual estimates.
#'
#' @param fittedModel fitted latent-trait or beta MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param M number of posterior predictive samples. As a maximum, the number of posterior samples in \code{fittedModel} is used.
#' @param expected if \code{TRUE}, the expected frequencies per person are returned (without additional sampling from a multinomial distribution)
#' @param nCPU number of CPUs used for parallel sampling. For large models and many participants, this requires a lot of memory.
#' @return a list of \code{M} matrices with individual frequencies (rows=participants, columns=MPT categories). For \code{M=1}, a single matrix is returned.
#' @export
#' @importFrom parallel clusterExport makeCluster stopCluster parLapply parApply
#' @importFrom  stats cor cov2cor density rmultinom
posteriorPredictive <- function(fittedModel, M=100, expected=FALSE, nCPU=4){

  mptInfo <- fittedModel$mptInfo
  # get information about model:
  dat <- mptInfo$dat
  tree <- mptInfo$MPT$Tree
  TreeNames <- unique(tree)
  # selection list for mapping of categories to trees:
  sel.cat <- lapply(TreeNames, function(tt) tree %in% tt)

  N <- nrow(mptInfo$data)
  S <- length(mptInfo$thetaUnique)
  numTrees <- length(TreeNames)
  chains <- length(fittedModel$runjags$mcmc)
  numItems <-   t(apply(mptInfo$data, 1,
                        function(x)  tapply(x, mptInfo$MPT$Tree, sum)))

  sample <- nrow(fittedModel$runjags$mcmc[[1]])
  max.samp <- min(sample, ceiling(M/chains))
  # fixed effects:
  sel.thetaFE <- grep("thetaFE", varnames(fittedModel$runjags$mcmc), fixed=TRUE)
  # standard random effects:
  sel.var <- setdiff(grep("theta", varnames(fittedModel$runjags$mcmc), fixed=TRUE), sel.thetaFE)

  n.thetaFE <- length(sel.thetaFE)

  par.ind <- par.thetaFE <- c()
  for(m in 1:chains){
    sel.samp <- sample(1:sample, max.samp)
    par.ind <- rbind(par.ind, as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, sel.var]))
    if(n.thetaFE > 0){
      par.thetaFE <- rbind(par.thetaFE, as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, sel.thetaFE]))
    }else{
      par.thetaFE <- NULL
    }
  }

  expectedFreq <- function(n, theta, thetaFE){
    sapply(mptInfo$MPT$Equation,USE.NAMES = FALSE,
           function(ff) {
             eval(parse(text=ff),
                  envir = list(n=n, theta=theta, thetaFE=thetaFE) )
           })
  }

  getPostPred <- function(tt){

    # single replication theta:
    theta <- matrix(tt[(n.thetaFE+1):length(tt)], S, N)
    if(n.thetaFE>0){
      thetaFE <- tt[1:n.thetaFE]
    }else{
      thetaFE <- NULL
    }

    freq.exp <- t(sapply(1:N, expectedFreq,
                         theta = theta, thetaFE = thetaFE))*numItems[,tree]

    if(!expected){
      # multinomial sampling:
      for(k in 1:length(TreeNames)){
        freq.exp[,sel.cat[[k]]] <- t(apply(freq.exp[,sel.cat[[k]]], 1,
                                           function(x)
                                             rmultinom(1, size=sum(x), prob=x/sum(x))))
      }
    }
    colnames(freq.exp) <- mptInfo$MPT$Category

    list(freq.exp)
  }

  if(nCPU >1){
    cl <- makeCluster(nCPU)
    clusterExport(cl, c("S","N","mptInfo","numItems","expected","expectedFreq",
                        "tree","TreeNames","n.thetaFE","sel.cat"), envir=environment())
    # loop across replications:
    freq.list <- parApply(cl, cbind(par.thetaFE, par.ind), 1, getPostPred)
    stopCluster(cl)
  }else{
    freq.list <- apply(cbind(par.thetaFE, par.ind), 1, getPostPred)
  }

  # remove strange list structure:
  freq.list <- lapply(freq.list, function(xx) xx[[1]])

  if(M == 1){
    freq.list[[1]]
  }else{
    freq.list[1:min(M, length(freq.list))]
  }
}
