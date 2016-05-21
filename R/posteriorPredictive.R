
################# posterior predicted frequencies ###########################


#' Get Posterior Predictive Samples
#'
#' Draw frequencies based on posterior distribution of individual estimates.
#'
#' @param fittedModel fitted latent-trait or beta MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param M number of posterior predictive samples. As a maximum, the number of posterior samples in \code{fittedModel} is used.
#' @param expected if \code{TRUE}, the expected frequencies per person are returned (without additional sampling from a multinomial distribution)
#' @param nCPU number of CPUs used for parallel sampling
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
  sel <- lapply(TreeNames, function(tt) tree %in% tt)

  N <- nrow(mptInfo$data)
  S <- length(mptInfo$thetaUnique)
  numTrees <- length(TreeNames)
  chains <- length(fittedModel$runjags$mcmc)
  numItems <-   t(apply(mptInfo$data, 1,
                        function(x)  tapply(x, mptInfo$MPT$Tree, sum)))

  sample <- nrow(fittedModel$runjags$mcmc[[1]])
  max.samp <- min(sample, ceiling(M/chains))
  sel.var <- grep("theta", varnames(fittedModel$runjags$mcmc), fixed=TRUE)
  par.ind <- c()
  for(m in 1:chains){
    sel.samp <- sample(1:sample, max.samp)
    par.ind <- rbind(par.ind, fittedModel$runjags$mcmc[[m]][sel.samp, sel.var])
  }

  cl <- makeCluster(nCPU)
  clusterExport(cl, c("S","N","mptInfo","numItems", "TreeNames"), envir=environment())

  # loop across replications:
  freq.list <- parApply(cl, par.ind, 1, function(tt){

    # single replication theta:
    theta <- matrix(tt, S, N)

    # loop across participants:
    freq.rep <- t(sapply(1:N, function(n){

      # single participant:
      freq <- rep(NA, nrow(mptInfo$MPT))
      names(freq) <- mptInfo$MPT$Category
      prob <- sapply(mptInfo$MPT$Equation,
                     function(ff) eval(parse(text=ff)))

      # EXPECTED frequencies:
      freq <- prob*numItems[n,tree]
      names(freq) <- mptInfo$MPT$Category
      return(freq)
    }
    ))

    # return list to avoid parsing to matrix (loss of dimensions)
    list(freq.rep)
  })

  # remove strange list structure:
  freq.list <- lapply(freq.list, function(xx) xx[[1]])

  # posterior PREDICTIVE sampling:
  if(!expected){
    # multinomial sampling:
    freq.list <- parLapply(cl, freq.list, function(fe){
      for(k in 1:length(TreeNames)){
        fe[,sel[[k]]] <- t(apply(fe[,sel[[k]]], 1,
                                 function(x) rmultinom(1, size=sum(x), prob=x/sum(x))))
      }
      fe
    })
  }
  stopCluster(cl)


  if(M == 1){
    freq.list[[1]]
  }else{
    freq.list[1:min(M, length(freq.list))]
  }
}
