
################# posterior predicted frequencies ###########################


#' Get Posterior Predictive Samples
#'
#' Draw frequencies based on posterior distribution of individual estimates.
#'
#' @param fittedModel fitted latent-trait or beta MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param M number of posterior predictive samples. As a maximum, the number of posterior samples in \code{fittedModel} is used.
#' @param expected if \code{TRUE}, the expected frequencies per person are returned (without additional sampling from a multinomial distribution)
#' @return a list of \code{M} matrices with individual frequencies (rows=participants, columns=MPT categories). For \code{M=1}, a single matrix is returned.
#' @export
posteriorPredictive <- function(fittedModel, M=100, expected=FALSE){

  # get information about model:
  dat <- fittedModel$mptInfo$dat
  tree <- fittedModel$mptInfo$MPT$Tree
  TreeNames <- unique(tree)
  # selection list for mapping of categories to trees:
  sel <- lapply(TreeNames, function(tt) tree %in% tt)

  N <- nrow(fittedModel$mptInfo$data)
  S <- length(fittedModel$mptInfo$thetaUnique)
  numTrees <- length(TreeNames)
  chains <- length(fittedModel$runjags$mcmc)
  numItems <-   t(apply(fittedModel$mptInfo$data, 1,
                        function(x)  tapply(x, fittedModel$mptInfo$MPT$Tree, sum)))

  sample <- nrow(fittedModel$runjags$mcmc[[1]])
  max.samp <- min(sample, ceiling(M/chains))
  sel.var <- grep("theta", varnames(fittedModel$runjags$mcmc), fixed=TRUE)
  par.ind <- c()
  for(m in 1:chains){
    sel.samp <- sample(1:sample, max.samp)
    par.ind <- rbind(par.ind, fittedModel$runjags$mcmc[[m]][sel.samp, sel.var])
  }

  # loop across replications:
  freq.list <- apply(par.ind, 1, function(tt){

    # single replication theta:
    theta <- matrix(tt, S, N)

    # loop across participants:
    freq.rep <- t(sapply(1:N, function(n){

      # single participant:
      freq <- rep(NA, nrow(fittedModel$mptInfo$MPT))
      names(freq) <- fittedModel$mptInfo$MPT$Category
      prob <- sapply(fittedModel$mptInfo$MPT$Equation,
                     function(ff) eval(parse(text=ff)))

      # EXPECTED frequencies:
      freq <- prob*numItems[n,tree]
      names(freq) <- fittedModel$mptInfo$MPT$Category
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
    freq.list <- lapply(freq.list, function(fe){
      for(k in 1:length(TreeNames)){
        fe[,sel[[k]]] <- t(apply(fe[,sel[[k]]], 1,
                                 function(x) rmultinom(1, size=sum(x), prob=x/sum(x))))
      }
      fe
    })
  }

  if(M == 1){
    freq.list[[1]]
  }else{
    freq.list[1:min(M, length(freq.list))]
  }
}
