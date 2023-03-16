
################# posterior predicted frequencies ###########################


#' Get Posterior Predictive Samples
#'
#' Draw predicted frequencies based on posterior distribution of (a) individual estimates (default) or (b) for a new participant (if \code{numItems} is provided; does not consider continuous or discrete predictors in traitMPT).
#'
#' @param fittedModel fitted latent-trait or beta MPT model (\code{\link{traitMPT}}, \code{\link{betaMPT}})
#' @param M number of posterior predictive samples. As a maximum, the number of posterior samples in \code{fittedModel} is used.
#' @param numItems optional: a vector with the number of items per MPT tree to sample predicted data for a new participant (first, a participant vector \eqn{\theta} is sampled from the hierarchical posterior; second, frequencies are generated).
#' @param expected if \code{TRUE}, the expected frequencies per person are returned (without additional sampling from a multinomial distribution)
#' @param nCPU number of CPUs used for parallel sampling. For large models and many participants, this requires considerable computer-memory resources (as a remedy, use \code{nCPU=1}).
#' @return by default, a list of \code{M} posterior-predictive samples (i.e., matrices) with individual frequencies (rows=participants, columns=MPT categories). For \code{M=1}, a single matrix is returned. If \code{numItems} is provided, a matrix with samples for a new participant is returned (rows=samples)
#' @export
#' @examples
#' \dontrun{
#' # add posterior predictive samples to fitted model
#' #     (facilitates plotting using ?plotFit)
#' fittedModel$postpred$freq.pred <-
#'   posteriorPredictive(fittedModel, M = 1000)
#' }
#' @importFrom parallel clusterExport makeCluster stopCluster parLapply parApply
#' @importFrom  stats cor cov2cor density rmultinom
posteriorPredictive <- function(fittedModel,
                                M = 100,
                                numItems = NULL,
                                expected = FALSE,
                                nCPU = 4) {
  mptInfo <- fittedModel$mptInfo
  # get information about model:
  tree <- mptInfo$MPT$Tree
  TreeNames <- unique(tree)
  # selection list for mapping of categories to trees:
  sel.cat <- lapply(TreeNames, function(tt) tree %in% tt)

  S <- length(mptInfo$thetaUnique)
  numTrees <- length(TreeNames)
  chains <- length(fittedModel$runjags$mcmc)
  sample <- nrow(fittedModel$runjags$mcmc[[1]])
  max.samp <- min(sample, ceiling(M / chains))
  # fixed effects:
  sel.thetaFE <- grep("thetaFE", varnames(fittedModel$runjags$mcmc), fixed = TRUE)
  n.thetaFE <- length(sel.thetaFE)


  ########### sample data for participants in data set!
  if (missing(numItems) || is.null(numItems)) {
    pred.new <- FALSE
    N <- nrow(mptInfo$data)
    treeLabels <- unique(mptInfo$MPT$Tree)
    numItems <- matrix(NA, nrow(mptInfo$data), length(treeLabels))
    colnames(numItems) <- treeLabels
    for (tl in treeLabels) {
      numItems[, tl] <- rowSums(mptInfo$data[, mptInfo$MPT$Tree == tl, drop = FALSE])
    }
    # numItems <-   t(apply(mptInfo$data[,mptInfo$MPT$Category,drop=FALSE], 1,
    #                       function(x)  tapply(x, mptInfo$MPT$Tree, sum)))

    ########### sample data for NEW participant!
  } else {
    pred.new <- TRUE
    if (length(numItems) != length(unique(fittedModel$mptInfo$MPT$Tree))) {
      stop("Length of 'numItems' does not match number of MPT trees.")
    }
    if (is.null(names(numItems))) {
      warning("numItems are ordered alphabetically:\n  ", paste(TreeNames, sep = ", "))
      names(numItems) <- TreeNames
    }
    numItems <- matrix(numItems, nrow = 1, dimnames = list(NULL, names(numItems)))
    N <- 1
  }

  ############# SAMPLING
  par.ind <- par.thetaFE <- c()
  for (m in 1:chains) {
    sel.samp <- sample(1:sample, max.samp)

    ######## standard random effects:
    if (!pred.new) {
      sel.var <- setdiff(
        grep("theta", varnames(fittedModel$runjags$mcmc), fixed = TRUE),
        sel.thetaFE
      )
      par.tmp <- as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, sel.var, drop = FALSE])
      if (ncol(par.tmp) == 0){
        stop("No MCMC samples for the individual-level MPT parameters (theta) were stored. \n",
             "Please re-fit model with the argument:  monitorIndividual = TRUE")
      }
    } else {
      par.tmp <- matrix(NA, max.samp, S,
        dimnames = list(NULL, paste0("theta[", 1:S, ",1]"))
      )
      # sample hierarchical values, then individuals
      if (inherits(fittedModel, "betaMPT")) {
        if (S == 1) {
          alpha <- as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, "alph", drop = FALSE])
          beta <- as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, "bet", drop = FALSE])
        } else {
          alpha <- as.matrix(fittedModel$runjags$mcmc[[m]][
            sel.samp, paste0("alph[", 1:S, "]"),
            drop = FALSE
          ])
          beta <- as.matrix(fittedModel$runjags$mcmc[[m]][
            sel.samp, paste0("bet[", 1:S, "]"),
            drop = FALSE
          ])
        }
        for (i in 1:S) {
          par.tmp[, i] <- rbeta(max.samp, alpha[, i], beta[, i])
        }
      } else if (inherits(fittedModel, "traitMPT")) {
        if (S == 1) {
          mu <- as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, "mu", drop = FALSE])
          sig <- as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, "sigma", drop = FALSE])
        } else {
          mu <- as.matrix(fittedModel$runjags$mcmc[[m]][
            sel.samp, paste0("mu[", 1:S, "]"),
            drop = FALSE
          ])
          sig <- as.matrix(fittedModel$runjags$mcmc[[m]][
            sel.samp, paste0("sigma[", 1:S, "]"),
            drop = FALSE
          ])
        }
        sel.rho <- grep("rho", varnames(fittedModel$runjags$mcmc))
        rho <- as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, sel.rho, drop = FALSE])
        for (mm in 1:max.samp) {
          Sig <- matrix(rho[mm, ], S, S) * (sig[mm, ] %*% t(sig[mm, ]))
          par.tmp[mm, ] <- pnorm(mvrnorm(1, mu[mm, ], Sig))
        }
      }
    }
    par.ind <- rbind(par.ind, par.tmp)
    if (n.thetaFE > 0) {
      par.thetaFE <- rbind(
        par.thetaFE,
        as.matrix(fittedModel$runjags$mcmc[[m]][sel.samp, sel.thetaFE])
      )
    } else {
      par.thetaFE <- NULL
    }
  }

  expectedFreq <- function(n, theta, thetaFE) {
    # p <- tapply(mptInfo$MPT$c*apply(theta[,n]^t(mptInfo$MPT$a) * theta[,n]^t(mptInfo$MPT$b),2,prod),
    #             mptInfo$MPT$map, sum)
    # names(p) <- mptInfo$MPT$cat.names
    # p[mptInfo$MPT$Category]
    # mptInfo$MPT$Category
    sapply(mptInfo$MPT$Equation,
      USE.NAMES = FALSE,
      function(ff) {
        eval(parse(text = ff),
          envir = list(n = n, theta = theta, thetaFE = thetaFE)
        )
      }
    )
  }

  getPostPred <- function(tt) {
    # single posterior sample for all individual parameters theta:
    theta <- matrix(tt[(n.thetaFE + 1):length(tt)], S, N, byrow = FALSE)
    if (n.thetaFE > 0) {
      thetaFE <- tt[1:n.thetaFE]
    } else {
      thetaFE <- NULL
    }

    freq.exp <- t(sapply(1:N, expectedFreq,
      theta = theta, thetaFE = thetaFE
    )) * numItems[, tree]

    if (!expected) {
      # multinomial sampling:
      rmultinom_stable <- function(x) {
        if (sum(x) > 0) {
          rand <- rmultinom(1, size = round(sum(x)), prob = x / sum(x))
        } else {
          rand <- matrix(rep(0, length(x)))
        }
        rand
      }
      for (k in 1:length(TreeNames)) {
        freq.exp[, sel.cat[[k]]] <- t(apply(
          freq.exp[, sel.cat[[k]], drop = FALSE], 1,
          rmultinom_stable
        ))
      }
    }
    colnames(freq.exp) <- mptInfo$MPT$Category

    list(freq.exp)
  }

  # reorder by name!
  var.names <- apply(expand.grid("theta[", 1:S, ",", 1:N, "]"), 1, paste0, collapse = "")
  par.ind <- par.ind[, gsub(" ", "", var.names, fixed = TRUE)]
  if (nCPU > 1) {
    cl <- makeCluster(nCPU)
    clusterExport(cl, c(
      "S", "N", "mptInfo", "numItems", "expected", "expectedFreq",
      "tree", "TreeNames", "n.thetaFE", "sel.cat"
    ), envir = environment())
    # loop across replications:
    freq.list <- parApply(cl, cbind(par.thetaFE, par.ind), 1, getPostPred)
    stopCluster(cl)
  } else {
    freq.list <- apply(cbind(par.thetaFE, par.ind), 1, getPostPred)
  }

  # remove strange list structure:
  freq.list <- lapply(freq.list, function(xx) xx[[1]])

  if (M == 1) {
    freq.list[[1]]
  } else if (pred.new) {
    do.call("rbind", freq.list)[1:M, , drop = FALSE]
  } else {
    freq.list[1:min(M, length(freq.list))]
  }
}
