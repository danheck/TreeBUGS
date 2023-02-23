#' Prior Predictive Samples
#'
#' Samples full data sets (i.e., individual response frequencies) or group-level
#' MPT parameters based on prior distribution for group-level parameters.
#'
#' @inheritParams betaMPT
#' @param prior a named list defining the priors. For the \link{traitMPT}, the
#'   default is \code{list(mu = "dnorm(0,1)", xi="dunif(0,10)", V=diag(S),
#'   df=S+1)}, where S is the number of free parameters. For the \link{betaMPT},
#'   the default is \code{list(alpha ="dgamma(1,.1)", beta = "dgamma(1,.1)")}.
#'   Note that the normal distribution \code{"dnorm(mu,prec)"} is parameterized
#'   as in JAGS by the mean and precision (= 1/variance).
#' @param numItems vector with the number of items per MPT tree (either named or
#'   assigned to alphabetically ordered tree labels)
#' @param N number of participants per replication
#' @param level either \code{"data"} (returns individual frequencies) or
#'   \code{"parameter"} (returns group-level MPT parameters; \code{M} and
#'   \code{numItems} are ignored)
#' @param M number of prior predictive samples (i.e., data sets with \code{N}
#'   participants).
#' @param nCPU number of CPUs used for parallel sampling. For large models and
#'   many participants, this may require a lot of memory.
#' @return a list of \code{M} matrices with individual frequencies
#'   (rows=participants, columns=MPT categories). A single matrix is returned if
#'   \code{M=1} or \code{level="parameter"}.
#'
#' @examples
#' eqnfile <- system.file("MPTmodels/2htm.eqn",
#'   package = "TreeBUGS"
#' )
#' ### beta-MPT:
#' prior <- list(
#'   alpha = "dgamma(1,.1)",
#'   beta = "dgamma(1,.1)"
#' )
#'
#' ### prior-predictive frequencies:
#' priorPredictive(prior, eqnfile,
#'   restrictions = list("g=.5", "Do=Dn"),
#'   numItems = c(50, 50), N = 10, M = 1, nCPU = 1
#' )
#'
#' ### prior samples of group-level parameters:
#' priorPredictive(prior, eqnfile,
#'   level = "parameter",
#'   restrictions = list("g=.5", "Do=Dn"),
#'   M = 5, nCPU = 1
#' )
#'
#' ### latent-trait MPT
#' priorPredictive(
#'   prior = list(
#'     mu = "dnorm(0,1)", xi = "dunif(0,10)",
#'     df = 3, V = diag(2)
#'   ),
#'   eqnfile, restrictions = list("g=.5"),
#'   numItems = c(50, 50), N = 10, M = 1, nCPU = 1
#' )
#'
#' @importFrom parallel clusterExport makeCluster stopCluster parLapply parApply
#' @importFrom  stats cor cov2cor density rmultinom
#' @export
priorPredictive <- function(
    prior,
    eqnfile,
    restrictions,
    numItems,
    level = "data",
    N = 1,
    M = 100,
    nCPU = 4
) {
  if (missing(restrictions)) restrictions <- NULL
  # 1. get MPT model
  mpt <- readEQN(eqnfile, restrictions = restrictions)
  merged <- thetaHandling(mpt, restrictions)
  S <- max(merged$SubPar$theta)
  thetaNames <- merged$SubPar[, 1:2]
  thetaUnique <- thetaNames[rownames(unique(thetaNames[2])), ]$Parameter

  # 2. sample group-level parameters
  phi <- sampleHyperprior(prior = prior, M = M, S = S, nCPU = nCPU)
  if ("alpha" %in% names(prior)) {
    colnames(phi$alpha) <- colnames(phi$beta) <- thetaUnique
  } else {
    colnames(phi$mu) <- colnames(phi$sigma) <- thetaUnique
    dimnames(phi$rho) <- list(thetaUnique, thetaUnique, NULL)
    phi$mean <- phi$sd <- NULL
  }

  if (level == "data") {
    # 3. sample participants
    freq.list <- list()
    getData <- function(mm) {
      if (!"rho" %in% names(phi)) {
        freq <- genBetaMPT(N, numItems, eqnfile, restrictions,
          alpha = phi$alpha[mm, ], beta = phi$beta[mm, ],
          warning = FALSE
        )$data
      } else {
        freq <- genTraitMPT(N, numItems, eqnfile, restrictions,
          mean = pnorm(phi$mu[mm, ]), sigma = phi$sigma[mm, ],
          rho = phi$rho[, , mm], warning = FALSE
        )$data
      }
      freq
    }
    if (nCPU > 1) {
      cl <- makeCluster(nCPU)
      clusterExport(cl, c("S", "N", "numItems", "eqnfile", "restrictions"),
        envir = environment()
      )
      # loop across replications:
      freq.list <- parLapply(cl, 1:M, getData)
      stopCluster(cl)
    } else {
      freq.list <- lapply(1:M, getData)
    }

    # remove strange list structure:
    if (M == 1) {
      res <- freq.list[[1]]
    } else {
      res <- freq.list[1:min(M, length(freq.list))]
    }
    attr(res, "true") <- phi
  } else {
    res <- phi
  }
  res
}
