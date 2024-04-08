#' Generate Data for Beta MPT Models
#'
#' Generating a data file with known parameter structure using the Beta-MPT.
#' Useful for simulations and robustness checks.
#'
#' @inheritParams betaMPT
#' @inheritParams genMPT
#' @param N number of participants
#' @param mean Named vector of true group means of individual MPT parameters. If
#'   the vector is not named, the internal order of parameters is used (can be
#'   obtained using \code{\link{readEQN}}).
#' @param sd named vector of group standard deviations of individual MPT
#'   parameters.
#' @param alpha Alternative specification of the group-level distribution using
#'   the shape parameters of the beta distribution (see \link{dbeta}).
#' @param beta see \code{alpha}
#'
#' @details Data are generated in a two-step procedure. First, person parameters
#' are sampled from the specified beta distributions for each paramter (either
#' based on mean/sd or based on alpha/beta). In a second step, response
#' frequencies are sampled for each person using \code{\link{genMPT}}.
#'
#' @return a list including the generated frequencies (\code{data}) and the
#'   true, underlying parameters (\code{parameters}) on the group and individual
#'   level.
#' @seealso \code{\link{genMPT}}
#' @references Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial
#'   processing tree models for addressing individual differences. Journal of
#'   Mathematical Psychology, 54, 167-183.
#'
#' @examples
#' # Example: Standard Two-High-Threshold Model (2HTM)
#' EQNfile <- system.file("MPTmodels/2htm.eqn", package = "TreeBUGS")
#' genDat <- genBetaMPT(
#'   N = 100,
#'   numItems = c(Target = 250, Lure = 250),
#'   eqnfile = EQNfile,
#'   mean = c(Do = .7, Dn = .5, g = .5),
#'   sd = c(Do = .1, Dn = .1, g = .05)
#' )
#' head(genDat$data, 3)
#' plotFreq(genDat$data, eqn = EQNfile)
#' @importFrom stats  rbeta
#' @export
genBetaMPT <- function(
    N,
    numItems,
    eqnfile,
    restrictions,
    mean = NULL,
    sd = NULL,
    alpha = NULL,
    beta = NULL,
    warning = TRUE
) {
  if (missing(restrictions)) {
    restrictions <- NULL
  }

  # REMINDER 2024-04: restrictions were ineffective in readEQN(); only processed by thetaHandling
  Tree <- readEQN(eqnfile)
  mergedTree <- mergeBranches(Tree)
  Tree.restr <- thetaHandling(mergedTree, restrictions)
  thetaNames <- Tree.restr$SubPar[, 1:2]
  thetaNames <- thetaNames[rownames(unique(thetaNames[2])), ]$Parameter
  treeLabels <- unique(mergedTree$Tree)
  S <- length(thetaNames)


  if (!is.null(mean) && !is.null(sd)) {
    mean <- checkNaming(S, thetaNames, mean, "mean",
      interval = c(0, 1), warning = warning
    )
    sd <- checkNaming(S, thetaNames, sd, "sd",
      interval = c(0, Inf), warning = warning
    )

    alpha <- ((1 - mean) / sd^2 - 1 / mean) * mean^2
    beta <- alpha * (1 / mean - 1)
    if (any(alpha <= 0) | any(beta <= 0)) {
      stop(
        "Check numerical values for mean and sd, result in negative alpha/beta\n",
        "parameters of beta-hyperprior distribution."
      )
    }
  } else if (!is.null(alpha) && !is.null(beta)) {
    alpha <- checkNaming(S, thetaNames, alpha, "alpha",
      interval = c(0, Inf), warning = warning
    )
    beta <- checkNaming(S, thetaNames, beta, "beta",
      interval = c(0, Inf), warning = warning
    )
  } else {
    stop("Either 'mean'/'sd' or 'alpha'/'beta' must be provided.")
  }

  # individual parameters, drawn from hierarchical distribution:
  theta <- c()
  for (s in 1:S) {
    if (!is.null(sd) && sd[s] == 0) {
      theta <- cbind(theta, rep(mean[s], N))
    } else {
      theta <- cbind(theta, rbeta(N, shape1 = alpha[s], shape2 = beta[s]))
    }
  }
  colnames(theta) <- thetaNames

  # response frequencies:
  freq <- genMPT(theta, numItems, eqnfile, restrictions, warning = warning)
  list(
    data = freq,
    parameters = list(theta = theta, mean = mean, sd = sd, alpha = alpha, beta = beta)
  )
}
