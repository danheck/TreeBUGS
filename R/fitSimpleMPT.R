#' C++ Sampler for Standard (Nonhierarchical) MPT Models
#'
#' Fast Gibbs sampler in C++ that is tailored to the standard fixed-effects MPT
#' model (i.e., fixed-effects, non-hierarchical MPT). Assumes independent
#' parameters per person if a matrix of frequencies per person is supplied.
#'
#' @inheritParams betaMPT
#' @inheritParams betaMPTcpp
#' @param alpha first shape parameter(s) for the beta prior-distribution of the
#'   MPT parameters \eqn{\theta_s} (can be a named vector to use a different
#'   prior for each MPT parameter)
#' @param beta second shape parameter(s)
#'
#' @details Beta distributions with fixed shape parameters \eqn{\alpha} and
#'   \eqn{\beta} are used. The default \eqn{\alpha=1} and \eqn{\beta=1} assumes
#'   uniform priors for all MPT parameters.
#' @author Daniel Heck
#'
#' @examples
#' \dontrun{
#' # fit nonhierarchical MPT model for aggregated data (see ?arnold2013):
#' EQNfile <- system.file("MPTmodels/2htsm.eqn", package = "TreeBUGS")
#' d.encoding <- subset(arnold2013, group == "encoding", select = -(1:4))
#' fit <- simpleMPT(EQNfile, colSums(d.encoding),
#'   restrictions = list("D1=D2=D3", "d1=d2", "a=g")
#' )
#' # convergence
#' plot(fit)
#' summary(fit)
#' }
#' @importFrom parallel parLapply stopCluster detectCores
#' @export
simpleMPT <- function(
    eqnfile,
    data,
    restrictions,
    n.iter = 2000,
    n.burnin = 500,
    n.thin = 3,
    n.chains = 3,
    ppp = 0,
    alpha = 1,
    beta = 1,
    parEstFile,
    posteriorFile,
    cores = 1
) {
  hyperprior <- list(alpha = alpha, beta = beta)
  if (!is.character(data) && is.null(dim(data))) {
    data <- matrix(data, nrow = 1, dimnames = list(NULL, names(data)))
  }

  fittedModel <- fitModelCpp("simpleMPT",
                             eqnfile = eqnfile,
    data = data, restrictions = restrictions,
    hyperprior = hyperprior,
    n.iter = n.iter,
    n.burnin = n.burnin, n.thin = n.thin,
    n.chains = n.chains, ppp = ppp,
    parEstFile = parEstFile,
    posteriorFile = posteriorFile,
    call = match.call(),
    cores = cores
  )

  fittedModel
}
