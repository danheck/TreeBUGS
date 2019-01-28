#' WAIC: Widely Applicable Information Criterion
#'
#' Implementation of the WAIC for model comparison.
#'
#' @inheritParams posteriorPredictive
#' @param n.adapt number of adaptation samples.
#' @param n.chains number of chains (no parallel computation).
#' @param n.iter number of iterations after burnin.
#' @param n.thin Thinning rate.
#' @param summarize whether to summarize the WAIC values which are computed separately
#'   for each observed node.
#'
#' @details
#' WAIC provides an approximation of predictive accuracy with respect to out-of-sample deviance.
#' The uncertainty of the WAIC for the given number of observed nodes
#' (i.e., number of free categories times the number of participants)
#' is quantified by the standard error of WAIC \code{"se_waic"} (cf. Vehtari et al., 2017).
#' In contrast, to assess whether the approximation uncertainty due to MCMC sampling
#' (not sample size) is sufficiently low, it is a good idea to fit each model twice and compute WAIC
#' again to assess the stability of the WAIC values.
#'
#' For more details, see Vehtari et al. (2017) and the following discussion about
#' the JAGS implementation (which is currently an experimental feature of JAGS 4.3.0):
#'
#' \url{https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/}
#'
#' @return
#' If \code{summarize=TRUE} (default), a vector containing the WAIC penalty term
#' \code{"p_waic"}, \code{"deviance"}, \code{"waic"}, and the corresponding
#' standard error \code{"se_waic"}.
#'
#' If \code{summarize=FALSE}, a list containing three vectors \code{p_waic},
#' \code{deviance}, and \code{waic} with separate values for each observed node
#'  (i.e., for all combinations of persons and free categories).
#'
#' @examples
#' \dontrun{
#'
#' #### WAIC for a latent-trait MPT model:
#' fit <- traitMPT(...)
#' WAIC(fit)
#'
#'
#' #### pairwise comparison of two models:
#' # (1) compute WAIC per model
#' waic1 <- WAIC(fit1, summarize = FALSE)
#' waic2 <- WAIC(fit2, summarize = FALSE)
#' # (2) WAIC differences (per observation!)
#' waic_diff <- waic1$waic - waic2$waic
#' # (3) standard error of the WAIC differences:
#' n_obs <- length(waic_diff)
#' c(diff = sum(waic_diff),
#'   se_diff = sqrt(n_obs) * sd(waic_diff))
#' }
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. Statistics and
#' Computing, 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4
#' @importFrom rjags jags.samples jags.model load.module
#' @export
WAIC <- function(fittedModel, n.adapt = 1000, n.chains = 3, n.iter = 10000, n.thin = 1,
                 summarize = TRUE){

  stopifnot(class(fittedModel) %in% c("traitMPT", "betaMPT"))
  load.module("dic")

  # use last MCMC samples as initial values
  mcmc <- fittedModel$runjags$mcmc
  M <- nrow(mcmc[[1]])
  init <- mcmc[M,]
  cc <- length(init)
  initvec <- c(init, init[sample(cc, max(0, n.chains - cc), replace = TRUE)])[1:n.chains]
  inits <- list()
  for(i in seq_along(initvec)){
    inits[[i]] <- list("mu" = initvec[[i]][grep("mu", names(initvec[[i]]))])
  }

  # extract data list from fitted runjags model:
  dat <- strsplit(fittedModel$runjags$data, "\\n")[[1]]
  datlist <- lapply(dat, function(x) eval(parse(text = x)))
  names(datlist) <- sapply(dat, function(x) sub(" .*", "", gsub("\\\"", "", x)))

  # construct new JAGS model
  mod <- jags.model(textConnection(fittedModel$runjags$model),
                    datlist[-length(datlist)], n.chains = n.chains,
                    n.adapt=n.adapt, inits = inits)

  # the WAIC feature is still experimental! (requires JAGS 4.3.0)
  # returns deviance + waic penalty for each observed node (= free category per person)
  # cf.: https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/
  samples <- jags.samples(mod, c("deviance", "WAIC"), type= "mean",
                          n.iter = n.iter, thin = n.thin)

  samples$p_waic <- samples$WAIC
  samples$waic <- samples$deviance + samples$p_waic
  attributes(samples$waic) <- attributes(samples$p_waic) <- attributes(samples$deviance)
  samples$WAIC <- NULL

  if (!summarize)
    return(samples)

  # standard error for WAIC (cf. Vehtari et al., 2017)
  n_obs <- length(samples$waic)  # = number of free categories times nubmer of persons
  se_waic <- sqrt(n_obs) * sd(samples$waic)

  tmp <- sapply(samples, sum)
  waic <- c("p_waic" = tmp[["p_waic"]],
            "deviance" = tmp[["deviance"]],
            "waic" = tmp[["waic"]],
            "se_waic" = se_waic)
  return(waic)
}

