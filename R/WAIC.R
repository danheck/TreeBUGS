#' Widely Applicable Information Criterion
#'
#' Experimental implementation of the WAIC.
#'
#' @inheritParams posteriorPredictive
#' @param n.adapt number of adaptation samples.
#' @param n.chains number of chains (no parallel computation).
#' @param n.iter number of iterations after burnin.
#' @param n.thin Thinning rate.
# ' @param  ... further arguments passed to \link[rjags]{jags.samples}.
#'
#' @details
#' Please note that the implementation uses an experimental feature of the JAGS version
#' 4.3.0. It might be a good idea to fit each model and compute WAIC twice to
#' assess the stability of the WAIC values. See the following discussion for
#' details:
#'
#' \url{https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/}
#'
#' @return
#' a vector with the WAIC penalty \code{"p_waic"}, deviance \code{"deviance"}, and the WAIC \code{"waic"}.
#'
#' @examples
#' \dontrun{
#' fit <- traitMPT(...)
#' WAIC(fit)
#' }
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. Statistics and
#' Computing, 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4
#' @importFrom rjags jags.samples jags.model load.module
#' @export
WAIC <- function(fittedModel, n.adapt = 1000, n.chains = 3, n.iter = 10000, n.thin = 1){

  load.module("dic")

  # use last MCMC samples as initial values
  mcmc <- fittedModel$runjags$mcmc
  M <- nrow(mcmc[[1]])
  init <- mcmc[M,]
  cc <- length(init)
  initvec <- c(init, init[sample(cc, n.chains - cc, replace = TRUE)])[1:n.chains]
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
  # cf.: https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/
  samples <- jags.samples(mod, c("deviance", "WAIC"), type="mean",
                          n.iter = n.iter, thin = n.thin)
  tmp <- sapply(samples, sum)

  waic <- c("p_waic" = tmp[["WAIC"]],
            "deviance" = tmp[["deviance"]],
            "waic" = sum(tmp))
  waic
}

