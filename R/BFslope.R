#' Bayes Factor for Slope Parameters in Latent-Trait MPT
#'
#' Uses the Savage-Dickey method to compute the Bayes factor that the slope
#' parameter of a continuous covariate is zero vs. positive/negative/unequal to zero.
#'
#' @param fittedModel a fitted latent-trait model with predictor variables that have
#'     been defined via \code{predStructure}.
#' @param parameter name of the slope parameter (e.g., \code{"slope_d_covariate"})
#' @param direction alternative hypothesis: whether slope is smaller or larger than zero (\code{"<"} or \code{">"})
#'     or unequal to zero (\code{"!="})
#' @param plot if \code{TRUE}, the prior and posterior densities and the ratio at slope=0 are plotted.
#' @param ... further arguments passed to \code{\link[logspline]{logspline}}, which is used to
#'     approximate the density of the posterior distribution.
#'
#' @examples
#' \dontrun{
#' # fit beta-MPT model for encoding condition (see ?arnold2013):
#' EQNfile <- system.file("MPTmodels/2htsm.eqn", package="TreeBUGS")
#'
#' d.encoding <- subset(arnold2013, group == "encoding", select = -(1:4))
#' d.cov <- subset(arnold2013, group == "encoding", select = "age")
#' fit <- traitMPT(EQNfile, d.encoding, n.thin=5,
#'                 restrictions=list("D1=D2=D3","d1=d2","a=g"),
#'                 covData = d.cov, predStructure = list("D1 ; age"))
#' # convergence
#' plot(fit, parameter = "slope", type = "default")
#' summary(fit)
#' BayesFactorSlope(fit, "slope_D1_age", direction = "<")
#' }
#' @export
#' @importFrom stats dcauchy
#' @importFrom logspline dlogspline logspline
######## Bayes factors:
## H0: slope beta = 0
## H1: slope beta < 0   (i.e., beta ~ Cauchy)
BayesFactorSlope <- function (fittedModel, parameter,
                              direction = "<", plot = TRUE, ...){
  if (any(fittedModel$mptInfo$hyperprior$IVprec != "dchisq(1)"))
    stop("BayesFactorSlope requires that default priors are used for the slope parameter!")

  if (length(parameter)!= 1 || !parameter %in% rownames(fittedModel$mcmc.summ))
    stop("'parameter' not in model or not of length=1.")
  tmp <- strsplit(parameter, "_")
  tmp[[1]][3] <- paste(tmp[[1]][-c(1:2)], collapse = "_")
  tmp <- tmp[[1]][1:3]
  if (length(tmp) != 3)
    stop("Parameter must be of the form 'slope_MPTparam_cov' .")
  if (tmp[[1]][1] != "slope")
    stop("Only valid for slope parameters.")
  cov <- tmp[3]
  s <- apply(fittedModel$mptInfo$covData, 2, sd)[cov]
  samples <- unlist(fittedModel$runjags$mcmc[,parameter]) * s   # slope parameters are not standardized wrt covariate!

  # approximation of posterior density
  if (direction == "<"){
    sel <- samples < 0
    xlim <- c(min(samples[sel], -.05), 0)
    posterior <- logspline(samples[sel], ubound = 0, ...)
  } else if (direction == ">"){
    sel <- samples > 0
    xlim <- c( 0, max(samples[sel], .05))
    # knots <- seq(0, 3, .5)
    posterior <- logspline(samples[sel], lbound = 0, ...)
  } else if (direction == "!="){
    sel <- rep(TRUE, length(samples))
    xlim <- c(min(samples[sel], -.05), max(samples[sel], .05))
    posterior <- logspline(samples[sel], ...)
  } else {
    stop("'direction' must be '>', '<', or '!=' ")
  }

  # posterior and prior density for beta=0:
  post0 <- dlogspline(0, posterior)
  prior0 <- dcauchy(0) * ifelse(direction == "!=", 1, 2)  # one-sided

  # BF in favor of effect:
  bf <- data.frame(post0/prior0, prior0/post0)
  colnames(bf) <- paste0("BF_", c(0, direction),  c(direction, 0))

  # illustration of Savage-Dickey method:
  dcauchy_trunc <- function(x){
    if (direction == ">"){
      dx <- 2*dcauchy(x)*ifelse(x > 0, 1, 0)
    } else if (direction == "<"){
      dx <- 2*dcauchy(x)*ifelse(x < 0, 1, 0)
    } else if (direction == "!="){
      dx <- dcauchy(x)
    }
    dx
  }
  if (plot){
    hist(samples[sel], col = adjustcolor("gray", alpha.f =.3), 100,
         freq = FALSE, xlim = xlim,
         main = "Bayes factor: Prior (red) vs. posterior (blue)",
         xlab = paste0("Standardized slope parameter: ", tmp[2], " ~ ", tmp[3]))
    plot(posterior, add = TRUE, col = 4, lwd = 2, n = 1000)
    curve(dcauchy_trunc, col=2, add=TRUE, n= 1001, lwd = 2)
    points(c(0,0), c(post0, prior0), col=c(4,2), pch=16, lwd = 2, cex = 2)
  }

  bf
}


