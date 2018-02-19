#' Bayes Factor for Slope Parameters in Latent-Trait MPT
#'
#' Uses the Savage-Dickey method to compute the Bayes factor that the slope
#' parameter of a continuous covariate in \code{\link{traitMPT}} is
#' zero vs. positive/negative/unequal to zero.
#'
#' @param fittedModel a fitted latent-trait model fitted with \code{\link{traitMPT}}
#'     with predictor variables that have been defined via \code{predStructure}.
#' @param parameter name of the slope parameter (e.g., \code{"slope_d_covariate"}).
#' @param direction alternative hypothesis: whether slope is smaller or larger
#'     than zero (\code{"<"} or \code{">"}) or unequal to zero (\code{"!="}).
#' @param approx how to approximate the posterior density of the slope parameter at zero:
#'     \code{approx="normal"} uses a normal approximation to all samples and \code{approx="logspline"}
#'     uses a nonparametric density estimate of the package \link[logspline]{logspline}.
#'     Usually, both methods provide similar results.
#' @param plot if \code{TRUE}, the prior and posterior densities and the ratio at slope=0 are plotted.
#' @param ... further arguments passed to \code{\link[logspline]{logspline}}, which is used to
#'     approximate the density of the posterior distribution.
#'
#' @details
#' The Bayes factor is computed with the Savage-Dickey method, which is defined as the ratio of the
#' density of the posterior and the density of the prior evalauted at \code{slope=0}.
#' Note that this method cannot be used with default JZS priors (\code{IVprec="dgamma(.5,.5)"})
#' if more than one predictor is added for an MPT parameter. As a remedy, a g-prior (normal distribution)
#' can be used on the slopes by setting the hyperprior parameter \eqn{g} to a
#' fixed constant when fitting the model: \code{traitMPT(..., IVprec = 1)}.
#'
#' @examples
#' \dontrun{
#' # latent-trait MPT model for the encoding condition (see ?arnold2013):
#' EQNfile <- system.file("MPTmodels/2htsm.eqn", package="TreeBUGS")
#' d.enc <- subset(arnold2013, group == "encoding")
#'
#' fit <- traitMPT(EQNfile, data = d.enc[, -(1:4)], n.thin=5,
#'                 restrictions=list("D1=D2=D3","d1=d2","a=g"),
#'                 covData = d.enc[,c("age", "pc")],
#'                 predStructure = list("D1 ; age"))
#' plot(fit, parameter = "slope", type = "default")
#' summary(fit)
#'
#' BayesFactorSlope(fit, "slope_D1_age", direction = "<")
#' }
#' @export
#' @importFrom stats dcauchy
#' @importFrom logspline dlogspline logspline
######## Bayes factors:
## H0: slope beta = 0
## H1: slope beta < 0   (i.e., beta ~ Cauchy)
BayesFactorSlope <- function (fittedModel, parameter,
                              direction = "!=", approx = "normal",
                              plot = TRUE, ...){
  approx <- match.arg(approx, c("normal", "logspline"))

  # hyperprior for "g" (numeric: g-prior; dgamma: JZS)
  IVprec <- fittedModel$mptInfo$hyperprior$IVprec
  if (length(IVprec) != 1)
    stop("Fitted model must use the same 'IVprec' for all slope parameters.")
  if (is.numeric(IVprec)){
    IVfamily <- "constant"
  } else {
    IVsplit <- strsplit(IVprec, "[\\(,\\)]")[[1]]
    IVfamily <- IVsplit[[1]]
    IVpars <- as.numeric(sapply(IVsplit[- 1],
                                function(g) eval(as.expression(g))))
  }

  parlabels <- rownames(fittedModel$mcmc.summ)
  if (length(parameter) != 1 || !parameter %in% parlabels)
    stop("'parameter' not in model or not of length=1.")

  if (substr(parameter, 1, 6) != "slope_")
    stop("Only valid for slope parameters.")

  thetas <- fittedModel$mptInfo$thetaUnique
  theta <- thetas[sapply(paste0("_", thetas, "_"), grepl, parameter)]
  if (length(theta) != 1)
    stop("Check that slope parameter includes the correct label of the MPT parameters:\n",
         "  theta:     ", paste(thetas, collapse = ", "), "  \n  parameter: ", parameter)
  cov <- substr(parameter, 7 + nchar(theta) + 1, 999)
  if(!cov %in% colnames(fittedModel$mptInfo$covData))
    stop("Covariate", cov, "not in covariate data.")

  if (IVfamily == "dgamma" && sum(grepl(paste0("slope_", theta), parlabels)) > 1)
    stop("The Savage-Dickey method provides incorrect Bayes factors if:\n",
         "     (A) one of the MPT parameters has more than one predictors AND \n",
         "     (B) the slope parameters have the default prior (known as JZS or Cauchy prior).\n",
         "     (cf. Heck (2018): Computing Bayes factors for regression models: \n",
         "                       A caveat on the Savage-Dickey density ratio)\n\n",
         "  As a remedy, refit the model \n",
         "      (1) with a maximum of one predictor per MPT parameter OR\n",
         "      (2) with standard-normal priors (g-prior) on the regression slopes: traitMPT(..., IVprec = 1)\n")

  # slope parameters are not standardized wrt covariate! => standardization
  s <- apply(fittedModel$mptInfo$covData, 2, sd)[cov]
  samples <- unlist(fittedModel$runjags$mcmc[,parameter]) * s

  # approximation of posterior density
  lbnd <- switch(direction, "<" = -Inf, ">" = 0, "!=" = -Inf, NA)
  ubnd <- switch(direction, "<" = 0, ">" = Inf, "!=" = Inf, NA)
  if (is.na(lbnd)) stop("'direction' must be one of: '>', '<', or '!=' ")
  sel <- samples > lbnd & samples < ubnd

  xlim <- switch(direction,
                 "<" = c(min(samples[sel], -.05), 0),
                 ">" = c( 0, max(samples[sel], .05)),
                 "!="= c(min(samples[sel], -.05), max(samples[sel], .05)))

  if (sum(sel) == 0){
    warning("No posterior samples are in line with the predicted direction of the effect!",
            "\n   Bayes factor can only be approximated with approx='normal'.")
    approx <- "normal"
  } else if (sum(sel) < 1000){
    warning("Less than 1000 posterior samples are in line with the predicted direction of the effect!",
            "\n    This might result in imprecise estimates of the Bayes factor")
  }

  # posterior and prior density for beta=0:
  post0 <- NA
  if (approx == "logspline"){
    try({
      posterior <- switch(direction,
                          "<" = logspline(samples[sel], ubound = 0, ...),
                          ">" = logspline(samples[sel], lbound = 0, ...),
                          "!=" = logspline(samples[sel], ...))
      post0 <- dlogspline(0, posterior)
    })
  } else if(approx == "normal"){
    mm <- mean(samples)
    ss <- sd(samples)
    posterior <- function(x)
      ifelse(x <= ubnd & x >= lbnd, 1, 0) * dnorm(x, mm, ss) /
      (pnorm(ubnd, mm, ss) - pnorm(lbnd, mm, ss))
    post0 <- posterior(0)
  }

  if (IVfamily == "dgamma"){
    if (IVpars[[1]] != .5)
      stop("First parameter in 'IVprec' must be equal to .5, that is: 'IVprec=dgamma(.5,.5*s^2)'.")
    scale <- sqrt(IVpars[[2]]*2)
    prior0 <- dcauchy(0, 0, scale) * ifelse(direction == "!=", 1, 2)  # one-sided

    # illustration of Savage-Dickey method:
    dprior <- function(x){
      if (direction == ">"){
        dx <- dcauchy(x, 0, scale)*ifelse(x > 0, 2, 0)
      } else if (direction == "<"){
        dx <- dcauchy(x, 0, scale)*ifelse(x < 0, 2, 0)
      } else if (direction == "!="){
        dx <- dcauchy(x, 0, scale)
      }
      dx
    }

  } else if (IVfamily == "constant") {
    prior0 <- dnorm(0, 0, 1/sqrt(IVprec)) * ifelse(direction == "!=", 1, 2)  # one-sided

    dprior <- function(x){
      if (direction == ">"){
        dx <- dnorm(x, 0, 1/sqrt(IVprec))*ifelse(x > 0, 2, 0)
      } else if (direction == "<"){
        dx <- dnorm(x, 0, 1/sqrt(IVprec))*ifelse(x < 0, 2, 0)
      } else if (direction == "!="){
        dx <- dnorm(x, 0, 1/sqrt(IVprec))
      }
      dx
    }
  } else {
    stop("Either a JZS prior 'IVprec=dgamma(.5,.5*s^2)' or g-prior 'IVprec=g' must be used.\n
         (where s and g are numeric constants).")
  }

  # BF in favor of effect:
  bf <- data.frame(post0/prior0, prior0/post0)
  colnames(bf) <- paste0("BF_", c(0, direction),  c(direction, 0))


  if (plot){
    hist(samples[sel], col = adjustcolor("gray", alpha.f =.3), 70,
         freq = FALSE, xlim = xlim,
         main = paste0("Bayes factor B_10=", round(bf[1,2], 3),
                       " (prior red; posterior blue)"),
         xlab = paste0("Standardized slope parameter: ", theta, " ~ ", cov))
    if(is.function(posterior))
      curve(posterior, add = TRUE, col = 4, lwd = 2, n = 1000)
    else
      try(plot(posterior, add = TRUE, col = 4, lwd = 2, n = 1000))
    curve(dprior, col=2, add=TRUE, n= 1001, lwd = 2)
    points(c(0,0), c(post0, prior0), col=c(4,2), pch=16, lwd = 2, cex = 2)
    # text(xlim[1]+.1, 3, col = 2,
    #      label = paste("Prior:\n", ifelse(IVfamily == "dgamma", "JZS (Cauchy)", "g-prior (normal)")))
  }

  bf
}


