printSummary <- function(
    x,
    model,
    ...
) {
  cat("Call: \n")
  print(x$call)
  cat("\n")

  if (model == "betaMPT") {
    cat("Group-level means of MPT parameters:\n")
  } else if (model == "traitMPT") {
    cat("Group-level medians of MPT parameters (probability scale):\n")
  }
  print(round(x$groupParameters$mean, x$round))

  if (model == "betaMPT") {
    cat("\nStandard deviation of parameters across individuals:\n")
    print(round(x$groupParameters$SD, x$round))
    cat("\nAlpha parameters of beta distributions:\n")
    print(round(x$groupParameters$alpha, x$round))
    cat("\nBeta parameters of beta distributions:\n")
    print(round(x$groupParameters$beta, x$round))
  } else if (model == "traitMPT") {
    cat("\nMean/Median of latent-trait values (probit-scale) across individuals:\n")
    print(round(x$groupParameters$mu, x$round))
    cat("\nStandard deviation of latent-trait values (probit scale) across individuals:\n")
    print(round(x$groupParameters$sigma, x$round))
  } else if (model == "simpleMPT" && !all(is.na(x$groupParameters$SD[, 1]))) {
    cat("\nStandard deviation of (fixed-effects) MPT parameters:\n")
    print(round(x$groupParameters$SD, x$round))
  }


  if (!is.null(x$groupParameters$thetaFE)) {
    cat("\nFixed effects MPT parameters (= identical for all subjects):\n")
    print(round(x$groupParameters$thetaFE, x$round))
  }

  if (model == "traitMPT") {
    if (any(abs(x$groupParameters$rho.matrix) > 1)) {
      cat("\n[note: latent-trait version without explicit correlation parameters]\n")
    } else if (nrow(x$groupParameters$rho.matrix) != 1) {
      cat("\nCorrelations of latent-trait values on probit scale:\n")
      print(round(x$groupParameters$rho, x$round))
      cat("\nCorrelations (posterior mean estimates) in matrix form:\n")
      print(round(x$groupParameters$rho.matrix, x$round))
    }
  }

  cat(
    "\n\n##############\n",
    "Model fit statistics (posterior predictive p-values):\n"
  )
  if (!is.null(x$fitStatistics$overall)) {
    print(round(x$fitStatistics$overall, x$round))
    cat("\nT1 per person:\n")
    print(round(x$fitStatistics$individual$T1.p, x$round))
  } else {
    cat("Use PPP(fittedModel) to get T1 and T2 posterior predictive checks.\n")
  }
  if (!is.null(x$dic)) {
    print(x$dic)
    cat("\n")
  }

  if (!is.null(x$transformedParameters)) {
    cat("\nTransformed parameters:\n")
    print(round(x$transformedParameters, x$round))
  }

  if (!is.null(x$groupParameters$slope)) {
    cat("\nSlope parameters for predictor variables:\n")
    print(round(x$groupParameters$slope, x$round))
  }

  if (!is.null(x$groupParameters$factor)) {
    cat("\nEffects of factors on latent scale (additive shift from overall mean):\n")
    print(round(x$groupParameters$factor, x$round))
    cat("\nFactor SD on latent scale:\n")
    print(round(x$groupParameters$factorSD, x$round))
  }

  if (!is.null(x$groupParameters$correlation) && !nrow(x$groupParameters$correlation) == 0) {
    cat(
      "\nSampled correlations of MPT parameters with covariates:\n",
      "\n (only quantifies the uncertainty with respect to parameter estimation",
      "\n  not with respect to sample size! See ?correlationPosterior):\n"
    )
    print(round(x$groupParameters$correlation, x$round))
  }
}



#' @export
print.summary.betaMPT <- function(x, ...) {
  printSummary(x, "betaMPT")
}

#' @export
print.summary.traitMPT <- function(x, ...) {
  printSummary(x, "traitMPT")
}

#' @export
print.summary.simpleMPT <- function(x, ...) {
  printSummary(x, "simpleMPT")
}

#' @export
summary.betaMPT <- function(object, round = 3, ...) {
  summ <- object$summary
  summ$call <- object$call
  summ$round <- round
  return(summ)
}

#' @export
summary.traitMPT <- function(object, round = 3, ...) {
  summ <- object$summary
  summ$call <- object$call
  summ$round <- round
  return(summ)
}

#' @export
summary.simpleMPT <- function(object, round = 3, ...) {
  summ <- object$summary
  summ$call <- object$call
  summ$round <- round
  return(summ)
}

#' @export
print.betaMPT <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("\n")
  print(round(cbind(
    "Group Mean" = x$summary$groupParameters$mean[, 1],
    "Group SD" = x$summary$groupParameters$SD[, 1]
  ), 4))

  cat("\nUse 'summary(fittedModel)' or 'plot(fittedModel)' to get a more detailed summary.")
}

#' @export
print.traitMPT <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("\n")
  print(round(cbind(
    "Group median (probability scale)" = x$summary$groupParameters$mean[, 1],
    "Group mean/median (latent-probit)" = x$summary$groupParameters$mu[, 1],
    "Group SD (latent-probit)" = x$summary$groupParameters$sigma[, 1]
  ), 4))

  cat("\nUse 'summary(fittedModel)' or 'plot(fittedModel)' to get a more detailed summary.")
}

#' @export
print.simpleMPT <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("\n")
  print(round(cbind(
    "Mean(MPT Parameters)" = x$summary$groupParameters$mean[, 1],
    "SD(MPT parameters)" = x$summary$groupParameters$SD[, 1]
  ), 4))

  cat("\nUse 'summary(fittedModel)' or 'plot(fittedModel)' to get a more detailed summary.")
}
