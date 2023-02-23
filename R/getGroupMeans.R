#' Get Mean Parameters per Group
#'
#' For hierarchical latent-trait MPT models with discrete predictor variables
#' as fitted with \code{traitMPT(..., predStructure = list("f"))}.
#'
#' @param traitMPT a fitted \code{\link{traitMPT}} model
#' @param factor whether to get group estimates for all combinations of factor
#'     levels (default) or only for specific factors (requires the names of the covariates in covData)
#' @param probit whether to use probit scale or probability scale
#' @param file filename to export results in .csv format (e.g., \code{file="fit_group.csv"})
#' @param mcmc if \code{TRUE}, the raw MCMC samples for the group means are returned
#'    as an \code{\link[coda]{mcmc.list}} object. This allows pairwise tests of group means
#'    (see \code{\link{transformedParameters}}).
#'
#' @examples
#' \dontrun{
#' # save group means (probability scale):
#' getGroupMeans(traitMPT, file = "groups.csv")
#' }
#' @seealso \code{\link{getParam}} for parameter estimates
#' @author Daniel Heck
#' @export
getGroupMeans <- function(traitMPT, factor = "all", probit = FALSE,
                          file = NULL, mcmc = FALSE) {
  if (is.null(traitMPT$mptInfo$predTable)) {
    stop("Model does not contain discrete predictors.")
  }

  uniqueNames <- traitMPT$mptInfo$thetaUnique
  facLevelNames <- traitMPT$mptInfo$predFactorLevels
  if (all(factor == "all")) {
    includeFactors <- names(facLevelNames)
  } else {
    if (any(!factor %in% names(facLevelNames))) {
      stop(
        "Factors defined by 'factor' are not factors in the fitted model!\n      Possible factors:  ",
        paste(names(facLevelNames)[!sapply(facLevelNames, is.null)], collapse = ", ")
      )
    }
    includeFactors <- factor
    facLevelNames <- facLevelNames[factor]
  }
  S <- length(uniqueNames)

  summaryMat <- matrix(NA, 0, 6)

  mcmc.list <- mcmc.list(lapply(
    traitMPT$runjags$mcmc,
    function(x) mcmc(matrix(NA, niter(x), 0))
  ))
  for (s in 1:S) {
    # parameter s of S
    select <- grep(paste0("mu", ifelse(S == 1, "", paste0("[", s, "]"))),
      varnames(traitMPT$runjags$mcmc),
      fixed = TRUE
    )
    muPosterior <- do.call("rbind", traitMPT$runjags$mcmc[, select])

    # select only relevant predictors
    predTable <- traitMPT$mptInfo$predTable
    predTable <- subset(
      predTable,
      predTable$prefix == "factor" &
        predTable$theta == s &
        predTable$Covariate %in% includeFactors
    )

    if (nrow(predTable) > 0) {
      factors <- unique(as.character(predTable$Covariate))
      combinations <- expand.grid(facLevelNames[factors])
      combNames <- sapply(combinations, as.character)

      # get posterior samples for differences due to factors
      facPosterior <- list()
      for (j in 1:length(factors)) {
        parnam <- as.character(subset(predTable, predTable$Covariate == factors[j])[, "covPar"])
        if (length(facLevelNames[[factors[j]]]) > 1) {
          grep.est <- setdiff(
            grep(parnam, varnames(traitMPT$runjags$mcmc)),
            grep("SD_", varnames(traitMPT$runjags$mcmc))
          )
        }
        facPosterior[[j]] <- do.call("rbind", traitMPT$runjags$mcmc[, grep.est])
        colnames(facPosterior[[j]]) <- facLevelNames[[factors[j]]]
      }

      # splitMat <- unique(traitMPT$mptInfo$covData[,factors, drop=FALSE])
      for (i in 1:nrow(combinations)) {
        samples <- muPosterior
        for (j in 1:ncol(combinations)) {
          samples <- samples + facPosterior[[j]][, combinations[i, j]]
        }
        if (mean(samples) > mean(muPosterior)) {
          pval <- mean(samples < muPosterior)
        } else {
          pval <- mean(samples > muPosterior)
        }
        if (!probit) {
          samples <- pnorm(samples)
        }

        label <- paste0(uniqueNames[[s]], "_", paste0(paste0(
          colnames(combinations),
          "[", combNames[i, ], "]"
        ), collapse = "_"))
        newMat <- rbind(c(
          Mean = mean(samples), SD = sd(samples),
          quantile(samples, c(.025, .5, .975)),
          "p(one-sided vs. overall)" = pval
        ))
        rownames(newMat) <- label
        summaryMat <- rbind(summaryMat, newMat)

        rownames(samples) <- rep(label, nrow(samples))
        for (c in 1:nchain(mcmc.list)) {
          mcmc.list[[c]] <- cbind(mcmc.list[[c]], t(samples)[, c, drop = FALSE])
        }
      }
    }
  }

  if (!is.null(file)) {
    write.csv(summaryMat, file = file)
  }

  mcmc.list <- mcmc.list(lapply(mcmc.list, mcmc))
  if (mcmc) {
    return(mcmc.list)
  } else {
    summaryMat
  }
}
