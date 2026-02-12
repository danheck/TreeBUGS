#' Plot Posterior Predictive Mean Frequencies
#'
#' Plots observed means/covariances of individual frequencies against the
#' means/covariances sampled from the posterior distribution (posterior
#' predictive distribution).
#'
#' @inheritParams posteriorPredictive
#' @param stat whether to plot mean frequencies (\code{"mean"}), covariances
#'   of frequencies (\code{"cov"}), standard deviations (\code{"sd"}), or
#'   correlations (\code{"cor"})
#' @param main main title for plot
#' @param ylab label for y-axis
#' @param col color for boxplots of predicted values
#' @param ... further arguments passed to \code{\link{boxplot}}
#'
#' @details If posterior predictive p-values were computed when fitting the
#'   model (e.g., by adding the argument \code{traitMPT(...,ppp=1000)} ), the
#'   stored posterior samples are re-used for plotting. Note that the last
#'   category in each MPT tree is dropped, because one category per multinomial
#'   distribution is fixed.
#'
#' @examples
#' \dontrun{
#' # add posterior predictive samples to fitted model (optional step)
#' fittedModel$postpred$freq.pred <-
#'   posteriorPredictive(fittedModel, M = 1000)
#'
#' # plot model fit
#' plotFit(fittedModel, stat = "mean")
#' }
#' @export
plotFit <- function(
    fittedModel,
    M = 1000,
    stat = c("mean", "cov", "cor", "sd"),
    main = NULL,
    ylab = NULL,
    col = "gray",
    ...
) {
  stat <- match.arg(stat)

  # get information about model:
  tree <- fittedModel$mptInfo$MPT$Tree
  cats <- fittedModel$mptInfo$MPT$Category
  dat <- fittedModel$mptInfo$dat[, cats]
  TreeNames <- unique(tree)

  # free categories (drop last category per tree):
  free_cats <- unlist(tapply(
    X = cats, INDEX = tree,
    FUN = function(cat) cat[-length(cat)]
  ))

  # get posterior predictive:
  if (is.null(fittedModel$postpred) | M != 1000) {
    freq.list <- posteriorPredictive(fittedModel, M = M)
  } else {
    freq.list <- fittedModel$postpred$freq.pred
  }

  if (stat == "mean") {
    # Plot mean frequencies:
    predicted <- t(sapply(freq.list, colMeans))[, free_cats, drop = FALSE]
    observed  <- colMeans(dat[, free_cats, drop = FALSE])
    if(is.null(ylab)) ylab <- "Mean frequency"
    if(is.null(main)) main <- "Observed (red) and predicted (gray) mean frequencies"
    x_labels <- free_cats


  } else if (stat == "cov") {

    nams <- outer(free_cats, free_cats, paste, sep = "-")
    sel_cov <- nams[lower.tri(nams, diag = TRUE)]


    # observed/predicted
    c.obs <- cov(dat[, free_cats])
    c.pred <- sapply(freq.list, function(xx) {
      cc <- cov(xx[, free_cats])
      cc[lower.tri(cc, diag = TRUE)]
    })

    predicted <- t(c.pred)
    observed  <- c.obs[lower.tri(c.obs, diag = TRUE)]
    if(is.null(ylab)) ylab <- "Covariance"
    if(is.null(main)) main <- "Observed (red) and predicted (gray) covariances"
    x_labels <- sel_cov
  } else if (stat == "cor") {
    nams <- outer(free_cats, free_cats, paste, sep = "-")
    sel_cov <- nams[lower.tri(nams, diag = FALSE)]

    c.obs <- cor(dat[, free_cats])
    c.pred <- sapply(freq.list, function(xx) {
      cc <- cor(xx[, free_cats])
      cc[lower.tri(cc, diag = FALSE)]
    })
    predicted <- t(c.pred)
    observed  <- c.obs[lower.tri(c.obs, diag = FALSE)]
    if(is.null(ylab)) ylab <- "Correlation"
    if(is.null(main)) main <- "Observed (red) and predicted (gray) correlations"
    x_labels <- sel_cov

  } else if (stat == "sd") {
    nams <- free_cats
    s.obs <- sqrt(diag(cov(dat[, free_cats])))
    s.pred <- sapply(freq.list, function(xx) {
      sqrt(diag(cov(xx[, free_cats])))
    })
    predicted <- t(s.pred)
    observed <- s.obs
    if(is.null(ylab)) ylab <- "Standard deviation"
    if(is.null(main)) main <- "Observed (red) and predicted (gray) standard deviations"
    x_labels <- nams
  }


  out <- list()

  out$boxplot <- boxplot(
    x = predicted
    , col = col
    , ylab = ylab
    , main = main
    , xaxt = "n"
    , ...
  )
  out$points <- points(x = seq_along(observed), y = observed, col = 2, pch = 17)
  out$axis   <- axis(side = 1, seq_along(x_labels), labels = x_labels, las = 2)
  abline(h = 0, lty = 1, col = "gray")




  if(stat %in% c("mean", "sd")) {
    # We add tree names if means or SDs are plotted
    xx <- by(seq_along(free_cats), tree[cats %in% free_cats], mean)
    axis(1, xx, TreeNames, tick = F, line = NA, mgp = c(3, 2.5, 0))
    abline(v = cumsum(table(tree) - 1)[1:(length(TreeNames) - 1)] + .5, col = "gray")
  } else {
    # If co-variances or correlations are plotted, we add separators between blocks of columns for each tree.
    p <- nrow(c.obs)

    if (stat == "cov") {
      block_sizes <- p:1
    } else if (stat == "cor") {
      block_sizes <- (p - 1):1
    }

    # positions after each completed column block
    ends <- cumsum(block_sizes)

    abline(v = ends[-length(ends)] + 0.5, col = "lightgray")
  }

  # invisibly return
  invisible(out)
}
