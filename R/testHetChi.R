#' Chi-Square Test of Heterogeneity
#'
#' Tests whether whether participants (items) are homogeneous under the
#' assumption of item (participant) homogeneity.
#'
#' @param freq matrix with observed frequencies (rows: persons/items; columns:
#'   categories). Can also be the path to a .csv file with frequencies
#'   (comma-separated; first line defines category labels)
#' @param tree a vector defining which columns of x belong to separate
#'   multinomial distributions (i.e., MPT trees). For instance, if \code{x} has
#'   five categories from two MPT trees: \code{tree=c(1,1,2,2,2)} or
#'   \code{tree=c("t1","t1","t2","t2","t2")}
#'
#' @details If an item/person has zero frequencies on all categories in an MPT
#' tree, these zeros are neglected when computing mean frequencies per column.
#' As an example, consider a simple recognition test with a fixed assignments of
#' words to the learn/test list. In such an experiment, all learned words will
#' result in hits or misses (i.e., the MPT tree of old items), whereas new words
#' are always false alarms/correct rejections and thus belong to the MPT tree of
#' new items (this is not necessarily the case if words are assigned randomly).
#'
#' Note that the test assumes independence of observations and item homogeneity
#' when testing participant heterogeneity. The latter assumption can be dropped
#' when using a permutation test (\code{\link{testHetPerm}}).
#' @seealso \code{\link{testHetPerm}}, \code{\link{plotFreq}}
#' @author Daniel W. Heck
#' @references Smith, J. B., & Batchelder, W. H. (2008). Assessing individual
#'   differences in categorical data. Psychonomic Bulletin & Review, 15,
#'   713-731. \doi{10.3758/PBR.15.4.713}
#'
#' @examples
#' # some made up frequencies:
#' freq <- matrix(
#'   c(
#'     13, 16, 11, 13,
#'     15, 21, 18, 13,
#'     21, 14, 16, 17,
#'     19, 20, 21, 18
#'   ),
#'   ncol = 4, byrow = TRUE
#' )
#' # for a product-binomial distribution:
#' # (categories 1 and 2 and categories 3 and 4 are binomials)
#' testHetChi(freq, tree = c(1, 1, 2, 2))
#' # => no significant deviation from homogeneity (low power!)
#' @export
testHetChi <- function(
    freq,
    tree
) {
  # gen2htm <- genBetaMPT(24, c(Target=100,Lure=100), htm,
  #                       mean=c(Do=.7, Dn=.7, g=.4),
  #                       sd=c(Do=0.05, Dn=0.05, g=0.05))
  # freq <- gen2htm$data #matrix(round(runif(items*subjs, .5, 14.49)),subjs, 10)
  # tree <- c(1,1,2,2) # c(rep(1,6), 2,2,3,3)

  if (is.character(freq)) {
    freq <- read.csv(file = freq)
  }

  freq <- as.matrix(freq)

  if (any(is.character(freq) | freq != round(freq))) {
    stop("The data ('freq') may only contain frequencies!")
  }
  if (missing(tree)) {
    warning(
      "It is assumed that all columns of 'freq' stem from one multinomial distribution",
      "\n  (i.e., from a single MPT tree)"
    )
    tree <- rep(1, ncol(freq))
  }

  if (length(tree) != ncol(freq)) {
    stop("Length of vector 'tree' must be identical to number of columns of 'freq'.")
  }

  # number of participants/items:
  N <- nrow(freq)
  # numer of items per person/tree:
  if (length(unique(tree)) == 1) {
    M <- matrix(rowSums(freq))
    colnames(M) <- unique(tree)
  } else {
    M <- t(apply(freq, 1, function(xx) tapply(xx, tree, sum)))
  }
  M[M == 0] <- NA
  # number of categories per tree:
  K <- tapply(freq[1, ], tree, length)

  # compute mean frequencies across proportions
  # (not raw frequencies: different M!)
  prop <- freq / M[, tree]
  mu <- colMeans(prop, na.rm = TRUE)

  freq.exp <- M[, tree] * matrix(mu, nrow(M), length(mu), byrow = TRUE)

  # for(t in seq_along(tree)){
  #   sel <- rowSums(freq[ , tree == tree[t] ]) == 0
  #   if(any(sel))
  #     freq.exp[sel , tree == tree[t] ] <- NA
  # }

  # freq.exp <- matrix(mu, N, sum(K), byrow = TRUE)
  chi <- sum((freq - freq.exp)^2 / freq.exp, na.rm = TRUE)

  df <- sum(K - 1) * (N - 1)
  list(
    chisq = chi, df = df,
    prob = pchisq(chi, df, lower.tail = FALSE)
  )
}
