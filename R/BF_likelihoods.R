
# full log-likelihood of simple MPT models
llMPT <- function(pars, mod, dataset = 1) {
  if (is.vector(pars)) pars <- t(pars)
  pars <- as.matrix(pars)
  ll <- loglikMPT(
    pars,
    unlist(mod$mptInfo$data[dataset, ]),
    mod$mptInfo$MPT$a, mod$mptInfo$MPT$b,
    mod$mptInfo$MPT$c, mod$mptInfo$MPT$map
  )
  const <- logMultinomCoefficient(mod, dataset = dataset)
  ll + const
}

# product-multinomial constant for density
logMultinomCoefficient <- function(mod, dataset = 1) {
  tree <- mod$mptInfo$MPT$tree.idx
  data <- mod$mptInfo$data[dataset, ]
  logCoef <- tapply(t(data), list(tree), function(n) {
    lgamma(sum(n) + 1) - sum(lgamma(n + 1))
  })
  sum(logCoef)
}


# density for product of beta distributions
dProductBeta <- function(x, shapes, log = TRUE) {
  x <- as.matrix(x)
  ll <- 0
  for (i in 1:nrow(shapes)) {
    ll <- ll + dbeta(x[, i], shapes[i, 1],
      shapes[i, 2],
      log = TRUE
    )
  }
  if (!log) ll <- exp(ll)
  ll
}


rProductBeta <- function(n, shapes) {
  apply(shapes, 1, function(s) rbeta(n, s[1], s[2]))
}
