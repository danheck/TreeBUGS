###############################################################################
### Posterior distribution for Pearson's correlation coefficient
###
### Code by Alexander Ly
### * https://github.com/AlexanderLyNL/jasp-desktop/blob/development/JASP-Engine/JASP/R/correlationbayesian.R#L1581
###
### * Ly, A., Marsman, M., Wagenmakers, E.-J. (2015).
###   Analytic Posteriors for Pearson’s Correlation Coefficient.
###   Manuscript submitted for publication. https://arxiv.org/abs/1510.01188
###
### * Alexander Ly, Udo Boehm, Andrew Heathcote, Brandon M. Turner, Birte Forstmann,
###   Maarten Marsman, and Dora Matzke (2016).
###   A flexible and efficient hierarchical Bayesian approach to the exploration of
###   individual differences in cognitive-model-based neuroscience. https://osf.io/evsyv/
###
###############################################################################

# require("hypergeo")
#' @importFrom hypergeo hypergeo

## Auxilary functions ------------------------------------------------------------
# 0. Prior specification
.scaledBeta <- function(rho, alpha, beta){
  result <- 1/2*dbeta((rho+1)/2, alpha, beta)
  return(result)
}


.priorRho <- function(rho, kappa=1) {
  .scaledBeta(rho, 1/kappa, 1/kappa)
}

# 1.0. Built-up for likelihood functions
.aFunction <- function(n, r, rho) {
  hyper.term <- Re(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=(1/2), z=(r*rho)^2))
  result <- (1-rho^2)^((n-1)/2)*hyper.term
  return(result)
}

.bFunction <- function(n, r, rho) {
  hyper.term <- Re(hypergeo::genhypergeo(U=c(n/2, n/2), L=(3/2), z=(r*rho)^2))
  log.term <- 2*(lgamma(n/2)-lgamma((n-1)/2))+((n-1)/2)*log(1-rho^2)
  result <- 2*r*rho*exp(log.term)*hyper.term
  return(result)
}

.hFunction <- function(n, r, rho) {
  result <- .aFunction(n, r, rho) + .bFunction(n, r, rho)
  return(result)
}

.jeffreysApproxH <- function(n, r, rho) {
  result <- ((1 - rho^(2))^(0.5*(n - 1)))/((1 - rho*r)^(n - 1 - 0.5))
  return(result)
}


#
# 2.1 Two-sided main Bayes factor ----------------------------------------------
.bf10Exact <- function(n, r, kappa=1) {
  # Ly et al 2015
  # This is the exact result with symmetric beta prior on rho
  # with parameter alpha. If kappa = 1 then uniform prior on rho
  #
  #
  if (n <= 2){
    return(1)
  } else if (any(is.na(r))){
    return(NaN)
  }
  # TODO: use which
  check.r <- abs(r) >= 1 # check whether |r| >= 1
  if (kappa >= 1 && n > 2 && check.r) {
    return(Inf)
  }
  #log.hyper.term <- log(hypergeo::hypergeo(((n-1)/2), ((n-1)/2), ((n+2/kappa)/2), r^2))
  log.hyper.term <- log(hypergeo::genhypergeo(U=c((n-1)/2, (n-1)/2), L=((n+2/kappa)/2), z=r^2))
  log.result <- log(2^(1-2/kappa))+0.5*log(pi)-lbeta(1/kappa, 1/kappa)+
    lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)+log.hyper.term
  real.result <- exp(Re(log.result))
  #return(realResult)
  return(real.result)
}

# 2.2 Two-sided secondairy Bayes factor
.bf10JeffreysIntegrate <- function(n, r, kappa=1) {
  # Jeffreys' test for whether a correlation is zero or not
  # Jeffreys (1961), pp. 289-292
  # This is the exact result, see EJ
  ##
  if (n <= 2){
    return(1)
  } else if ( any(is.na(r)) ){
    return(NaN)
  }

  # TODO: use which
  if (n > 2 && abs(r)==1) {
    return(Inf)
  }
  hyper.term <- Re(hypergeo::genhypergeo(U=c((2*n-3)/4, (2*n-1)/4), L=(n+2/kappa)/2, z=r^2))
  log.term <- lgamma((n+2/kappa-1)/2)-lgamma((n+2/kappa)/2)-lbeta(1/kappa, 1/kappa)
  result <- sqrt(pi)*2^(1-2/kappa)*exp(log.term)*hyper.term
  return(result)
}


# 4.1 Two-sided
.posteriorRho <- function(n, r, rho, kappa=1){
  if (!is.na(r) && !r==0){
    return(1/.bf10Exact(n, r, kappa)*.hFunction(n, r, rho)*.priorRho(rho, kappa))
  } else if (!is.na(r) && r==0){
    return(1/.bf10JeffreysIntegrate(n, r, kappa)*.jeffreysApproxH(n, r, rho)*.priorRho(rho, kappa))
  }
}
