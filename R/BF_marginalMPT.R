
#' Marginal Likelihood for Simple MPT
#'
#' Computes the marginal likelihood for simple (fixed-effects, nonhierarchical) MPT models.
#'
#' @inheritParams simpleMPT
#' @param method either \code{"importance"} (importance sampling using a mixture of uniform and beta-aproximation of the posterior) or \code{"prior"} (brute force Monte Carlo sampling from prior)
#' @param posterior number of posterior samples used to approximate importance-sampling densities (i.e., beta distributions)
#' @param mix mixture proportion of the uniform distribution for the importance-sampling density
#' @param scale how much should posterior-beta approximations be downscaled to get fatter importance-sampling density
#' @param sample number of paramter samples
#' @param batches number of batches to compute standard error
#' @param show whether to show progress
#' @param nCPU number of CPUs used
#'
#' @details
#' Currently, this is only implemented for a single data set!
#'
#' If \code{method = "prior"}, a brute-force Monte Carlo method is used and parameters are directly sampled from the prior.Then, the likelihood is evaluated for these samples and averaged (fast, but inefficient).
#'
#' Alternatively, an importance sampler is used if \code{method = "importance"}, and the posterior distributions of the MPT parameters are approximated by independent beta distributions. Then each parameter \eqn{s} is sampled from the importance density:
#'
#' \eqn{mix*U(0,1) + (1-mix)*Beta(scale*a_s, scale*b_s)}
#'
#' @examples
#' \dontrun{
#' # 2-High-Threshold Model
#' eqn <- "## 2HTM ##
#'    Target  Hit  d
#'    Target  Hit  (1-d)*g
#'    Target  Miss (1-d)*(1-g)
#'    Lure    FA   (1-d)*g
#'    Lure    CR   (1-d)*(1-g)
#'    Lure    CR   d"
#' data <- c(Hit = 46, Miss = 14,
#'           FA = 14, CR = 46)
#'
#' # informative priors for guessing
#' aa <- c(d = 1, g = 2)
#' bb <- c(d = 1, g = 2)
#' curve(dbeta(x, aa["g"], bb["g"]))
#'
#' # compute marginal likelihood
#' htm <- marginalMPT(eqn, data,
#'                    alpha = aa, beta = bb,
#'                    posterior = 200,
#'                    sample = 500, nCPU = 1)
#' # second model: g=.50
#' htm.g50 <- marginalMPT(eqn, data, list("g=.5"),
#'                        alpha = aa, beta = bb,
#'                        posterior = 200,
#'                        sample = 500, nCPU = 1)
#'
#' # Bayes factor
#' # (per batch to get estimation error)
#' bf <- htm.g50$p.per.batch / htm$p.per.batch
#' mean(bf)                 # BF
#' sd(bf)/sqrt(length(bf))  # standard error
#' }
#'
#' @seealso \code{\link{BayesFactorMPT}}
#' @references
#' Vandekerckhove, J. S., Matzke, D., & Wagenmakers, E. (2015). Model comparison and the principle of parsimony. In Oxford Handbook of Computational and Mathematical Psychology (pp. 300-319). New York, NY: Oxford University Press.
#' @export
marginalMPT <- function (eqnfile,
                         data,
                         restrictions,
                         alpha = 1,
                         beta = 1,
                         method = "importance",
                         posterior = 500,
                         mix = .01,
                         scale = .9,
                         sample = 50000,
                         batches = 10,
                         show = TRUE,
                         nCPU = 4){
  t0 <- Sys.time()
  if (mix<0 | mix>1 | scale<0 | scale>1)
    stop ("The tuning parameters 'mix' and 'scale' must be in the interval [0,1].")
  if (is.character(data)) data <- readData(data)
  if (!is.vector(data) && nrow(data) > 1)
    stop ("Computation of marginal likelihoods not suported for",
          "multiple data sets/participants.")

  ############### 1. fit MPT / get MPT structure
  tmp <- capture.output(
    mod <- simpleMPT(eqnfile = eqnfile, data = data, restrictions = restrictions,
                     n.iter = ifelse(method=="importance",posterior+200,5),
                     n.burnin = ifelse(method=="importance", 200, 2),
                     n.thin = 1, n.chains = 1, alpha = alpha, beta = beta))
  if (method == "importance"){
    betapar <- approximatePosterior(mod)*scale
  }
  S <- length(mod$mptInfo$thetaUnique)

  if (show) cat("Sampling parameters to estimate marginal likelihood...\n")
  marginal <- rep(NA, batches)
  cl <- makeCluster(nCPU)
  samp <- matrix(NA, sample, batches)
  const <- logMultinomCoefficient(mod)
  for (b in 1:batches){
    if (show) cat(b,"")

    if (method == "importance"){
      # sample from mixture: mix*U(0,1) + (1-mix)*Beta(a,b)
      rmix <- function (i, sample, betapar, mix) ifelse(runif(sample)<mix,
                                                       runif(sample),
                                                       rbeta(sample, betapar[i,1], betapar[i,2]))
      xx <- parSapply(cl, 1:nrow(betapar), rmix,
                      sample=sample, betapar=betapar, mix=mix)
      # sampling density:
      dmix <- function (x, betapar, mix)
        sum(log(mix + (1-mix)*dbeta(x, betapar[,1], betapar[,2])))
      gx <- parApply(cl, xx, 1, dmix, betapar=betapar, mix=mix)
      # prior:
      px <- rep(0, sample)
      for(s in 1:S)
        px <- px + dbeta(xx[,s], log = TRUE,
                         shape1 = mod$mptInfo$hyperprior$alpha[s],
                         shape2 = mod$mptInfo$hyperprior$alpha[s])
    } else if (method == "prior"){
      xx <- matrix(NA, sample, S)
      for (s in 1:S)
        xx[,s] <- rbeta(sample,
                        shape1 = mod$mptInfo$hyperprior$alpha[s],
                        shape2 = mod$mptInfo$hyperprior$beta[s])
    }

    # likelihood:
    # fx <- parapply(cl, xx, 1, )
    fx <- #apply(xx, 1,
      c(loglikMPT(xx,
                  h=unlist(mod$mptInfo$data),
                  a = mod$mptInfo$MPT$a,
                  b = mod$mptInfo$MPT$b,
                  c = mod$mptInfo$MPT$c,
                  map = mod$mptInfo$MPT$map))

    # importance weights wx <- px-gx
    if (method == "importance"){
      samp[,b] <- exp(fx + px - gx + const)
    } else if (method == "prior"){
      samp[,b] <- exp(fx + const)
    }
  }
  stopCluster(cl)
  t1 <- Sys.time()

  ############### 4. Aproximate integral + batch SE
  p <- c(mean = mean(c(samp)), SE=sd(c(samp))/sqrt(sample*batches))
  p.batch <- colMeans(samp)
  p.batch <- c("p" = mean(p.batch),
               "SE_p" = sd(p.batch)/sqrt(batches),
               "log_p" = mean(log(p.batch)),
               "SE_log_p" = sd(log(p.batch))/sqrt(batches))

  res <- list("p" = p,
              "p.batch" = p.batch,
              "p.per.batch" = colMeans(samp),
              "prior" = do.call("cbind", mod$mptInfo$hyperprior),
              "data" = mod$mptInfo$data,
              "sampler" = list("method" = method,
                               "sample" = sample,
                               "mix" = mix,
                               "scale" = scale,
                               "time" = t1-t0))
  class(res) <- "marginalMPT"
  if (show)
    cat("\nFinished in ", format(t1-t0))
  return (res)
}


