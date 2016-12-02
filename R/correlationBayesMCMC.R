#' Posterior Distribution for Correlations
#'
#' Adjusts the posterior distribution of correlations for the sampling error of a population correlation according to the sample size (i.e., the number of participants; Ly, Marsman, & Wagenmakers, 2016).
#' @param fittedModel a fitted \link{betaMPT} or \link{traitMPT} model with covariates (added during fitting by the argument \code{covData})
#' @param r optional: a vector of posterior correlations (instead of \code{fittedModel})
#' @param N only if \code{r} is used: the number of participants the correlation is based on
#' @param kappa parameter for the prior of the correlation, that is, a scaled beta distribution: Beta(1/kappa, 1/kappa). The default \code{kappa=1} defines a uniform distribution on [-1,1], whereas \code{kappa<1} defines a unimodal prior centered around zero.
#' @param precision precision on the interval [-1,1] to approximate the posterior density
#' @param plot whether to plot (a) the unadjusted posterior correlations (gray histogram) and (b) the corrected posterior (black line with red credibility intervals)
#' @param nCPU number of CPUs used for parallel computation of posterior distribution
#' @param M number of subsamples from the fitted model
#' @param ci credibility interval
#' @details
#' This function (1) uses all posterior samples of a correlation to (2) derive the posterior of the correlation corrected for sampling error and (3) averages these densities across the posterior samples. Thereby, the method accounts for estimation uncertainty of the MPT model (due to the use of the posterior samples) and also for sampling error of the population correlation due to sample size (cf. Ly, Boehm, Heathcote, Turner, Forstmann, Marsman, & Matzke, 2016).
#' @references
#' Ly, A., Marsman, M., Wagenmakers, E.-J. (2015). Analytic Posteriors for Pearson's Correlation Coefficient. Manuscript submitted for publication. https://arxiv.org/abs/1510.01188
#'
#' Ly, A., Boehm, U., Heathcote, A., Turner, B. M. , Forstmann, B., Marsman, M., and Matzke, D. (2016). A flexible and efficient hierarchical Bayesian approach to the exploration of individual differences in cognitive-model-based neuroscience. https://osf.io/evsyv/
#' @author Daniel W. Heck, Alexander Ly
#' @importFrom parallel clusterEvalQ
#' @examples
#' # test effect of number of participants:
#' cors <- rbeta(50, 100, 70)
#' correlationPosterior(r=cors, N=10, nCPU=1)
#' correlationPosterior(r=cors, N=100, nCPU=1)
#' @export
correlationPosterior <- function(fittedModel,
                                 r, N,
                                 kappa=1, ci=.95,
                                 M = 1000, precision=.005, plot=TRUE, nCPU=4){

  rho <- seq(-1, 1, by = precision)
  if(!missing(fittedModel) && !is.null(fittedModel)){
    N <- nrow(fittedModel$mptInfo$data)
    chains <- length(fittedModel$runjags$mcmc)
    M.fit <- nrow(fittedModel$runjags$mcmc[[1]])

    sel.idx <- grep("cor_", varnames(fittedModel$runjags$mcmc))
    # if(class(fittedModel) == "betaMPT")
    #   sel.idx <- c(sel.idx, grep("rho", varnames(fittedModel$runjags$mcmc)))
    r <- do.call("rbind",
                      fittedModel$runjags$mcmc[sample(M.fit, min(M.fit, ceiling(M/chains))),
                                               sel.idx])
  }else{
    if(missing(N) || is.na(N) || N != round(N))
      stop("Number of participants 'N' missing or not an integer!")
    if(missing(r) || is.na(r) || is.null(r))
      stop("Correlation posterior samples 'r' missing!")

    r <- as.matrix(r)
    if(is.null(colnames(r))) colnames(r) <- paste0("corr",1:ncol(r))
  }

  singleCorrelation <- function(r.samples){
    rowMeans(sapply(r.samples,
                    function(ss)  .posteriorRho(N, ss, rho=rho, kappa=kappa)))
  }

  if(nCPU>1){
    cl <- makeCluster(nCPU)
    tmp <- clusterEvalQ(cl, library(hypergeo))
    clusterExport(cl, c("kappa", "N","precision","rho",
                        ".scaledBeta", ".priorRho",".aFunction",".bFunction", ".hFunction",
                        ".jeffreysApproxH", ".bf10Exact", ".bf10JeffreysIntegrate",".posteriorRho" ),
                  envir = environment())
    r.post <- t(parApply(cl, r, 2, singleCorrelation))
    stopCluster(cl)
  }else{
    r.post <- t(apply(r, 2, singleCorrelation))
  }

  colnames(r.post) <- rho
  attr(r.post, "kappa") <- kappa

  idx <- apply(r.post*precision, 1, function(fx)
    c(max(which(cumsum(fx) < (1-ci)/2)),
      max(which(cumsum(fx) < .50)),
      min(which(cumsum(fx) > (1+ci)/2))))
  summ <- t(apply(idx,2,function(i) rho[i]))
  colnames(summ) <- c(paste0(100*(1-ci)/2,"%") , "50%", paste0(100*(1+ci)/2,"%") )

  if(plot){
    mfrow <- par()$mfrow
    if(ncol(r)==2)
      par(mfrow = c(1,2))
    if(ncol(r)>=4)
      par(mfrow = c(2,2))
    for(i in 1:ncol(r)){
      hist(r[,i],min(60, max(25,round(nrow(r)/3))),
           freq=FALSE, xlim=c(-1,1), col="gray",border="gray",
           ylab="Density",
           xlab="Correlation", main=rownames(r.post)[i])
      lines(rho, r.post[i,])
      abline(v=summ[i,c(1,3)], col=2)
    }
    par(mfrow = mfrow)
  }

  summ
}
