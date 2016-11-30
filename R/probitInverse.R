
#' Probit-Inverse of Group-Level Normal Distribution
#'
#' Transform latent group-level normal distribution (latent-trait MPT) into mean and SD on probability scale.
#' @param mu latent-probit mean of normal distribution
#' @param sigma latent-probit SD of normal distribution
#' @param fittedModel optional: fitted \link{traitMPT} model. If provided, the bivariate inverse-probit transform is applied to all MCMC samples (and \code{mu} and \code{sigma} are ignored).
#' @return implied mean and SD on probability scale
#' @examples
#' ####### compare bivariate vs. univariate transformation
#' probitInverse(mu=.8, sigma=c(.25,.5,.75,1))
#' pnorm(.8)
#'
#' # full distribution
#' prob <- pnorm(rnorm(10000, .8, .7))
#' hist(prob, 80, col="gray", xlim=0:1)
#'
#' \dontrun{
#' # transformation for fitted model
#' mean_sd <- probitInverse(fittedModel=fit)
#' summarizeMCMC(mean_sd)
#' }
#' @export
#' @importFrom coda varnames
#' @importFrom stats integrate
probitInverse <- function(mu, sigma, fittedModel=NULL){

  probitInverseVec <- Vectorize(
    function(mu, sigma){
      mp <- vp <- NA
      try({
        mp <- integrate(function(x)
          pnorm(x)*dnorm(x, mu, sigma),
          # -Inf, Inf)$value
          mu-5*sigma, mu+5*sigma)$value
        vp <- integrate(function(x)
          (pnorm(x)-mp)^2*dnorm(x, mu, sigma),
          # -Inf, Inf)$value
          mu-5*sigma, mu+5*sigma)$value
      })
      if(!is.na(vp) && vp<0){
        vp <- NA
      }else if(is.na(vp)){
        mp <- NA
      }
      return(c(mean=mp, sd = sqrt(vp)))
    }, c("mu","sigma"))

  if(missing(fittedModel) || is.null(fittedModel)){
    res <- t(probitInverseVec(mu, sigma))
    if(any(is.na(res))) {
      cat("Transformation resulted in NAs for:\n")
      print(cbind(mu=mu, sigma=sigma, res)[apply(is.na(res),1,any),])
    }
    return(res)
  }else{
    if(!class(fittedModel) == "traitMPT")
      stop("'fittedModel' must be a latent-trait MPT model.")

    samp <- fittedModel$runjags$mcmc[,c()]
    thetaNames <- fittedModel$mptInfo$thetaUnique
    sel.mu <- grep("mu", varnames(fittedModel$runjags$mcmc))
    sel.sig <- grep("sigma", varnames(fittedModel$runjags$mcmc))
    s.mu <- fittedModel$runjags$mcmc[,sel.mu]
    s.sig <- fittedModel$runjags$mcmc[,sel.sig]
    # cl <- makeCluster(nCPU)
    # clusterExport(cl, c("probitInverseVec","s.mu","s.sig", "thetaNames"))
    for(cc in 1:length(s.mu)){
      for(s in 1:ncol(s.mu[[cc]])){
        res <- probitInverseVec(mu = s.mu[[cc]][,s], sigma=s.sig[[cc]][,s])
        rownames(res) <- paste0(c("mean_","sd_"), thetaNames[s])
        samp[[cc]] <- cbind(samp[[cc]], t(res))
      }
      samp[[cc]] <- mcmc(samp[[cc]])
      attr(samp[[cc]], "mcpar") <-  attr(fittedModel$runjags$mcmc[[cc]], "mcpar")
    }
    return(as.mcmc.list(samp))
  }
}
