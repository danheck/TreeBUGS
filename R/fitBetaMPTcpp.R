
#' C++ Sampler for Hierarchical Beta-MPT Model
#'
#' Fast Gibbs sampler in C++ that is tailored to the beta-MPT model.
#'
#' @inheritParams betaMPT
#' @param shape shape parameter(s) of Gamma-hyperdistribution for the hierarchical beta-parameters \eqn{\alpha_s} and \eqn{\beta_s} (can be a named vector to provide different hyperpriors for each parameter)
#' @param rate rate parameter(s) of Gamma-hyperdistribution
#' @author Daniel Heck
#' @importFrom parallel parLapply stopCluster detectCores
#' @export
betaMPTcpp <- function(eqnfile, data, restrictions,
                       covData, corProbit=FALSE,
                       n.iter=20000, n.burnin = 2000,
                       n.thin = 5,  n.chains=3, ppp = 0,
                       shape = 1, rate = .1, parEstFile){

  hyperprior <- list(shape = shape, rate = rate)

  fittedModel <- fitModelCpp("betaMPT",
                             eqnfile=eqnfile,
                             data=data,
                             restrictions=restrictions,
                             covData=covData,
                             corProbit=corProbit,
                             hyperprior=hyperprior,
                             n.iter=n.iter,
                             n.burnin=n.burnin, n.thin=n.thin,
                             n.chains=n.chains,
                             ppp = ppp,
                             parEstFile=parEstFile,
                             call = match.call())
  fittedModel
}
