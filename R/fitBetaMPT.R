#' Fit a Hierarchical Beta-MPT Model
#'
#' Fits a Beta-MPT model (Smith & Batchelder, 2010) based on a standard MPT model file (.eqn) and individual data table (.csv).
#'
#' @param eqnfile The (full path to the) file that specifies the MPT model (standard .eqn syntax). Note that category labels must start with a letter (different to multiTree) and match the column names of \code{data}
#' @param data The (full path to the) csv file with the data (comma separated; category labels in first row). Alternatively: a data frame or matrix (rows=individuals, columns = individual category frequencies, category labels as column names)
#' @param restrictions  Specifies which parameters should be (a) constant (e.g., \code{"a=b=.5"}) or (b) constrained to be identical (e.g., \code{"Do=Dn"}) or (c) treated as fixed effects (i.e., identical for all participants; \code{"a=b=FE"}). Either given as the path to a text file with restrictions per row or as a list of restrictions, e.g., \code{list("D1=D2","g=0.5")}
#' @param covData Data that contains covariates, for which correlations with individual MPT parameters will be sampled. Either the path to a .csv file (comma-separated: rows=individuals in the same order as \code{data}; first row must contain covariate labels); or alternatively: a data frame or matrix (rows=individuals, columns = variables; covariate labels as column names). Note that in \code{betaMPT}, correlatios are computed for discrete variables that are coded numerically (in \code{traitMPT}, this can be suppressed by using \code{predType="f"})
#' @param transformedParameters list with parameter transformations that should be computed based on the posterior samples (e.g., for testing parameter differences: \code{list("diffD=Do-Dn")})
#' @param modelfilename Name that the modelfile that is made by the function to work with JAGS should get.
#'        Default is to write this information to the tempdir as required by CRAN standards.
#' @param corProbit whether to use probit-transformed MPT parameters to compute correlations (probit-values of \code{+Inf} are truncated to \code{max(5,max(probit))}; similarly for \code{-Inf}). Default for beta-MPT: MPT parameters are used on the probability scale [0,1].
#' @param alpha Hyperprior for the shape parameters \eqn{\alpha} of the group-level beta distributions (in JAGS syntax). Default: Gamma distributions for \eqn{\alpha} and \eqn{\beta} with shape 1 and rate .1. To use uniform priors on the interval [.01,5000] as proposed by Smith and Batchelder (2008), use \code{alpha = "dunif(.01,5000)"} and \code{beta = "dunif(.01,5000)"}. Note that a vector can be used to specify separate hyperpriors for each MPT parameter (the order of parameters is determined by the names of the vector or by the default order as shown in \code{\link{readEQN}} with \code{paramOrder = TRUE}).
#' @param beta Hyperprior for \eqn{\beta} of group-level distributions, see \code{alpha}
#' @param parEstFile Name of the file to with the estimates should be stored (e.g., "parEstFile.txt")
#' @param n.iter Number of iterations per chain (including burnin samples). See \code{\link[runjags]{run.jags}} for details.
#' @param n.adapt number of adaption samples to adjust MCMC sampler in JAGS. The sampler will be more efficient if it is tuned well.
#' @param n.burnin Number of samples for burnin (samples will not be stored and removed from n.iter)
#' @param n.thin Thinning rate.
#' @param n.chains number of MCMC chains (sampled in parallel).
#' @param dic whether to compute DIC using \code{\link[runjags]{extract}}, which requires additional sampling. Can also be computed and added after fitting the model by \code{fittedModel$dic <- extract(fittedModel$runjags, "dic")}
#' @param ppp number of samples to compute  posterior predictive p-value (see \code{\link{posteriorPredictive}})
#' @param autojags if provided (as an empty list or with arguments passed to \link[runjags]{autoextend.jags}), JAGS runs repeatedly until the MCMC chains converges . E.g., use \code{list(max.time="30m")} to restrict sampling to 30 minutes (similarly for hours, days, and weeks)
#' @param ... Arguments to be passed to the JAGS sampling function (i.e., to \code{\link[runjags]{run.jags}}.
#'
#' @details Note that, in the Beta-MPT model, correlations of individual MPT parameters with covariates are sampled. Hence, the covariates do not affect the estimation of the actual Beta-MPT parameters. Therefore, the correlation of covariates with the individual MPT parameters can equivalently be performed after fitting the model using the sampled posterior parameter values stored in \code{betaMPT$mcmc}
#'
#' @return a list of the class \code{betaMPT} with the objects:
#' \itemize{
#'  \item \code{summary}: MPT tailored summary. Use \code{summary(fittedModel)}
#'  \item \code{mptInfo}: info about MPT model (eqn and data file etc.)
#'  \item \code{runjags}: the object returned from the MCMC sampler. Note that the object \code{fittedModel$runjags} is an \link[runjags]{runjags} object, whereas \code{fittedModel$runjags$mcmc} is a \code{mcmc.list} as used by the coda package (\link[coda]{mcmc})
#' }
#' @author Daniel Heck, Nina R. Arnold, Denis Arnold,
#' @references Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial processing tree models for addressing individual differences. Journal of Mathematical Psychology, 54, 167-183.
#' @export

betaMPT <- function(eqnfile, data, restrictions,
                    covData, #covStructure,
                    transformedParameters,
                    corProbit=FALSE,
                    alpha = "dgamma(1,.1)",
                    beta = "dgamma(1,.1)",
                    # alpha = "dunif(0,1)",
                    # beta = "dnorm(0,1)T(0,)",

                    # MCMC stuff:
                    n.iter=20000, n.adapt=2000,
                    n.burnin=2000, n.thin=5,
                    n.chains=3, dic =FALSE,
                    ppp = 0,

                    # File Handling stuff:
                    modelfilename, parEstFile,
                    autojags = NULL,   ...){

  hyperprior <- list(alpha=alpha, beta=beta)

  fittedModel <- fitModel(type="betaMPT", eqnfile=eqnfile,
                          data=data,restrictions=restrictions,
                          covData=covData,#covStructure=covStructure,
                          transformedParameters=transformedParameters,
                          corProbit=corProbit, hyperprior=hyperprior,
                          n.iter=n.iter, n.adapt = n.adapt,
                          n.burnin=n.burnin, n.thin=n.thin,
                          n.chains=n.chains, dic =dic,  ppp = ppp,
                          modelfilename=modelfilename,
                          parEstFile=parEstFile,
                          autojags=autojags,
                          call = match.call(),
                          ...)

  return(fittedModel)
}
