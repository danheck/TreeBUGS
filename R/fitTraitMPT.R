#' Fit a Hierarchical Latent-Trait MPT Model
#'
#' Fits a latent-trait MPT model (Klauer, 2010) based on a standard MPT model file
#' (.eqn) and individual data table (.csv).
#'
#' @inheritParams betaMPT
#' @param predStructure  Defines which variables in \code{covData} are included
#'     as predictors for which MPT parameters. Either the path to the file that
#'     specifies the assigment of MPT parameters to covariates (that is, each row
#'     assigns one or more MPT parameters to one or more covariates, separated
#'     by a semicolon, e.g., \code{Do g; age extraversion}). Can also be provided
#'     as a list, e.g., \code{list("Do Dn ; age", "g ; extraversion"}).
#'     Note that no correlations of MPT parameters and predictors are computed.
#' @param predType a character vector specifying the type of continuous or
#'     discrete predictors in each column of \code{covData}:
#'     \code{"c"} = continuous covariate (which are centered to have a mean of zero);
#'     \code{"f"} = discrete predictor, fixed effect (default for character/factor variables);
#'     \code{"r"} = discrete predictor, random effect.
#' @param mu hyperprior for group means of probit-transformed parameters in JAGS syntax.
#'     Default is a standard normal distribution, which implies a uniform
#'     distribution on the MPT probability parameters. A named vector can be used to
#'     specify separate hyperpriors for each MPT parameter (the order of parameters
#'     is determined by the names of the vector or by the default order as shown
#'     in \code{\link{readEQN}} with \code{paramOrder = TRUE}).
#' @param xi hyperprior for scaling parameters of the group-level parameter variances.
#'     Default is a uniform distribution on the interval [0,10].
#'     Similarly as for \code{mu}, a vector of different priors can be used.
#'     Less informative priors can be used (e.g., \code{"dunif(0,100)")}) but
#'     might result in reduced stability.
#' @param V  S x S matrix used as a hyperprior for the inverse-Wishart hyperprior
#'     parameters with as many rows and columns as there are core MPT parameters.
#'     Default is a diagonal matrix.
#' @param df degrees of freedom for the inverse-Wishart hyperprior for the individual
#'     parameters. Minimum is S+1, where S gives the number of core MPT parameters.
#' @param IVprec hyperprior on the precision (i.e., the inverse of the variance)
#'     of the slope parameters for the z-standardized continuous predictors.
#'     For ease of interpretation, TreeBUGS reports unstandardized regression coefficients.
#'     See details below.
#'
#' @section Regression Extensions:
#' Continuous (discrete) predictors are added on the latent-probit scale via:
#' \deqn{\theta = \Phi(\mu + X \beta +\delta ),}
#' where X is a design matrix with centered (!) continuous covariates and recoded factor
#' variables (using an orthogonal contrast coding scheme; cf. Rouder et al., 2012).
#' TreeBUGS reports unstandardized regression coefficients that correspond to the
#' scale/SD of the predictor variables. However, since the regression is on the latent
#' probit scale, the coefficients are not standardized with respect to the 'depend variable'
#' as in the standard linear regression as in Rouder & Morey (2012).
#'
#' For continuous predictors, the default prior \code{IVprec = "dchisq(1)"} implies
#' a Cauchy prior on each standardized \eqn{\beta} (similar to the Jeffreys-Zellner-Siow prior
#' with scale parameter \eqn{s=1}; for details, see: Rouder et. al, 2012; Rouder & Morey, 2012).
#' If small effects are expected, smaller scale values \eqn{s} can be used by changing the default to
#' \code{IVprec = 'dgamma(1/2,(s^2)/2)'}.
#' To use a standard-normal priors on the standardized slopes, use \code{IVprec = 'dcat(1)'}.
#'
#' @section Uncorrelated Latent-Trait Values:
#' The standard latent-trait MPT model assumes a multivariate normal distribution
#' of the latent-trait values, where the covariance matrix follows a scaled-inverse
#' Wishart distribution. As an alternative, the parameters can be assumed to be
#' independent (this is equivalent to a diagonal covariance matrix).
#' If the assumption of uncorrelated parameters is justified, such a simplified model
#' has less parameters and is more parsimonious, which in turn might result in more
#' robust estimation and more precise parameter estimates.
#'
#' This alternative method can be fitted in TreeBUGS (but not all of the features
#' of TreeBUGS might be compatible with this alternative model structure).
#' To fit the model, the scale matrix \code{V} is set to \code{NA}
#' (V is only relevant for the multivariate Wishart prior) and the prior
#' on \code{xi} is changed: \code{traitMPT(..., V=NA, xi="dnorm(0,1)")}.
#' The model assumes that the latent-trait values \eqn{\delta[i]} (=random-intercepts)
#' are decomposed by the scaling parameter \eqn{\xi} and the raw deviation \eqn{\epsilon[i]}
#' (cf. Gelman, 2006):
#'      \deqn{\delta[i] = \xi * \epsilon[i]}
#'      \deqn{\epsilon[i] ~ Normal(0,\sigma^2)}
#'      \deqn{\sigma^2 ~ Inverse_Chisquare(df)}
#' Note that the default prior for \eqn{\xi} should be changed to \code{xi="dnorm(0,1)"}, which
#' results in a half-Cauchy prior (Gelman, 2006).
#'
#' @return a list of the class \code{traitMPT} with the objects:
#' \itemize{
#'  \item \code{summary}: MPT tailored summary. Use \code{summary(fittedModel)}
#'  \item \code{mptInfo}: info about MPT model (eqn and data file etc.)
#'  \item \code{mcmc}: the object returned from the MCMC sampler.
#'                     Note that the object \code{fittedModel$mcmc} is an
#'                     \link[runjags]{runjags} object, whereas
#'                     \code{fittedModel$mcmc$mcmc} is an \code{mcmc.list} as used by
#'                     the coda package (\link[coda]{mcmc})
#' }
#' @author Daniel Heck, Denis Arnold, Nina R. Arnold
#' @references
#' Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models
#' (comment on article by Browne and Draper). Bayesian Analysis, 1(3), 515-534.
#'
#' Klauer, K. C. (2010).
#' Hierarchical multinomial processing tree models: A latent-trait approach.
#' Psychometrika, 75, 70-98.
#'
#' Matzke, D., Dolan, C. V., Batchelder, W. H., & Wagenmakers, E.-J. (2015).
#' Bayesian estimation of multinomial processing tree models with heterogeneity in participants and items.
#' Psychometrika, 80, 205-235.
#'
#' Rouder, J. N., Morey, R. D., Speckman, P. L., & Province, J. M. (2012).
#' Default Bayes factors for ANOVA designs.
#' Journal of Mathematical Psychology, 56, 356-374.
#'
#' Rouder, J. N., & Morey, R. D. (2012).
#' Default Bayes Factors for Model Selection in Regression.
#' Multivariate Behavioral Research, 47, 877-903.
#'
#' @examples
#' \dontrun{
#' # fit beta-MPT model for encoding condition (see ?arnold2013):
#' EQNfile <- system.file("MPTmodels/2htsm.eqn", package="TreeBUGS")
#' d.encoding <- subset(arnold2013, group == "encoding", select = -(1:4))
#' fit <- traitMPT(EQNfile, d.encoding, n.thin=5,
#'                 restrictions=list("D1=D2=D3","d1=d2","a=g"))
#' # convergence
#' plot(fit, parameter = "mean", type = "default")
#' summary(fit)
#' }
#' @export
traitMPT <- function(eqnfile, data, restrictions, covData, predStructure,
                     predType,  # one of: c("c," f", "r")
                     transformedParameters, corProbit=TRUE,

                     # hyperpriors:
                     mu = "dnorm(0,1)", xi = "dunif(0,10)",
                     V, df, IVprec = "dchisq(1)",  # change to "dcat(1)" to set beta ~ dnorm(0,1)

                     # MCMC stuff:
                     n.iter=20000, n.adapt = 2000, n.burnin=2000,  n.thin=5,
                     n.chains=3, dic =FALSE, ppp = 0,

                     # File Handling stuff:
                     modelfilename, parEstFile,posteriorFile,
                     autojags=NULL,   ...){

  if(missing(V)) V <- NULL
  if(missing(df)) df <- NULL
  hyperprior <- list(mu=mu, xi=xi, V=V, df=df, IVprec=IVprec)

  fitModel(type="traitMPT", eqnfile=eqnfile,
           data=data,restrictions=restrictions,
           covData=covData,
           predStructure=predStructure,
           predType=predType,    # c("c," f", "r")
           transformedParameters=transformedParameters,
           corProbit=corProbit,
           hyperprior=hyperprior,
           n.iter=n.iter,
           n.adapt = n.adapt,
           n.burnin=n.burnin,
           n.thin=n.thin,
           n.chains=n.chains,
           dic =dic,
           ppp = ppp,
           modelfilename=modelfilename,
           parEstFile=parEstFile,
           posteriorFile=posteriorFile,
           autojags=autojags,
           call = match.call(),
           ...)
}

