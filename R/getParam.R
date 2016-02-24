#' Get Parameter Posterior Estimates
#'
#' Returns posterior statistics (e.g., mean, median) for the parameters of a hierarchical MPT model.
#'
#' @param fittedModel a fitted latent-trait MPT model (see \code{\link{traitMPT}}) or beta MPT model (see \code{\link{betaMPT}})
#' @param parameter which parameter(s) should be returned? (see below for details)
#' @param posterior whether to show the posterior \code{"mean"}, \code{"median"}, or \code{"sd"}
#' @details
#' In the latent-trait MPT, the following parameters are being estimated:
#' \itemize{
#' \item \code{"mean"} (group means on probability scale)
#' \item \code{"mu"} (group means on probit scale)
#' \item \code{"sigma"} (SD on probit scale)
#' \item \code{"rho"} (correlations on probit scale)
#' \item \code{"theta"} (individual MPT parameters)
#' }
#'
#' In the beta MPT, the following parameters are being estimated:
#' \itemize{
#' \item \code{"mean"} (group means on probability scale)
#' \item \code{"sd"} (SD on probability scale)
#' \item \code{"alph"},\code{"bet"} (group parameters of beta distribution)
#' \item \code{"theta"} (individual MPT parameters)
#' }
#'
#' Note that this function is only a wrapper to conveniently access the information stored in \code{fittedModel$mcmc$BUGSoutput$mean}
#' @seealso \code{\link{traitMPT}}, \code{\link{betaMPT}}
#' @author Daniel Heck
#' @export
getParam <- function(fittedModel, parameter="mean", posterior="mean"){
  if(! class(fittedModel) %in% c("betaMPT", "traitMPT"))
    stop("Only for hierarchical MPT models (see ?traitMPT & ?betaMPT).")

   fittedModel$mcmc$BUGSoutput[[posterior]][[parameter]]

}

