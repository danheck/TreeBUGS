#' Get Parameter Posterior Estimates
#'
#' Returns posterior statistics (e.g., mean, median) for the parameters of a hierarchical MPT model.
#'
#' @param fittedModel a fitted latent-trait MPT model (see \code{\link{traitMPT}}) or beta MPT model (see \code{\link{betaMPT}})
#' @param parameter which parameter(s) should be returned? (see below for details)
#' @param posterior whether to show the posterior \code{"Mean"}, \code{"50\%"} (median), or \code{"SD"}
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
#' Note that this function is only a wrapper to conveniently access the information stored in \code{fittedModel$mcmc.summ}
#' @seealso \code{\link{traitMPT}}, \code{\link{betaMPT}}
#' @author Daniel Heck
#' @export
getParam <- function(fittedModel, parameter="mean", posterior="Mean"){
  if(! class(fittedModel) %in% c("betaMPT", "traitMPT"))
    stop("Only for hierarchical MPT models (see ?traitMPT & ?betaMPT).")

  thetaUnique <- fittedModel$mptInfo$thetaUnique
  S <- length(thetaUnique)
  summ <- fittedModel$mcmc.summ # summary(fittedModel$mcmc)
  allnam <- rownames(summ)
  select <- setdiff(grep(parameter,allnam) , grep(".pred",allnam))
  if(length(select) == 0)
    stop("parameter not found.")
  par <-  summ[select, posterior]
  if(length(par) == S){
    names(par) <- paste0(names(par), "_", thetaUnique)
  }else if(parameter == "theta"){
    par <- matrix(par, ncol=S, byrow=TRUE)
    colnames(par) <- thetaUnique
  }else if(parameter == "rho"){
    par <- getRhoMatrix (thetaUnique, par)
  }

  par
}

