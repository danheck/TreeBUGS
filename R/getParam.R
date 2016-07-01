#' Get Parameter Posterior Statistics
#'
#' Returns posterior statistics (e.g., mean, median) for the parameters of a hierarchical MPT model.
#'
#' @param fittedModel a fitted latent-trait MPT model (see \code{\link{traitMPT}}) or beta MPT model (see \code{\link{betaMPT}})
#' @param parameter which parameter(s) should be returned? (see below for details)
#' @param statistic whether to get the posterior \code{"mean"}, \code{"median"}, \code{"sd"}, or \code{"summary"} (includes mean, SD, and 95\% credibility interval)
#' @param file filename to export results in .csv format (e.g., \code{file="est_param.csv"})
#'
#' @details This function is a convenient way to get the information stored in \code{fittedModel$mcmc.summ}.
#'
#' The latent-trait MPT includes the following parameters:
#' \itemize{
#' \item \code{"mean"} (group means on probability scale)
#' \item \code{"mu"} (group means on probit scale)
#' \item \code{"sigma"} (SD on probit scale)
#' \item \code{"rho"} (correlations on probit scale)
#' \item \code{"theta"} (individual MPT parameters)
#' }
#'
#' The beta MPT includes the following parameters:
#' \itemize{
#' \item \code{"mean"} (group means on probability scale)
#' \item \code{"sd"} (SD on probability scale)
#' \item \code{"alph"},\code{"bet"} (group parameters of beta distribution)
#' \item \code{"theta"} (individual MPT parameters)
#' }
#'
#' @examples
#' \dontrun{
#' # mean estimates per person:
#' getParam(fittedModel, parameter = "theta")
#'
#' # save summary of individual estimates:
#' getParam(fittedModel, parameter = "theta",
#'          statistic = "summary", file= "ind_summ.csv")}
#' @author Daniel Heck
#' @seealso \code{\link{getGroupMeans}} mean group estimates
#' @export
getParam <- function(fittedModel,
                     parameter="mean",
                     statistic="mean",
                     file = NULL){

  if(! class(fittedModel) %in% c("betaMPT", "traitMPT"))
    stop("Only for hierarchical MPT models (see ?traitMPT & ?betaMPT).")

  thetaUnique <- fittedModel$mptInfo$thetaUnique
  S <- length(thetaUnique)
  summ <- fittedModel$mcmc.summ # summary(fittedModel$mcmc)
  allnam <- rownames(summ)
  select <- setdiff(grep(parameter,allnam) , grep(".pred",allnam))
  if(length(select) == 0)
    stop("parameter not found.")

  label <- c("Mean", "SD", "2.5%", "97.5%")
  sel.stat <- switch(statistic,
                     "mean" = "Mean",
                     "sd" = "SD",
                     "median" = "50%",
                     "summary" = label,
                     stop("statistic not supported."))
  par <-  summ[select, sel.stat, drop = FALSE]

  if(statistic != "summary"){
    if(length(par) == S){
      names(par) <- paste0(names(par), "_", thetaUnique)
    }else if(parameter == "theta"){
      par <- matrix(par, ncol=S, byrow=TRUE)
      colnames(par) <- thetaUnique
    }else if(parameter == "rho"){
      par <- getRhoMatrix (thetaUnique, par)
    }
  }else{
    if(length(select) == S){
      rownames(par) <- paste0(rownames(par), "_", thetaUnique)
    }else if(parameter == "theta"){
      par <- matrix(t(par), ncol=S*4, byrow=TRUE)
      colnames(par) <- paste0( rep(thetaUnique, each=4), "_" , label)
    }
  }

  if(!is.null(file)){
    if(is.null(dim(par))) par <- t(par)
    write.csv(par, file = file, row.names = FALSE)
  }

  par
}

