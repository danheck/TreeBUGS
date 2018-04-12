#' Get Posterior Samples from Fitted MPT Model
#'
#' Extracts MCMC posterior samples as an \code{coda::mcmc.list} and relabels the MCMC variables.
#'
#' @inheritParams getParam
#' @inheritParams plotParam
#' @param names whether and how to rename the variables in the MCMC output:
#'     \code{par} (internal parameter labels such as \code{mu[1]}),
#'     \code{label} (MPT label from EQN file such as \code{"d"}), or
#'     \code{par_name} (concatenation of both such as \code{"mu[1]_d"}).
#'
#' @importFrom coda varnames niter nvar
#' @examples
#' \dontrun{
#' getSamples(fittedModel, "mu", select = c("d", "g"))
#' }
#' @export
getSamples <- function (fittedModel, parameter = "mean", select = "all",
                        names = "par_label"){

  parnames <- fittedModel$mptInfo$thetaUnique
  if (missing(select) || select == "all")
    select <- parnames
  else
    if (!all(select %in% parnames))
      stop("Check arguments: Not all parameters in 'select' are included in the MPT model!\n",
           "Parameters are: ", paste(parnames, collapse=", "))

  S <- length(parnames)
  var <- ""
  # matches <- grep(parameter, varnames(fittedModel$runjags$mcmc), value = TRUE)
  idx <-  match(select, parnames)
  if (S > 1 || (parameter == "theta" && S == 1))
    var <- paste0("[" , idx, "]")
  else if (S > 1 && parameter == "theta")
    var <- paste0("[", outer(idx, seq(nrow(fittedModel$mptInfo$data)),
                             "paste", sep = ","), "]")
  else if (S > 1 && parameter == "rho")
    var <-  paste0("[", outer(idx, idx, "paste", sep = ","), "]")
  mcmc <- fittedModel$runjags$mcmc[,paste0(parameter, var), drop = FALSE]

  if (parameter != "rho")
    mcmc <- rename_mcmc(mcmc, names, select)
  mcmc
}

rename_mcmc <- function(mcmc, names, parnames){
  if (nvar(mcmc) == length(parnames)){
    if (names == "par_label")
      coda::varnames(mcmc) <- paste0(varnames(mcmc), "_", parnames)
    else if (names == "label")
      coda::varnames(mcmc) <- parnames
  }
  mcmc
}
