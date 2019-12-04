
#' Between-Subject Comparison of Parameters
#'
#' Computes differencesor other statistics of MPT parameters for two hierarchical MPT models fitted separately to between-subjects data
#'
#' @param model1 fitted hierarchical MPT model for first between-subjects condition
#' @param model2 fitted hierarchical MPT model for second between-subjects condition
#' @param par1 label of parameter from first model for which statistic should be computed
#' @param par2 label of parameter from second model. Default: The same parameter as in the first model
#' @param stat one or more functions of the parameters using \code{"x"} and \code{"y"} as placeholders for the parameters from the first and second model, respectively. Default: Compute (A) the difference between parameters and (B) a Bayesian p-value (by counting how often x<y).
#' @param plot whether to plot the convergence of the difference in parameters
#'
#' @return a list of the class \code{betweenMPT} with the values:
#' \itemize{
#'  \item \code{summary}: Summary for parameter difference
#'  \item \code{mptInfo1}, \code{mptInfo2}: info about MPT models (eqn and data file etc.)
#'  \item \code{mcmc}: the MCMC samples of the differences in parameters
#' }
#' @author Daniel Heck
#' @export
#' @importFrom coda mcmc.list mcmc
betweenSubjectMPT <- function(model1, model2,
                              par1, par2=par1,
                              stat = c("x-y","x<y"),
                              # level = "group",
                              plot = FALSE){

  ################ check input

  if(! inherits(model1, c("betaMPT","traitMPT", "simpleMPT")) ||
     ! inherits(model2, c("betaMPT","traitMPT", "simpleMPT")))
    stop("The two models 'model1' and 'model2' must be (hierarchical) MPT models",
         "\n  of the class (both traitMPT / betaMPT / simpleMPT).")
  if(missing(par1) || length(par1) != 1 || !is.character(par1))
    stop("'par1' must be a single par1 label!")
  if(length(par2) != 1 || !is.character(par2))
    stop("'par2' must be a single par1 label!")
  if( any(model1$mptInfo$thetaUnique !=  model2$mptInfo$thetaUnique))
    warning("Hierarchical MPT models have different sets of parameters.")
  if(! par1 %in%  model1$mptInfo$thetaUnique |
     ! par2 %in%  model2$mptInfo$thetaUnique)
    stop("MPT parameters not found in models. Must be in the set:\n  ",
         paste0(union(model1$mptInfo$thetaUnique,
                      model2$mptInfo$thetaUnique), collapse = "; "))


  ################ get MCMC samples of both models
  idx1 <- match(par1, model1$mptInfo$thetaUnique)
  idx2 <- match(par2, model2$mptInfo$thetaUnique)
  if (length(model1$mptInfo$thetaUnique) > 1)
    nam1 <- paste0("mean[", idx1,"]")
  else
    nam1 <- "mean"
  if (length(model2$mptInfo$thetaUnique) > 1)
    nam2 <- paste0("mean[", idx2,"]")
  else
    nam2 <- "mean"

  pp1 <- model1$runjags$mcmc[,nam1, drop=FALSE]
  pp2 <- model2$runjags$mcmc[,nam2, drop=FALSE]

  namDiff <- gsub("x",  paste0(par1,".m1"), stat, fixed = TRUE)
  namDiff <- gsub("y",  paste0(par2,".m2"), namDiff, fixed = TRUE)

  ################ check MCMC iterations
  if(length(pp1) != length(pp2)){
    warning("Different number of chains for the two models. Some chains are dropped!")
    chains <- min(length(pp1), length(pp2))
    pp1 <- pp1[1:chains]
    pp2 <- pp2[1:chains]
  }

  if(nrow(pp1[[1]]) != nrow(pp2[[1]])){
    warning("Different numbers of iterations for both models. Some iterations are dropped!")
    nn <- min(nrow(pp1[[1]]), nrow(pp2[[1]]))
    for(mm in 1:length(pp1)){
      round(seq(1, nrow(pp1[[mm]]), length.out = nn))
      pp1[[mm]] <- pp1[[mm]][round(seq(1, nrow(pp1[[mm]]), length.out = nn)),,drop = FALSE]
      pp2[[mm]] <- pp2[[mm]][round(seq(1, nrow(pp2[[mm]]), length.out = nn)),,drop = FALSE]
    }
  }

  ################ compute differences
  res <- mcmc.list()
  for(mm in 1:length(pp1)){
    res[[mm]] <- mcmc(sapply(stat, function(ss)
      eval(parse(text = ss),
           envir = list(x = pp1[[mm]],
                        y = pp2[[mm]]))))
    colnames(res[[mm]]) <- namDiff
  }
  if(plot)
    plot(res)

  summ <- summarizeMCMC(res)
  out <- list(summary = summ,
              mptInfo1 = model1$mptInfo,
              mptInfo2 = model2$mptInfo,
              mcmc = res)
  class(out) <- "betweenMPT"
  out
}


#' @export
print.betweenMPT <- function(x, round = 3, ...){
  print(round(x$summary,round))
}




############### OLD CODE (differences for individuals, meaningless)


# }else{
#   if(N1 != N2)
#     stop("Sample sizes of participants do not match for both models.",
#          "\n  Use 'level=\"group\" ' to compare the parameter on the mean level.")
#   nam1 <- paste0("theta[",idx,",",1:N1,"]")
#   nam2 <- paste0("theta[",idx,",",1:N2,"]")
# }




# if(level != "group"){
#   namDiff <- paste0(namDiff, "[",1:N1, "]")
# }
