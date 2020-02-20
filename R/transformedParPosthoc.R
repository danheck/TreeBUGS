
#' Get Transformed Parameters
#'
#' Computes transformations of MPT parameters based on the MCMC posterior samples
#' (e.g., differences of parameters).
#'
#' @param fittedModel either a fitted latent-trait or beta MPT model
#'     (\code{\link{traitMPT}}, \code{\link{betaMPT}}) or an \code{\link[coda]{mcmc.list}}.

#' @param transformedParameters list with parameter transformations that should
#'     be computed based on the posterior samples (e.g., for testing parameter
#'     differences: \code{list("diffD=Do-Dn")}).
#' @param level whether to compute transformations of \code{"group"} or
#'     \code{"individual"} estimates
#' @param nCPU number of CPU cores across which the MCMC chains are distributed
#'
#' @return an \link[coda]{mcmc.list} of posterior samples for the transformed parameters
#' @examples
#' \dontrun{
#' tt <- transformedParameters(fittedModel,
#'                             list("diff=a-b","p=a>b"),
#'                             level="individual")
#' summary(tt)
#' }
#' @export
transformedParameters <- function (fittedModel, transformedParameters,
                                   level = "group", nCPU = 4){

  if (inherits(fittedModel, c("mcmc", "mcmc.list"))){
    mcmc <- fittedModel
    if (inherits(fittedModel, "mcmc"))
      mcmc <- mcmc.list(mcmc)
    S <- N <- P <- parsed <- level <- NA  # only for cpu > 1

    # mm <- fittedModel[[1]]
    getTrans <- function(mm){

      new.mcmc <- matrix(NA, nrow(mm), length(transformedParameters))
      newnames <- vapply(unlist(transformedParameters), function(tp)
        strsplit(tp, "=", fixed = TRUE)[[1]][1], "a")
      colnames(new.mcmc) <- gsub(" ", "", newnames, fixed = TRUE)

      for(i in seq_along(transformedParameters)){
        select <- which(sapply(colnames(mm), grepl, x = transformedParameters[[i]], fixed = TRUE))
        for(ii in 1:nrow(mm)){
          parsed <- transformedParameters[[i]]
          for(ss in 1:length(select)){
            parsed <- gsub(names(select[ss]), mm[ii,select[ss]], parsed, fixed = TRUE)
          }
          new.mcmc[ii, i] <- eval(parse(text = parsed))
        }
      }
      as.mcmc(new.mcmc)
    }
    tmp <- tryCatch(getTrans(mcmc[[1]][1,,drop =FALSE]),
                    error = function(e) stop("Check definition of transformations! \n", e))

  } else {
    mcmc <- fittedModel$runjags$mcmc
    parsed <- getTransformed(fittedModel$mptInfo$thetaNames,
                             transformedParameters=transformedParameters,
                             mergeString=FALSE)
    S <- length(fittedModel$mptInfo$thetaUnique)
    P <- length(parsed$transformedParameters)
    N <- nrow(fittedModel$mptInfo$data)

    samp <- list()
    #### for a single MCMC chain:
    getTrans <- function (mm){
      if (level == "group"){
        if (S == 1){
          sel.mean <- "mean"
        } else {
          sel.mean <- paste0("mean[",1:S,"]")
        }
        samp <- mm[,rep(1,P), drop=FALSE]
        colnames(samp) <- parsed$transformedParameters
        for(i in 1:P){
          samp[,i] <- apply(mm[,sel.mean,drop=FALSE], 1,
                            function(xx) eval(parse(text=parsed$modelstring[i]),
                                              list(mean=unlist(xx))))
        }
      }else{
        samp <- mm[,rep(1,P*N),drop=FALSE]
        for(i in 1:P){
          # if(S==1){
          #   sel.theta <- paste0("theta[")
          # }else{
          sel.theta <- paste0("theta[",1:S,",")
          # }
          idx <- (i-1)*N+1:N
          colnames(samp)[idx] <- paste0(parsed$transformedParameters[i],"[", 1:N,"]")
          for(n in 1:N){
            samp[,idx[n]] <- apply(mm[,paste0(sel.theta,n,"]"),drop=FALSE], 1,
                                   function(xx) eval(parse(text=parsed$modelstring[1]),
                                                     list(mean=unlist(xx))))
          }
        }
      }
      samp
    }
  }

  if (nCPU > 1){
    cl <- makeCluster(nCPU)
    clusterExport(cl, c("S","N","P","parsed",
                        "transformedParameters","level"), envir=environment())
    # loop across MCMC chains:
    samp <- parLapply(cl, mcmc, getTrans)
    stopCluster(cl)
  } else {
    samp <- lapply(mcmc, getTrans)
  }

  as.mcmc.list(samp)
}

