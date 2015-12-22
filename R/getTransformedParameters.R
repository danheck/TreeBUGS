

# get quantities of interest / transformed parameters

# returns the vector of parameters that needs to be sampled and the appropriate model string for JAGS

# model: either "betaMPT" or "traitMPT"
# thetaNames: matrix which assigns parameter labels (first column) to the parameter index (second column)
# transformedParameters: either string with file location or list with transformed parameters
getTransformed <- function (model,
                            thetaNames,
                            transformedParameters=NULL){

  if(is.null(transformedParameters)){
    return(list(transformedParameters=NULL,
         modelstring="\n### No tranformed parameters specified # \n"))
  }else{
    if(is.character(transformedParameters)){
      # read file:

    }else if(!is.list(transformedParameters)){
      warning("The argument 'transformedParameters' must either be a list of parameter transformations or the path to such a file.")
    }
    splitEqual <- sapply(transformedParameters,strsplit, split="=", fixed=TRUE)
    pars <- sapply(splitEqual, function(x) x[1])

    S <- length(pars)
    selCriticalName <- pars %in% c(thetaNames,
                                   "mu", "sd","mu", "sigma",
                                   "beta","alpha", "rho","theta","xi")
    if(any(selCriticalName)){
      error <- paste0("Use different label for transformed parameters:\n  ",
                      paste(pars[selCriticalName], collapse=", "))
      stop(error)
      }

    if(length(unique(pars)) != S){
      stop("The argument 'transformedParameters' does not specifcy unique names for the transformed parameters")
    }

    modelstring <- "### Transformed Parameters (on group level) ###\n"
    for(i in 1:S){
      replacedString <- splitEqual[[i]][2]
      for(k in 1:nrow(thetaNames)){
        replacedString <- gsub(pattern = thetaNames[k,1],
                                 replacement = paste0(
                                   "XXXXXXXXXXXXXX[",
                                   thetaNames[k,2],"]"),
                                 x = replacedString)
      }
      # test whether transformed parameters are proper function: (not working at the moment)
      test <- try(eval(parse(text = replacedString),
               list("XXXXXXXXXXXXXX"=runif(nrow(thetaNames)))), silent=TRUE
        )
      if(class(test) == "try-error"){
        error <- paste0("Check transformedParameter: ", pars[i],
                        ".\n  Function may contain an invalid equation or unknown model parameters.",
                        "\n  Currently, it is defined as: \n  ",
                        gsub("XXXXXXXXXXXXXX","mean", replacedString) )
        warning(error)
      }
      modelstring <- paste(modelstring, "\n",
                           pars[i], "<-", replacedString)
    }
    modelstring <- paste(modelstring, "\n")
    modelstring <- gsub("XXXXXXXXXXXXXX","mean", modelstring)
    return(list(transformedParameters=pars, modelstring=modelstring))
  }
}
