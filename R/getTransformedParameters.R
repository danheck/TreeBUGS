

# get quantities of interest / transformed parameters
# returns the vector of parameters that needs to be sampled and the appropriate model string for JAGS
# thetaNames matrix which assigns parameter labels (first column) to the parameter index (second column)
# transformedParameters: either string with file location or list with transformed parameters
getTransformed <- function (thetaNames, transformedParameters=NULL){

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

    P <- length(pars)

    if(length(unique(pars)) != P){
      stop("The argument 'transformedParameters' does not specifcy unique names for the transformed parameters")
    }

    modelstring <- "### Transformed Parameters ###\n"
    for(i in 1:P){
      replacedString <- splitEqual[[i]][2]
      for(k in 1:nrow(thetaNames)){
        replacedString <- gsub(pattern = thetaNames[k,1],
                                 replacement = paste0("mnb[",thetaNames[k,2],"]"),
                                 x = replacedString)
      }
      # test whether transformed parameters are proper function: (not working at the moment)
#       try(eval(
#         substitute(eval(replacedString),env=list(mnb=runif(1000)))
#         function(mnb=runif(1000)) replacedString
#         ) )
      modelstring <- paste(modelstring, "\n",
                           pars[i], "<-", replacedString)
    }
    modelstring <- paste(modelstring, "\n")
    return(list(transformedParameters=pars, modelstring=modelstring))
  }
}
