
# write summary to file
writeSummary <- function(fittedModel, parEstFile=NULL){

  if(!(missing(parEstFile) || is.null(parEstFile))){
    sink(file = parEstFile,type = "o")
    try({
      print(summary(fittedModel))

#       cat("\n\n#################################\n#### Group Parameter Estimates\n")
#       print(fittedModel$summary$groupParameters)

      cat("\n\n#################################\n#### Individual Parameter Estimates\n")
      printIndividualPar(fittedModel$summary$individParameters)

      if(!is.null(fittedModel$summary$transformedParameters)){
        cat("\n\n#################################\n#### Transformed Parameters (Group level)\n")
        print(fittedModel$summary$transformedParameters)
      }


      cat("\n\n#################################\n#### Model information\n")
      print(fittedModel$mptInfo)
    })
    sink()
  }
}


# print array of individual estimates
printIndividualPar <- function(array){
  dd <- dim(array)
  par <- dimnames(array)[[1]]
  for(i in 1:dd[1]){
    cat("Parameter ", par[i], "\n")
    print(array[i,,])
  }
}
