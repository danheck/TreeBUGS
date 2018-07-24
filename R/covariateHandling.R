


################### covariate handling of: covData, covStructure   ###################

# test files (Daniel):
# covData <- "../first_tests/covDatFac.csv"
# covData <- "../first_tests/covDat.csv"
# debug(TreeBUGS:::covHandling)



# covHandling: returns a nice table that assigns parameters to covariates (including indices)
#
# thetanames: data frame with assignment of parameters to indices / restrictions etc.
# predType: character vector with values "c"/"f"/"r" for continuous/fixed/random effects
covHandling <- function(covData,
                        covStructure=NULL,
                        N,
                        thetaNames,
                        predType=NULL,
                        defaultExclude=NULL,
                        T1group=NULL){

  # get a clean list
  if(! (is.null(covData) | is.list(covStructure) & length(covStructure) == 0 |
        is.null(covStructure) & "ALL_COVARIATES" %in% defaultExclude)){

    if(!("ALL_COVARIATES" %in% defaultExclude)){     # for correlation structure
      covData <- covData[,predType == "c" &          # select only continuous variables
                           !colnames(covData)%in% defaultExclude, drop=FALSE]   # and exclude predictor variables
    }
    if(ncol(covData) != 0){
      covNames <- colnames(covData)
      if(is.character(covStructure)){
        covStructure <- as.list(read.csv(covStructure, header=F,stringsAsFactors=F, sep="}")$V1)
      }


      if(is.null(covStructure)){
        # default: all covariates for all parameters included
        covStructure <- lapply(unique(thetaNames$Parameter),
                               function(tt) paste(tt, ";", paste(covNames, collapse=" ")))
      }


      # recode factors from character values => integer with number of factor level as index
      covTmp <- covRecodeFactor(covData, predType)
      covData <- covTmp$covData
      predFactorLevels <- covTmp$predFactorLevels
      names(predFactorLevels) <- colnames(covData)

      # set up table and iterate across all combinations of parameters and covariates
      covTable <- data.frame()
      for(i in 1:length(covStructure)){
        covStructure[[i]] <- gsub("\\,", " ", covStructure[[i]])
        sss <- strsplit(covStructure[[i]], ";")[[1]]
        if(length(sss)  != 2)
          stop("Check predStructure (in each argument, exactly one semicolon is required to separate parameters (left hand) and predictors (right hand)):\n",
               covStructure[[i]])
        pars <- strsplit(sss[1], " +")[[1]]
        covs <- strsplit(sss[2], " +")[[1]]
        pars <- pars[pars != ""]
        covs <- covs[covs != ""]

        # parameters:
        for(pp in 1:length(pars)){
          # covariates:
          for(cc in 1:length(covs)){
            covTable <- rbind(covTable, data.frame(Parameter = pars[pp], Covariate = covs[cc]))
          }
        }
      }

      # replace the constrained parameters by free parameters and remove redundant rows
      parSel <- match(covTable$Parameter,  thetaNames$Parameter)
      if(any(is.na(parSel)))
        stop("Check parameter names in predStructure. Problematic right now:\n  ",
             covTable$Parameter[is.na(parSel)])

             # ifelse( any(defaultExclude != "ALL_COVARIATES"),
             #         "\n  (note that for correlations in covStructure, only continuous variables are allowed)",""))
      covTable$theta <- thetaNames$theta[parSel]
      covTable$Parameter <- thetaNames$Parameter[parSel]
      covTable$covIdx <- (1:ncol(covData))[match(covTable$Covariate,  covNames)]
      covTable <- covTable[!duplicated(paste(covTable$theta, covTable$covIdx)),]

      covTable$predType <- predType[covTable$covIdx]


      # get parameter labels for JAGS
      covTable$prefix <- ifelse(covTable$predType=="c", "slope", "factor")
      covTable$covPar <- apply(covTable[,c("prefix", "Parameter", "Covariate")],
                               1, paste, collapse="_")
    }else{
      covTable <- NULL
      predFactorLevels <- NULL
    }
  }else{  #(is.null(covData) || ncol(covData) == 0 || length(covData) == 0)
    # no covariates
    covTable <- NULL
    predFactorLevels <- NULL
  }

  return(list(covTable = covTable,
              covData=covData,
              predFactorLevels = predFactorLevels))
}





