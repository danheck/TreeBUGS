# add model string that computes correlation of MPT model parameters with external covariate
# does not change estimation of the model!

covStringCorrelation <- function(covTable, corProbit=FALSE){

  checkCorrelations <-  any(covTable$predType != "c")
  if(checkCorrelations)
    stop("To compute correlations, only continuous covariates are allowed.")

  modelString <- "\n## Covariate Handling: Correlations ##\n"

  if(!is.null(covTable)){

    # use probit-transformed values (e.g., in traitMPT model)
    corProbitString <- ifelse(corProbit, "theta.probit", "theta")

    pars <- unique(covTable$Parameter)
    modelString <- paste0(modelString,
                          "for(i in 1:S){",
                          ifelse(corProbit,
                                 "\n  for(n in 1:subjs){\n    theta.probit[i,n] <- probit(theta[i,n])\n  }" ,
                                 ""),

                          "\nthetaSD[i] <- sd(",  corProbitString,
                          "[i,])
  }
                          ")
    cnt <- 0
    for(pp in 1:length(pars)){
      sel <- covTable$Parameter == pars[pp]
      thetaIdx <- covTable$theta[sel][1]
      covs <- covTable$Covariate[sel]
      for(cc in 1:length(covs)){
        covIdx <-  covTable$covIdx[sel][cc]

        ### requires commputation of correlation using only sum / sd
        #  (par-mean(par))*(cov-mean(cov))/(sd(par) * sd(cov))

        modelString <- paste0(modelString,
                              "cor_", pars[pp],"_",covs[cc],
                              # " <- mean( (theta[",thetaIdx, ",]-mean(theta[",thetaIdx,
                              " <- mean( (",corProbitString ,"[",thetaIdx, ",]- mean(",
                              corProbitString,"[",thetaIdx,
                              ",]))*(covData[,",covIdx, "]-mean(covData[,",covIdx,
                              "])) )/covSD[",covIdx,
                              "] / thetaSD[",thetaIdx,"]\n")
        # nice, not supported: cor(theta[",thetaIdx, ",], covData[,",covIdx, "])
        covTable$covPar[cnt<-cnt+1] <-  paste0("cor_", pars[pp],"_",covs[cc])
      }
    }
  }


  # lines in JAGS:
  # additionally monitored variable: covPars <- paste0("cor_", sapply(covList, function(ll, ll$Par) ))
  ###################  ###################   ###################

  return(list(modelString = modelString,
              covPars = covTable$covPar))
}
