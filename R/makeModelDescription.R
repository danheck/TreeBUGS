# #' Make the model description file needed by JAGS / WinBUGS etc

# #' @param model either "betaMPT" or "traitMPT"
# #' @param filename The name that the model description file should get
# #' @param mergedTree Output of mergeBranches()
# #' @param S The number of thetas
# #' @param hyperprior named list, either with entries "alpha" and "beta" for the "betaMPT", or with "mu" and "xi" for the latent-trait MPT
# #' @param parString  additional lines added to sample transformed parameters (e.g., parameter differences for testing)
# #' @author Daniel Heck
# #' @export
makeModelFile <-function(model, # either "betaMPT" or "traitMPT"
                         filename, mergedTree, S,
                         hyperprior,       # list (either with alpha+beta or with mu+xi)
                         corString=NULL,   # model string to compute correlations
                         predString=NULL,  # model string to include predictors (in traitMPT)
                         parString="",
                         groupMatT1=NULL,    # a G x 2 matrix that contains the grouping indices for T1 per group (column 1:2 = from:to)
                         fixedPar=NULL
){

  treeNames <- as.character(sort(unique(mergedTree$Tree)))
  NOT=length(treeNames)
  count=1
  # number of categories per tree:
  ncatPerTree <- as.vector(by(mergedTree$Category, mergedTree$Tree, length))


  #  ###################################### MODEL ###################
  cat(ifelse(model=="traitMPT",
             "####### Hierarchical latent-trait MPT with TreeBGUS #####\n\n",
             "####### Hierarchical beta MPT with TreeBGUS #####\n\n"),
      "data{\nfor(s in 1:S){\n  zeros[s] <- 0\n}\n}\n\n",     # for zeros-trick and dmnorm(zeros, Tau.prec)
      "model{\n\n",
      "for (n in 1: subjs){\n\n### MPT model equations:\n",
      file=filename)

  for(i in 1:NOT){
    for(j in 1:ncatPerTree[i]){
      cat(treeNames[i],"[n,",j,"] <- ",mergedTree$Equation[count],
          "\n",sep="",file=filename,append=T)
      count=count+1
    }
    cat("\n",file=filename,append=T)
  }

  cat("\n### Multinomial likelihood:\n",file=filename,append=T)

  for(i in 1:NOT){
    cat("response.",treeNames[i],"[n,1:",ncatPerTree[i],"] ~ dmulti(",
        treeNames[i],"[n,1:",ncatPerTree[i],"],items.",treeNames[i],
        "[n])\n",sep="",file=filename,append=T)
  }

  cat("}\n",file=filename,append=T)

  ######################################### MODEL SPECIFIC HYPERPRIOR PART ###########

  hyperprior <- switch(model,
                       "betaMPT" = makeBetaHyperprior(S =S ,
                                                      alpha = hyperprior$alpha,
                                                      beta = hyperprior$beta),
                       "traitMPT" = makeTraitHyperprior(S = S,
                                                        predString = predString,
                                                        mu = hyperprior$mu,
                                                        xi = hyperprior$xi,
                                                        wishart = !anyNA(hyperprior$V))
  )

  cat("\n\n### Hierarchical structure:",
      hyperprior, file=filename, append=TRUE)
  if(! is.null(fixedPar)){
    S.fixed <- length(unique(fixedPar$theta))
    cat("\nfor(i in 1:", S.fixed,"){\n",
        "      thetaFE[i] ~ dunif(0,1)\n}\n", file=filename, append=TRUE)
  }

  ######################################### END OF MODEL SPECIFIC HYPERPRIOR PART #####

  if( !is.null(corString)){
    cat(corString, file=filename, append=T)
  }

  ### Transformed parameters:
  cat(parString, file=filename, append = TRUE)

  cat("}\n",file=filename,append=T)
}



################### Beta-MPT specific hyperprior part
makeBetaHyperprior <- function(S, alpha = "dunif(1,5000)", beta = "dunif(1,5000)"){

  if(class(alpha) != "character" || !length(alpha) %in% c(1,S)){
    stop("Hyperprior for 'alpha' must be a character vector of length 1
     if the same prior should be used for all MPT parameters (default)
     or a vector of the same length as the number of parameters (=", S, "; to
     check the order see ?readEQN).")
  }
  if(class(beta) != "character" || !length(beta) %in% c(1,S)){
    stop("Hyperprior for 'beta' must be a character vector of length 1
     if the same prior should be used for all MPT parameters (default)
     or a vector of the same length as the number of parameters (=", S, "; to
     check the order see ?readEQN).")
  }

  modelString <- paste0("

for(s in 1:S){
  for(n in 1:subjs) {
    theta[s,n] ~ dbeta(alph[s], bet[s])
  }
}

for(s in 1:S){
  mean[s] <- alph[s]/(alph[s]+bet[s])
  sd[s] <- sqrt(alph[s]*bet[s]/(pow(alph[s]+bet[s],2)*(alph[s]+bet[s]+1)))
}\n\n")

  for(s in 1:S){
    if(alpha[s] == "zero" | beta[s] == "zero"){
      modelString <- paste0(modelString,
                            "alph[", s, "] ~ dunif(.01,5000)\n",
                            "bet[", s, "] ~ dunif(.01,5000)\n",
                            "zeros[",s,"] ~ dpois(phi[",s,"])\n",
                            "phi[",s,"] <- -log(1/pow(alph[",s,"]+bet[",s,"],5/2))\n")
    }else{
      modelString <- paste0(modelString,
                            "alph[", s, "] ~ ", alpha[s], "\n",
                            "bet[",  s, "] ~ ", beta[s],  "\n")
    }
  }
  modelString <- paste0(modelString, "\n")

  return(modelString)
}


################### latent-trait-MPT specific hyperprior part
makeTraitHyperprior <- function(S, predString, mu = "dnorm(0,1)",
                                xi = "dunif(0,100)", wishart = TRUE){

  if(class(mu) != "character" || !length(mu) %in% c(1,S)){
    stop("Hyperprior for 'mu' must be a character vector of length 1
     if the same prior should be used for all MPT parameters (default)
     or a vector of the same length as the number of parameters (=", S, "; to
     check the order see ?readEQN).")
  }
  if(class(xi) != "character" || !length(xi) %in% c(1,S)){
    stop("Hyperprior for 'xi' must be a character vector of length 1
     if the same prior should be used for all MPT parameters (default)
     or a vector of the same length as the number of parameters (=", S, "; to
     check the order see ?readEQN).")
  }

  if (wishart){

    modelString <- paste0(predString, "

    # hyperpriors
    for(i in 1:subjs) {
      delta.part.raw[1:S,i] ~ dmnorm(zeros,T.prec.part[1:S,1:S])
    }

    ",
                          ######################## special case if S=1:
                          ifelse(S > 1,"
    T.prec.part[1:S,1:S] ~ dwish(V, df)",
                                 "
    T.prec.part[1,1] ~ dchisq(df)"),

                          "
    Sigma.raw[1:S,1:S] <- inverse(T.prec.part[,])
    for(s in 1:S){
      mean[s] <- phi(mu[s])
      for(q in 1:S){
        Sigma[s,q] <- Sigma.raw[q,s]*xi[s]*xi[q]
      }
    }

    for(s in 1:S){
      for(q in 1:S){
        # Off-diagonal elements of S (correlations not affected by xi)
        rho[s,q] <- Sigma[s,q]/sqrt(Sigma[s,s]*Sigma[q,q])
      }
      # Diagonal elements of S (rescale sigma)
      sigma[s] <- sqrt(Sigma[s,s])
    }")
  } else {
    modelString <- paste0(predString, "

    # hyperpriors
    for(s in 1:S) {
      for(i in 1:subjs) {
        delta.part.raw[s,i] ~ dnorm(0,tau[s])
      }

      mean[s] <- phi(mu[s])
      tau[s] ~ dgamma(.5, df / 2)  # = chi-square
      sigma[s] <- abs(xi[s]) / sqrt(tau[s])
      for(s2 in 1:S){
        rho[s,s2] <- -99
      }
    }")
  }

  paste0(modelString, "\n\n",
         paste0("\nmu[", 1:S, "] ~ ", mu, collapse = ""),
         paste0("\nxi[", 1:S, "] ~ ", xi, collapse = ""),"\n\n")
}






############################## Posterior predictive checks:
#
# 	cat("\n############# T1: posterior predictive check of mean frequencies\n", file=filename,append=T)
#
#   cat("for(n in 1:subjs) {\n",file=filename,append=T)
#   for(i in 1:NOT){
#     # sample predicted: frequencies from posterior
#     cat("response.",treeNames[i],".pred[n,1:",ncatPerTree[i],"]~dmulti(",
#         treeNames[i],"[n,1:",ncatPerTree[i],"],items.",treeNames[i],
#         "[n])\n",sep="",file=filename,append=T)
#
#     # compute expected frequencies:
#     cat("n.expected.",treeNames[i],"[n,1:",ncatPerTree[i],"] <- ",
#         treeNames[i],"[n,1:",ncatPerTree[i],"] * items.",treeNames[i],
#         "[n]\n",sep="",file=filename,append=T)
#
#     # get means of expected and predicted frequencies:
#
#   }
#
#   cat("\n\n### T1statistic for individual data:\n",
#       "T1ind.obs[n] <- 0", file=filename,append=T )
#   for(i in 1:NOT){
#     cat("+ sum( (response.",treeNames[i],"[n,] - n.expected.",treeNames[i],"[n,])^2/n.expected.",treeNames[i],"[n,])",
#         sep="", file=filename,append=T)
#   }
#
#   # T1 for predicted means:
#   cat("\nT1ind.pred[n] <- 0", file=filename,append=T )
#   for(i in 1:NOT){
#     cat("+ sum( (response.",treeNames[i],".pred[n,] - n.expected.",treeNames[i],
#         "[n,])^2/n.expected.",treeNames[i],"[n,])", sep="", file=filename,append=T)
#   }
#
#   cat("\np.T1ind[n] <- T1ind.pred[n] > T1ind.obs[n]\n",
#       "}\n",file=filename,append=T)
#
#
#
#   cat("\n\n### T1 statistic for aggregated data:\n", file=filename,append=T)
#   for(i in 1:NOT){
#     cat("for(k in 1:",ncatPerTree[i],") {\n",
#         "n.expected.",treeNames[i],".mean[k] <- mean(n.expected.",treeNames[i],"[,k])\n",
#         "response.",treeNames[i],".pred.mean[k] <- mean(response.",treeNames[i],".pred[,k])\n",
#         sep="",file=filename,append=T)
#
#     # T1 per group
# #     for(s in 1:S){
# #       pred.group.mean[s,k] <- mean(pred[study.idx[s,1]:study.idx[s,2],k])
# #     }
#     if(!is.null(groupMatT1)){
#       cat("for(g in 1:",nrow(groupMatT1),") {\n",
#
#           "group.n.exp.",treeNames[i],".mean[g,k] <- mean(n.expected.",treeNames[i],  # expected per group
#           "[groupMatT1[g,1:NgroupT1[g]],k])\n",
#           "group.resp.",treeNames[i],".pred.mean[g,k] <- mean(response.",treeNames[i],  # sampled per group
#           ".pred[groupMatT1[g,1:NgroupT1[g]],k])\n",
#
#           "}\n", sep="",file=filename,append=T)
#     }
#     cat("}\n", file=filename,append=T )
#   }
#   if(!is.null(groupMatT1)){
#     # T1 for observed means:
#     cat("\n###### T1 (per group)\n",
#         "for(g in 1:",nrow(groupMatT1),") {\n",
#         "T1.group.obs[g] <- 0", file=filename,append=T )
#     for(i in 1:NOT){
#       cat("+ sum( (group.resp.",treeNames[i],".mean[g,] - group.n.exp.",treeNames[i],".mean[g,])^2/",
#           "group.n.exp.",treeNames[i],".mean[g,])",
#           sep="", file=filename,append=T)
#     }
#
#     # T1 for predicted means:
#     cat("\nT1.group.pred[g] <- 0", file=filename,append=T )
#     for(i in 1:NOT){
#       cat("+ sum( (group.resp.",treeNames[i],".pred.mean[g,] - group.n.exp.",treeNames[i],".mean[g,])^2/",
#           "group.n.exp.",treeNames[i],".mean[g,])",
#           sep="", file=filename,append=T)
#     }
#     cat("\n## T1 (per group) comparison:\n",
#         "p.T1.group[g] <- T1.group.pred[g] > T1.group.obs[g]\n}\n",  file=filename,append=T)
#   }
#
#   # T1 for observed means:
#   cat("T1.obs <- 0", file=filename,append=T )
#   for(i in 1:NOT){
#     cat("+ sum( (response.",treeNames[i],".mean - n.expected.",treeNames[i],".mean)^2/n.expected.",treeNames[i],".mean)",
#         sep="", file=filename,append=T)
#   }
#
#   # T1 for predicted means:
#   cat("\nT1.pred <- 0", file=filename,append=T )
#   for(i in 1:NOT){
#     cat("+ sum( (response.",treeNames[i],".pred.mean - n.expected.",treeNames[i],".mean)^2/n.expected.",treeNames[i],".mean)", sep="", file=filename,append=T)
#   }
#   cat("\np.T1 <- T1.pred > T1.obs\n",  file=filename,append=T)



# ############# T2: posterior predictive check of covariance structure
#
# # covariance of expected frequencies: math!
#
# for(i in 1:N){
#   # multinomial variance and covariance:
#   for(k in 1:4){
#     for(l in 1:4){
#       cov.pers[i,k,l] <- ifelse(k==l, J[i,1]*p.g[i,k]*(1-p.g[i,k]),-J[i,1]*p.g[i,k]*p.g[i,l])
#     }
#     for(l in 5:16){
#       cov.pers[i,k,l] <- 0
#       cov.pers[i,l,k] <- 0
#     }
#   }
#   for(k in 5:8){
#     for(l in 5:8){
#       cov.pers[i,k,l] <- ifelse(k==l, J[i,1]*p.k[i,k-4]*(1-p.k[i,k-4]),-J[i,1]*p.k[i,k-4]*p.k[i,l-4])
#     }
#     for(l in 9:16){
#       cov.pers[i,k,l] <- 0
#       cov.pers[i,l,k] <- 0
#     }
#   }
#   for(k in 9:16){
#     for(l in 9:16){
#       cov.pers[i,k,l] <- ifelse(k==l, J[i,1]*p.r[i,k-8]*(1-p.r[i,k-8]),-J[i,1]*p.r[i,k-8]*p.r[i,l-8])
#     }
#   }
# }
#
# for(k in 1:16){
#   for(l in 1:16){
#     cov.exp[k,l] <- inprod(n.exp[,k]-mean(n.exp[,k]), n.exp[,l]-mean(n.exp[,l]) )/N + (N-1)/N*mean(cov.pers[,k,l])
#     cov.pred[k,l] <- inprod(pred[,k]-mean(pred[,k]), pred[,l]-mean(pred[,l]) )/N
#
#     dev.obs[k,l] <-  (cov.obs[k,l]-cov.exp[k,l])^2/sqrt(cov.exp[k,k]*cov.exp[l,l])
#     dev.pred[k,l] <- (cov.pred[k,l]-cov.exp[k,l])^2/sqrt(cov.exp[k,k]*cov.exp[l,l])
#   }
# }
# T2.obs <- sum(dev.obs)
# T2.pred <- sum(dev.pred)
# p.T2 <- T2.pred>T2.obs



