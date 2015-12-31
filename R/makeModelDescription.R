# #' Make the model description file needed by JAGS / WinBUGS etc

# #' @param model either "betaMPT" or "traitMPT"
# #' @param filename The name that the model description file should get
# #' @param mergedTree Output of mergeBranches()
# #' @param S The number of thetas
# #' @param hyperprior named list, either with entries "alpha" and "beta" for the "betaMPT", or with "mu" and "xi" for the latent-trait MPT
# #' @param sampler Which sampler should be used? Default is "JAGS", further options are "OpenBUGS" and "WinBugs"
# #' @param parString  additional lines added to sample transformed parameters (e.g., parameter differences for testing)
# #' @author Daniel Heck, Nina R. Arnold, Denis Arnold
# #' @export
makeModelFile <-function(model, # either "betaMPT" or "traitMPT"
                         filename,
                         mergedTree,
                         S,
                         hyperprior, # list (either with alpha+beta or with mu+xi)
                         covString=NULL,
                         sampler="JAGS",
                         parString=""){

  treeNames <- as.character(sort(unique(mergedTree$Tree)))
  NOT=length(treeNames)
  count=1
  # number of categories per tree:
  ncatPerTree <- as.vector(by(mergedTree$Category, mergedTree$Tree, length))
#   ncatPerTree <- rep(NA, NOT)
#   for(i in 1:NOT) {
#     ncatPerTree[i] <- length(grep(treeNames[i],mergedTree$Tree))
#   }

  if(sampler%in%c("jags","JAGS")){

    ################################### DATA #######################
    cat("data\n",file=filename)
    cat("{\n",file=filename,append=T)
    # cat("for(s in 1:",S,"){\n",sep="",file=filename,append=T)
    # cat("zero[s] <- 0\n",file=filename,append=T)
    # cat("}\n ",file=filename,append=T)


    ####### for posterior predictive cehck: mean frequencies
    for(i in 1:NOT){
      cat("for(k in 1:",ncatPerTree[i],") {\n",sep="",file=filename,append=T)
        cat("response.",treeNames[i],".mean[k] <- mean(response.",treeNames[i],"[,k])\n",sep="",file=filename,append=T)
      cat("}\n", file=filename,append=T)
    }
    cat("}\n ",file=filename,append=T)

    ###################################### MODEL ###################
    cat("model\n",file=filename,append=T)
  }else{
    cat("model\n",file=filename)
  }

	cat("{\n",file=filename,append=T)
	cat("for (n in 1: subjs){\n",file=filename,append=T)


	for(i in 1:NOT){
		for(j in 1:ncatPerTree[i]){
		cat(treeNames[i],"[n,",j,"]<-",mergedTree$Equation[count],"\n",sep="",file=filename,append=T)
		count=count+1
		}
	}

	for(i in 1:NOT){
		cat("response.",treeNames[i],"[n,1:",ncatPerTree[i],"]~dmulti(",
			treeNames[i],"[n,1:",ncatPerTree[i],"],items.",treeNames[i],
			"[n])\n",sep="",file=filename,append=T)
	}

	cat("}\n",file=filename,append=T)

	######################################### MODEL SPECIFIC HYPERPRIOR PART #################################

	hyperprior <- switch(model,
	                     "betaMPT" = makeBetaHyperprior(S,
	                                                    alpha = hyperprior$alpha,
	                                                    beta = hyperprior$beta,
	                                                    sampler = sampler),
	                     "traitMPT" = makeTraitHyperprior(S,
	                                                     mu = hyperprior$mu,
	                                                     xi = hyperprior$xi,
	                                                     sampler = sampler)
	)

	cat(hyperprior, file=filename, append=T)

	######################################### END OF MODEL SPECIFIC HYPERPRIOR PART ##########################


	if(model == "betaMPT" & !is.null(covString)){
	  cat(covString, file=filename, append=T)
	}


	############################## Posterior predictive checks:
	if(sampler%in%c("jags","JAGS")){

	cat("\n############# T1: posterior predictive check of mean frequencies\n", file=filename,append=T)

  cat("for(n in 1:subjs) {\n",file=filename,append=T)
  for(i in 1:NOT){
    # sample predicted: frequencies from posterior
    cat("response.",treeNames[i],".pred[n,1:",ncatPerTree[i],"]~dmulti(",
        treeNames[i],"[n,1:",ncatPerTree[i],"],items.",treeNames[i],
        "[n])\n",sep="",file=filename,append=T)

    # compute expected frequencies:
    cat("n.expected.",treeNames[i],"[n,1:",ncatPerTree[i],"] <- ",
        treeNames[i],"[n,1:",ncatPerTree[i],"] * items.",treeNames[i],
        "[n]\n",sep="",file=filename,append=T)

    # get means of expected and predicted frequencies:

  }

  cat("\n\n### T1statistic for individual data:\n", file=filename,append=T)
  cat("T1ind.obs[n] <- 0", file=filename,append=T )
  for(i in 1:NOT){
    cat("+ sum( (response.",treeNames[i],"[n,] - n.expected.",treeNames[i],"[n,])^2/n.expected.",treeNames[i],"[n,])",
        sep="", file=filename,append=T)
  }

  # T1 for predicted means:
  cat("\nT1ind.pred[n] <- 0", file=filename,append=T )
  for(i in 1:NOT){
    cat("+ sum( (response.",treeNames[i],".pred[n,] - n.expected.",treeNames[i],"[n,])^2/n.expected.",treeNames[i],"[n,])", sep="", file=filename,append=T)
  }
  # 	T1.obs <- sum( (freq.mean-n.mean)^2/n.mean)
  # 	T1.pred <- sum( (pred.mean-n.mean)^2/n.mean)
  cat("\np.T1ind[n] <- T1ind.pred[n] > T1ind.obs[n]\n",  file=filename,append=T)


  ### individual analysis:
#   T1ind.obs[n] <- 0+
#     sum( (response.A[n,] - n.expected.A[n,])^2/n.expected.A[n,])+
#     sum( (response.B[n,] - n.expected.B[n,])^2/n.expected.B[n,])+
#     sum( (response.N[n,] - n.expected.N[n,])^2/n.expected.N[n,])
#
#   T1ind.pred[n] <- 0+
#     sum( (response.A.pred[n,] - n.expected.A[n,])^2/n.expected.A[n,])+
#     sum( (response.B.pred[n,] - n.expected.B[n,])^2/n.expected.B[n,])+
#     sum( (response.N.pred[n,] - n.expected.N[n,])^2/n.expected.N[n,])
#
#   p.T1ind[n] <- T1ind.pred[n] > T1ind.obs[n]

  cat("}\n",file=filename,append=T)



  cat("\n\n### T1 statistic for aggregated data:\n", file=filename,append=T)
  for(i in 1:NOT){
    cat("for(k in 1:",ncatPerTree[i],") {\n",sep="",file=filename,append=T)
    cat("n.expected.",treeNames[i],".mean[k] <- mean(n.expected.",treeNames[i],"[,k])\n",
        sep="",file=filename,append=T)
    cat("response.",treeNames[i],".pred.mean[k] <- mean(response.",treeNames[i],".pred[,k])\n",
        sep="",file=filename,append=T)
    cat("}\n", file=filename,append=T)
  }

  # T1 for observed means:
  cat("T1.obs <- 0", file=filename,append=T )
  for(i in 1:NOT){
    cat("+ sum( (response.",treeNames[i],".mean - n.expected.",treeNames[i],".mean)^2/n.expected.",treeNames[i],".mean)",
        sep="", file=filename,append=T)
  }

  # T1 for predicted means:
  cat("\nT1.pred <- 0", file=filename,append=T )
  for(i in 1:NOT){
    cat("+ sum( (response.",treeNames[i],".pred.mean - n.expected.",treeNames[i],".mean)^2/n.expected.",treeNames[i],".mean)", sep="", file=filename,append=T)
  }
# 	T1.obs <- sum( (freq.mean-n.mean)^2/n.mean)
# 	T1.pred <- sum( (pred.mean-n.mean)^2/n.mean)
	  cat("\np.T1 <- T1.pred > T1.obs\n",  file=filename,append=T)

}
	### Transformed parameters:
  cat(parString, file=filename, append = TRUE)

  cat("}\n",file=filename,append=T)
}





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








################### Beta-MPT specific hyperprior part
makeBetaHyperprior <-function(S,
                              alpha = "dunif(1,5000)",
                              beta = "dunif(1,5000)",
                              sampler = "JAGS"){
  modelString <- paste0("

for(s in 1:S){
  for(n in 1:subjs) {
    theta[s,n] ~dbeta(alph[s],bet[s])
  }
}

for(s in 1:S){",

ifelse(sampler%in%c("openbugs","OpenBUGS","winbugs","WinBUGS"),
       "\nzero[s] <- 0", "\n"),
"
  zero[s] ~dpois(phi[s])
  phi[s] <- -log(1/pow(alph[s]+bet[s],5/2))
  alph[s] ~",alpha,"
  bet[s] ~",beta,"
  mean[s] <- alph[s]/(alph[s]+bet[s])
  sd[s] <- sqrt(alph[s]*bet[s]/(pow(alph[s]+bet[s],2)*(alph[s]+bet[s]+1)))
}

")

  return(modelString)
}


################### Beta-MPT specific hyperprior part
makeTraitHyperprior <-function(S,
                               mu = "dnorm(0,1)",
                               xi = "dunif(0,100)",
                               sampler = "JAGS"){
  modelString <- paste0("

for(i in 1:subjs) {
  # probit transformation
  for(s in 1:S){
    theta[s,i] <- phi(mu[s] + xi[s]*delta.part.raw[s,i])
  }

  # hyperpriors
  delta.part.raw[1:S,i] ~ dmnorm(mu.delta.raw[1:S],T.prec.part[1:S,1:S])
}

T.prec.part[1:S,1:S] ~ dwish(V, df)
Sigma.Tau.raw[1:S,1:S] <- inverse(T.prec.part[,])


for(s in 1:S){
  mu.delta.raw[s] <- 0
  mu[s] ~ ", mu, "
  mean[s] <- phi(mu[s])
  xi[s] ~ ", xi, "
  for(q in 1:S){
    # Off-diagonal elements of S
    rho[s,q] <- Sigma.Tau.raw[s,q]/sqrt(Sigma.Tau.raw[s,s]*Sigma.Tau.raw[q,q])
  }
  # Diagonal elements of S
  sigma[s] <- xi[s]*sqrt(Sigma.Tau.raw[s,s])
}
")

  ###################################################### Hierarchical model ####################
#   for(p in 1:P){
#     theta.probit[i,p] <- mu[p] + xi[p]*delta.part.raw[i,p]
#   }
  # delta.part.raw[i,1:P] ~ dmnorm(mu.delta.raw[1:P],T.prec.part[1:P,1:P])


# T.prec.part[1:P,1:P] ~ dwish(V, df)
# Sigma.Tau.raw[1:P,1:P] <- inverse(T.prec.part[,])
#
# for(p in 1:P){
#   mu.delta.raw[p] <- 0
#   mu[p] ~ dnorm(0,1)
#   xi[p] ~ dunif (0, 100)
#
#
#   for(q in 1:P){
#     # Off-diagonal elements of S
#     rho[p,q] <- Sigma.Tau.raw[p,q]/sqrt(Sigma.Tau.raw[p,p]*Sigma.Tau.raw[q,q])
#   }
#   # Diagonal elements of S
#   sigma[p] <- xi[p]*sqrt(Sigma.Tau.raw[p,p])
# }

  return(modelString)
}
