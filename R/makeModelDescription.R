#' Make the model description file needed by WinBUGS.
#' @param filename The name that the model description file should get
#' @param DataFormulas Output of mergeBranches()
#' @param nthetas The number of thetas
#' @param alpha Value for alpha
#' @param beta Value for betha
#' @param sampler Which sampler should be used? Default is "JAGS", further options are "OpenBUGS" and "WinBugs"
#' @author Nina R. Arnold, Denis Arnold, Daniel Heck
#' @export


makeModelDescription<-function(filename,
                               DataFormulas,
                               nthetas,
                               alpha="dunif(1,5000)",
                               beta="dunif(1,5000)",
                               sampler="JAGS"){

  Trees=unique(DataFormulas$Trees)
  NOT=length(Trees)
  count=1
  # number of categories per tree:
  ncatPerTree <- rep(NA, NOT)
  for(i in 1:NOT) {
    ncatPerTree[i] <- length(grep(Trees[i],DataFormulas$Trees))
  }

  if(sampler%in%c("jags","JAGS")){

    ################################### DATA #######################
    cat("data\n",file=filename)
    cat("{\n",file=filename,append=T)
    cat("for(s in 1:",nthetas,"){\n",sep="",file=filename,append=T)
    cat("zero[s] <- 0\n",file=filename,append=T)
    cat("}\n ",file=filename,append=T)


    ####### for posterior predictive cehck: mean frequencies
    for(i in 1:NOT){
      cat("for(k in 1:",ncatPerTree[i],") {\n",sep="",file=filename,append=T)
        cat("response.",Trees[i],".mean[k] <- mean(response.",Trees[i],"[,k])\n",sep="",file=filename,append=T)
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
		cat(Trees[i],"[n,",j,"]<-",DataFormulas$Formulas[count],"\n",sep="",file=filename,append=T)
		count=count+1
		}
	}

	for(i in 1:NOT){
		cat("response.",Trees[i],"[n,1:",ncatPerTree[i],"]~dmulti(",
			Trees[i],"[n,1:",ncatPerTree[i],"],items.",Trees[i],
			"[n])\n",sep="",file=filename,append=T)
	}

	cat("}\n",file=filename,append=T)

	cat("for(s in 1:",nthetas,"){\n",sep="",file=filename,append=T)
	cat("for(n in 1:subjs) {\n",file=filename,append=T)
	cat("theta[s,n] ~dbeta(alph[s],bet[s])\n }\n }\n",file=filename,append=T)

	cat("for(s in 1:",nthetas,"){\n",sep="",file=filename,append=T)
	if(sampler%in%c("openbugs","OpenBUGS","winbugs","WinBUGS")){
	cat("zero[s] <- 0\n",file=filename,append=T)
	}
	cat("zero[s] ~dpois(phi[s])\n",file=filename,append=T)
	cat("phi[s] <- -log(1/pow(alph[s]+bet[s],5/2))\n",file=filename,append=T)

	cat("alph[s] ~",alpha,"\n",sep="",file=filename,append=T)
	cat("bet[s] ~",beta,"\n",sep="",file=filename,append=T)

	cat("mnb[s] <- alph[s]/(alph[s]+bet[s]) varp[s]<-sqrt(alph[s]*bet[s]/(pow(alph[s]+bet[s],2)*(alph[s]+bet[s]+1)))\n }\n \n",file=filename,append=T)



	if(sampler%in%c("jags","JAGS")){
	############################## Posterior predictive checks:

	cat("\n############# T1: posterior predictive check of mean frequencies\n", file=filename,append=T)

  cat("for(n in 1:subjs) {\n",file=filename,append=T)
  for(i in 1:NOT){
    # sample predicted: frequencies from posterior
    cat("response.",Trees[i],".pred[n,1:",ncatPerTree[i],"]~dmulti(",
        Trees[i],"[n,1:",ncatPerTree[i],"],items.",Trees[i],
        "[n])\n",sep="",file=filename,append=T)

    # compute expected frequencies:
    cat("n.expected.",Trees[i],"[n,1:",ncatPerTree[i],"] <- ",
        Trees[i],"[n,1:",ncatPerTree[i],"] * items.",Trees[i],
        "[n]\n",sep="",file=filename,append=T)

    # get means of expected and predicted frequencies:

  }
  cat("}\n",file=filename,append=T)

  for(i in 1:NOT){
    cat("for(k in 1:",ncatPerTree[i],") {\n",sep="",file=filename,append=T)
    cat("n.expected.",Trees[i],".mean[k] <- mean(n.expected.",Trees[i],"[,k])\n",
        sep="",file=filename,append=T)
    cat("response.",Trees[i],".pred.mean[k] <- mean(response.",Trees[i],".pred[,k])\n",
        sep="",file=filename,append=T)
    cat("}\n", file=filename,append=T)
  }

  # T1 for observed means:
  cat("T1.obs <- 0", file=filename,append=T )
  for(i in 1:NOT){
    cat("+ sum( (response.",Trees[i],".mean - n.expected.",Trees[i],".mean)^2/n.expected.",Trees[i],".mean)",
        sep="", file=filename,append=T)
  }

  # T1 for predicted means:
  cat("\nT1.pred <- 0", file=filename,append=T )
  for(i in 1:NOT){
    cat("+ sum( (response.",Trees[i],".pred.mean - n.expected.",Trees[i],".mean)^2/n.expected.",Trees[i],".mean)", sep="", file=filename,append=T)
  }
# 	T1.obs <- sum( (freq.mean-n.mean)^2/n.mean)
# 	T1.pred <- sum( (pred.mean-n.mean)^2/n.mean)
	  cat("\np.T1 <- T1.pred > T1.obs\n",  file=filename,append=T)

}
	### Transformed parameters:


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
