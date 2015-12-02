#' Make the model description file needed by WinBUGS.
#' @param filename The name that the model description file should get
#' @param DataFormulas Output of mergeBranches()
#' @param nthetas The number of thetas
#' @param alpha Value for alpha
#' @param beta Value for betha
#' @param sampler Which sampler should be used? Default is "JAGS", further options are "OpenBUGS" and "WinBugs"
#' @author Nina R. Arnold, Denis Arnold
#' @export


makeModelDescription<-function(filename,DataFormulas,nthetas,alpha="dunif(1,5000)",beta="dunif(1,5000)",sampler="JAGS"){

  if(sampler%in%c("jags","JAGS")){

    cat("data\n",file=filename)
    cat("{\n",file=filename,append=T)
    cat("for(s in 1:",nthetas,"){\n",sep="",file=filename,append=T)
    cat("zero[s] <- 0\n",file=filename,append=T)
    cat("}\n }\n",file=filename,append=T)
    cat("model\n",file=filename,append=T)
  }else{
    cat("model\n",file=filename)
  }

	cat("{\n",file=filename,append=T)
	cat("for (n in 1: subjs){\n",file=filename,append=T)

	Trees=unique(DataFormulas$Trees)
	NOT=length(Trees)
	count=1

	for(i in 1:NOT){
		for(j in 1:length(grep(Trees[i],DataFormulas$Trees))){
		cat(Trees[i],"[n,",j,"]<-",DataFormulas$Formulas[count],"\n",sep="",file=filename,append=T)
		count=count+1
		}
	}

	for(i in 1:NOT){
		cat("response.",Trees[i],"[n,1:",length(grep(Trees[i],DataFormulas$Trees)),"]~dmulti(",
			Trees[i],"[n,1:",length(grep(Trees[i],DataFormulas$Trees)),"],items.",Trees[i],
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

	cat("mnb[s] <- alph[s]/(alph[s]+bet[s]) varp[s]<-sqrt(alph[s]*bet[s]/(pow(alph[s]+bet[s],2)*(alph[s]+bet[s]+1)))\n }\n }\n",file=filename,append=T)


}
