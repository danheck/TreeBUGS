#'  Calling the Sampler (JAGS or OpenBUGS or WinBUGS)
#'  @param DataFormulas Ouptut of mergeBranches()
#'  @param data data frame with frequencies
#'  @param modelfile Full path to the model description file. (Made with makeModelDescription())
#'  @param numberOfParameters How many parameters does the model have
# not necessary:  @param parameters List with parameters that should be returned from the sampler.
#'  @param n.iter Number of iterations.
#'  @param n.burnin Burnin period.
#'  @param n.update Update paramater for JAGS
#'  @param n.thin Thinning rate.
#'  @param n.chains number of MCMC chains
#'  @param sampler The Sampler to be used. Options are JAGS, OpenBUGS or WinBUGS.
#'         Note that you need to install the sampler on your computer.
#' @param autojags whether to run JAGS until the MCMC chains converge (see \link{autojags}). Can take a lot of time for large models.
# not necessary:  @param savetable name of the file to which the results should be saved on the hard drive.
#'  @param ... Arguments to be passed to other methods.
#'  @author Nina R. Arnold, Denis Arnold, Daniel W. Heck
#'  @export
callingBetaMPT<-function(DataFormulas,
                         data,
                         modelfile,
                         numberOfParameters=0,
                         # not necessary:parameters = list("theta", "alph", "bet", "mnb", "varp"),
                         n.iter=100000,
                         n.burnin=NULL,
                         n.update= 10,
                         n.thin=2,
                         n.chains=3,
                         sampler="JAGS",
                         autojags=TRUE,
                         # savetable = NULL,
                         ...){

  if(is.na(n.burnin)){n.burnin=n.iter/2}
  parameters <- list("theta", "alph", "bet", "mnb", "varp")

  responses <- data
  subjs     = nrow(responses)

  ### prepare data for WinBUGS

  Trees=unique(DataFormulas$Trees)
  NresponsesTree=vector("numeric",length=length(Trees))

  for(i in 1:length(Trees)){
    NresponsesTree[i]=length(which(DataFormulas$Trees==Trees[i]))
  }

  index=0
  for(i in 1:length(Trees)){
    assign(paste("items",Trees[i],sep="."),rowSums(responses[(index+1):(index+NresponsesTree[i])]))
    assign(paste("response",Trees[i],sep="."),matrix(as.vector(t(responses[(index+1):(index+NresponsesTree[i])])),ncol=NresponsesTree[i],nrow=subjs,byrow=TRUE))
    index=index+NresponsesTree[i]
  }

  #make data

  data=list()
  index=1
  for(i in 1:length(Trees)){
    data[[index]]=paste("items",Trees[i],sep=".")
    index=index+1
  }
  for(i in 1:length(Trees)){
    data[[index]]=paste("response",Trees[i],sep=".")
    index=index+1
  }
  data[[index]]="subjs"


  # call Sampler

  if(sampler%in%c("jags","JAGS")){
    W=diag(numberOfParameters)
    data[[index+1]]="W"
    parametervector=c(unlist(parameters), "T1.obs","T1.pred","p.T1")
    samples = jags.parallel(data,
                            inits=NULL,
                            parameters.to.save=parametervector,
                            model.file = modelfile,
                            n.iter=n.iter,
                            n.burnin=n.burnin,
                            n.chains=n.chains,
                            DIC=T,...)
    if(autojags){
      recompile(samples)
      samples.upd <- autojags(samples, n.update = n.update)
      samples=samples.upd
    }

  }else{
    if(sampler%in%c("openbugs","OpenBUGS","winbugs","WinBUGS")){

      samples = bugs(data,
                     inits=NULL,
                     parameters,
                     model.file = modelfile,
                     n.chains=n.chains,
                     n.iter=n.iter,
                     n.burnin=n.burnin,
                     n.thin=n.thin,
                     DIC=T,
                     codaPkg=F,
                     debug=F,
                     program = sampler,
                     ...)
      # get coda samples (better to process afterwards)
      # tmp <- capture.output({samples <- read.bugs(out)})
    }else{
      print(paste("Unknown sampler:",sampler))
    }
  }

  return(samples)
}
