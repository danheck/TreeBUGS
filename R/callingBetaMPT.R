#'  Calling the Sampler (JAGS or OpenBUGS or WinBUGS)
#'  @param DataFormulas Ouptut of mergeBranches()
#'  @param datafile Full path to cvs file with the data.
#'  @param modelfile Full path to the model description file. (Made with makeModelDescription())
#'  @param numberOfParameters How many parameters does the model have
#'  @param parameters List with parameters that should be returned from the sampler.
#'  @param parestfile Name of the file to with the estimates should be stored.
#'  @param n.iter Number of iterations.
#'  @param n.burnin Burnin period.
#'  @param n.update Update paramater for JAGS
#'  @param n.thin Thinning rate.
#'  @param sampler The Sampler to be used. Options are JAGS, OpenBUGS or WinBUGS.
#'         Note that you need to install the sampler on your computer.
#'  @param printresults Whether a summary of the results should be printed.
#'  @param savetable Whether the results should be saved to the hard drive.
#'  @param ... Arguments to be passed to other methods.
#'  @author Nina R. Arnold, Denis Arnold
#'  @export

callingBetaMPT<-function(DataFormulas,
                         datafile,
                         modelfile,
                         numberOfParameters=0,
                         parameters = list("theta", "alph", "bet", "mnb", "varp"),
                         parestfile = "ResultsBetaMPT.txt",
                         n.iter=100000,
                         n.burnin=NULL,
                         n.update= 10,
                         n.thin=2,
                         sampler="JAGS",
                         printresults = FALSE,
                         savetable = FALSE,
                         ...){

  if(is.na(n.burnin)){n.burnin=n.iter/2}

  ### read data

  responses = read.csv(datafile, header=TRUE)
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
  for(i in 1:length(Trees)){data[[index]]=paste("items",Trees[i],sep=".");index=index+1}
  for(i in 1:length(Trees)){data[[index]]=paste("response",Trees[i],sep=".");index=index+1}
  data[[index]]="subjs"


  # call Sampler

  if(sampler%in%c("jags","JAGS")){
    W=diag(numberOfParameters)
    data[[index+1]]="W"
    parametervector=unlist(parameters)
    samples = jags.parallel(data, inits=NULL, parameters.to.save=parametervector, model.file = modelfile, n.chains=3, DIC=T)
    recompile(samples)
    samples.upd <- autojags(samples, n.update = n.update)
    samples=samples.upd

  }else{
  if(sampler%in%c("openbugs","OpenBUGS","winbugs","WinBUGS")){

  samples = bugs(data,
                 inits=NULL,
                 parameters,
                 model.file = modelfile,
                 n.chains=3,
                 n.iter=n.iter,
                 n.burnin=n.burnin,
                 n.thin=n.thin,
                 DIC=T,
                 codaPkg=F,
                 debug=F,
                 program = sampler,
                 ...)
  }else{
    print(paste("Unknown sampler:",sampler))
  }
  }
  # write results

  if(savetable){
  write.table(samples$summary[,], file=parestfile, sep ="\t",na="NA",dec=".",row.names=T,col.names=T,quote=F)
  }

  # print results

  if(printresults){
  print(samples)
  }

  return(samples)
}
