#' Userinterface to use BUGS
#' @param eqnfile The (full path to the) file that specifies the Multitree.
#' @param subjdata The (full path to the) cvs file with the data.
#' @param restrictions  The (full path to the) file that specifies which parameters should be constants and which should be equal.
#' @param modelfilename Name that the modelfile that is made by the function to work with WinBUGS should get.
#'        Default is to write this information to the tempdir as required by CRAN standards.
#' @param alpha Model parameter.
#' @param beta Model parameter.
#' @param parameters List of parameters that WinBUGS should return.
#' @param parestfile Name of the file to which the results should be stored.
#' @param n.iter Number of iterations.
#' @param n.burnin Burnin period.
#' @param n.thin Thinning rate.
#' @param sampler Which sampler should be used? Default is "JAGS", further options are "OpenBUGS" and "WinBugs"
#' @param ... Arguments to be passed to other methods.
#' @author Nina R. Arnold, Denis Arnold
#' @export

mpt2BetaMPT<-function(eqnfile,
                   subjdata,
                   restrictions,
                   modelfilename=NULL,
                   alpha="dunif(1,5000)",
                   beta="dunif(1,5000)",
                   parameters=list("theta", "alph", "bet", "mnb", "varp"),
                   parestfile,
                   n.iter=100000,
                   n.burnin=50000,
                   n.thin=2,
                   sampler="JAGS",
                   ...){

  if(is.null(modelfilename)){
    modelfilename=tempfile(pattern="MODELFILE",fileext=".txt")
  }
  Tree=readMultiTree(eqnfile)
  SubjData=readSubjectData(subjdata,unique(Tree$Answers))
  Tree=mergeBranches(Tree,names(SubjData))
  tHoutput=thetaHandling(Tree,restrictions)
  SubPar=tHoutput[[1]]
  Tree=tHoutput[[2]]
  makeModelDescription(modelfilename,Tree,max(SubPar$theta),alpha,beta,sampler)
  callingBetaMPT(Tree,
                 subjdata,
                 modelfile=modelfilename,
                 numberOfParameters=max(SubPar$theta),
                 parameters,
                 parestfile,
                 n.iter=n.iter,
                 n.burnin=n.burnin,
                 n.thin=n.thin,
                 sampler=sampler,
                 ...)

}
