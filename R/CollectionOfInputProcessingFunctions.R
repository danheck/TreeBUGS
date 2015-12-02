#	A collection of functions to read multiTree files and generate WinBUGS control files.
#	(c) 2015 Nina Arnold and Denis Arnold

#' Read a multiTree file
#' @param file The (full path to the) file that specifies the Multitree
#' @author Nina R. Arnold, Denis Arnold
#' @export



readMultiTree<-function(file){

	multiTreeDefinition = read.csv(file,header=F,stringsAsFactors=F)	#read file

	# Variables
	numberOfBranches=-1
	checkNOB=F

	# Checking whether everything is like we expect it.
	options(warn=-1)

	try(if((grep("^[0-9]+$",multiTreeDefinition$V1[1])==1)){numberOfBranches=as.numeric(multiTreeDefinition$V1[1])},silent=T)

	try(if((dim(multiTreeDefinition)[1]==(as.numeric(multiTreeDefinition$V1[1])+1))){checkNOB=T},silent=T)

	if(numberOfBranches==-1){
		cat("Wrong format of the definition file: Only number of branches expected in first row!\n")
		return(-1)
	}
	if(checkNOB==F){
		cat("Wrong format of the definition file: Number of branches not as specified in the first row!\n")
		return(-1)
	}

	# Get Trees, Answers and Formulas out
	Elements=multiTreeDefinition[-1,]
	Elements=gsub("[\ ]+","\ ",Elements)
	Elements=unlist(strsplit(Elements,split="\ "))


	# Check whether nothing went wrong
	if((length(Elements)/numberOfBranches)!=3){
		cat("Getting the elements failed! Is everything seperated by \"\ \"? Are any \"\ \" in the names or formulas?")
		return(-1)
	}

	# Make a nice data.frame to return
	Trees=Elements[seq(1,length(Elements),3)]
	Answers=Elements[seq(2,length(Elements),3)]
	Formulas=Elements[seq(3,length(Elements),3)]

	TreeData=data.frame(Trees=Trees,Answers=Answers,Formulas=Formulas,stringsAsFactors=F)

	return(TreeData)
}

#' Unique branches for each Tree with the same answer and sum the corresponing formulas
#' @param TreeData Data returned from readMultiTree()
#' @param DataNames Vector containing the names in order of the data file
#' @author Nina R. Arnold, Denis Arnold
#' @export

mergeBranches<-function(TreeData, DataNames){ # Unique branches for each Tree with one answer and sum the corresponing formulas

	Trees=unique(TreeData$Trees)
	TreesUniqueAnswers=list()

	for(i in 1:length(Trees)){
		TreesUniqueAnswers[[i]]=unique(TreeData[TreeData$Trees==Trees[i],]$Answers)
	}

	NewTrees=data.frame(Trees="character",Answers="character",Formulas="character",stringsAsFactors=F)
	NewTreesFormula=TreesUniqueAnswers

	for(i in 1:length(Trees)){
		for(j in 1:length(TreesUniqueAnswers[[i]])){
			NewTreesFormula[[i]][j]=paste(TreeData$Formulas[TreeData$Trees==Trees[i]&TreeData$Answers==TreesUniqueAnswers[[i]][j]],collapse="+")
		}
	}

	for(i in 1:length(Trees)){
		for(j in 1:length(TreesUniqueAnswers[[i]])){
			NewTrees=rbind(NewTrees,data.frame(Trees=Trees[i],Answers=TreesUniqueAnswers[[i]][j],Formulas=NewTreesFormula[[i]][j],stringsAsFactors=F))
		}
	}

	MergedTree=NewTrees[-1,]
	MergedTree=MergedTree[order(MergedTree$Answers)[order(DataNames)],]

	return(MergedTree)
}

#' Extract all parameter from the formulas of a given model
#' @param TreeData Data returned from readMultiTree()
#' @author Nina R. Arnold, Denis Arnold
#' @export

getParameter<-function(TreeData){

	Parameter=unique(unlist(strsplit(TreeData$Formulas,split="\\*|\\(|\\)|\\-|\\+")))
	r=c(which(nchar(Parameter)==0),grep("^[0-9]+$|^[0-9]+\\.[0-9]+",Parameter))
	Parameter=Parameter[-r]
	Parameter=c(sort(Parameter[grep("[A-Z]",Parameter)]),sort(Parameter[-grep("[A-Z]",Parameter)]))
	return(Parameter)

}


#' Read the subject data from file
#' @param file is the data stored in a .csv file
#' @param answers is the unique of the $Answers from TreeData which is $Answers after mergeBrachnes()
#' @author Nina R. Arnold, Denis Arnold
#' @export

readSubjectData<-function(file,answers){

	Data=read.csv(file)

	if(sum(names(Data)%in%answers)!=length(answers)){
		if(dim(Data)[2]!=length(answers)){
			cat("Number of answers differs in eqn and csv file")
			return(-1)
		}
		else{
			cat("At least one name of the answers differs in eqn and csv file")
			return(-1)
		}
	}

	return(Data)
}


#'  Set constants, replace parameters for thetas
#' @param MergedTree is the output of mergeBranches
#' @param file The (full path to the) file that specifies which parameters should be constants and which should be equal
#' @author Nina R. Arnold, Denis Arnold
#' @export

thetaHandling<-function(MergedTree,file){

  Parameter=getParameter(MergedTree)
  partrack=Parameter

  # Handle the constants

  SubPar=read.csv(file,header=F,stringsAsFactors=F)
  numberSubs=grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1)
  replaceConst=data.frame(const="character",sub="character",stringsAsFactors=F)
  if(length(numberSubs)>0){
    for(i in 1:length(numberSubs)){
      X=unlist(strsplit(SubPar$V1[numberSubs[i]],"="))
      number=grep("^[0-9]+$|^[0-9]+\\.[0-9]+$",X)
      if(length(number)>1){
        cat("Variable set to multiple constants")
        return(-1)
      }
      else{
        for(i in 1:length(X)){
          if(i!=number){
            replaceConst=rbind(replaceConst,data.frame(const=X[i],sub=X[number]))
          }
        }

      }
    }
    replaceConst = replaceConst[-1,]
    for(i in 1:dim(replaceConst)[1]){
      MergedTree$Formulas=gsub(paste("(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",replaceConst$const[i],"(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)",sep=""),paste("\\1",replaceConst$sub[i],"\\2",sep=""),MergedTree$Formulas,perl=T)
    }
  }

  partrack=partrack[-which(partrack%in%replaceConst$const)]

  # Handle parameters set equal
  if(length(grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1))>0){
    SubPar=SubPar[-grep("=[0-9]+$|=[0-9]+\\.[0-9]+$",SubPar$V1),]
  }
  SubPar=strsplit(SubPar,split="=")

  theta=numeric()

  for(i in 1:length(SubPar)){
    for(j in 1:length(SubPar[[i]])){
      theta=c(theta,i)
    }

  }
  SubPar=data.frame(Par=unlist(SubPar),theta=theta,stringsAsFactors=F)
  SubPar$sub=paste("theta[",SubPar$theta,",n]",sep="")

  partrack=partrack[-which(partrack%in%SubPar$Par)]
  if(length(partrack)>0){
    for(i in 1:length(partrack)){
      SubPar=rbind(SubPar,data.frame(Par=partrack[i],theta=max(SubPar$theta+1),sub=paste("theta[",max(SubPar$theta)+1,",n]",sep="")))
    }
  }

  for(j in 1:length(SubPar$Par)){
    MergedTree$Formulas=gsub(paste("(^|\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\))",SubPar$Par[j],"(\\*|\\(|\\)|\\-|\\+|\\+\\(|\\*\\(|\\-\\(|\\+\\)|\\*\\)|\\-\\)|$)",sep=""),paste("\\1",SubPar$sub[j],"\\2",sep=""),MergedTree$Formulas,perl=T)
  }

  output=list(SubPar,MergedTree)

  return(output)

}
