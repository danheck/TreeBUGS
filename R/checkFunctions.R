

### CHECK FUNCTIONS FOR INPUT


checkParEstFile <- function(parEstFile){
  if(!missing(parEstFile) && !is.null(parEstFile)){
    if(!is.character(parEstFile)){
      stop("'parEstFile' must be a character string pointing\n to the output file in an existing directory.")
    }
  }
  NULL
}


checkModelfilename <- function(modelfilename){

  if(missing(modelfilename) || is.null(modelfilename)){
    modelfilename <- tempfile(pattern = "MODELFILE",fileext = ".txt")
  }else if(!is.character(modelfilename)){
    stop("'parEstFile' must be a character string pointing\n to the mode file in an existing directory.")
  }

  modelfilename
}



readData  <- function(data){
  if(is.matrix(data) | is.data.frame(data)){
    data <- as.data.frame(data)
  }else{
    data = read.csv(data, header=TRUE, sep=",")
  }
  if(any(is.na(data))){
    stop("Missings in the data file!")
  }
  colnames(data) <- gsub(" ", "", colnames(data), fixed=TRUE)
  data
}




check.hyperprior <- function(par, thetaUnique, label="parameter"){
  if(length(par) == length(thetaUnique) && !is.null(names(par))){
    if(any(thetaUnique != sort(names(par))))
      stop("Names of the hyperprior vector '", label, "' do not match model parameters.",
           "\n  Use read.EQN(.., paramOrder=TRUE) to get the correct parameter labels.")
    par <- par[thetaUnique]
  }else if(length(par) == 1){
    par <- rep(par, length(thetaUnique))
    names(par) <- thetaUnique
  }

  par
}
