
#' Permutation Test of Heterogeneity
#'
#' Tests whether whether participants (items) are homogeneous without assuming item (participant) homogeneity.
#'
#' @param data matrix or data frame with three columns: person code/index, item label, response category. Can also be the path to a .csv file with frequencies (comma-separated; first line defines category labels)
#' @param tree a list that defines which categories belong to the same multinomial distribution (i.e., the the same MPT tree). For instance: \code{tree = list(tree.old = c("hit","cr"), tree.new = c("fa","miss"))}. Category labels must match the values of the third column of \code{data}
#' @param source whether to test for \code{"person"} or \code{"item"} homogeneity
# @param stat which statistic to use (either \code{"var"} or \code{"chisq"})
#' @param rep number of permutations to be sampled
#' @param nCPU number of CPUs used for parallel Monte Carlo sampling of permutations
#' @details
#' If an item/person has zero frequencies on all categories in an MPT tree, these zeros are neglected when computing mean frequencies per column. As an example, consider a simple recognition test with a fixed assignments of words to the learn/test list. In such an experiment, all learned words will result in hits or misses (i.e., the MPT tree of old items), whereas new words are always false alarms/correct rejections and thus belong to the MPT tree of new items (this is not necessarily the case if words are assigned randomly).
#'
#' Note that the test does still assume independence of observations. However, it does not require item homogeneity when testing participant heterogeneity (in contrast to the chi-square test: \code{\link{testHetChi}}).
#' @seealso \code{\link{testHetChi}}, \code{\link{plotFreq}}
#' @author Daniel W. Heck
#' @references Smith, J. B., & Batchelder, W. H. (2008). Assessing individual differences in categorical data. Psychonomic Bulletin & Review, 15, 713-731. \doi{10.3758/PBR.15.4.713}
#' @examples
#' # generate homogeneous data
#' # (N=15 participants, M=30 items)
#' data <- data.frame(id = rep(1:15, each=30),
#'                    item = rep(1:30, 15))
#' data$cat <- sample(c("h","cr","m","fa"),15*30,
#'                    replace = TRUE,
#'                    prob = c(.7,.3,.4,.6))
#' head(data)
#' tree <- list(old = c("h","m"),
#'              new = c("fa", "cr"))
#'
#' # test participant homogeneity:
#' tmp <- testHetPerm(data, tree, rep=200, nCPU=1)
#' tmp[2:3]
#' @export
testHetPerm <- function(data,
                        tree,
                        source = "person",
                        # stat = "var",
                        rep = 1000,
                        nCPU = 4){

  # some made up data for a simple memory recognition experiment:
  # data <- data.frame(id =   c(1,1,1,1,1,1,1,3,3,3,3,3,3,3,3),
  #                    item = c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5),
  #                    cat = c("h","cr","h","m","fa",
  #                            "h","cr","fa","fa","m",
  #                            "m","cr","h","fa","h"))
  # tree <- list(old = c("h","m"),
  #              new = c("fa", "cr"))

  if(is.character(data)){
    data <- read.csv(file = data)
  }
  data <- as.data.frame(data)

  if(ncol(data) != 3){
    stop("'data' must have three columns (person code, item label, response category).")
  }

  ################### NOT ELEGANT
  if(source == "item")
    data <- data[,c(2,1,3)]

  persons <- as.vector(unlist(data[,1]))
  items <- as.vector(unlist(data[,2]))
  cats <- as.vector(unlist(data[,3]))
  # # number of participants/items:
  N <- length(unique(persons))

  if(missing(tree)){
    warning("It is assumed that all columns of 'freq' stem from a single multinomial distribution",
            "\n  (i.e., from a single MPT tree)")
    tree <- list(single = unique(cats))
  }
  if(is.null(names(tree)))
    names(tree) <- paste0("t.", 1:length(tree))

  # number of categories per tree:
  K <- sapply(tree, length)
  if(sum(K) != length(unique(cats)))
    stop("Number of category labels in 'data' does not match number of labels in 'tree'.")


  ################# observed test statistic

  tree.labels <- rep(names(tree),sapply(tree, length))
  names(tree.labels) <- unlist(tree)
  chi.obs <- testHetChi(table(persons, cats),
                        tree = tree.labels[sort(unique(cats))])$chisq
  tree.vec <- tree.labels[data[,3]]


  ################# sampled test statistic

  single.rep <- function(i, data, tree.vec, tree.labels){
    # for each person: permute categorical responses:
    cat.perm <- ave(data[,3],             # categorical responses
                    data[,2], tree.vec,   # group by (A) items (B) MPT tree
                    FUN = sample)         # permutation
    testHetChi(table(data[,1], cat.perm),
               tree = tree.labels[sort(unique(data[,3]))])$chisq
  }

  # test:
  # cbind(data, replicate(3, ave(data[,3], data[,2], tree.vec, FUN = sample)))

  if(nCPU == 1){
    chi.samp <- sapply(1:rep, single.rep,
                       data=data, tree.vec=tree.vec,
                       tree.labels=tree.labels)
  }else{
    cl <- makeCluster(nCPU)
    chi.samp <- parSapply(cl, 1:rep, single.rep,
                          data=data, tree.vec=tree.vec,
                          tree.labels=tree.labels)
    stopCluster(cl)
  }

  list(chi.samp = chi.samp,
       chi.obs = chi.obs,
       prob = mean(chi.samp > chi.obs))
}



# # The R code below tests for subject variability
# reps<-10000 #number of permutations
# subjs<-29 #number of participants in data set
# items<-50 #number of items in data set
# cate<-4 #number of response categories
# permutedata <-matrix(0,subjs,items)
# contingency<-matrix(0,subjs,cate)
# x <- matrix(round(runif(items*subjs, -.49, 4.49)),subjs,items)
#
#
# #initiate variable to hold contingency tables from permuted data sets
# for (i in 1:reps){
#   for (ii in 1:items) { #randomly permutes each row of data matrix D
#     permutedata[,ii]<-x[sample(subjs,subjs),ii]
#   }
#   for (cc in 1:cate) { #converts data matrix to Nx2 contingency table
#     contingency[cc,]<-apply(permutedata==cc,1,sum)
#     #calculated permuted contingency table
#   }
#   #any test of interest can then be run on the permuted contingency table
#   var(colSums(contingency))
# }
#
#
# contingency=zeros(subjs, cate);
# %initiate variable to hold contingency tables from permuted data sets
# for i=1:reps
# for ii=1:subjs %randomly permutes each row of data matrix D
# permuteddata(ii,:)=datamatrix(ii,randperm(items));
# end
# for cc=1:cate %converts data matrix to Nx2 contingency table
# contingency(:,cc)=sum(permuteddata==cc)';
#  end
#  %the statistical test of interest
# end
