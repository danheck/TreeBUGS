#' Generate EQN Files for Within-Subject Designs
#'
#' Replicates an MPT model multiple times with different tree, category, and
#' parameter labels for within-subject factorial designs.
#'
#' @inheritParams betaMPT
#' @param labels a character vector defining the labels that are added to the
#'   parameters in each within-subject condition
#' @param constant optional: a character vector defining which parameters are
#'   constrained to be constant across within-conditions
#' @param save optional: path to an EQN output file. By default, the model is
#'   return as a string character
#'
#' @examples
#' # Example: Standard Two-High-Threshold Model (2HTM)
#' EQNfile <- system.file("MPTmodels/2htm.eqn",
#'   package = "TreeBUGS"
#' )
#' withinSubjectEQN(EQNfile, c("high", "low"), constant = c("g"))
#' @export
withinSubjectEQN <- function(
    eqnfile,
    labels,
    constant,
    save
) {
  tree <- readEQN(eqnfile)
  param <- colnames(readEQN(eqnfile, parse = TRUE)$a)
  if (!missing(constant)) {
    param <- setdiff(param, constant)
  }

  if (length(unique(labels)) != length(labels)) {
    stop("The within-subject 'labels' must be unique!")
  }
  tree.list <- list(tree)[rep(1, length(labels))]

  for (w in 1:length(labels)) {
    tree.list[[w]]$Tree <- paste0(labels[w], "_", substr(tree$Tree, 3, 999))
    tree.list[[w]]$Category <- paste0(labels[w], "_", tree$Category)
    for (b in 1:nrow(tree)) {
      for (p in 1:length(param)) {
        # look behind mechanism: check for dots in parameter labels via  (?!\\.)
        # https://stackoverflow.com/questions/23094532/java-regular-expression-word-without-ending-with-dot
        # requires perl=TRUE

        tree.list[[w]]$Equation[b] <- gsub(paste0("\\b", param[p], "\\b(?!\\.)"),
          paste0(param[p], "_", labels[w]),
          tree.list[[w]]$Equation[b],
          perl = TRUE
        )

        # tree.list[[w]]$Equation[b] <- ifelse(tree.list[[w]]$Equation[b] == param[p],
        #                                      paste0(param[p],"_",labels[w]),
        #                                      tree.list[[w]]$Equation[b])
        # tree.list[[w]]$Equation[b] <- gsub(paste0("1-",param[p]),
        #                                    paste0("1-",param[p],"_",labels[w]),
        #                                    tree.list[[w]]$Equation[b],fixed=TRUE)
        # tree.list[[w]]$Equation[b] <- gsub(paste0("+",param[p]),
        #                                    paste0("+",param[p],"_",labels[w]),
        #                                    tree.list[[w]]$Equation[b],fixed=TRUE)
        # tree.list[[w]]$Equation[b] <- gsub(paste0("*",param[p]),
        #                                    paste0("*",param[p],"_",labels[w]),
        #                                    tree.list[[w]]$Equation[b],fixed=TRUE)
      }
    }
  }

  res <- do.call("rbind", tree.list)
  if (!missing(save)) {
    write.table(res, file = save, quote = FALSE, row.names = FALSE, sep = "     ")
  }
  res
}
