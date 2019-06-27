#' Get biomarkers from average data differences matrix (per type)
#'
#' Use this function to find either positive or negative biomarkers across many
#' performance classification group matchings based on a given threshold between
#' 0 and 1. The logic behind the biomarker selection is that if there is at
#' least one value in a column of
#' the \code{diff.mat} matrix that surpasses the threshold given, then the
#' corresponding node (name of the column) is return as a biomarker. This means
#' that for a single node, if at least one value that represents an average data
#' difference (for example, the average activity state difference) between any
#' of the given classification group comparisons (below) the threshold (negative
#' threshold), then a \emph{positive} (\emph{negative}) biomarker is reported.
#'
#' @param diff.mat a matrix whose rows are vectors of average node data
#' differences between two groups of models based on some kind of classification
#' (e.g. number of TP predictions) and whose names are set in the \code{rownames}
#' attribute of the matrix (usually denoting the diffferent classification
#' groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
#' that predicted 2 TP synergies, if the classification is done by number of TP
#' predictions). The columns represent the network's node names.
#' @param threshold numeric. A number in the [0,1] interval, above which (or
#' below its negative value) a biomarker will be registered in the returned result.
#' Values closer to 1 translate to a more strict threshold and thus less
#' biomarkers are found.
#' @param type character. Accepted values are \emph{positive} or \emph{negative}.
#'
#' @return a character vector that includes the node names that were found
#' either as \emph{positive} or \emph{negative}.
#'
#' @family biomarker functions
#'
#' @export
get_biomarkers_per_type = function(diff.mat, threshold, type) {
  stopifnot(threshold >= 0 & threshold <= 1)
  stopifnot(type == "positive" | type == "negative")

  dimen = dim(diff.mat)
  rows = dimen[1]
  nodes.num = dimen[2]

  biomarkers = character(0)
  for(node.index in 1:nodes.num) {
    node.name = colnames(diff.mat)[node.index]
    for (row.index in 1:rows) {
      if (type == "positive") {
        if (diff.mat[row.index, node.index] > threshold) {
          biomarkers = c(biomarkers, node.name)
          break
        }
      } else { # negative
        if (diff.mat[row.index, node.index] < -threshold) {
          biomarkers = c(biomarkers, node.name)
          break
        }
      }
    }
  }

  return(biomarkers)
}

#' Get total biomarkers from average data differences matrix
#'
#' Use this function to find all biomarkers across many
#' performance classification group matchings based on a given threshold between
#' 0 and 1. The logic behind the biomarker selection is that if there is at
#' least one value in a column of
#' the \code{diff.mat} matrix that surpasses the threshold given, then the
#' corresponding node (name of the column) is return as a biomarker. This means
#' that for a single node, if at least one value that represents an average data
#' difference (for example, the average activity state difference) between any
#' of the given classification group comparisons (below) the threshold (negative
#' threshold), then a \emph{positive} (\emph{negative}) biomarker is reported.
#'
#' @section Details:
#' This function uses the \code{\link{[emba](get_biomarkers_per_type)}} function
#' to get the biomarkers (nodes) of both types (positive and negative) from the
#' average data differences matrix. If a node though is found to surpass the
#' significance threshold
#' level given \emph{both negatively and positively}, we will keep it as a biomarker
#' in the category which corresponds to the \strong{comparison of the highest
#' classification groups}. For example, if the data comes from a model performance
#' classification based on the MCC score and in the comparison of the MCC classes
#' (1,3) the node of interest had an average difference of −0.89 (a negative
#' biomarker) while for the comparison of the (3,4) MCC classes it had a value
#' of 0.91 (a positive biomarker), then we will keep that node \emph{only as a
#' positive biomarker}. The logic behind this is that
#' the 'higher' performance-wise are the classification groups that we compare,
#' the more sure we are that the average data difference corresponds to a
#' \emph{better indicator} for the type of the biomarker found.
#'
#' @param diff.mat a matrix whose rows are vectors of average node data
#' differences between two groups of models based on some kind of classification
#' (e.g. number of TP predictions) and whose names are set in the \code{rownames}
#' attribute of the matrix (usually denoting the diffferent classification
#' groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
#' that predicted 2 TP synergies, if the classification is done by number of TP
#' predictions). The columns represent the network's node names.
#' @param threshold numeric. A number in the [0,1] interval, above which (or
#' below its negative value) a biomarker will be registered in the returned result.
#' Values closer to 1 translate to a more strict threshold and thus less
#' biomarkers are found.
#'
#' @return a list with two elements:
#'  \itemize{
#'    \item \code{biomarkers.pos}: a character vector that includes the node
#'    names that were found either as \emph{positive}
#'    \item \code{biomarkers.neg}: a character vector that includes the node
#'    names that were found either as \emph{negative}
#' }
#'
#' @family biomarker functions
#'
#' @importFrom usefun is_empty
#' @export
get_biomarkers = function(diff.mat, threshold) {
  stopifnot(threshold >= 0 & threshold <= 1)

  biomarkers.pos = get_biomarkers_per_type(diff.mat, threshold, type = "positive")
  biomarkers.neg = get_biomarkers_per_type(diff.mat, threshold, type = "negative")
  common.biomarkers = intersect(biomarkers.pos, biomarkers.neg)

  if (!is_empty(common.biomarkers)) {
    # in case of common biomarkers, remove them
    biomarkers.pos = biomarkers.pos[!biomarkers.pos %in% common.biomarkers]
    biomarkers.neg = biomarkers.neg[!biomarkers.neg %in% common.biomarkers]

    # find the proper category of the biomarkers and add them there
    for (biomarker in common.biomarkers) {
      logical.vector = diff.mat[, biomarker] > threshold |
        diff.mat[, biomarker] < (-threshold)
      comparison.index = max(which(logical.vector == TRUE))
      if (diff.mat[comparison.index, biomarker] > threshold)
        biomarkers.pos = append(biomarkers.pos, biomarker)
      else
        biomarkers.neg = append(biomarkers.neg, biomarker)
    }
  }

  # check: no common biomarkers
  stopifnot(is_empty(intersect(biomarkers.pos, biomarkers.neg)))

  res.list = list()
  res.list$biomarkers.pos = biomarkers.pos
  res.list$biomarkers.neg = biomarkers.neg

  return(res.list)
}
