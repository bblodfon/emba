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
#' corresponding node (name of the column) is returned as a biomarker. This means
#' that for a single node, if at least one value that represents an average data
#' difference (for example, the average activity state difference) between any
#' of the given classification group comparisons is above (below) the threshold
#' (negative threshold), then a \emph{positive} (\emph{negative}) biomarker is
#' reported.
#'
#' @section Details:
#' This function uses the \code{\link{get_biomarkers_per_type}} function
#' to get the biomarkers (nodes) of both types (positive and negative) from the
#' average data differences matrix. If a node though is found to surpass the
#' significance threshold
#' level given \emph{both negatively and positively}, we will keep it as a biomarker
#' in the category which corresponds to the \strong{comparison of the highest
#' classification groups}. For example, if the data comes from a model performance
#' classification based on the MCC score and in the comparison of the MCC classes
#' (1,3) the node of interest had an average difference of -0.89 (a negative
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
#'    names of the \emph{positive} biomarkers
#'    \item \code{biomarkers.neg}: a character vector that includes the node
#'    names of the \emph{negative} biomarkers
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

#' Get synergy biomarkers from dir
#'
#' This function reads the synergy biomarker files inside the given directory and merges
#' the results into a \code{data.frame} which it returns. This functions should
#' be used when the synergy biomarker results are in seperate files inside the
#' directory given (see \code{biomarkers.dir} parameter).
#'
#' @param predicted.synergies a character vector of the synergies (drug
#' combination names) that were predicted by \strong{at least one} of the models
#' in the dataset.
#' @param biomarkers.dir string. It specifies the full path name of the
#' directory which holds the biomarker files. The biomarker files must be
#' formatted as: \emph{\%drug.comb\%_biomarkers_active} or
#' \emph{\%drug.comb\%_biomarkers_inhibited}, where \%drug.comb\% is an element
#' of the \code{predicted.synergies} vector.
#' @param models.dir string. A directory with \emph{.gitsbe} files/models. It's
#' needed in order to call \code{\link{get_node_names}}.
#' @param node.names a character vector which has the names of the nodes. If it's
#' not NULL, then it will be used instead of the \code{models.dir} parameter.
#' The \code{node.names} should include all the nodes that are reported as
#' biomarkers in the biomarker files inside the \code{biomarkers.dir} directory.
#' Default value: NULL.
#'
#' @return a data.frame, whose columns represent the network nodes and the
#' rows the predicted synergies. Possible values for each \emph{synergy-node}
#' element are either \emph{1} (\emph{active state} biomarker), \emph{-1}
#' (\emph{inhibited state} biomarker) or \emph{0} (not a biomarker).
#'
#' @importFrom utils read.table
#' @export
get_synergy_biomarkers_from_dir =
  function(predicted.synergies, biomarkers.dir, models.dir, node.names = NULL) {
    stopifnot(!is.null(models.dir) || !is.null(node.names))

    # get the node names
    if (is.null(node.names))
      node.names = get_node_names(models.dir)

    # initialize res data.frame
    res = as.data.frame(matrix(0, ncol = length(node.names),
                               nrow = length(predicted.synergies)))
    colnames(res) = node.names
    rownames(res) = predicted.synergies

    for (drug.comb in predicted.synergies) {
      # insert the active biomarkers
      active.biomarkers.file =
        paste0(biomarkers.dir, drug.comb, "_biomarkers_active")

      if (file.exists(active.biomarkers.file)
          && file.size(active.biomarkers.file) != 0) {
        biomarkers.active =
          read.table(active.biomarkers.file, stringsAsFactors = FALSE)
        biomarkers.active.names = biomarkers.active[,1]
        # biomarkers.active.values = biomarkers.active[,2]

        res[drug.comb, biomarkers.active.names] = 1
      }

      # insert the inhibited biomarkers
      inhibited.biomarkers.file =
        paste0(biomarkers.dir, drug.comb, "_biomarkers_inhibited")

      if (file.exists(inhibited.biomarkers.file)
          && file.size(inhibited.biomarkers.file) != 0) {
        biomarkers.inhibited =
          read.table(inhibited.biomarkers.file, stringsAsFactors = FALSE)
        biomarkers.inhibited.names = biomarkers.inhibited[,1]
        # biomarkers.inhibited.values = biomarkers.inhibited[,2]

        res[drug.comb, biomarkers.inhibited.names] = -1
      }
    }

    return(res)
  }

#' Get synergy biomarkers per cell line
#'
#' Use this function to get the synergy biomarkers for each cell line.
#' The biomarkers must be stored in a single file inside each given cell line-specific
#' directory.
#'
#' @param biomarkers.dirs a character vector of the biomarker directories, in the
#' form of \emph{\{path\}/cell_line_name/\{dir\}}. The cell line name directory
#' should be different for each element of the vector as we use it to fill in the
#' \code{rownames} of each cell line-specific \code{data.frame} object.
#' Inside each \emph{\{dir\}}, we read the synergy biomarkers from a file (if it
#' exists and is non-empty) with the name \emph{biomarkers_per_synergy}. This file
#' has as first row the node names (columns) while every next row starts with the row name
#' (drug combination name) followed by a series of numbers from the ternary set
#' \{1,-1,0\}, denoting thus which nodes where found as active biomarkers for that
#' synergy, inhibited or not at all as biomarkers.
#'
#' @return a list of cell line-specific data frames (each element
#' from the list takes its name from the respective cell line).
#' Each cell-line specific \code{data.frame} object has as rows the
#' \strong{true positive predicted synergies} for that particular cell line
#' and columns the network nodes (should be the same for all cell lines).
#' Possible values for each \emph{synergy-node}
#' element in each cell line-specific \code{data.frame} are either \emph{1}
#' (\emph{active state} biomarker), \emph{-1}
#' (\emph{inhibited state} biomarker) or \emph{0} (not a biomarker).
#'
#' @importFrom utils read.table
#' @importFrom usefun get_parent_dir
#' @export
get_synergy_biomarkers_per_cell_line = function(biomarkers.dirs) {
  biomarkers.per.synergy = list()

  for (biomarker.dir in biomarkers.dirs) {
    biomarkers.file = paste0(biomarker.dir, "/biomarkers_per_synergy")
    cell.line = basename(get_parent_dir(biomarker.dir))

    if (file.exists(biomarkers.file) && file.size(biomarkers.file) != 0) {
      biomarkers.per.synergy[[cell.line]] =
        read.table(file = biomarkers.file, stringsAsFactors = FALSE,
                   check.names = FALSE)
    } else {
      biomarkers.per.synergy[[cell.line]] = data.frame()
    }
  }

  return(biomarkers.per.synergy)
}

#' Get performance biomarkers per cell line
#'
#' Use this function to get the performance biomarkers from the respective
#' files inside the given list of directories.
#'
#' @param biomarkers.dirs a character vector of the biomarker directories, in the
#' form of \emph{\{path\}/cell_line_name/\{dir\}}. The cell line name directory
#' should be different for each element of the vector as we use it to fill in the
#' \code{rownames} of the result \code{data.frame} object. Inside each \emph{\{dir\}}
#' (the directory name does not matter, but 'biomarkers' is a good choice),
#' we read the biomarkers from two files (if they exist and are non-empty):
#' \emph{biomarkers_active} and \emph{biomarkers_inhibited}, which have the
#' active and inhibited performance biomarkers for each cell line (these files
#' have a list of node names/biomarkers, one in each line).
#' @param node.names a character vector of the node names used in the analysis.
#' The biomarker names taken from the files inside the given directories must be
#' a subset of this vector.
#'
#' @return a data.frame, whose columns represent the network nodes and the
#' rows the cell lines. Possible values for each \emph{cell line-node}
#' element are either \emph{1} (\emph{active state} biomarker), \emph{-1}
#' (\emph{inhibited state} biomarker) or \emph{0} (not a biomarker).
#'
#' @importFrom utils read.table
#' @importFrom usefun get_parent_dir add_row_to_ternary_df
#' @export
get_perf_biomarkers_per_cell_line = function(biomarkers.dirs, node.names) {

  # initialize res data.frame
  res = as.data.frame(matrix(data = NA, nrow = 0, ncol = length(node.names)))
  colnames(res) = node.names

  for (biomarker.dir in biomarkers.dirs) {
    # active biomarkers
    active.biomarkers.file = paste0(biomarker.dir, "/biomarkers_active")
    if (file.exists(active.biomarkers.file)
        && file.size(active.biomarkers.file) != 0) {
      active.biomarkers = read.table(active.biomarkers.file, stringsAsFactors = FALSE)$V1
    } else {
      active.biomarkers = c()
    }

    # inhibited biomarkers
    inhibited.biomarkers.file = paste0(biomarker.dir, "/biomarkers_inhibited")
    if (file.exists(inhibited.biomarkers.file)
        && file.size(inhibited.biomarkers.file) != 0) {
      inhibited.biomarkers = read.table(inhibited.biomarkers.file, stringsAsFactors = FALSE)$V1
    } else {
      inhibited.biomarkers = c()
    }

    # cell line
    cell.line = basename(get_parent_dir(biomarker.dir))
    res = add_row_to_ternary_df(df = res, values.pos = active.biomarkers,
                                values.neg = inhibited.biomarkers,
                                row.name = cell.line)
  }

  return(res)
}
