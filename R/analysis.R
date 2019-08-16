#' Subset the model predictions to the (true) observed synergies
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param observed.synergies a character vector with elements the names of the
#' drug combinations that were found as synergistic
#'
#' @return a \code{data.frame} object with rows the models
#' and columns the drug combinations that were found/observed as \strong{synergistic}
#' (\emph{positive results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#'
#' @export
get_observed_model_predictions = function(model.predictions, observed.synergies) {
  drug.combinations.tested = colnames(model.predictions)
  return(model.predictions[,sapply(drug.combinations.tested, function(drug.comb) {
    is_comb_element_of(drug.comb, observed.synergies)
  })])
}

#' Subset the model predictions to the (false) non-observed synergies
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param observed.synergies a character vector with elements the names of the
#' drug combinations that were found as synergistic
#'
#' @return a \code{data.frame} object with rows the models
#' and columns the drug combinations that were found/observed as \strong{non-synergistic}
#' (\emph{negative results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#'
#' @export
get_unobserved_model_predictions = function(model.predictions, observed.synergies) {
  drug.combinations.tested = colnames(model.predictions)
  return(model.predictions[,sapply(drug.combinations.tested, function(drug.comb) {
    !is_comb_element_of(drug.comb, observed.synergies)
  })])
}

#' Find the number of predictive models for every synergy subset
#'
#' Use this function to find for each possible subset of drug combinations out
#' of a given list of synergies, the number of models that predicted it given
#' the models' predictions. So, if for example the set of synergies is this one:
#' \{'A-B','C-D','E-F'\}, we want to know how many models predicted none of them,
#' just the single subsets (e.g. the \{'A-B'\}),
#' the two-element subsets (e.g. the \{'A-B','C-D'\}) and all 3 of them.
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models).
#' @param synergies a character vector with elements the synergistic drug
#' combinations. Note that these synergies should be a subset of the column
#' names of the \code{model.predictions} data.frame.
#'
#' @return an integer vector with elements the number of models the predicted
#' each synergy subset. The \emph{names} attribute has the names of each
#' synergistic drug combination subset, which are the drug combinations comma
#' separated (e.g. 'A-B,C-D').
#'
#' @importFrom rje powerSet
#' @export
get_synergy_subset_stats = function(model.predictions, synergies) {
  # check: the predicted synergies is a subset of the observed (positive) ones
  stopifnot(all(synergies %in% colnames(model.predictions)))

  # get the powerset
  synergies.powerset = powerSet(synergies)

  # order by subset size and give the list better names
  synergies.powerset = synergies.powerset[
    order(sapply(synergies.powerset, length))
  ]
  names(synergies.powerset) =
    sapply(synergies.powerset, function(drug.comb.set) {
      paste(drug.comb.set, collapse = ",")
    })

  synergy.subset.stats =
    sapply(synergies.powerset, function(drug.comb.set) {
      count_models_that_predict_synergies(drug.comb.set, model.predictions)
    })

  return(synergy.subset.stats)
}

#' Count models that predict list of synergies
#'
#' Use this function to find the number of models that predict a given list of
#' drug combinations (usually the ones found as synergies).
#'
#' @param drug.comb.list a list with (synergistic) drug combinations as elements
#' (each drug combination is a string in the form \emph{A-B} - no spaces between
#' the names and the hyphen '-')
#' @param model.predictions a \code{data.frame} object with rows as the models
#' and columns the drug combinations tested. Possible values for each
#' \emph{model-drug combination element} are either \emph{0} (no synergy
#' predicted), \emph{1} (synergy was predicted) or \emph{NA}
#'
#' @return the number of models that predict the given drug combination set
#' (have a value of 1 in the respective columns of the \code{model.predictions}
#' data.frame). If the given set is empty, we return the number of models that
#' predicted no synergies at all (after the \emph{NA} values are discarded, the
#' row in the \code{model.predictions} data.frame is all zeros)
#'
#' @export
count_models_that_predict_synergies =
  function(drug.comb.list, model.predictions) {
    synergy.vector = unlist(drug.comb.list)
    if (length(synergy.vector) == 0) {
      count = sum(apply(model.predictions, 1, function(x) {
        all(x == 0, na.rm = T)
      }))
    } else if (length(synergy.vector) == 1) {
      count = sum(model.predictions[, synergy.vector], na.rm = T)
    } else {
      count = sum(apply(model.predictions[, synergy.vector], 1, function(x) {
        all(x == 1)
      }), na.rm = T)
    }

    return(count)
}

#' Calculate the Matthews correlation coefficient for each model
#'
#' @param observed.model.predictions \code{data.frame} object with rows the models
#' and columns the drug combinations that were found as \strong{synergistic}
#' (\emph{positive results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param unobserved.model.predictions \code{data.frame} object with rows the models
#' and columns the drug combinations that were found as \strong{non-synergistic}
#' (\emph{negative results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param number.of.drug.comb.tested numeric. The total number of drug
#' combinations tested, which should be equal to the sum of the columns of the
#' \code{observed.model.predictions} and the \code{unobserved.model.predictions}.
#'
#' @return a numeric vector of MCC values, each value being in the [-1,1]
#' interval or \emph{NaN}. The \emph{names} attribute holds the models' names.
#'
#' @family confusion matrix calculation functions
#'
#' @export
calculate_models_mcc = function(observed.model.predictions, unobserved.model.predictions,
                                number.of.drug.comb.tested) {
    # Count the true positives (TP)
    models.synergies.tp = calculate_models_synergies_tp(observed.model.predictions)

    # Count the false negatives (FN)
    models.synergies.fn = calculate_models_synergies_fn(observed.model.predictions)

    # Count the False Positives (FP)
    models.synergies.fp = calculate_models_synergies_fp(unobserved.model.predictions)

    # Count the True Negatives (TN)
    models.synergies.tn = calculate_models_synergies_tn(unobserved.model.predictions)

    # P = TP + FN (Positives)
    positives = ncol(observed.model.predictions)
    models.synergies.p = models.synergies.tp + models.synergies.fn

    # N = FP + TN (Negatives)
    negatives = ncol(unobserved.model.predictions)
    models.synergies.n = models.synergies.fp + models.synergies.tn

    # checks
    stopifnot(models.synergies.p == positives)
    stopifnot(models.synergies.n == negatives)
    stopifnot(positives + negatives == number.of.drug.comb.tested)

    # Calculate Matthews Correlation Coefficient (MCC)
    models.mcc = calculate_mcc(models.synergies.tp, models.synergies.tn,
                               models.synergies.fp, models.synergies.fn,
                               positives, negatives)
    return(models.mcc)
}


#' Count the predictions of the observed synergies per model (TP)
#'
#' Since the given \code{observed.model.predictions} data.frame has only the
#' positive results, this function returns the total number of 1's in each row.
#'
#' @param observed.model.predictions \code{data.frame} object with rows the models
#' and columns the drug combinations that were found/observed as \strong{synergistic}
#' (\emph{positive results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#'
#' @return an integer vector with elements the number of true positive predictions
#' per model. The model names are given in the \emph{names} attribute (same order
#' as in the \emph{rownames} attribute of the observed.model.predictions
#' \code{data.frame}).
#'
#' @family confusion matrix calculation functions
#'
#' @export
calculate_models_synergies_tp = function(observed.model.predictions) {
  return(apply(observed.model.predictions, 1, sum, na.rm = T))
}

#' Count the non-synergies of the observed synergies per model (FN)
#'
#' Since the given \code{observed.model.predictions} data.frame has only the
#' positive results, this function returns the total number of 0's \emph{and}
#' NA's in each row.
#'
#' @param observed.model.predictions \code{data.frame} object with rows the models
#' and columns the drug combinations that were found/observed as \strong{synergistic}
#' (\emph{negative results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#'
#' @return an integer vector with elements the number of false negative predictions
#' per model. The model names are given in the \emph{names} attribute (same order
#' as in the \emph{rownames} attribute of the observed.model.predictions
#' \code{data.frame}).
#'
#' @family confusion matrix calculation functions
#'
#' @export
calculate_models_synergies_fn = function(observed.model.predictions) {
  return(apply(observed.model.predictions, 1, function(x) {
    sum( x == 0 | is.na(x) )
  }))
}

#' Count the predictions of the non-synergistic drug combinations per model (FP)
#'
#' Since the given \code{unobserved.model.predictions} data.frame has only the
#' negative results, this function returns the total number of 1's in each row.
#'
#' @param unobserved.model.predictions \code{data.frame} object with rows the models
#' and columns the drug combinations that were found/observed as \strong{non-synergistic}
#' (\emph{negative results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#'
#' @return an integer vector with elements the number of false positive predictions
#' per model. The model names are given in the \emph{names} attribute (same order
#' as in the \emph{rownames} attribute of the unobserved.model.predictions
#' \code{data.frame}).
#'
#' @family confusion matrix calculation functions
#'
#' @export
calculate_models_synergies_fp = function(unobserved.model.predictions) {
  return(apply(unobserved.model.predictions, 1, sum, na.rm = T))
}

#' Count the non-synergies of the non-synergistic drug combinations per model (TN)
#'
#' Since the given \code{unobserved.model.predictions} data.frame has only the
#' negative results, this function returns the total number of 0's \emph{and}
#' NA's in each row.
#'
#' @param unobserved.model.predictions \code{data.frame} object with rows the models
#' and columns the drug combinations that were found/observed as \strong{non-synergistic}
#' (\emph{negative results}). Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#'
#' @return an integer vector with elements the number of true negative predictions
#' per model. The model names are given in the \emph{names} attribute (same order
#' as in the \emph{rownames} attribute of the unobserved.model.predictions
#' \code{data.frame}).
#'
#' @family confusion matrix calculation functions
#'
#' @export
calculate_models_synergies_tn = function(unobserved.model.predictions) {
  return(apply(unobserved.model.predictions, 1, function(x) {
    sum( x == 0 | is.na(x))
  }))
}

#' Calculate the Matthews correlation coefficient vector
#'
#' Use this function to calculate the MCC values given vectors of \emph{TP} (true
#' positives), \emph{FP} (false positives), \emph{TN} (true negatives), \emph{FN}
#' (false negatives), \emph{P} (positives) and \emph{N} (negatives). Note that
#' the input vectors have to be of the same size and have one-to-one value
#' correspondence for the output MCC vector values to make sense.
#'
#' @param tp numeric vector of TPs
#' @param tn numeric vector of TNs
#' @param fp numeric vector of FPs
#' @param fn numeric vector of FNs
#' @param p numeric vector of Ps (p = tp + fn)
#' @param n numeric vector of Ns (n = tn + fp)
#'
#' @return a numeric vector of MCC values, each value being in the [-1,1]
#' interval or \emph{NaN}.
#'
#' @family confusion matrix calculation functions
#'
#' @export
calculate_mcc = function(tp, tn, fp, fn, p, n) {
  return(
    (tp * tn - fp * fn) / sqrt((tp + fp) * p * n * (tn + fn))
  )
}

#' Update biomarker files for a specific synergy
#'
#' This function gets the (previously-found or 'old') synergy biomarkers from their
#' respective files and if any of these files are empty (no 'old' biomarkers
#' found) or non-existent, the 'new' biomarkers (given as input vector parameters) are
#' automatically saved. When the 'new' biomarkers \strong{share common nodes}
#' with the 'old' biomarkers, there exist 3
#' possible ways to combine the results, given by the \code{method} parameter.
#' If \strong{no common nodes} exist, no matter the \code{method} selected,
#' the 'new' biomarkers are added to the 'old' ones.
#'
#' @param biomarkers.dir string. It specifies the full path name of the
#' directory which holds the biomarker files for the synergistic drug combination
#' specified in the parameter \code{drug.comb}. The biomarker files must be
#' formatted as: \emph{\%drug.comb\%_biomarkers_active} or
#' \emph{\%drug.comb\%_biomarkers_inhibited}, where \%drug.comb\% is the value
#' of the \code{drug.comb} parameter.
#' @param drug.comb string. The drug combination (e.g. "A-B") that will be used
#' to identify the related biomarker files.
#' @param biomarkers.active.new a numeric vector whose \emph{names} attribute
#' includes the node names of the (newly found) \emph{active biomarkers} for the specified
#' synergy. The values of the vector are the average activity difference of
#' each node, derived from a comparison between 2 different groups of models.
#' @param biomarkers.inhibited.new a numeric vector whose \emph{names} attribute
#' includes the node names of the (newly found) \emph{inhibited biomarkers} for the specified
#' synergy. The values of the vector are the average activity difference of
#' each node, derived from a comparison between 2 different groups of models.
#' @param method string. It specifies the method to use to update the biomarker
#' files when there are \emph{common nodes} between the 'old' and 'new' biomarkers:
#' \enumerate{
#'   \item \code{replace}(DEFAULT): we discard the 'old' biomarkers and keep
#'   only the 'new' ones
#'   \item \code{prune.to.common}: we keep only the common biomarkers
#'   \item \code{extend}: we add to the 'old' set of biomarkers the extra ones
#'   from the 'new' set that are not non-common to the 'old' ones, extending
#'   thus the 'old' biomarker set
#' }
#'
#' @importFrom usefun save_vector_to_file add_vector_to_df
#' @importFrom utils read.table
#' @export
update_biomarker_files =
  function(biomarkers.dir, drug.comb, biomarkers.active.new,
           biomarkers.inhibited.new, method = "replace") {
    # check method
    stopifnot(method %in% c("replace", "prune.to.common", "extend"))

    # update the active biomarkers
    active.biomarkers.file =
      paste0(biomarkers.dir, drug.comb, "_biomarkers_active")

    if (file.size(active.biomarkers.file) == 0
        || !file.exists(active.biomarkers.file)) {
      save_vector_to_file(vector = biomarkers.active.new,
                          file = active.biomarkers.file,
                          with.row.names = TRUE)
    } else {
      biomarkers.active.prev =
        read.table(active.biomarkers.file, stringsAsFactors = FALSE)
      biomarkers.active.prev.names = biomarkers.active.prev[,1]
      biomarkers.active.new.names = names(biomarkers.active.new)

      biomarkers.active.common = intersect(biomarkers.active.prev.names,
                                           biomarkers.active.new.names)

      if (length(biomarkers.active.common) == 0) {
        biomarkers.active = add_vector_to_df(biomarkers.active.prev,
                                             biomarkers.active.new)
        biomarkers.active = transform(biomarkers.active,
                                      V2 = as.numeric(biomarkers.active$V2))
        save_vector_to_file(vector = biomarkers.active,
                            file = active.biomarkers.file)
      } else {
        if (method == "replace") {
          biomarkers.active = biomarkers.active.new
          save_vector_to_file(vector = biomarkers.active,
                              file = active.biomarkers.file,
                              with.row.names = TRUE)
        } else if (method == "extend") {
          # find the extra 'new' biomarkers and add them to the 'old' ones
          biomarkers.active.to.add = biomarkers.active.new[
            !(biomarkers.active.new.names %in% biomarkers.active.prev.names)
            ]
          biomarkers.active = add_vector_to_df(biomarkers.active.prev,
                                               biomarkers.active.to.add)
          biomarkers.active = transform(biomarkers.active,
                                        V2 = as.numeric(biomarkers.active$V2))
          save_vector_to_file(vector = biomarkers.active,
                              file = active.biomarkers.file)
        } else if (method == "prune.to.common") {
          biomarkers.active = biomarkers.active.new[biomarkers.active.common]
          save_vector_to_file(vector = biomarkers.active,
                              file = active.biomarkers.file,
                              with.row.names = TRUE)
        }
      }
    }

    # update the inhibited biomarkers
    inhibited.biomarkers.file =
      paste0(biomarkers.dir, drug.comb, "_biomarkers_inhibited")

    if (file.size(inhibited.biomarkers.file) == 0) {
      save_vector_to_file(vector = biomarkers.inhibited.new,
                          file = inhibited.biomarkers.file,
                          with.row.names = TRUE)
    } else {
      biomarkers.inhibited.prev =
        read.table(inhibited.biomarkers.file, stringsAsFactors = FALSE)
      biomarkers.inhibited.prev.names = biomarkers.inhibited.prev[,1]
      biomarkers.inhibited.new.names = names(biomarkers.inhibited.new)

      biomarkers.inhibited.common = intersect(biomarkers.inhibited.prev.names,
                                              biomarkers.inhibited.new.names)

      if (length(biomarkers.inhibited.common) == 0) {
        biomarkers.inhibited = add_vector_to_df(biomarkers.inhibited.prev,
                                                biomarkers.inhibited.new)
        biomarkers.inhibited = transform(biomarkers.inhibited,
                                         V2 = as.numeric(biomarkers.inhibited$V2))
        save_vector_to_file(vector = biomarkers.inhibited,
                            file = inhibited.biomarkers.file)
      } else {
        if (method == "replace") {
          biomarkers.inhibited = biomarkers.inhibited.new
          save_vector_to_file(vector = biomarkers.inhibited,
                              file = inhibited.biomarkers.file,
                              with.row.names = TRUE)
        } else if (method == "extend") {
          # find the extra 'new' biomarkers and add them to the 'old' ones
          biomarkers.inhibited.to.add = biomarkers.inhibited.new[
            !(biomarkers.inhibited.new.names %in% biomarkers.inhibited.prev.names)
            ]
          biomarkers.inhibited = add_vector_to_df(biomarkers.inhibited.prev,
                                                  biomarkers.inhibited.to.add)
          biomarkers.inhibited = transform(biomarkers.inhibited,
                                           V2 = as.numeric(biomarkers.inhibited$V2))
          save_vector_to_file(vector = biomarkers.inhibited,
                              file = inhibited.biomarkers.file)
        } else if (method == "prune.to.common") {
          biomarkers.inhibited = biomarkers.inhibited.new[biomarkers.inhibited.common]
          save_vector_to_file(vector = biomarkers.inhibited,
                              file = inhibited.biomarkers.file,
                              with.row.names = TRUE)
        }
      }
    }
  }

# merge the results of the performance (active and inhibited) biomarkers
# from each cell line to a common `data.frame` object
merge_perf_biomarkers =
  function(node.names, cell.lines, biomarkers.perf.active, biomarkers.perf.inhibited) {
    # initialize res data.frame
    res = as.data.frame(matrix(0, ncol = length(node.names),
                                  nrow = length(cell.lines)))
    colnames(res) = node.names
    rownames(res) = cell.lines

    for (cell.line in cell.lines) {
      biomarkers.perf.active.vec = unlist(biomarkers.perf.active[cell.line])
      biomarkers.perf.inhibited.vec = unlist(biomarkers.perf.inhibited[cell.line])

      for (biomarker.active in biomarkers.perf.active.vec) {
        res[cell.line, biomarker.active] = 1
      }

      for (biomarker.inhibited in biomarkers.perf.inhibited.vec) {
        res[cell.line, biomarker.inhibited] = -1
      }
    }

    return(res)
}

# merge the observed synergies from each cell line to a common
# `data.frame` object
merge_observed_synergies =
  function(drug.combinations.tested, cell.lines, observed.synergies.per.cell.line) {
    # initialize res data.frame
    res = as.data.frame(matrix(0, ncol = length(drug.combinations.tested),
                                  nrow = length(cell.lines)))
    colnames(res) = drug.combinations.tested
    rownames(res) = cell.lines

    for (cell.line in cell.lines) {
      observed.synergies = observed.synergies.per.cell.line[cell.line]
      for (synergy in observed.synergies) {
        res[cell.line, synergy] = 1
      }
    }

    return(res)
}

# `df.list` is a list of cell line data frames with rows the true positive
# predicted synergies for each cell line and columns the network nodes.
# `observed.synergies.res` is a data.frame with rows the cell lines and columns
# the observed synergies
# `predicted.synergies.vector` is a subset of the `observed.synergies.res` column
# names
# Returns a list of data frames for each synergy with rows the cell lines and
# columns the network nodes, so a re-arrangement (and grouping) of the `df.list`
# object
#' @importFrom dplyr select
arrange_by_synergy =
  function(df.list, observed.synergies.res, predicted.synergies.vector, node.names) {
    stopifnot(all(predicted.synergies.vector %in%
                  colnames(observed.synergies.res)))

    # prune the observed synergies to predicted and order by increasing number
    # of cell line predictions
    observed.synergies.res =
      select(observed.synergies.res, predicted.synergies.vector)
    observed.synergies.res =
      observed.synergies.res[order(colSums(observed.synergies.res))]

    # check that `node.names` is in the correct order
    for (df in df.list) {
      stopifnot(all(node.names == colnames(df)))
    }

    # initialize `res` list (correct dimensions, all values zero)
    res = list()
    for (synergy in colnames(observed.synergies.res)) {
      row.df = select(observed.synergies.res, synergy)

      synergy.res = as.data.frame(matrix(0, ncol = length(node.names),
                                            nrow = sum(row.df)))
      colnames(synergy.res) = node.names
      rownames(synergy.res) = rownames(row.df)[row.df == 1]
      res[[synergy]] = synergy.res
    }

    # fill in `res` list
    for (cell.line in names(df.list)) {
      df = df.list[[cell.line]]
      for (synergy in rownames(df)) {
        res[[synergy]][cell.line,] = df[synergy, ]
      }
    }

    return(res)
}

#' Get synergy comparison sets
#'
#' This helper function identifies pairs of (\emph{set}, \emph{subset}) for each
#' synergy (implicity given throught the \code{synergy.subset.stats} object) where
#' each respective \emph{subset} misses just one synergy from the larger \emph{set}.
#'
#' @param synergy.subset.stats integer vector with values the amount of models
#' that predicted each synergy subset, defined as a comma-seperated string of
#' drug combinations in the \emph{names} attribute of the vector. It can be the
#' result of using the function \code{\link{get_synergy_subset_stats}}.
#'
#' @return \code{data.frame} object with 3 columns. For each row, the 1st column defines a
#' \emph{single synergy} of interest (e.g. drug combination "A-B"), the 2nd a
#' \emph{synergy set} that includes the single one (e.g. the set "F-G,A-B,C-D")
#' and the 3rd the \emph{synergy subset} of the \emph{set} that does not include
#' the single synergy of the first column (e.g. "F-G,C-D").
#'
#' @importFrom usefun outersect
#' @export
get_synergy_comparison_sets = function(synergy.subset.stats) {
  # keep only the synergy sets where we have at least one model predicting them
  synergy.sets = synergy.subset.stats[synergy.subset.stats > 0]

  # remove the zero set (models that predicted none of the synergy sets)
  synergy.sets = synergy.sets[!names(synergy.sets) == ""]

  # find the single drug combinations
  synergy.set.sizes = numeric(0)
  for (set in names(synergy.sets)) {
    synergy.set.size = length(unlist(strsplit(set, split = ",")))
    synergy.set.sizes = c(synergy.set.sizes, synergy.set.size)
  }

  # get the maximum size of a synergy set
  max.size = max(synergy.set.sizes)

  synergies = NULL
  sets      = NULL
  subsets   = NULL

  for (size in 1:(max.size - 1)) {
    small.size = size
    large.size = size + 1

    small.synergy.sets = names(synergy.sets[synergy.set.sizes == small.size])
    large.synergy.sets = names(synergy.sets[synergy.set.sizes == large.size])

    for (small.synergy.set in small.synergy.sets) {
      small.synergies = unlist(strsplit(small.synergy.set, split = ","))
      for (large.synergy.set in large.synergy.sets) {
        large.synergies = unlist(strsplit(large.synergy.set, split = ","))
        if (all(small.synergies %in% large.synergies) &
            synergy.sets[small.synergy.set] > synergy.sets[large.synergy.set]) {
          synergy.to.test = outersect(small.synergies, large.synergies)

          synergies = append(synergies, synergy.to.test)
          sets      = append(sets, large.synergy.set)
          subsets   = append(subsets, small.synergy.set)
        }
      }
    }
  }

  res.df = data.frame(synergies, sets, subsets, stringsAsFactors = FALSE)
  res.df = res.df[order(res.df$synergies), ]
  rownames(res.df) = 1:length(res.df$synergies)

  return(res.df)
}

