#' Biomarker analysis based on TP model classification
#'
#' Use this function to perform a full biomarker analysis on an ensemble boolean model
#' dataset where the model classification is based on the number of \emph{true
#' positive} (TP) predictions. This analysis enables the discovery of \emph{performance
#' biomarkers}, nodes whose activity and/or boolean model parameterization (link
#' operator) affects the prediction performance of the models (as measured by
#' the number of TPs).
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models).
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names whereas the column names
#' specify the name of the network nodes (gene, proteins, etc.).
#' Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node). Note that the
#' rows (models) have to be in the same order as in the \code{model.predictions}
#' parameter.
#' @param models.link.operator a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names whereas the column names specify
#' the name of the network nodes (gene, proteins, etc.). Possible values for each
#' \emph{model-node element} are either \emph{0} (\strong{AND NOT} link operator),
#' \emph{1} (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted
#' by both activating and inhibiting regulators (no link operator). Default value:
#' NULL (no analysis on the models parameterization regarding the mutation of the
#' boolean equation link operator will be done).
#' @param observed.synergies a character vector with elements the names of the
#' drug combinations that were found as synergistic. This should be a subset of
#' the tested drug combinations, that is the column names of the \code{model.predictions}
#' parameter.
#' @param threshold numeric. A number in the [0,1] interval, above which (or
#' below its negative value) a biomarker will be registered in the returned result.
#' Values closer to 1 translate to a more strict threshold and thus less
#' biomarkers are found.
#' @param calculate.subsets.stats logical. If \emph{TRUE}, then the results will
#' include a vector of integers, representing the number of models that predicted
#' every subset of the given \code{observed.synergies} (where at least one model
#' predicts every synergy in the subset). The default value is \emph{FALSE}, since
#' the powerset of the predicted \code{observed.synergies} can be very large to compute.
#'
#' @return a list with various elements:
#' \itemize{
#'   \item \code{observed.model.predictions}: the part of the \code{model.predictions}
#'   data that includes the \code{observed.synergies}.
#'   \item \code{unobserved.model.predictions}: the complementary part of the
#'   \code{model.predictions} data that does not include the \code{observed.synergies}
#'   \item \code{predicted.synergies}: a character vector of the synergies (drug
#'   combination names) that were predicted by \strong{at least one} of the models
#'   in the dataset.
#'   \item \code{synergy.subset.stats}: an integer vector with elements the number
#'   of models the predicted each \strong{observed synergy subset} if the
#'   \emph{calculate.subsets.stats} option is enabled.
#'   \item \code{models.synergies.tp}: an integer vector of true positive (TP)
#'   values, one for each model.
#'   \item \code{diff.tp.mat}: a matrix whose rows are \strong{vectors of
#'   average node activity state differences} between two groups of models where
#'   the classification was based on the number of true positive predictions.
#'   Rows represent the different classification group matchings, e.g. (1,2) means the
#'   models that predicted 1 TP synergy vs the models that predicted 2 TP
#'   synergies and the columns represent the network's node names.
#'   Values are in the [-1,1] interval.
#'   \item \code{biomarkers.tp.active}: a character vector whose elements are
#'   the names of the \emph{active state} biomarkers. These nodes appear as more
#'   active in the better performance models.
#'   \item \code{biomarkers.tp.inhibited}: a character vector whose elements are
#'   the names of the \emph{inhibited state} biomarkers. These nodes appear as more
#'   inhibited in the better performance models.
#'   \item \code{diff.link.tp.mat}: a matrix whose rows are \strong{vectors of
#'   average node link operator differences} between two groups of models where
#'   the classification was based on the number of true positive predictions.
#'   Rows represent the different classification group matchings, e.g. (1,2) means the
#'   models that predicted 1 TP synergy vs the models that predicted 2 TP
#'   synergies and the columns represent the network's node names.
#'   Values are in the [-1,1] interval.
#'   \item \code{biomarkers.tp.or}: a character vector whose elements are
#'   the names of the \emph{OR} link operator biomarkers. These nodes have
#'   mostly the \emph{OR} link operator in their respective boolean equations
#'   in the better performance models.
#'   \item \code{biomarkers.tp.and}: a character vector whose elements are
#'   the names of the \emph{AND} link operator biomarkers. These nodes have
#'   mostly the \emph{AND} link operator in their respective boolean equations
#'   in the better performance models.
#' }
#'
#' @family general analysis functions
#' @export
biomarker_tp_analysis =
  function(model.predictions, models.stable.state, models.link.operator = NULL,
           observed.synergies, threshold, calculate.subsets.stats = FALSE) {
  # check input
  stopifnot(threshold >= 0 & threshold <= 1)
  models = rownames(model.predictions)
  stopifnot(all(models == rownames(models.stable.state)))

  # Split model.predictions to positive (observed) and negative (non-observed) results
  observed.model.predictions =
    get_observed_model_predictions(model.predictions, observed.synergies)
  unobserved.model.predictions =
    get_unobserved_model_predictions(model.predictions, observed.synergies)

  # check
  stopifnot(ncol(observed.model.predictions)
            + ncol(unobserved.model.predictions) == ncol(model.predictions))

  # get the predicted synergies (at least one model should predict it)
  predicted.synergies = names(which(colSums(observed.model.predictions, na.rm = TRUE) > 0))

  # check: the predicted synergies is a subset of the observed ones
  stopifnot(all(predicted.synergies %in% observed.synergies))

  # Find the number of predictive models for every synergy subset
  if (calculate.subsets.stats) {
    synergy.subset.stats = get_synergy_subset_stats(observed.model.predictions, predicted.synergies)
  }

  # Count the predictions of the observed synergies per model (TP)
  models.synergies.tp = calculate_models_synergies_tp(observed.model.predictions)

  # Make all possible classification group matchings and get the
  # average state differences
  diff.state.tp.mat = get_avg_activity_diff_mat_based_on_tp_predictions(
    models, models.synergies.tp, models.stable.state
  )

  # find the active and inhibited biomarkers based on the TP classification groups
  biomarkers.state.list = get_biomarkers(diff.state.tp.mat, threshold)

  # return all necessary data as elements of a list
  res.list = list()
  res.list$observed.model.predictions = observed.model.predictions
  res.list$unobserved.model.predictions = unobserved.model.predictions
  res.list$predicted.synergies = predicted.synergies
  if (calculate.subsets.stats) {
    res.list$synergy.subset.stats = synergy.subset.stats
  }
  res.list$models.synergies.tp = models.synergies.tp
  res.list$diff.state.tp.mat = diff.state.tp.mat
  res.list$biomarkers.tp.active = biomarkers.state.list$biomarkers.pos
  res.list$biomarkers.tp.inhibited = biomarkers.state.list$biomarkers.neg

  if (!is.null(models.link.operator)) {
    # check
    stopifnot(all(models == rownames(models.link.operator)))

    # Make all possible classification group matchings and get the average
    # link operator differences
    diff.link.tp.mat = get_avg_link_operator_diff_mat_based_on_tp_predictions(
      models, models.synergies.tp, models.link.operator
    )
    # find the 'OR' and 'AND' biomarkers based on the TP classification groups
    biomarkers.link.list = get_biomarkers(diff.link.tp.mat, threshold)

    res.list$diff.link.tp.mat = diff.link.tp.mat
    res.list$biomarkers.tp.or = biomarkers.link.list$biomarkers.pos
    res.list$biomarkers.tp.and = biomarkers.link.list$biomarkers.neg
  }

  return(res.list)
}

#' Biomarker analysis based on MCC model classification
#'
#' Use this function to perform a full biomarker analysis on an ensemble boolean model
#' dataset where the model classification is based on the \emph{Matthews correlation
#' coefficient score (MCC)}. This analysis enables the discovery of \emph{performance
#' biomarkers}, nodes whose activity and/or boolean model parameterization (link
#' operator) affects the prediction performance of the models (as measured by
#' the MCC score).
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models).
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names whereas the column names
#' specify the name of the network nodes (gene, proteins, etc.).
#' Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node). Note that the
#' rows (models) have to be in the same order as in the \code{model.predictions}
#' parameter.
#' @param models.link.operator a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names whereas the column names specify
#' the name of the network nodes (gene, proteins, etc.). Possible values for each
#' \emph{model-node element} are either \emph{0} (\strong{AND NOT} link operator),
#' \emph{1} (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted
#' by both activating and inhibiting regulators (no link operator). Default value:
#' NULL (no analysis on the models parameterization regarding the mutation of the
#' boolean equation link operator will be done).
#' @param observed.synergies a character vector with elements the names of the
#' drug combinations that were found as synergistic. This should be a subset of
#' the tested drug combinations, that is the column names of the \code{model.predictions}
#' parameter.
#' @param threshold numeric. A number in the [0,1] interval, above which (or
#' below its negative value) a biomarker will be registered in the returned result.
#' Values closer to 1 translate to a more strict threshold and thus less
#' biomarkers are found.
#' @param num.of.mcc.classes numeric. A positive integer larger than 2 that
#' signifies the number of mcc classes (groups) that we should split the models
#' MCC values (excluding the 'NaN' values).
#' @param include.NaN.mcc.class logical. Should the models that have NaN MCC value
#' (e.g. TP+FP = 0, models that predicted no synergies at all) be classified together
#' in one class - the 'NaN MCC Class' - and compared with the other model classes
#' in the analysis? If \emph{TRUE} (default), then the number of total MCC classes
#' will be \emph{num.of.mcc.classes + 1}.
#' @param calculate.subsets.stats logical. If \emph{TRUE}, then the results will
#' include a vector of integers, representing the number of models that predicted
#' every subset of the given \code{observed.synergies} (where at least one model
#' predicts every synergy in the subset). The default value is \emph{FALSE}, since
#' the powerset of the predicted \code{observed.synergies} can be very large to compute.
#'
#' @return a list with various elements:
#' \itemize{
#'   \item \code{observed.model.predictions}: the part of the \code{model.predictions}
#'   data that includes the \code{observed.synergies}.
#'   \item \code{unobserved.model.predictions}: the complementary part of the
#'   \code{model.predictions} data that does not include the \code{observed.synergies}
#'   \item \code{predicted.synergies}: a character vector of the synergies (drug
#'   combination names) that were predicted by \strong{at least one} of the models
#'   in the dataset.
#'   \item \code{synergy.subset.stats}: an integer vector with elements the number
#'   of models the predicted each \strong{observed synergy subset} if the
#'   \emph{calculate.subsets.stats} option is enabled.
#'   \item \code{models.mcc}: a numeric vector of MCC values (NaN's can be
#'   included), one for each model.
#'   \item \code{diff.state.mcc.mat}: a matrix whose rows are \strong{vectors of
#'   average node activity state differences} between two groups of models where
#'   the classification was based on the \emph{MCC score} of each model and was
#'   found using an optimal univariate k-means clustering method
#'   (\code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}).
#'   Rows represent the different classification group matchings, e.g. (1,2)
#'   means the models that were classified into the first MCC class vs the models
#'   that were classified in the 2nd class (higher is better). The columns
#'   represent the network's node names. Values are in the [-1,1] interval.
#'   \item \code{biomarkers.mcc.active}: a character vector whose elements are
#'   the names of the \emph{active state} biomarkers. These nodes appear more
#'   active in the better performance models.
#'   \item \code{biomarkers.mcc.inhibited}: a character vector whose elements are
#'   the names of the \emph{inhibited state} biomarkers. These nodes appear more
#'   inhibited in the better performance models.
#'   \item \code{diff.link.mcc.mat}: a matrix whose rows are \strong{vectors of
#'   average node link operator differences} between two groups of models where
#'   the classification was based on the \emph{MCC score} of each model and was
#'   found using an optimal univariate k-means clustering method
#'   (\code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}).
#'   Rows represent the different classification group matchings, e.g. (1,2)
#'   means the models that were classified into the first MCC class vs the models
#'   that were classified in the 2nd class (higher is better).
#'   The columns represent the network's node names. Values are in the [-1,1] interval.
#'   \item \code{biomarkers.mcc.or}: a character vector whose elements are
#'   the names of the \emph{OR} link operator biomarkers. These nodes have
#'   mostly the \emph{OR} link operator in their respective boolean equations
#'   in the better performance models.
#'   \item \code{biomarkers.mcc.and}: a character vector whose elements are
#'   the names of the \emph{AND} link operator biomarkers. These nodes have
#'   mostly the \emph{AND} link operator in their respective boolean equations
#'   in the better performance models.
#' }
#'
#' @family general analysis functions
#' @export
biomarker_mcc_analysis = function(model.predictions, models.stable.state,
  models.link.operator = NULL, observed.synergies, threshold, num.of.mcc.classes,
  include.NaN.mcc.class = TRUE, calculate.subsets.stats = FALSE) {

  # check input
  stopifnot(threshold >= 0 & threshold <= 1)
  stopifnot(num.of.mcc.classes >= 2)
  number.of.drug.comb.tested = ncol(model.predictions)
  models = rownames(model.predictions)

  stopifnot(all(models == rownames(models.stable.state)))

  # Split model.predictions to positive (observed) and negative (non-observed) results
  observed.model.predictions =
    get_observed_model_predictions(model.predictions, observed.synergies)
  unobserved.model.predictions =
    get_unobserved_model_predictions(model.predictions, observed.synergies)

  # check
  stopifnot(ncol(observed.model.predictions)
            + ncol(unobserved.model.predictions) == ncol(model.predictions))

  # get the predicted synergies (at least one model should predict it)
  predicted.synergies = names(which(colSums(observed.model.predictions, na.rm = TRUE) > 0))

  # check: the predicted synergies is a subset of the observed ones
  stopifnot(all(predicted.synergies %in% observed.synergies))

  # Find the number of predictive models for every synergy subset
  if (calculate.subsets.stats) {
    synergy.subset.stats = get_synergy_subset_stats(observed.model.predictions, predicted.synergies)
  }

  # Calculate Matthews Correlation Coefficient (MCC) for every model
  models.mcc = calculate_models_mcc(observed.model.predictions,
                                    unobserved.model.predictions,
                                    number.of.drug.comb.tested)

  # Make all possible classification group matchings and get the
  # average state differences
  diff.state.mcc.mat = get_avg_activity_diff_mat_based_on_mcc_clustering(
    models.mcc, models.stable.state, num.of.mcc.classes, include.NaN.mcc.class
  )

  # find the active and inhibited biomarkers based on the MCC classification groups
  biomarkers.state.list = get_biomarkers(diff.state.mcc.mat, threshold)

  # return all necessary data as elements of a list
  res.list = list()
  res.list$observed.model.predictions = observed.model.predictions
  res.list$unobserved.model.predictions = unobserved.model.predictions
  res.list$predicted.synergies = predicted.synergies
  if (calculate.subsets.stats) {
    res.list$synergy.subset.stats = synergy.subset.stats
  }
  res.list$models.mcc = models.mcc
  res.list$diff.state.mcc.mat = diff.state.mcc.mat
  res.list$biomarkers.mcc.active = biomarkers.state.list$biomarkers.pos
  res.list$biomarkers.mcc.inhibited = biomarkers.state.list$biomarkers.neg

  if (!is.null(models.link.operator)) {
    # check
    stopifnot(all(models == rownames(models.link.operator)))

    # Make all possible classification group matchings and get the average
    # link operator differences
    diff.link.mcc.mat = get_avg_link_operator_diff_mat_based_on_mcc_clustering(
      models.mcc, models.link.operator, num.of.mcc.classes, include.NaN.mcc.class
    )
    # find the 'OR' and 'AND' biomarkers based on the TP classification groups
    biomarkers.link.list = get_biomarkers(diff.link.mcc.mat, threshold)

    res.list$diff.link.mcc.mat = diff.link.mcc.mat
    res.list$biomarkers.mcc.or = biomarkers.link.list$biomarkers.pos
    res.list$biomarkers.mcc.and = biomarkers.link.list$biomarkers.neg
  }

  return(res.list)
}

#' Biomarker analysis per synergy predicted
#'
#' Use this function to discover \emph{synergy biomarkers}, i.e. nodes whose
#' activity and/or boolean equation parameterization (link operator) affect the
#' manifestation of synergies in the models. Models are classified based on
#' whether they predict or not each of the predicted synergies.
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models).
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names whereas the column names
#' specify the name of the network nodes (gene, proteins, etc.).
#' Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node). Note that the
#' rows (models) have to be in the same order as in the \code{model.predictions}
#' parameter.
#' @param models.link.operator a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names whereas the column names specify
#' the name of the network nodes (gene, proteins, etc.). Possible values for each
#' \emph{model-node element} are either \emph{0} (\strong{AND NOT} link operator),
#' \emph{1} (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted
#' by both activating and inhibiting regulators (no link operator). Default value:
#' NULL (no analysis on the models parameterization regarding the mutation of the
#' boolean equation link operator will be done).
#' @param observed.synergies a character vector with elements the names of the
#' drug combinations that were found as synergistic. This should be a subset of
#' the tested drug combinations, that is the column names of the \code{model.predictions}
#' parameter.
#' @param threshold numeric. A number in the [0,1] interval, above which (or
#' below its negative value) a biomarker will be registered in the returned result.
#' Values closer to 1 translate to a more strict threshold and thus less
#' biomarkers are found.
#' @param calculate.subsets.stats logical. If \emph{TRUE}, then the results will
#' include a vector of integers, representing the number of models that predicted
#' every subset of the given \code{observed.synergies} (where at least one model
#' predicts every synergy in the subset). The default value is \emph{FALSE}, since
#' the powerset of the predicted \code{observed.synergies} can be very large to compute.
#'
#' @return a list with various elements:
#' \itemize{
#'   \item \code{observed.model.predictions}: the part of the \code{model.predictions}
#'   data that includes the \code{observed.synergies}.
#'   \item \code{unobserved.model.predictions}: the complementary part of the
#'   \code{model.predictions} data that does not include the \code{observed.synergies}
#'   \item \code{predicted.synergies}: a character vector of the synergies (drug
#'   combination names) that were predicted by \strong{at least one} of the models
#'   in the dataset.
#'   \item \code{synergy.subset.stats}: an integer vector with elements the number
#'   of models the predicted each \strong{observed synergy subset} if the
#'   \emph{calculate.subsets.stats} option is enabled.
#'   \item \code{diff.state.synergies.mat}: a matrix whose rows are
#'   \strong{vectors of average node activity state differences} between two
#'   groups of models where the classification for each individual row was based
#'   on the prediction or not of a specific synergistic drug combination. The
#'   row names are the predicted synergies, one per row, while the columns
#'   represent the network's node names. Values are in the [-1,1] interval.
#'   \item \code{activity.biomarkers}: a \code{data.frame} object with rows
#'   the \code{predicted synergies} and columns the nodes (column names of the
#'   \code{models.stable.states} matrix). Possible values for each
#'   \emph{synergy-node} element are either \emph{1} (\emph{active state}
#'   biomarker), \emph{-1} (\emph{inhibited state} biomarker) or \emph{0} (not
#'   a biomarker) for the given \code{threshold} value.
#'   \item \code{diff.link.synergies.mat}: a matrix whose rows are
#'   \strong{vectors of average node link operator differences} between two
#'   groups of models where the classification for each individual row was
#'   based on the prediction or not of a specific synergistic drug combination.
#'   The row names are the predicted synergies, one per row, while the columns
#'   represent the network's node names. Values are in the [-1,1] interval.
#'   \item \code{link.operator.biomarkers}: a \code{data.frame} object with rows
#'   the \code{predicted synergies} and columns the nodes (column names of the
#'   \code{models.link.operator} matrix). Possible values for each
#'   \emph{synergy-node} element are either \emph{1} (\emph{OR} link operator
#'   biomarker), \emph{-1} (\emph{AND} link operator biomarker) or \emph{0} (not
#'   a biomarker) for the given \code{threshold} value.
#' }
#'
#' @family general analysis functions
#'
#' @importFrom usefun get_ternary_class_id
#' @export
biomarker_synergy_analysis =
  function(model.predictions, models.stable.state, models.link.operator = NULL,
           observed.synergies, threshold, calculate.subsets.stats = FALSE) {
    # check input
    stopifnot(threshold >= 0 & threshold <= 1)
    models = rownames(model.predictions)
    stopifnot(all(models == rownames(models.stable.state)))

    # Split model.predictions to positive (observed) and negative (non-observed) results
    observed.model.predictions =
      get_observed_model_predictions(model.predictions, observed.synergies)
    unobserved.model.predictions =
      get_unobserved_model_predictions(model.predictions, observed.synergies)

    # check
    stopifnot(ncol(observed.model.predictions)
              + ncol(unobserved.model.predictions) == ncol(model.predictions))

    # get the predicted synergies (at least one model should predict it)
    predicted.synergies = names(which(colSums(observed.model.predictions, na.rm = TRUE) > 0))

    # check: the predicted synergies is a subset of the observed ones
    stopifnot(all(predicted.synergies %in% observed.synergies))

    # Find the number of predictive models for every synergy subset
    if (calculate.subsets.stats) {
      synergy.subset.stats = get_synergy_subset_stats(observed.model.predictions, predicted.synergies)
    }

    # get the average activity state differences for each predicted synergy
    diff.state.synergies.mat =
      get_avg_activity_diff_mat_based_on_specific_synergy_prediction(
        model.predictions, models.stable.state, predicted.synergies
      )

    # find the active and inhibited biomarkers for each predicted synergy
    activity.biomarkers = as.data.frame(
      apply(diff.state.synergies.mat, c(1,2), get_ternary_class_id, threshold))

    # return all necessary data as elements of a list
    res.list = list()
    res.list$observed.model.predictions = observed.model.predictions
    res.list$unobserved.model.predictions = unobserved.model.predictions
    res.list$predicted.synergies = predicted.synergies
    if (calculate.subsets.stats) {
      res.list$synergy.subset.stats = synergy.subset.stats
    }
    res.list$diff.state.synergies.mat = diff.state.synergies.mat
    res.list$activity.biomarkers = activity.biomarkers

    if (!is.null(models.link.operator)) {
      # check
      stopifnot(all(models == rownames(models.link.operator)))

      # get the average link operator differences for each predicted synergy
      diff.link.synergies.mat =
        get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction(
          model.predictions, models.link.operator, predicted.synergies
        )

      # find the 'OR' and 'AND' biomarkers for each predicted synergy
      link.operator.biomarkers = as.data.frame(
        apply(diff.link.synergies.mat, c(1,2), get_ternary_class_id, threshold))

      res.list$diff.link.synergies.mat = diff.link.synergies.mat
      res.list$link.operator.biomarkers = link.operator.biomarkers
    }

    return(res.list)
  }
