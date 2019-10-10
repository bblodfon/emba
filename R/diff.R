#' Get average activity difference matrix based on the number of true positives
#'
#' This function finds all the TP values of the models given (e.g. 0,1,2,3),
#' generates every pairwise combination (e.g. the group matchings: (0,1), (1,3),
#' etc.) and uses the \code{\link{get_avg_activity_diff_based_on_tp_predictions}}
#' function on each generated classification group matching, comparing thus all
#' groups of models with different true positive (TP) values.
#'
#' @param models character vector. The model names.
#' @param models.synergies.tp an integer vector of TP values. The \emph{names}
#' attribute holds the models' names and have to be in the same order as in the
#' \code{models} parameter.
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names (same order as in the \code{models}
#' parameter) whereas the column names specify the name of the network nodes
#' (gene, proteins, etc.). Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node).
#'
#' @return a matrix whose rows are \strong{vectors of
#' average node activity state differences} between two groups of models where
#' the classification was based on the number of true positive predictions.
#' Rows represent the diffferent classification group matchings, e.g. (1,2) means the
#' models that predicted 1 TP synergy vs the models that predicted 2 TP
#' synergies and the columns represent the network's node names.
#' Values are in the [-1,1] interval.
#'
#' @family average data difference functions
#'
#' @importFrom utils combn
#' @export
get_avg_activity_diff_mat_based_on_tp_predictions =
  function(models, models.synergies.tp, models.stable.state) {
    # check
    stopifnot(all(names(models.synergies.tp) == rownames(models.stable.state)))

    tp.values = sort(unique(models.synergies.tp))
    tp.values.comb = t(combn(tp.values, 2))

    diff.tp.mat = apply(tp.values.comb, 1, function(comb) {
      return(get_avg_activity_diff_based_on_tp_predictions(
        models, models.synergies.tp, models.stable.state,
        num.low = comb[1], num.high = comb[2]))
    })

    tp.comb.names = apply(tp.values.comb, 1, function(row) {
      return(paste0("(", paste(row, collapse = ","), ")"))
    })
    colnames(diff.tp.mat) = tp.comb.names
    diff.tp.mat = t(diff.tp.mat)

    return(diff.tp.mat)
  }

#' Get average link operator difference matrix based on the number of true positives
#'
#' This function uses the \code{\link{get_avg_activity_diff_mat_based_on_tp_predictions}}
#' function with the parameter \code{models.link.operator} as input in the place of
#' \code{models.stable.state}, since the two matrices representing the two inputs
#' have the same data format (rows represent models, columns represent nodes,
#' and each value is a number in the [0,1] interval).
#'
#' @param models character vector. The model names.
#' @param models.synergies.tp an integer vector of TP values. The \emph{names}
#' attribute holds the models' names and have to be in the same order as in the
#' \code{models} parameter.
#' @param models.link.operator matrix (nxm) with n models and m nodes. The row names of the matrix
#' specify the models' names whereas the column names specify the name of the
#' network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
#' element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
#' (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
#' both activating and inhibiting regulators (no link operator).
#'
#' @return a matrix whose rows are \strong{vectors of average node link operator
#' differences} between two groups of models based on some kind of classification
#' (e.g. number of TP predictions) and whose names are set in the \code{rownames}
#' attribute of the matrix (usually denoting the diffferent classification
#' groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
#' that predicted 2 TP synergies, if the classification is done by number of TP
#' predictions). The columns represent the network's node names. Values are in
#' the [-1,1] interval.
#'
#' @section Details:
#' So, if a node has a value close to -1 it means that on average,
#' this node's boolean equation has the \strong{AND NOT} link operator in the
#' 'good' models compared to the 'bad' ones while a value closer to 1 means that
#' the node's boolean equation has mostly the \strong{OR NOT} link operator
#' in the 'good' models. A value closer to 0 indicates that the link operator in
#' the node's boolean equation is \strong{not so much different} between the
#' 'good' and 'bad' models and so it won't not be a node of interest when
#' searching for indicators of better performance (higher number of true positives)
#' in the parameterization of the good models (the boolean equations). A value
#' exactly equal to 0 can also mean that this node didn't not have a link operator
#' in its boolean equation, again making it a non-important indicator of difference
#' in model performance.
#'
#' @family average data difference functions
#' @export
get_avg_link_operator_diff_mat_based_on_tp_predictions =
  function(models, models.synergies.tp, models.link.operator) {
    get_avg_activity_diff_mat_based_on_tp_predictions(
      models, models.synergies.tp, models.stable.state = models.link.operator
    )
  }

#' Get the average activity difference based on the number of true positives
#'
#' This function splits the models to 'good' and 'bad' based on the number of true
#' positive predictions: \emph{num.high} TPs (good) vs \emph{num.low} TPs (bad).
#' Then, for each network node, it finds the node's average activity in each of
#' the two classes (a value in the [0,1] interval) and then subtracts the
#' 'bad' average activity value from the good' one.
#'
#' @param models character vector. The model names.
#' @param models.synergies.tp an integer vector of TP values. The \emph{names}
#' attribute holds the models' names and have to be in the same order as in the
#' \code{models} parameter.
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names (same order as in the \code{models}
#' parameter) whereas the column names specify the name of the network nodes
#' (gene, proteins, etc.). Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node).
#' @param num.low integer. The number of true positives representing the 'bad'
#' model class.
#' @param num.high integer. The number of true positives representing the 'good'
#' model class. This number has to be strictly higher than \code{num.low}.
#'
#' @return a numeric vector with values in the [-1,1] interval (minimum and maximum
#' possible average difference) and with the \emph{names} attribute representing the name
#' of the nodes.
#'
#' @section Details:
#' So, if a node has a value close to -1 it means that on average,
#' this node is more \strong{inhibited} in the 'good' models compared to the
#' 'bad' ones while a value closer to 1 means that the node is more \strong{activated}
#' in the 'good' models. A value closer to 0 indicates that the activity of that
#' node is \strong{not so much different} between the 'good' and 'bad' models and
#' so it won't not be a node of interest when searching for indicators of better
#' performance (higher number of true positives) in the good models.
#'
#' @family average data difference functions
#' @export
get_avg_activity_diff_based_on_tp_predictions =
  function(models, models.synergies.tp, models.stable.state, num.low, num.high) {
    if (num.low >= num.high) {
      stop("`num.low` needs to be smaller than `num.high`")
    }

    good.models = models[models.synergies.tp == num.high]
    bad.models  = models[models.synergies.tp == num.low]

    # `good.models` != `bad.models` (disjoing sets of models)
    stopifnot(!(good.models %in% bad.models))

    # check: no models in some category
    stopifnot(length(good.models) > 0, length(bad.models) > 0)

    if (length(good.models) == 1) {
      warning( "Only 1 \'good\' model in TP class: ", num.high, " - very biased analysis\n")
      good.avg.activity = models.stable.state[good.models, ]
    } else {
      good.avg.activity = apply(models.stable.state[good.models, ], 2, mean)
    }

    if (length(bad.models) == 1) {
      warning(paste0("Only 1 \'bad\' model in TP class: ", num.low, " - very biased analysis\n"))
      bad.avg.activity = models.stable.state[bad.models, ]
    } else {
      bad.avg.activity = apply(models.stable.state[bad.models, ], 2, mean)
    }

    return(good.avg.activity - bad.avg.activity)
  }

#' Get average activity difference matrix based on MCC clustering
#'
#' This function splits the Matthews correlation coefficient (MCC) scores
#' of the models to specific groups using the \pkg{Ckmeans.1d.dp}
#' package for the clustering (groups are denoted by ids, e.g. NaN,1,2,3, etc.
#' where a larger id corresponds to a group of models with higher MCC scores)
#' and for each pairwise
#' combination of group id matchings (e.g. (0,1), (1,3), etc.), it uses the
#' \code{\link{get_avg_activity_diff_based_on_mcc_clustering}}
#' function, comparing thus all groups of models that belong to different
#' MCC classes.
#'
#' @param models.mcc a numeric vector of Matthews Correlation Coefficient (MCC)
#' scores, one for each model. The \emph{names} attribute holds the models' names.
#' Can be the result of using the function \code{\link{calculate_models_mcc}}.
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names (same order as in the \code{models.mcc}
#' parameter) whereas the column names specify the name of the network nodes
#' (gene, proteins, etc.). Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node).
#' @param num.of.mcc.classes numeric. A positive integer larger than 2 that
#' signifies the number of mcc classes (groups) that we should split the models
#' MCC values (excluding the 'NaN' values).
#' @param include.NaN.mcc.class logical. Should the models that have NaN MCC value
#' (e.g. TP+FP = 0, models that predicted no synergies at all) be classified together
#' in one class - the 'NaN MCC Class' - and compared with the other model classes
#' in the analysis? If \emph{TRUE}, then the number of total MCC classes will be
#' \emph{num.of.mcc.classes + 1}.
#'
#' @return a matrix whose rows are \strong{vectors of
#' average node activity state differences} between two groups of models where
#' the classification was based on the models' MCC values.
#' Rows represent the diffferent classification group matchings, e.g. (1,2) means the
#' models that belonged to the 1st group of MCC values vs the models that
#' belonged to the 2nd group. The columns represent the network's node names.
#' Values are in the [-1,1] interval.
#'
#' @family average data difference functions
#'
#' @importFrom Ckmeans.1d.dp Ckmeans.1d.dp
#' @export
get_avg_activity_diff_mat_based_on_mcc_clustering =
  function(models.mcc, models.stable.state, num.of.mcc.classes, include.NaN.mcc.class) {
    stopifnot(num.of.mcc.classes >= 2)
    mcc.class.ids = 1:num.of.mcc.classes

    models.mcc.no.nan = models.mcc[!is.nan(models.mcc)]
    models.mcc.no.nan.sorted = sort(models.mcc.no.nan)

    # find the clusters
    res = Ckmeans.1d.dp(x = models.mcc.no.nan.sorted, k = num.of.mcc.classes)
    models.cluster.ids = res$cluster

    # add NaN class (if applicable)
    if (include.NaN.mcc.class & (sum(is.nan(models.mcc)) > 0)) {
      mcc.class.ids = append(mcc.class.ids, values = NaN, after = 0)
    }

    mcc.class.id.comb = t(combn(1:length(mcc.class.ids), 2))

    diff.mcc.mat = apply(mcc.class.id.comb, 1, function(comb) {
      return(get_avg_activity_diff_based_on_mcc_clustering(
        models.mcc, models.stable.state,
        mcc.class.ids, models.cluster.ids,
        class.id.low = comb[1], class.id.high = comb[2]))
    })

    mcc.classes.comb.names = apply(mcc.class.id.comb, 1, function(comb) {
      return(paste0("(", mcc.class.ids[comb[1]], ",", mcc.class.ids[comb[2]], ")"))
    })

    colnames(diff.mcc.mat) = mcc.classes.comb.names
    diff.mcc.mat = t(diff.mcc.mat)

    return(diff.mcc.mat)
  }

#' Get average link operator difference matrix based on MCC clustering
#'
#' This function uses the \code{\link{get_avg_activity_diff_mat_based_on_mcc_clustering}}
#' function with the parameter \code{models.link.operator} as input in the place of
#' \code{models.stable.state}, since the two matrices representing the two inputs
#' have the same data format (rows represent models, columns represent nodes,
#' and each value is a number in the [0,1] interval).
#'
#' @param models.mcc a numeric vector of Matthews Correlation Coefficient (MCC)
#' scores, one for each model. The \emph{names} attribute holds the models' names.
#' Can be the result of using the function \code{\link{calculate_models_mcc}}.
#' @param models.link.operator matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names (same order as in the
#' \code{models.mcc} parameter) whereas the column names specify the name of the
#' network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
#' element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
#' (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
#' both activating and inhibiting regulators (no link operator).
#' @param num.of.mcc.classes numeric. A positive integer larger than 2 that
#' signifies the number of mcc classes (groups) that we should split the models
#' MCC values (excluding the 'NaN' values).
#' @param include.NaN.mcc.class logical. Should the models that have NaN MCC value
#' (e.g. TP+FP = 0, models that predicted no synergies at all) be classified together
#' in one class - the 'NaN MCC Class' - and compared with the other model classes
#' in the analysis? If \emph{TRUE}, then the number of total MCC classes will be
#' \emph{num.of.mcc.classes + 1}.
#'
#' @return a matrix whose rows are \strong{vectors of average node link operator
#' differences} between two groups of models where
#' the classification was based on the models' MCC values.
#' Rows represent the diffferent classification group matchings, e.g. (1,2) means the
#' models that belonged to the 1st group of MCC values vs the models that
#' belonged to the 2nd group. The columns represent the network's node names.
#' Values are in the [-1,1] interval.
#'
#' @section Details:
#' So, if a node has a value close to -1 it means that on average,
#' this node's boolean equation has the \strong{AND NOT} link operator in the
#' 'good' models compared to the 'bad' ones while a value closer to 1 means that
#' the node's boolean equation has mostly the \strong{OR NOT} link operator
#' in the 'good' models. A value closer to 0 indicates that the link operator in
#' the node's boolean equation is \strong{not so much different} between the
#' 'good' and 'bad' models and so it won't not be a node of interest when
#' searching for indicators of better performance (higher average MCC value)
#' in the parameterization of the good models (the boolean equations). A value
#' exactly equal to 0 can also mean that this node didn't not have a link operator
#' in its boolean equation, again making it a non-important indicator of difference
#' in model performance.
#'
#' @family average data difference functions
#' @export
get_avg_link_operator_diff_mat_based_on_mcc_clustering =
  function(models.mcc, models.link.operator, num.of.mcc.classes, include.NaN.mcc.class) {
    get_avg_activity_diff_mat_based_on_mcc_clustering(
      models.mcc, models.stable.state = models.link.operator,
      num.of.mcc.classes, include.NaN.mcc.class
    )
  }

#' Get the average activity difference based on MCC clustering
#'
#' This function splits the models to 'good' and 'bad' based on an MCC value
#' clustering method: \emph{class.id.high} denotes the group id with the higher MCC
#' values (good model group) vs \emph{class.id.low} which denotes the group id with
#' the lower MCC values (bad model group). Then, for each network node, the function
#' finds the node's average activity in each of the two classes (a value in
#' the [0,1] interval) and then subtracts the bad class average activity value from
#' the good one.
#'
#' @param models.mcc a numeric vector of Matthews Correlation Coefficient (MCC)
#' scores, one for each model. The \emph{names} attribute holds the models' names.
#' Can be the result of using the function \code{\link{calculate_models_mcc}}.
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names (same order as in the \code{models.mcc}
#' parameter) whereas the column names specify the name of the network nodes
#' (gene, proteins, etc.). Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node).
#' @param mcc.class.ids a numeric vector of group/class ids starting from NaN if
#' models with NaN MCC score are included or 1 otherwise. E.g. c(1,2,3), where
#' we have 3 MCC classes and no NaN values.
#' @param models.cluster.ids a numeric vector of cluster ids assigned to each
#' model. It is the result of using \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}
#' with input the sorted vector of the models' MCC values with no NaNs included.
#' @param class.id.low integer. This number specifies the MCC class id of the
#' 'bad' models.
#' @param class.id.high integer. This number specifies the MCC class id of the
#' 'good' models and needs to be strictly higher than \code{class.id.low}.
#'
#' @return a numeric vector with values in the [-1,1] interval (minimum and maximum
#' possible average difference) and with the \emph{names} attribute representing the name
#' of the nodes.
#'
#' @section Details:
#' So, if a node has a value close to -1 it means that on average,
#' this node is more \strong{inhibited} in the 'good' models compared to the
#' 'bad' ones while a value closer to 1 means that the node is more \strong{activated}
#' in the 'good' models. A value closer to 0 indicates that the activity of that
#' node is \strong{not so much different} between the 'good' and 'bad' models and
#' so it won't not be a node of interest when searching for indicators of better
#' performance (higher MCC score/class) in the good models.
#'
#' @family average data difference functions
#' @export
get_avg_activity_diff_based_on_mcc_clustering =
  function(models.mcc, models.stable.state, mcc.class.ids, models.cluster.ids,
           class.id.low, class.id.high) {
    if (class.id.low >= class.id.high) {
      stop("`class.id.low` needs to be smaller than `class.id.high`")
    }

    models.mcc.no.nan = models.mcc[!is.nan(models.mcc)]
    models.mcc.no.nan.sorted = sort(models.mcc.no.nan)

    bad.class.id  = mcc.class.ids[class.id.low]
    good.class.id = mcc.class.ids[class.id.high]

    # find the 'good' models
    good.models = get_models_based_on_mcc_class_id(
      good.class.id, models.cluster.ids, models.mcc.no.nan.sorted
    )

    # find the 'bad' models
    if (is.nan(bad.class.id)) {
      # the `NaN` MCC scored models (can only be 'bad' ones)
      bad.models = names(models.mcc)[is.nan(models.mcc)]
    } else {
      bad.models =
        get_models_based_on_mcc_class_id(
          bad.class.id, models.cluster.ids, models.mcc.no.nan.sorted
        )
    }

    # `good.models` != `bad.models` (disjoing sets of models)
    stopifnot(!(good.models %in% bad.models))

    # check: no models in some category
    stopifnot(length(good.models) > 0, length(bad.models) > 0)

    if (length(good.models) == 1) {
      warning(paste0("Only 1 \'good\' model in MCC class: ", class.id.high, " - very biased analysis\n"))
      good.avg.activity = models.stable.state[good.models, ]
    } else {
      good.avg.activity = apply(models.stable.state[good.models, ], 2, mean)
    }

    if (length(bad.models) == 1) {
      warning(paste0("Only 1 \'bad\' model in MCC class: ", class.id.low, " - very biased analysis\n"))
      bad.avg.activity = models.stable.state[bad.models, ]
    } else {
      bad.avg.activity = apply(models.stable.state[bad.models, ], 2, mean)
    }

    return(good.avg.activity - bad.avg.activity)
  }

#' Get models based on the MCC class id
#'
#' This helper function finds all the models that belong to a specific MCC
#' cluster, i.e. their MCC values belong to the same cluster id.
#'
#' @param class.id an integer specifying the class id.
#' @param models.cluster.ids a numeric vector of cluster ids assigned to each
#' model. It is the result of using \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}
#' with input the sorted vector of the models' MCC values with no NaNs included.
#' @param models.mcc a numeric sorted vector of Matthews Correlation Coefficient (MCC)
#' scores, one for each model (no NaNs included).
#' The \emph{names} attribute holds the models' names.
#'
#' @return a character vector of model names
get_models_based_on_mcc_class_id =
  function(class.id, models.cluster.ids, models.mcc) {
    return(names(models.mcc[models.cluster.ids == class.id]))
  }

#' Get average activity difference matrix based on specific synergy prediction
#'
#' This function uses the \code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}}
#' function on a vector of drug combinations that were observed as synergistic
#' (e.g. by experiments) but also found as such by at least one of the models
#' (these drug combinations are the \emph{predicted synergies}).
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names
#' whereas the column names specify the name of the network nodes
#' (gene, proteins, etc.). Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node).
#' @param predicted.synergies a character vector of the synergies (drug
#' combination names) that were predicted by \strong{at least one} of the models
#' in the dataset. It must be a subset of the column names (the drug combinations)
#' of the \code{model.predictions} object.
#'
#' @return a matrix whose rows are \strong{vectors of
#' average node activity state differences} between two groups of models where
#' the classification for each individual row was based on the prediction or not
#' of a specific synergistic drug combination.
#' The row names are the predicted synergies, one per row, while the columns
#' represent the network's node names. Values are in the [-1,1] interval.
#'
#' @family average data difference functions
#' @export
get_avg_activity_diff_mat_based_on_specific_synergy_prediction =
  function(model.predictions, models.stable.state, predicted.synergies) {
    stopifnot(all(predicted.synergies %in% colnames(model.predictions)))

    diff.synergies.mat =
      sapply(predicted.synergies, function(drug.comb) {
        get_avg_activity_diff_based_on_specific_synergy_prediction(
          model.predictions, models.stable.state, drug.comb)
      })
    diff.synergies.mat = t(diff.synergies.mat)

    return(diff.synergies.mat)
  }

#' Get average link operator difference matrix based on specific synergy prediction
#'
#' This function uses the \code{\link{get_avg_activity_diff_mat_based_on_specific_synergy_prediction}}
#' function with the parameter \code{models.link.operator} as input in the place of
#' \code{models.stable.state}, since the two matrices representing the two inputs
#' have the same data format (rows represent models, columns represent nodes,
#' and each value is a number in the [0,1] interval).
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param models.link.operator matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names (same order as in the
#' \code{models.mcc} parameter) whereas the column names specify the name of the
#' network nodes (gene, proteins, etc.). Possible values for each \emph{model-node
#' element} are either \emph{0} (\strong{AND NOT} link operator), \emph{1}
#' (\strong{OR NOT} link operator) or \emph{0.5} if the node is not targeted by
#' both activating and inhibiting regulators (no link operator).
#' @param predicted.synergies a character vector of the synergies (drug
#' combination names) that were predicted by \strong{at least one} of the models
#' in the dataset. It must be a subset of the column names (the drug combinations)
#' of the \code{model.predictions} object.
#'
#' @return a matrix whose rows are \strong{vectors of average node link operator
#' differences} between two groups of models where the classification for each
#' individual row was based on the prediction or not of a specific synergistic
#' drug combination.
#' The row names are the predicted synergies, one per row, while the columns
#' represent the network's node names. Values are in the [-1,1] interval.
#' the classification was based on the models' MCC values.
#'
#' @section Details:
#' So, if a node has a value close to -1 it means that on average,
#' this node's boolean equation has the \strong{AND NOT} link operator in the
#' models that predicted the specified synergy while a value closer to 1 means that
#' the node's boolean equation has mostly the \strong{OR NOT} link operator
#' in these models. A value closer to 0 indicates that the link operator in
#' the node's boolean equation is \strong{not so much different} between the
#' models that predicted the synergy and those that did not and so it won't not
#' be a node of interest when searching for \emph{synergy biomarkers} - nodes
#' whose parameterization (value of the link operator) affects the manifestation
#' of synergy. A value exactly equal to 0 can also mean that this node didn't
#' not have a link operator in its boolean equation (making it thus a non-important
#' node with regard to the parameterization).
#'
#' @family average data difference functions
#' @export
get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction =
  function(model.predictions, models.link.operator, predicted.synergies) {
    get_avg_activity_diff_mat_based_on_specific_synergy_prediction(
      model.predictions, models.stable.state = models.link.operator,
      predicted.synergies)
  }

#' Get average activity difference based on specific synergy prediction
#'
#' Given a specific drug combination, this function splits the models to
#' good (those that predicted that particular combination, i.e. found it
#' as synergistic - a value of \emph{1} in the \code{model.predictions}) and
#' bad (those that found it as non-synergistic - a value of \emph{0} in the
#' \code{model.predictions}). The models whose predicted value for that synergy is marked as
#' \emph{NA} are excluded from the analysis. Then, for each network node, the
#' function finds the node's average activity in each of the two model groups (a
#' value in the [0,1] interval) and then subtracts the bad group's average
#' activity value from the good one.
#'
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names
#' whereas the column names specify the name of the network nodes
#' (gene, proteins, etc.). Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node).
#' @param drug.comb string. The drug combination which will be used to split
#' the models. It must be included in the column names of the \code{model.predictions}
#' object.
#'
#' @return a numeric vector with values in the [-1,1] interval (minimum and maximum
#' possible average difference) and with the \emph{names} attribute representing
#' the name of the nodes.
#'
#' @section Details:
#' So, if a node has a value close to -1 it means that on average,
#' this node is more \strong{inhibited} in the models that predicted the specific
#' drug combination given, whereas a value closer to 1 means that the node is more
#' \strong{activated} in these models. A value closer to 0 indicates that the activity of that
#' node is \strong{not so much different} between the models that predicted the synergy and
#' those that did not and so it won't not be a node of interest when searching
#' for \emph{synergy biomarkers} - nodes whose activity is important for the
#' manifestation of the synergy.
#'
#' @family average data difference functions
#'
#' @importFrom usefun is_empty
#' @export
get_avg_activity_diff_based_on_specific_synergy_prediction =
  function(model.predictions, models.stable.state, drug.comb) {
    stopifnot(drug.comb %in% colnames(model.predictions))

    good.models = rownames(model.predictions)[
      model.predictions[, drug.comb] == 1 & !is.na(model.predictions[, drug.comb])
      ]
    bad.models  = rownames(model.predictions)[
      model.predictions[, drug.comb] == 0 & !is.na(model.predictions[, drug.comb])
      ]
    # na.models = rownames(model.predictions)[is.na(model.predictions[, drug.comb])]

    # check: no empty list of either good or bad models
    stopifnot(!is_empty(bad.models))
    stopifnot(!is_empty(good.models))

    if (length(good.models) == 1) {
      good.avg.activity = models.stable.state[good.models, ]
    } else {
      good.avg.activity = apply(models.stable.state[good.models, ], 2, mean)
    }

    if (length(bad.models) == 1) {
      bad.avg.activity = models.stable.state[bad.models, ]
    } else {
      bad.avg.activity = apply(models.stable.state[bad.models, ], 2, mean)
    }

    return(good.avg.activity - bad.avg.activity)
  }

#' Get the average activity difference based on the comparison of two synergy sets
#'
#' This function splits the models to 'good' and 'bad' based on the predictions
#' of two different synergy sets, one of them being a subset of the other.
#' The 'good' models are those that predict the \code{synergy.set.str}
#' (e.g. "A-B,A-C,B-C") while the 'bad' models are those that predict the
#' \code{synergy.subset.str} (e.g. "A-B,B-C"). Then, for each network node,
#' the function finds the node's average activity in each of the two classes
#' (a value in the [0,1] interval) and then subtracts the bad class average
#' activity value from the good one.
#'
#' @param synergy.set.str a string of drug combinations, comma-seperated. The
#' number of the specified combinations must be larger than the ones defined
#' in the \code{synergy.subset.str} parameter. They also must be included in the
#' tested drug combinations, i.e. the columns of the \code{model.predictions}
#' parameter.
#' @param synergy.subset.str a string of drug combinations, comma-seperated.
#' There must be at least one combination defined and all of them should also
#' be included in the \code{synergy.set.str} parameter.
#' @param model.predictions a \code{data.frame} object with rows the models and
#' columns the drug combinations. Possible values for each \emph{model-drug combination
#' element} are either \emph{0} (no synergy predicted), \emph{1} (synergy was
#' predicted) or \emph{NA} (couldn't find stable states in either the drug
#' combination inhibited model or in any of the two single-drug inhibited models)
#' @param models.stable.state a matrix (nxm) with n models and m nodes. The row
#' names of the matrix specify the models' names
#' whereas the column names specify the name of the network nodes
#' (gene, proteins, etc.). Possible values for each \emph{model-node element}
#' are either \emph{0} (inactive node) or \emph{1} (active node).
#'
#' @return a numeric vector with values in the [-1,1] interval (minimum and
#' maximum possible average difference) and with the names attribute
#' representing the name of the nodes.
#'
#' @section Details:
#' So, if a node has a value close to -1 it means that on average,
#' this node is more \strong{inhibited} in the models that predicted the extra
#' synergy(-ies) that are included in the \code{synergy.set.str} but not in the
#' \code{synergy.subset.str}, whereas a value closer to 1 means that the node is
#' more \strong{activated} in these models. These nodes are \strong{potential
#' biomarkers} because their activity state can influence the prediction
#' performance of a model and make it predict the extra synergy(-ies).
#' A value closer to 0 indicates that the activity of that
#' node is \strong{not so much different} between the models that predicted the
#' synergy set and those that predicted it's subset, so it won't not be a node
#' of interest when searching for potential biomarkers for the extra synergy(-ies).
#'
#' @family average data difference functions
#'
#' @importFrom usefun outersect is_empty
#' @export
get_avg_activity_diff_based_on_synergy_set_cmp =
  function(synergy.set.str, synergy.subset.str, model.predictions,
           models.stable.state) {

    synergy.set = unlist(strsplit(synergy.set.str, split = ","))
    synergy.subset = unlist(strsplit(synergy.subset.str, split = ","))

    # some checks
    stopifnot(length(synergy.subset) > 0,
              length(synergy.set) > length(synergy.subset))
    stopifnot(all(synergy.subset %in% synergy.set))
    stopifnot(all(synergy.set %in% colnames(model.predictions)))

    # find models that predict the `synergy.set`
    if (length(synergy.set) == 1) {
      models.synergy.set = rownames(model.predictions)[
        model.predictions[, synergy.set] == 1 &
          !is.na(model.predictions[, synergy.set])]
    } else {
      models.synergy.set = rownames(model.predictions)[
        apply(model.predictions[, synergy.set], 1,
              function(x) all(x == 1 & !is.na(x)))]
    }

    # find models that predict the `synergy.subset`
    if (length(synergy.subset) == 1) {
      models.synergy.subset = rownames(model.predictions)[
        model.predictions[, synergy.subset] == 1 &
          !is.na(model.predictions[, synergy.subset])]
    } else {
      models.synergy.subset = rownames(model.predictions)[
        apply(model.predictions[, synergy.subset], 1,
              function(x) all(x == 1 & !is.na(x)))]
    }

    common.models = intersect(models.synergy.set, models.synergy.subset)
    good.models = common.models
    bad.models  = outersect(models.synergy.set, models.synergy.subset)

    # check: no good model inside the bad model list
    stopifnot(all(!(good.models %in% bad.models)))

    # check: no empty list of either good or bad models
    stopifnot(!is_empty(bad.models))
    stopifnot(!is_empty(good.models))

    if (length(good.models) == 1) {
      good.avg.activity = models.stable.state[good.models, ]
    } else {
      good.avg.activity = apply(models.stable.state[good.models, ], 2, mean)
    }

    if (length(bad.models) == 1) {
      bad.avg.activity = models.stable.state[bad.models, ]
    } else {
      bad.avg.activity = apply(models.stable.state[bad.models, ], 2, mean)
    }

    return(good.avg.activity - bad.avg.activity)
  }
