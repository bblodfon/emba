#' Get average activity difference matrix based on the number of true positives
#'
#' This function uses the \code{\link[emba]{get_avg_activity_diff_based_on_tp_predictions}}
#' function on each row vector of the given \code{models.stable.state} matrix.
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
#' @return a matrix whose rows are \strong{vectors of average node activity
#' state differences} between two groups of models based on some kind of classification
#' (e.g. number of TP predictions) and whose names are set in the \code{rownames}
#' attribute of the matrix (usually denoting the diffferent classification
#' groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
#' that predicted 2 TP synergies, if the classification is done by number of TP
#' predictions). The columns represent the network's node names. Values are in
#' the [-1,1] interval.
#'
#' @family average data difference functions
#'
#' @importFrom utils combn
#' @export
get_avg_activity_diff_mat_based_on_tp_predictions =
  function(models, models.synergies.tp, models.stable.state) {
    # check
    stopifnot(all(names(models.synergies.tp) == rownames(models.stable.state)))

    # TODO check if not all values are returned, does it affect results???
    tp.values = unique(models.synergies.tp)
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
#' This function uses the \code{\link[emba]{get_avg_activity_diff_mat_based_on_tp_predictions}}
#' function with the parameter \code{model.equations} as input in the place of
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
      models, models.synergies.tp, models.link.operator
    )
  }

#' Get the average activity difference based on the number of true positives
#'
#' This function splits the models to good and bad based on the number of true
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
#' average difference) and with the \emph{names} attribute representing the name
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

    # small number of models in some category: TP analysis not good :)
    stopifnot(length(good.models) > 1)
    stopifnot(length(bad.models) > 1)

    good.avg.activity = apply(models.stable.state[good.models, ], 2, mean)
    bad.avg.activity = apply(models.stable.state[bad.models, ], 2, mean)

    return(good.avg.activity - bad.avg.activity)
  }

#'`class.id.low` is the `mcc.interval` id for the 'bad' models
#' `class.id.high` is the `mcc.interval` id for the 'good' models
#'
#' @family average activity difference functions
get_avg_activity_diff_based_on_mcc_classification =
  function(models, models.mcc, mcc.intervals, models.stable.state,
           class.id.low, class.id.high) {
    if (class.id.low >= class.id.high) {
      stop("`class.id.low` needs to be smaller than `class.id.high`")
    }

    mcc.interval.low = mcc.intervals[class.id.low, ]
    mcc.interval.high = mcc.intervals[class.id.high, ]

    # find the 'good' models
    max.value = max(mcc.intervals, na.rm = TRUE)
    if (mcc.interval.high[2] == max.value) {
      good.models =
        get_models_based_on_mcc_interval(models, models.mcc, mcc.interval.high,
                                         include.high.value = TRUE)
    } else {
      good.models =
        get_models_based_on_mcc_interval(models, models.mcc, mcc.interval.high)
    }

    # find the 'bad' models
    if (is.na(mcc.interval.low[1])) {
      # the `NaN` MCC scored models (can only be 'bad' ones)
      bad.models = models[is.na(models.mcc)]
    } else {
      bad.models =
        get_models_based_on_mcc_interval(models, models.mcc, mcc.interval.low)
    }

    # `good.models` != `bad.models` (disjoing sets of models)
    stopifnot(!(good.models %in% bad.models))
    # small number of models in some category: need to redefine the MCC intervals
    stopifnot(length(good.models) > 1)
    stopifnot(length(bad.models) > 1)

    good.avg.activity = apply(models.stable.state[good.models, ], 2, mean)
    bad.avg.activity = apply(models.stable.state[bad.models, ], 2, mean)

    return(good.avg.activity - bad.avg.activity)
  }

#' Title mcc
#'
#' @family average activity difference functions
get_avg_activity_diff_based_on_mcc_clustering =
  function(models.mcc, models.mcc.no.nan.sorted, models.stable.state,
           mcc.class.ids, models.cluster.ids, class.id.low, class.id.high) {
    if (class.id.low >= class.id.high) {
      stop("`class.id.low` needs to be smaller than `class.id.high`")
    }

    bad.class.id  = mcc.class.ids[class.id.low]
    good.class.id = mcc.class.ids[class.id.high]

    # find the 'good' models
    good.models = get_models_based_on_mcc_class_id(
      good.class.id, models.cluster.ids, models.mcc.no.nan.sorted
    )

    # find the 'bad' models
    if(is.nan(bad.class.id)) {
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
    # small number of models in some category: need to redefine the MCC intervals
    stopifnot(length(good.models) > 1)
    stopifnot(length(bad.models) > 1)

    good.avg.activity = apply(models.stable.state[good.models, ], 2, mean)
    bad.avg.activity = apply(models.stable.state[bad.models, ], 2, mean)

    return(good.avg.activity - bad.avg.activity)
  }

# `models.cluster.ids` is a vector specifying the class id of the MCC score
# as defined in the `models.mcc` (one-to-one)
get_models_based_on_mcc_class_id =
  function(class.id, models.cluster.ids, models.mcc) {
    return(names(models.mcc[models.cluster.ids == class.id]))
  }

#' @importFrom usefun is_between
get_models_based_on_mcc_interval =
  function(models, models.mcc, mcc.interval, include.high.value = FALSE) {
    res = sapply(models.mcc, is_between, low.thres = mcc.interval[1],
                 high.thres = mcc.interval[2], include.high.value)
    # exclude the NA values
    res.pruned = res[!is.na(res)]
    models.pruned = models[!is.na(res)]

    return(models.pruned[res.pruned])
  }

#' Example use: `drug.comb` = "AK-PD"
#'
#' @family average activity difference functions
#'
#' @importFrom usefun is_empty
get_avg_activity_diff_based_on_specific_synergy_prediction =
  function(model.predictions, models.stable.state, drug.comb) {
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

#' To get meaningful results, one set must be a subset of the other
#' Example use:
#' synergy.set.str = "A-B,A-D,B-D,P-S"
#' synergy.subset.str = "A-B,B-D,P-S"

#' @family average activity difference functions
#'
#' @importFrom usefun outersect is_empty
get_avg_activity_diff_based_on_diff_synergy_prediction =
  function(synergy.set.str, synergy.subset.str, model.predictions,
           models.stable.state) {

    synergy.set = unlist(strsplit(synergy.set.str, split = ","))
    synergy.subset = unlist(strsplit(synergy.subset.str, split = ","))

    # some checks
    stopifnot(length(synergy.subset) > 0,
              length(synergy.set) > length(synergy.subset))
    stopifnot(all(synergy.subset %in% synergy.set))

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

get_mcc_intervals = function(mcc.values, interval.size) {
  min.mcc = min(mcc.values, na.rm = TRUE)
  max.mcc = max(mcc.values, na.rm = TRUE)
  mcc.points = seq(-1.0, 1.0, interval.size)
  mcc.points.pruned = mcc.points[min.mcc < (mcc.points + interval.size) &
                                   mcc.points < (max.mcc + interval.size)]

  mcc.intervals =
    matrix(numeric(), nrow = length(mcc.points.pruned) - 1, ncol = 2)
  for (i in 1:nrow(mcc.intervals)) {
    mcc.intervals[i, 1] = mcc.points.pruned[i]
    mcc.intervals[i, 2] = mcc.points.pruned[i + 1]
  }

  return(mcc.intervals)
}

get_mcc_classes = function(mcc.intervals) {
  number.of.intervals = nrow(mcc.intervals)
  mcc.classes = character(number.of.intervals)

  for (i in 1:number.of.intervals) {
    mcc.interval = mcc.intervals[i,]
    low.value = mcc.interval[1]
    high.value = mcc.interval[2]
    if (is.na(low.value)) mcc.classes[i] = "NaN"
    else if (i != number.of.intervals)
      mcc.classes[i] = paste0("[", low.value, ", " , high.value, ")")
    else
      mcc.classes[i] = paste0("[", low.value, ", " , high.value, "]")
  }

  return(mcc.classes)
}
