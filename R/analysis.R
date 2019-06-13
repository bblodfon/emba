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

#' Get the average activity difference based on the number of true positives
#'
#' This function splits the models to good and bad based on the number of true
#' positive predictions: \emph{num.high} TPs (good) vs \emph{num.low} TPs (bad).
#' Then, for each network node, finds the average activity in each of the two
#' classes and then subtracts the 'bad' average activity value from the good
#'
#' @param models s
#' @param models.synergies.tp a
#' @param models.stable.state ert
#' @param num.low integer. The number of true positives representing the 'bad'
#' model class
#' @param num.high integer. The number of true positives representing the 'good'
#' model class. This number has to be strictly higher than \code{num.low}
#'
#' @return
#'
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

# `class.id.low` is the `mcc.interval` id for the 'bad' models
# `class.id.high` is the `mcc.interval` id for the 'good' models
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
#'
#' @importFrom usefun outersect is_empty
get_avg_activity_diff_based_on_specific_synergy_prediction =
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

calculate_models_mcc =
  function(observed.model.predictions, unobserved.model.predictions,
           models.synergies.tp, number.of.drug.comb.tested) {
    # Count the false negatives (FN)
    models.synergies.fn = apply(observed.model.predictions, 1, function(x) {
      sum( x == 0 | is.na(x) )
    })

    # P = TP + FN (Positives)
    positives = ncol(observed.model.predictions)
    models.synergies.p = models.synergies.tp + models.synergies.fn

    # Count the predictions of the non-observed synergies per model (FP)
    models.synergies.fp =
      calculate_models_synergies_fp(unobserved.model.predictions)

    # Count the True Negatives (TN)
    models.synergies.tn = apply(unobserved.model.predictions, 1, function(x) {
      sum( x == 0 | is.na(x))
    })

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

calculate_models_synergies_fp = function(unobserved.model.predictions) {
  # Count the predictions of the non-observed synergies per model (FP)
  models.synergies.fp = apply(unobserved.model.predictions, 1, sum, na.rm = T)
  return(models.synergies.fp)
}

# inputs are vectors of same size and one-to-one value correspondence
calculate_mcc = function(tp, tn, fp, fn, p, n) {
  return(
    (tp * tn - fp * fn) / sqrt((tp + fp) * p * n * (tn + fn))
  )
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

# `diff.res` is a 2-dim matrix (rows = classification group matchings,
# columns = nodes)
# `type` = active or inhibited
# If there is at least one value in a column that surpasses the threshold given,
# the corresponding node is return as a biomarker
get_biomarkers = function(diff.res, threshold, type) {
  dimen = dim(diff.res)
  rows = dimen[1]
  nodes.num = dimen[2]

  biomarkers = character(0)
  for(node.index in 1:nodes.num) {
    node.name = colnames(diff.res)[node.index]
    for (row.index in 1:rows) {
      if (type == "active") {
        if (diff.res[row.index, node.index] > threshold) {
          biomarkers = c(biomarkers, node.name)
          break
        }
      } else { # inhibited
        if (diff.res[row.index, node.index] < -threshold) {
          biomarkers = c(biomarkers, node.name)
          break
        }
      }
    }
  }

  return(biomarkers)
}

# adds one more row to the `biomarkers.synergy.res` data.frame with the
# performance-related biomarkers
add_performance_biomarkers =
  function(biomarkers.synergy.res, biomarkers.active, biomarkers.inhibited) {
    # initialize `row` data.frame
    node.names = colnames(biomarkers.synergy.res)
    row = as.data.frame(matrix(0, ncol = length(node.names), nrow = 1))
    colnames(row) = node.names
    rownames(row) = "PERF"

    # add biomarkers
    row[colnames(row) %in% biomarkers.active] = 1
    row[colnames(row) %in% biomarkers.inhibited] = -1

    res = rbind(row, biomarkers.synergy.res)
    return(res)
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

#' helper function to check which synergy sets to compare (the small synergy set
#' misses just one synergy from the larger set)
#'
#' @importFrom usefun outersect
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

