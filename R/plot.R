#' Plot histogram of the MCC classes
#'
#' This function is a wrapper of the \code{\link[Ckmeans.1d.dp]{ahist}} function
#' for plotting nicely the distribution of the MCC models' values.
#'
#' @param models.mcc.no.nan.sorted a numeric sorted vector of Matthews
#' Correlation Coefficient (MCC) scores, one for each model (no NaNs included).
#' The \emph{names} attribute holds the models' names.
#' @param models.cluster.ids a numeric vector of cluster ids assigned to each
#' model. It is the result of using \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}
#' with input the sorted vector of the models' MCC values with no NaNs included
#' (\code{models.mcc.no.nan.sorted}).
#' @param num.of.mcc.classes numeric. A positive integer (>2) that signifies the
#' number of mcc classes (groups) that we should split the models MCC values
#' (excluding the 'NaN' values).
#' @param mcc.class.ids a numeric vector ranging from from 1 to
#' \code{num.of.mcc.classes}.
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' @importFrom Ckmeans.1d.dp ahist
#' @export
plot_mcc_classes_hist = function(models.mcc.no.nan.sorted, models.cluster.ids,
                                 num.of.mcc.classes, mcc.class.ids) {
  min.x.value = round(min(models.mcc.no.nan.sorted) - 0.1, digits = 1)
  max.x.value = round(max(models.mcc.no.nan.sorted) + 0.1, digits = 1)
  rainbow.colors = rainbow(num.of.mcc.classes)

  ahist(models.mcc.no.nan.sorted, k = num.of.mcc.classes,
        main = "Model MCC-Classification", xlab = "MCC value",
        sub = paste("n =", length(models.mcc.no.nan.sorted), "models, k =",
                    num.of.mcc.classes, "classes"),
        col = rainbow.colors, col.stick = rainbow.colors[models.cluster.ids],
        xlim = c(min.x.value, max.x.value))

  legend("topright", legend = mcc.class.ids, title = "MCC Classes",
         col = rainbow.colors, lty = 1, lwd = 10)
}

#' Bar plot of model stats
#'
#' Use this function to produce a bar plot when the input is the result of using
#' the \link[base]{table} function to a numeric vector
#'
#' @param models.stats table object, the result of using \link[base]{table} on
#' a (numeric) vector. Usually it represents some models statistics summary -
#' counts for each TP prediction value for example.
#' @param cell.line string. The name of the cell line to be used in the title
#' of the produced plot. Default value: NULL (the cell line name will not be
#' added to the title)
#' @param title string. The title of the plot
#' @param xlab string. The title of the x-axis
#' @param ylab string. The title of the y-axis
#' @param cont.values logical. If TRUE, the values of the x-axis will be trimmed
#' to 3 digits after the decimal point. Default value: FALSE.
#' @param threshold integer. Values from the \code{model.stats} that are \emph{less
#' or equal} to the threshold will be pruned. Use it when there too many
#' categories and the graph appears too dense. Default value: 0
#'
#' @importFrom graphics barplot axis par
#' @export
make_barplot_on_models_stats =
  function(models.stats, cell.line, title, xlab, ylab, cont.values = FALSE, threshold = 0) {
    if (!is.null(cell.line))
      cell.line.text = paste0(" (", cell.line, ")")
    else
      cell.line.text = ""

    # Find is there is just one `NaN` category
    there.is.one.NaN.category = FALSE
    nan.index = which(is.nan(as.numeric(names(models.stats))))
    if (length(nan.index) == 1) {
      there.is.one.NaN.category = TRUE
      nan.value = models.stats[nan.index]
    }

    # If there is just one `NaN` category, put it first in the `models.stats`
    # for presentation purposes in the barplot
    if (there.is.one.NaN.category) {
      models.stats = c(nan.value, models.stats[names(models.stats) != "NaN"])
    }

    # prune some bars :)
    models.stats = models.stats[models.stats > threshold]

    # If number of `NaN` values are lower then the `threshold` and
    # as such will be pruned, there will be no `NaN` bar in the plot
    if (there.is.one.NaN.category && nan.value <= threshold)
      there.is.one.NaN.category = FALSE

    x.axis.values =
      get_x_axis_values(models.stats, there.is.one.NaN.category, cont.values)
    y.axis.values = pretty(models.stats)

    bp = barplot(models.stats, col = "cornflowerblue",
                 names.arg = x.axis.values, yaxt = "n",
                 ylim = c(0, max(y.axis.values) + 500),
                 main = paste0(title, cell.line.text),
                 xlab = xlab, ylab = ylab)
    axis(2, at = y.axis.values, las = 1)

    add_numbers_above_the_bars(models.stats, bp, color = "red")

    # If there is just one `NaN` category, label it in the plot
    if (there.is.one.NaN.category) {
      text(x = bp[1], y = nan.value/2, labels = names(nan.value),
           col = "yellow", srt = 90, font = 2)
    }
}

#' Get the refined x-axis values
#'
#' This function returns the x-axis values that are going to be used by
#' \link[emba]{make_barplot_on_models_stats} to render the bar plot.
#'
#' @param models.stats table object, the result of using \link[base]{table} on
#' a (numeric) vector. Usually it represents some models statistics summary -
#' counts for each TP prediction value for example.
#' @param there.is.one.NaN.category logical. Is there one \emph{NaN} category?
#' (check is done before on the \emph{names} attribute of the \code{models.stats})
#' @param cont.values logical. If TRUE, the values of the x-axis will be trimmed
#' to 3 digits after the decimal point. Otherwise, they will be returned as they
#' are.
get_x_axis_values =
  function(models.stats, there.is.one.NaN.category, cont.values) {
    if (there.is.one.NaN.category) {
      # replace `NaN` value with empty space at the beginning of the x axis
      x.values = c(" ", names(models.stats)[names(models.stats) != "NaN"])
    } else x.values = names(models.stats)

    if (cont.values) {
      return(round(as.numeric(x.values), digits = 3))
    } else {
      return(x.values)
    }
}

#' Bar plot of observed synergy subsets
#'
#' Use this function to easily make a barplot that shows the amount of models
#' that predicted each synergy subset out of the set of all observed synergies.
#'
#' @param synergy.subset.stats integer vector with values the amount of models
#' that predicted each synergy subset, defined as a comma-separated string of
#' drug combinations in the \emph{names} attribute of the vector
#' @param threshold.for.subset.removal integer. Use it to discard elements of
#' the \code{synergy.subset.stats} vector that are strictly less than the
#' specified threshold
#' @param bottom.margin integer used to vertically fit in the names of the drug
#' combinations in the x-axis (specified in inches). The best \code{bottom.margin}
#' value depends on the \emph{maximum size} of a synergy subset as defined in the
#' \code{names} attribute of the \code{synergy.subset.stats}.
#' Some rules of thumb are:
#' size = 1 => bottom.margin = 4,
#' size = 2 => bottom.margin = 6,
#' size = 3 => bottom.margin = 9,
#' size = 4 => bottom.margin = 12, etc.
#' @param cell.line string. The name of the cell line to be used in the title
#' of the produced plot. Default value: NULL (the cell line name will not be
#' added to the title).
#'
#' @export
make_barplot_on_synergy_subset_stats =
  function(synergy.subset.stats, threshold.for.subset.removal, bottom.margin,
           cell.line = NULL) {
    if (!is.null(cell.line))
      cell.line.text = paste0(" (", cell.line, ")")
    else
      cell.line.text = ""

    synergy.subset.stats = synergy.subset.stats[
      !synergy.subset.stats < threshold.for.subset.removal
    ]

    par(mar = c(bottom.margin, 4, 4, 2)) # c(bottom, left, top, right)

    y.axis.values = pretty(synergy.subset.stats)
    bp = barplot(synergy.subset.stats, col = "green", space = 0.5, las = 2,
                 main = paste0("Model Synergy Predictions per Observed Synergy",
                                " Subset", cell.line.text),
                 ylab = "Number of models", yaxt = "n",
                 ylim = c(0, max(y.axis.values) + 500))
    axis(2, at = y.axis.values, las = 1)

    add_numbers_above_the_bars(synergy.subset.stats, bp, color = "red")
}

#' Add numbers horizontally above the bars of a barplot
#'
#' @param stats a numeric vector
#' @param bp the result of \strong{\code{barplot}} command, usually a numeric
#' vector or matrix
#' @param color string. The color for the numbers
#'
#' @importFrom graphics text
add_numbers_above_the_bars = function(stats, bp, color) {
  for (i in 1:length(stats)) {
    text(x = bp[i], y = stats[i], labels = stats[i],
         col = color, pos = 3)
  }
}

#' Plot the graph of average state differences (visNetwork)
#'
#' This function uses the \code{\link[visNetwork]{visNetwork}} package to plot a
#' network of nodes. The nodes are positioned according to the specified coordinates
#' given by the \code{layout} parameter and the colors are derived using the
#' \code{diff} values and the \code{\link{get_node_colors}} function. The color
#' of each node indicates how much more inhibited or active that node is, when
#' comparing the average model classified in the 'good' category vs the average
#' 'bad' one.
#'
#' @param net igraph graph object (to be translated to a \emph{visNetwork} object)
#' @param diff numeric vector. Every value is in the [-1,1] interval and
#' represents the average activity difference of each node. The node names have
#' to be specified in the \emph{names} attribute of the given vector. For example,
#' \code{diff} could be the result of using the function
#' \code{\link{get_avg_activity_diff_based_on_specific_synergy_prediction}}.
#' @param layout a (nx2) numeric matrix of x-y coordinates (2 columns) for each
#' of the nodes (n) in the \code{net} igraph object
#' @param title string. The title of the visNetwork plot
#'
#' @family network plotting functions
#'
#' @importFrom magrittr %>%
#' @importFrom visNetwork toVisNetworkData visNetwork visLegend
#' @export
plot_avg_state_diff_graph_vis = function(net, diff, layout, title) {
  data = toVisNetworkData(net)
  nodes = data$nodes
  edges = data$edges

  # colors for nodes (to be interpolated) matching one-to-one the diff values
  col = c("tomato", "grey", "gold")
  nodes$color = get_node_colors(net, diff, col)

  # set visualization graph attributes
  nodes$size = 30
  nodes$physics = FALSE
  nodes$shape = "dot"
  scale.factor = 60
  nodes$x = layout[,2] * scale.factor
  nodes$y = layout[,1] * scale.factor

  edges$smooth = FALSE
  edges$physics = FALSE
  edges$arrows = "to"

  # set legend properties
  legend.nodes = data.frame(
    label = c("More inhibited","No difference", "More activated"), color = col
  )

  # plot the network
  visNetwork(nodes, edges, main = title, width = "100%") %>%
    visLegend(addNodes = legend.nodes, useGroups = FALSE,
              main = "Good model activity state", zoom = FALSE)
}

#' Plot the graph of average state differences (igraph)
#'
#' This function uses the \code{\link[igraph]{igraph}} package to plot a network
#' of nodes. The nodes are positioned according to the specified coordinates
#' given by the \code{layout} parameter and the colors are derived using the
#' \code{diff} values and the \code{\link{get_node_colors}} function. The color
#' of each node indicates how much more inhibited or active that node is, when
#' comparing the average model classified in the 'good' category vs the average
#' 'bad' one.
#'
#' @param net igraph graph object
#' @param diff numeric vector. Every value is in the [-1,1] interval and
#' represents the average activity difference of each node. The node names have
#' to be specified in the \emph{names} attribute of the given vector. For example,
#' \code{diff} could be the result of using the function
#' \code{\link{get_avg_activity_diff_based_on_tp_predictions}}.
#' @param layout a (nx2) numeric matrix of x-y coordinates (2 columns) for each
#' of the nodes (n) in the \code{net} igraph object
#' @param title string. The title of the igraph plot
#'
#' @family network plotting functions
#'
#' @seealso
#' \code{\link{get_node_colors}}
#'
#' @importFrom igraph plot.igraph V
#' @importFrom graphics legend
#' @export
plot_avg_state_diff_graph = function(net, diff, layout, title) {
  # colors for nodes (to be interpolated) matching one-to-one the diff values
  col = c("tomato", "grey", "gold")
  V(net)$color = get_node_colors(net, diff, col)

  # plot the network
  par(mar = c(0, 0, 1, 0)) # c(bottom, left, top, right)
  plot.igraph(net, asp = 0, layout = layout, main = title)
  legend(x = -1.1, y = -0.7, pch = 21, col = "#777777",
        legend = c("More inhibited", "No difference", "More activated"),
        title = expression(bold("Good model activity state")),
        pt.bg = col, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)
}

#' Plot the graph of average link operator differences (igraph)
#'
#' This function uses the \code{\link[igraph]{igraph}} package to plot a network
#' of nodes. The nodes are positioned according to the specified coordinates
#' given by the \code{layout} parameter and the colors are derived using the
#' \code{diff} values and the \code{\link{get_node_colors}} function. The color
#' of each node indicates if the node's boolean function has on average the
#' \emph{AND NOT} or the \emph{OR NOT} link operator when comparing the average
#' model classified in the 'good' category vs the average bad' one. A non-colored
#' node (white) will indicate nodes that do not have the link operator in their
#' respective boolean equation (where they function as the target).
#'
#' @param net igraph graph object
#' @param diff numeric vector. Every value is in the [-1,1] interval and
#' represents the average link operator value difference of each node. The node
#' names have to be specified in the \emph{names} attribute of the given vector.
#' For example, \code{diff} could be the result of using the function
#' \code{\link{get_avg_link_operator_diff_mat_based_on_tp_predictions}} and
#' getting one vector row from the output matrix.
#' A value closer to -1 means that the 'good' models have more of the \emph{AND NOT}
#' link operator in their respective boolean equations while a value closer to 1
#' means that the 'good' models have more of the \emph{OR NOT} link operator.

#' @param layout a (nx2) numeric matrix of x-y coordinates (2 columns) for each
#' of the nodes (n) in the \code{net} igraph object
#' @param title string. The title of the igraph plot
#'
#' @family network plotting functions
#'
#' @seealso
#' \code{\link{get_node_colors}}
#'
#' @importFrom igraph plot.igraph V
#' @importFrom graphics legend
#' @export
plot_avg_link_operator_diff_graph = function(net, diff, layout, title) {
  # colors for nodes (to be interpolated) matching one-to-one the diff values
  col = c("tomato", "grey", "gold")
  col.2 = c("tomato", "gold", "grey", "white")
  V(net)$color = get_node_colors(net, diff, col)

  # plot the network
  par(mar = c(0, 0, 1, 0)) # c(bottom, left, top, right)
  plot.igraph(net, asp = 0, layout = layout, main = title)
  legend(x = -1.1, y = -0.6, pch = 21, col = "#777777",
         legend = c("AND NOT", "OR NOT", "No difference", "no link operator"),
         title = expression(bold("Good model link operator")),
         pt.bg = col.2, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)
}

#' Get the node colors
#'
#' This function splits the [-1,1] interval into \strong{2000} smaller
#' ones and matches each value of the \code{diff} vector to a specific hex color
#' code, using a spline interpolation of the colors as defined in the \code{col}
#' parameter.
#'
#' @param net an igraph graph object with the node names defined in \code{V(net)$name}
#' @param diff numeric vector. Every value is in the [-1,1] interval and
#' represents the average activity difference of each node. The node names have
#' to be specified in the \emph{names} attribute of the given \code{diff} vector
#' and have to be the same as in \code{V(net)$name}.
#' @param col a character vector of colors to do the color interpolation in the
#' [-1,1] interval. Usually a two-element vector specifying the colors matching
#' the start and end of the interval (-1 and 1 respectively) or a three-element
#' vector specifying the colors matching the values -1, 0 and 1 (can be more of
#' course, you get the idea).
#'
#' @return a character vector of hex color codes where the \emph{names} attribute
#' corresponds to the nodes of the given igraph object. Will be used to fill in
#' the \code{V(net)$color} property of the \code{net} object. If there are nodes
#' that are part of the network object \code{net} but not present in the \code{diff}
#' vector, then a \emph{NA} value will be given for the color of these nodes.
#'
#' @importFrom igraph V
#' @importFrom grDevices colorRampPalette
#' @export
get_node_colors = function(net, diff, col) {
  # 2000 equal-sized intervals for values between [-1,1]
  # for significance up to the 3rd decimal
  num.of.intervals = 2000

  # make the color of each node match the corresponding diff value
  color.palette = colorRampPalette(col, interpolate = "spline")
  color.values = color.palette(num.of.intervals)
  # check the colors
  # plot(x = 1:2000, y = 1:2000, cex = 10, pch = 20, col = color.values)

  diff.extra = c(diff, -1, 1)

  interval.ids.extra = as.numeric(
    cut(diff.extra, breaks = num.of.intervals, include.lowest = TRUE)
  )

  # remove the last two values
  interval.ids = interval.ids.extra[1:(length(interval.ids.extra) - 2)]
  diff.colors = color.values[interval.ids]
  names(diff.colors) = names(diff)

  # check that nodes for whom you found the colors are part of the net object nodes
  node.names = V(net)$name
  stopifnot(all(names(diff.colors) %in% node.names))

  # re-order based on the net object's node sequence
  # For the net object nodes that are not in the 'colored' nodes, we get NA values
  return(diff.colors[node.names])
}

#' Plot the graphs from an average state differences matrix
#'
#' This function presents a convenient way to use many times the
#' \code{\link{plot_avg_state_diff_graph}} function.
#'
#' @param net igraph graph object
#' @param diff.mat a matrix whose rows are \strong{vectors of average node activity
#' state differences} between two groups of models based on some kind of classification
#' (e.g. number of TP predictions) and whose names are set in the \code{rownames}
#' attribute of the matrix (usually denoting the different classification
#' groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
#' that predicted 2 TP synergies, if the classification is done by number of TP
#' predictions). The columns represent the network's node names.
#' @param layout a (nx2) numeric matrix of x-y coordinates (2 columns) for each
#' of the nodes (n) in the \code{net} igraph object
#'
#' @family network plotting functions
#'
#' @export
plot_avg_state_diff_graphs = function(net, diff.mat, layout) {
  for (row.index in 1:nrow(diff.mat)) {
    plot_avg_state_diff_graph(net, diff.mat[row.index, ], layout = layout,
                              title = rownames(diff.mat)[row.index])
  }
}

#' Plot the graphs from an average link operator differences matrix
#'
#' This function presents a convenient way to use many times the
#' \code{\link{plot_avg_link_operator_diff_graph}} function.
#'
#' @param net igraph graph object
#' @param diff.mat a matrix whose rows are \strong{vectors of average node link
#' operator differences} between two groups of models based on some kind of
#' classification (e.g. number of TP predictions) and whose names are set in the \code{rownames}
#' attribute of the matrix (usually denoting the different classification
#' groups, e.g. (1,2) means the models that predicted 1 TP synergy vs the models
#' that predicted 2 TP synergies, if the classification is done by number of TP
#' predictions). The columns represent the network's node names.
#' @param layout a (nx2) numeric matrix of x-y coordinates (2 columns) for each
#' of the nodes (n) in the \code{net} igraph object
#'
#' @family network plotting functions
#'
#' @export
plot_avg_link_operator_diff_graphs = function(net, diff.mat, layout) {
  for (row.index in 1:nrow(diff.mat)) {
    plot_avg_link_operator_diff_graph(net, diff.mat[row.index, ], layout = layout,
                                      title = rownames(diff.mat)[row.index])
  }
}

