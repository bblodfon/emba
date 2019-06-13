#' Plot histogram of the MCC classes
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' @importFrom Ckmeans.1d.dp ahist
plot_mcc_classes_hist = function(models.mcc.no.nan.sorted, models.cluster.ids,
                                 num.of.classes, mcc.class.ids) {
  min.x.value = round(min(models.mcc.no.nan.sorted) - 0.1, digits = 1)
  max.x.value = round(max(models.mcc.no.nan.sorted) + 0.1, digits = 1)
  rainbow.colors = rainbow(num.of.classes)

  ahist(models.mcc.no.nan.sorted, k = num.of.classes,
        main = "Model MCC-Classification", xlab = "MCC value",
        sub = paste("n =", length(models.mcc.no.nan.sorted), "models, k =",
                    num.of.classes, "classes"),
        col = rainbow.colors, col.stick = rainbow.colors[models.cluster.ids],
        xlim = c(min.x.value, max.x.value))

  legend("topright", legend = mcc.class.ids, title = "MCC Classes",
         col = rainbow.colors, lty = 1, lwd = 10)
}

#' Bar plot of model stats
#'
#' Use this function to produce a bar plot when the input is a \link[base]{table}
#' command to a numeric vector
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
#'
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
#' that predicted each synergy subset, defined as a comma-seperated string of
#' drug combinations in the \emph{names} attribute of the vector
#' @param threshold.for.subset.removal integer. Use it to discard elements of
#' the \code{synergy.subset.stats} vector that are stricly less than the
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

#' plot network using the `visNetwork` library
#'
#' A description
#'
#' @param net a net
#' @param diff a diff
#' @param layout layout
#' @param title a title
#'
#' @importFrom magrittr %>%
#' @importFrom visNetwork toVisNetworkData visNetwork visLegend
#' @export
plot_network_vis = function(net, diff, layout, title) {
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

#' plot network using the `igraph` library
#'
#' @importFrom igraph plot.igraph V
#' @importFrom graphics legend
plot_network = function(net, diff, layout, title) {
  # colors for nodes (to be interpolated) matching one-to-one the diff values
  col = c("tomato", "grey", "gold")
  V(net)$color = get_node_colors(net, diff, col)

  # plot the network
  par(mar = c(0, 0, 1, 0)) # c(bottom, left, top, right)
  plot.igraph(net, asp = 0, layout = layout, main = title)
  legend(x = -1.1, y = -0.7, pch = 21, col = "#777777",
        legend = c("More inhibited","No difference", "More activated"),
        title = expression(bold("Good model activity state")),
        pt.bg = col, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)
}

#' `net` is an igraph network object with node labels as: `V(net)$name`
#'
#' @importFrom igraph V
#' @importFrom grDevices colorRampPalette
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

  # re-order based on the net object's node sequence
  node.names = V(net)$name
  return(diff.colors[node.names])
}

# `diff.df` is a data.frame whose rows are classification group comparisons
# with their name set in the `rownames()` and the columns are the network's
# node names. This function is used to use `plot.network` multiple times on
# different `diff` vectors
plot_diff_df = function(net, diff.df, layout) {
  for (row.index in 1:nrow(diff.df)) {
    plot_network(net, diff.df[row.index, ], layout = layout,
                 title = rownames(diff.df)[row.index])
  }
}

