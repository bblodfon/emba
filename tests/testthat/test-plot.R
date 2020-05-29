context("Testing 'plot_mcc_classes_hist'")
test_that("plot object is created and returned", {
  models.mcc = c(-0.04, -0.17, 0.15, -0.24, -0.02 , 0.27, -0.42 , 0.38)
  models.cluster.ids = c(2,2,3,1,2,3,1,3)
  num.of.mcc.classes = 3

  expect_true(is.list(plot_mcc_classes_hist(models.mcc, models.cluster.ids, num.of.mcc.classes)))
  dev.off()
})

context("Testing 'make_barplot_on_models_stats'")
test_that("plot object is created and returned", {
  models.tp = c(rep(1,100), rep(2,423), rep(3,231), rep(NaN,531))

  expect_true(is.matrix(make_barplot_on_models_stats(models.stats = table(models.tp, useNA = "ifany"),
    title = "True Positives Distribution across models",
    xlab = "Number of TP values", ylab = "Number of models")))
  expect_true(is.matrix(make_barplot_on_models_stats(models.stats = table(models.tp),
    title = "True Positives Distribution across models", cell.line = "AGS",
    xlab = "Number of TP values", ylab = "Number of models", cont.values = TRUE)))
  dev.off()
})

context("Testing 'make_barplot_on_synergy_subset_stats'")
test_that("plot object is created and returned", {
  synergy.subset.stats = c(1,4,3,2)
  names(synergy.subset.stats) = c("A-B", "B-C", "C-A", "C-D")
  expect_true(is.matrix(make_barplot_on_synergy_subset_stats(
    synergy.subset.stats, threshold.for.subset.removal = 0,
    bottom.margin = 4, ylim.add = 0.5, cell.line = "AGS")))
  dev.off()
})

context("Testing 'plot_avg_state_diff_graph_vis'")
test_that("visNetwork object is created and returned", {
  topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
  net = construct_network(topology.file)
  diff = c(-0.95,-0.05,0.46,0.39,-0.04,0.72,-0.12,-0.51,-0.86,-0.80)
  names(diff) = c("A","C","B","D","W","I","E","J","F","K")
  res = plot_avg_state_diff_graph_vis(net, diff, nodes.size = 7.5, title = "TEST")

  # test node attributes
  expect_equal(res$x$nodes$id, names(diff))
  expect_equal(res$x$nodes$size, rep(7.5, 10))
  expect_equal(res$x$nodes$physics, rep(FALSE, 10))
  expect_equal(res$x$nodes$shape, rep("dot", 10))
  expect_equal(res$x$nodes$color, c("#F86854", "#BEBABF", "#CBD18C", "#C7CF98",
    "#BEBBBF", "#DFD654", "#BEB6C0", "#CE97A7", "#EE736A", "#E77A77"))

  # test edges attributes
  expect_equal(nrow(res$x$edges), 11)  # 11 edges
  expect_equal(res$x$edges$smooth, rep(FALSE, 11))
  expect_equal(res$x$edges$physics, rep(FALSE, 11))
  expect_equal(res$x$edges$arrows, rep("to", 11))

  # test plot attributes
  expect_equal(res$x$main$text, "TEST")
  expect_equal(res$x$width, "100%")
  expect_equal(as.vector(res$x$legend$nodes$label), c("More inhibited","No difference", "More activated"))
  expect_equal(as.vector(res$x$legend$nodes$color), c("tomato", "grey", "gold"))
  expect_null(res$x$groups)
})

context("Testing 'plot_avg_state_diff_graph'")
test_that("an igraph plot object is created and returned", {
  topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
  net = construct_network(topology.file)
  diff = c(-0.95,-0.05,0.46,0.39,-0.04,0.72,-0.12,-0.51,-0.86,-0.80)
  names(diff) = c("A","C","B","D","W","I","E","J","F","K")

  expect_true(is.list(plot_avg_state_diff_graph(net, diff, title = "TEST")))
  dev.off()
})

context("Testing 'plot_avg_link_operator_diff_graph'")
test_that("an igraph plot object is created and returned", {
  topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
  net = construct_network(topology.file)
  diff = c(-0.95,-0.05,0.46,0.39,-0.04,0.72,-0.12,-0.51)
  names(diff) = c("A","C","B","D","W","I","E","J")

  expect_true(is.list(plot_avg_link_operator_diff_graph(net, diff, title = "TEST")))
  dev.off()
})

context("Testing 'plot_avg_{state/link_operator}_diff_graphs'")
test_that("an igraph plot object is created and returned", {
  topology.file = system.file("extdata", "example.sif", package = "emba", mustWork = TRUE)
  net = construct_network(topology.file)
  diff = c(-0.95,-0.05,0.46,0.39,-0.04,0.72,-0.12,-0.51,-0.86,-0.80)

  diff.mat = matrix(data = diff, nrow = 1)
  colnames(diff.mat) = c("A","C","B","D","W","I","E","J","F","K")
  rownames(diff.mat) = "row1"

  diff.mat.2 = matrix(data = diff[1:7], nrow = 1)
  colnames(diff.mat.2) = colnames(diff.mat)[1:7]
  rownames(diff.mat.2) = "row1"

  expect_true(plot_avg_state_diff_graphs(net, diff.mat))
  expect_true(plot_avg_link_operator_diff_graphs(net, diff.mat.2))
  dev.off()
})
