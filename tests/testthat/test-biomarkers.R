# common data
m = matrix(0, 5, 5)
set.seed(1)
diff.mat = apply(m, c(1, 2), function(x) runif(n = 1, min = -1, max = 1))
colnames(diff.mat) = c("a", "b", "c", "d", "e")

context("Testing 'get_biomarkers_per_type'")
test_that("it returns proper results", {
  res1 = get_biomarkers_per_type(diff.mat, threshold = 0.87, type = "positive")
  res2 = get_biomarkers_per_type(diff.mat, threshold = 0.7, type = "negative")
  res3 = get_biomarkers_per_type(diff.mat, threshold = 0.99, type = "positive")

  expect_equal(res1, c("b", "d"))
  expect_equal(res2, c("b", "e"))
  expect_equal(res3, character(0))
})

test_that("it does correct input checks", {
  # threshold in [0,1]
  expect_error(get_biomarkers_per_type(diff.mat, threshold = -0.01))
  expect_error(get_biomarkers_per_type(diff.mat, threshold = 1.5))

  # biomarker type in {positive, negative}
  expect_error(get_biomarkers_per_type(diff.mat, threshold = 0.6, type = "inactive"))
  expect_error(get_biomarkers_per_type(diff.mat, threshold = 0.6, type = "wrong_type"))
})

context("Testing 'get_biomarkers'")
test_that("it returns proper results", {
  expect_error(get_biomarkers(diff.mat, threshold = -0.1))

  res1 = get_biomarkers(diff.mat, threshold = 0)
  res2 = get_biomarkers(diff.mat, threshold = 0.5)
  res3 = get_biomarkers(diff.mat, threshold = 0.6)
  res4 = get_biomarkers(diff.mat, threshold = 0.7)
  res5 = get_biomarkers(diff.mat, threshold = 0.8)
  res6 = get_biomarkers(diff.mat, threshold = 0.9)
  res7 = get_biomarkers(diff.mat, threshold = 1)

  expect_equal(res1$biomarkers.pos, c('c', 'd'))
  expect_equal(res1$biomarkers.neg, c('a', 'b', 'e'))
  expect_equal(res2$biomarkers.pos, c('d', 'c'))
  expect_equal(res2$biomarkers.neg, c('a', 'b', 'e'))
  expect_equal(res3$biomarkers.pos, c('a', 'd'))
  expect_equal(res3$biomarkers.neg, c('c', 'b', 'e'))
  expect_equal(res4$biomarkers.pos, c('a', 'd'))
  expect_equal(res4$biomarkers.neg, c('b', 'e'))
  expect_equal(res5$biomarkers.pos, c('a', 'd', 'e'))
  expect_equal(res5$biomarkers.neg, c('b'))
  expect_equal(res6$biomarkers.pos, c('d'))
  expect_equal(res6$biomarkers.neg, character())
  expect_equal(res7$biomarkers.pos, character())
  expect_equal(res7$biomarkers.neg, character())
})

context("Testing 'get_synergy_biomarkers_from_dir'")
test_that("it returns proper results", {
  biomarkers.dir = system.file("extdata", "biomarkers", package = "emba", mustWork = TRUE)

  # check: both `node.names` and `models.dir` are not specified
  expect_error(get_synergy_biomarkers_from_dir(biomarkers.dir = "", predicted.synergies = ""))

  res1 = get_synergy_biomarkers_from_dir(biomarkers.dir, predicted.synergies = c("A-B"),
    node.names = c('A1','A2','B3','A4','A3','B1'))
  res2 = get_synergy_biomarkers_from_dir(biomarkers.dir, predicted.synergies = c("A-B"),
    node.names = c('A1','B3','A2','A4','A3','B1','Noooo'))
  res3 = get_synergy_biomarkers_from_dir(biomarkers.dir, predicted.synergies = c("A-B"),
    node.names = c('A1','B3','A2','A4','Noooo','A3','B1'))
  res4 = get_synergy_biomarkers_from_dir(biomarkers.dir,
    predicted.synergies = c('A-B', 'C-D'), node.names = c('Noooo'))

  expect_equal(dim(res1), c(1,6))
  expect_equal(dim(res2), c(1,7))
  expect_equal(dim(res3), c(1,7))
  expect_equal(dim(res4), c(2,7))

  expect_equal(colnames(res1), c('A1','A2','B3','A4','A3','B1'))
  expect_equal(colnames(res2), c('A1','B3','A2','A4','A3','B1','Noooo'))
  expect_equal(colnames(res3), c('A1','B3','A2','A4','Noooo','A3','B1'))
  expect_equal(colnames(res4), c('Noooo','A1','A2','B1','A3','A4','B3'))

  expect_equal(unname(unlist(res1[1,])), c(1,1,-1,-1,-1,1))
  expect_equal(unname(unlist(res2[1,])), c(1,-1,1,-1,-1,1,0))
  expect_equal(unname(unlist(res3[1,])), c(1,-1,1,-1,0,-1,1))
  expect_equal(unname(unlist(res4[1,])), c(0,1,1,1,-1,-1,-1))
  expect_equal(unname(unlist(res4[2,])), rep(0,7))
})

context("Testing 'get_synergy_biomarkers_per_cell_line'")
test_that("it returns proper results", {
  biomarkers.dir = system.file("extdata", "AGS", "bio", package = "emba", mustWork = TRUE)

  res_list = get_synergy_biomarkers_per_cell_line(biomarkers.dirs = c(biomarkers.dir, tempdir()))

  expect_type(res_list, "list")
  expect_equal(length(res_list), 2)
  expect_equal(names(res_list)[1], "AGS")
  expect_equal(dim(res_list$AGS), c(2,5))
  expect_equal(colnames(res_list$AGS), paste0("x", 1:5))
  expect_equal(rownames(res_list$AGS), c("A-B", "C-D"))
  expect_equal(res_list$AGS["A-B", "x3"], 1)
  expect_equal(res_list$AGS["C-D", "x5"], -1)
})

context("Testing 'get_perf_biomarkers_per_cell_line'")
test_that("it returns proper results", {
  biomarkers.dir = system.file("extdata", "AGS", "bio", package = "emba", mustWork = TRUE)

  res = get_perf_biomarkers_per_cell_line(biomarkers.dirs = c(biomarkers.dir, tempdir()),
    node.names = paste0("x", 1:20))

  expect_s3_class(res, "data.frame")
  expect_equal(dim(res), c(2,20))
  expect_equal(colnames(res), paste0("x", 1:20))
  expect_equal(rownames(res)[2], "AGS")
  expect_equal(unname(unlist(res["AGS", c("x2","x4","x15")])), c(1,1,1))
  expect_equal(unname(unlist(res["AGS", c("x7","x11")])), c(-1,-1))
  expect_equal(unname(unlist(res[1,])), rep(0,20))
})
