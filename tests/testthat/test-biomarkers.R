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
