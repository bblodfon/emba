context("Testing 'get_biomarkers_per_type'")

test_that("it returns proper results", {
  m = matrix(0, 5, 5)
  set.seed(1)
  diff.mat = apply(m, c(1, 2), function(x) runif(n = 1, min = -1, max = 1))
  colnames(diff.mat) = c("a", "b", "c", "d", "e")

  res.1 = get_biomarkers_per_type(diff.mat, threshold = 0.87, type = "positive")
  expected.res.1 = c("b", "d")

  res.2 = get_biomarkers_per_type(diff.mat, threshold = 0.7, type = "negative")
  expected.res.2 = c("b", "e")

  res.3 = get_biomarkers_per_type(diff.mat, threshold = 0.99, type = "positive")
  expected.res.3 = character(0)

  expect_equal(res.1, expected.res.1)
  expect_equal(res.2, expected.res.2)
  expect_equal(res.3, expected.res.3)
})

test_that("it does correct input checks", {
  m = matrix(0, 5, 5)
  set.seed(1)
  diff.mat = apply(m, c(1, 2), function(x) runif(n = 1, min = -1, max = 1))
  colnames(diff.mat) = c("a", "b", "c", "d", "e")

  # threshold in [0,1]
  expect_error(get_biomarkers_per_type(diff.mat, threshold = -0.01))
  expect_error(get_biomarkers_per_type(diff.mat, threshold = 1.5))

  # biomarker type in {positive, negative}
  expect_error(get_biomarkers_per_type(diff.mat, threshold = 0.6, type = "inactive"))
  expect_error(get_biomarkers_per_type(diff.mat, threshold = 0.6, type = "wrong_type"))
})
