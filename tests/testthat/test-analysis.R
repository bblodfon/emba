context("Testing 'get.biomarkers'")

test_that("it returns proper results", {
  m = matrix(0, 5, 5)
  set.seed(1)
  diff.res = apply(m, c(1, 2), function(x) runif(n = 1, min = -1, max = 1))
  colnames(diff.res) = c("a", "b", "c", "d", "e")

  res.1 = get_biomarkers(diff.res, threshold = 0.87, type = "active")
  expected.res.1 = c("b", "d")

  res.2 = get_biomarkers(diff.res, threshold = 0.7, type = "inhibited")
  expected.res.2 = c("b", "e")

  res.3 = get_biomarkers(diff.res, threshold = 0.99, type = "active")
  expected.res.3 = character(0)

  expect_equal(res.1, expected.res.1)
  expect_equal(res.2, expected.res.2)
  expect_equal(res.3, expected.res.3)
})
