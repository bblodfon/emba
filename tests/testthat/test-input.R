context("Testing 'get_alt_drugname'")

test_that("it returns proper results", {
  expect_equal(get_alt_drugname("A-B"), "B-A")
  expect_equal(get_alt_drugname("BB-AA"), "AA-BB")
})

context("Testing 'is_comb_element_of'")

test_that("it returns proper results", {
  expect_true(is_comb_element_of("A-B", c("E-F", "A-B")))
  expect_true(is_comb_element_of("B-A", c("E-F", "A-B")))

  expect_false(is_comb_element_of("A-B", c("E-F", "A-D")))
  expect_false(is_comb_element_of("A-B", c()))
})

context("Testing 'get_model_predictions'")

test_that("it returns proper results", {
  model.predictions.file = system.file("extdata", "model_predictions", package = "emba", mustWork = TRUE)
  model.predictions = get_model_predictions(model.predictions.file)

  expect_equal(rownames(model.predictions)[1], "topology.sif_run_961__G19_M5.gitsbe")
  expect_equal(rownames(model.predictions)[2], "topology.sif_run_2402__G19_M9.gitsbe")

  expect_equal(colnames(model.predictions)[1], "5Z-AK")
  expect_equal(colnames(model.predictions)[2], "5Z-BI")

  expect_equal(model.predictions[1,39], 1)
  expect_equal(model.predictions[1,42], 0)
  expect_true(is.na(model.predictions[1,45]))
})
