context("Testing 'count_models_that_predict_synergies'")

test_that("it returns correct results", {
  RNGkind(sample.kind = "Rejection")
  set.seed(0)

  df = matrix(data = sample(c(0,1,NA), size = 10000, replace = TRUE), nrow = 1000, ncol = 10)
  colnames(df) = c('d-e','r-u','i-k','g-o','w-x','n-s','f-l','b-m','c-y','z-a')

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = character(0), model.predictions = df), 25)

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('g-o'), model.predictions = df), 360)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('b-m'), model.predictions = df), 312)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('r-u'), model.predictions = df), 329)

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('g-o', 'b-m'), model.predictions = df), 119)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('b-m', 'g-o'), model.predictions = df), 119)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('r-u','w-x','n-s','c-y'), model.predictions = df), 13)

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('r-u','w-x','n-s','c-y','b-m'), model.predictions = df), 5)
  expect_equal(count_models_that_predict_synergies(colnames(df), model.predictions = df), 0)
})
