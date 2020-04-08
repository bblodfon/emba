context("Testing 'count_models_that_predict_synergies'")

test_that("it returns correct results", {
  df = model_predictions_df # from R/sysdata.rda

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
