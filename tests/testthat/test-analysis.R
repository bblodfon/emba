library(dplyr)

context("Testing 'count_models_that_predict_synergies'")
test_that("it returns correct results", {
  df = model_predictions_df # from R/sysdata.rda

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = character(0),
    model.predictions = df), 25)

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('g-o'),
    model.predictions = df), 360)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('b-m'),
    model.predictions = df), 312)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('r-u'),
    model.predictions = df), 329)

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('g-o', 'b-m'),
    model.predictions = df), 119)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('b-m', 'g-o'),
    model.predictions = df), 119)
  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('r-u','w-x','n-s','c-y'),
    model.predictions = df), 13)

  expect_equal(count_models_that_predict_synergies(drug.comb.vec = c('r-u','w-x','n-s','c-y','b-m'),
    model.predictions = df), 5)
  expect_equal(count_models_that_predict_synergies(colnames(df), model.predictions = df), 0)
})

context("Testing 'calculate_mcc'")
test_that("it returns proper results", {
  tp = c(0, 2, 0, 1, 1)
  fp = c(2, 2, 0, 0, 1)
  tn = c(145, 145, 145, 1, 0)
  fn = c(5, 4, 3, 0, 0)

  expect_equal(calculate_mcc(tp, tn, fp, fn), c(-0.0213, 0.3889, 0, 1, 0), tolerance = 0.0001)
})

context("Testing 'calculate_models_mcc'")
test_that("it returns proper results", {
  df = model_predictions_df %>% as.data.frame() # from R/sysdata.rda
  dff = df[1:5,]
  rownames(dff) = c("model1", "model2", "model3", "model4", "model5")

  obs = dff %>% select(`d-e`, `r-u`)
  unobs = dff %>% select(`i-k`, `g-o`)

  res = calculate_models_mcc(observed.model.predictions = obs,
    unobserved.model.predictions = unobs,
    number.of.drug.comb.tested = 2+2)
  expect_equal(unname(res), c(0.5773503, 0, -0.5773503, 0, 0.5773503), tolerance = 0.0001)
  expect_equal(names(res), c("model1", "model2", "model3", "model4", "model5"))
})

context("Testing 'get_observed_model_predictions' and 'get_unobserved_model_predictions'")
test_that("it returns proper results", {
  df = model_predictions_df %>% as.data.frame() # from R/sysdata.rda
  dff = df[1:5,]
  rownames(dff) = c("model1", "model2", "model3", "model4", "model5")

  res = get_observed_model_predictions(model.predictions = dff,
    observed.synergies = c("d-e", "r-u"))
  expect_equal(dim(res), c(5, 2))

  res2 = get_unobserved_model_predictions(model.predictions = dff,
    observed.synergies = c("i-k", "g-o", "w-x", "n-s", "f-l", "b-m", "c-y", "z-a"))
  expect_equal(dim(res2), c(5, 2))
  expect_identical(colnames(res), colnames(res2))
})

context("Testing 'get_synergy_subset_stats' and 'get_synergy_comparison_sets'")
test_that("it returns proper results", {
  df = model_predictions_df %>% as.data.frame() # from R/sysdata.rda

  res = get_synergy_subset_stats(model.predictions = df,
    synergies = c("r-u", "i-k", "g-o", "f-l"))

  expect_equal(unname(res[1]), 25)
  expect_equal(unname(res["r-u"]), 329)
  expect_equal(unname(res["r-u,i-k"]), 111)
  expect_equal(unname(res["r-u,f-l"]), 111)
  expect_equal(unname(res["g-o,f-l"]), 120)
  expect_equal(unname(res["r-u,i-k,g-o"]), 44)
  expect_equal(unname(res["r-u,i-k,g-o,f-l"]), 13)

  res3_df = get_synergy_comparison_sets(synergy.subset.stats = res)
  expect_equal(dim(res3_df), c(28, 3))

  res4 = res3_df %>% filter(sets == "i-k,g-o")
  expect_equal(dim(res4), c(2, 3))
  expect_equal(res4 %>% pull(synergies), c("g-o", "i-k"))
  expect_equal(res4 %>% pull(sets), rep("i-k,g-o", 2))
  expect_equal(res4 %>% pull(subsets), c("i-k", "g-o"))
})

context("Testing 'update_biomarker_files'")
test_that("it returns proper results", {
  biomarkers.dir = system.file("extdata", "biomarkers", package = "emba", mustWork = TRUE)
  biomarkers.active.prev = read.table(paste0(biomarkers.dir, "/A-B_biomarkers_active"),
    stringsAsFactors = FALSE)
  biomarkers.inhibited.prev = read.table(paste0(biomarkers.dir, "/A-B_biomarkers_inhibited"),
    stringsAsFactors = FALSE)

  biomarkers.active.new = c(0.8, 0.9)
  names(biomarkers.active.new) = c("A2", "B2")

  biomarkers.inhibited.new = c(-0.8, -0.99)
  names(biomarkers.inhibited.new) = c("B3", "B4")

  # for not skipping the test
  expect_equal(unname(biomarkers.active.new), c(0.8, 0.9))

  update_biomarker_files(biomarkers.dir = biomarkers.dir, drug.comb = "A-B",
    biomarkers.active.new = biomarkers.active.new,
    biomarkers.inhibited.new = biomarkers.inhibited.new, method = "replace")
  update_biomarker_files(biomarkers.dir = biomarkers.dir, drug.comb = "A-B",
    biomarkers.active.new = biomarkers.active.new,
    biomarkers.inhibited.new = biomarkers.inhibited.new, method = "extend")
  update_biomarker_files(biomarkers.dir = biomarkers.dir, drug.comb = "A-B",
    biomarkers.active.new = biomarkers.active.new[1],
    biomarkers.inhibited.new = biomarkers.inhibited.new[1], method = "prune.to.common")
  # no common biomarkers
  update_biomarker_files(biomarkers.dir = biomarkers.dir, drug.comb = "A-B",
    biomarkers.active.new = biomarkers.inhibited.new[1],
    biomarkers.inhibited.new = biomarkers.active.new[1])
  # remove files
  file.remove(paste0(biomarkers.dir, "/A-B_biomarkers_active"))
  file.remove(paste0(biomarkers.dir, "/A-B_biomarkers_inhibited"))
  update_biomarker_files(biomarkers.dir = biomarkers.dir, drug.comb = "A-B",
    biomarkers.active.new = biomarkers.active.new,
    biomarkers.inhibited.new = biomarkers.inhibited.new)

  # restore biomarkers files
  usefun::save_vector_to_file(vector = biomarkers.active.prev,
    file = paste0(biomarkers.dir, "/A-B_biomarkers_active"))
  usefun::save_vector_to_file(vector = biomarkers.inhibited.prev,
    file = paste0(biomarkers.dir, "/A-B_biomarkers_inhibited"))
})

