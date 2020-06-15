## dff => model.predictions
df = model_predictions_df %>% as.data.frame() # from R/sysdata.rda
dff = df[c(567, 736, 810, 1000, 1000, 567, 736, 736),]
rownames(dff) = c("test", "test2", "test3", "test4", "test5", "test6", "test7", "test8")

## models_ss_cont => models.stable.states (from R/sysdata.rda)
models_ss = models_ss_cont[rownames(dff),] # keep the order as in the dff (from R/sysdata.rda)

## models_lo => model.link.operator
models_lo = round(models_ss)
models_lo["test2", c("MAP2K6", "NLK")] = 0.5

context("Test general functions")
test_that("proper input checking is done", {
  # threshold between 0 and 1
  expect_error(biomarker_tp_analysis(threshold = 3))
  expect_error(biomarker_mcc_analysis(threshold = -0.1))
  expect_error(biomarker_synergy_analysis(threshold = 2))

  # all observed synergies should be in `model.predictions` columns (and at least one should be given)
  expect_error(biomarker_tp_analysis(model.predictions = dff,
    models.stable.state = models_ss, observed.synergies = c("r-u", "x-r"),
    threshold = 0.7))
  expect_error(biomarker_mcc_analysis(model.predictions = dff,
    models.stable.state = models_ss, observed.synergies = c(), threshold = 0.7))
  expect_error(biomarker_synergy_analysis(model.predictions = dff,
    models.stable.state = models_ss, observed.synergies = c("a-a", "l-l"), threshold = 0.7))

  # same model name order
  models_ss_2 = models_ss
  rownames(models_ss_2) = rev(rownames(models_ss_2))
  expect_error(biomarker_tp_analysis(model.predictions = dff,
    models.stable.state = models_ss_2, observed.synergies = c("r-u"), threshold = 0.7))
  expect_error(biomarker_mcc_analysis(model.predictions = dff,
    models.stable.state = models_ss_2, observed.synergies = c("r-u"), threshold = 0.7))
  expect_error(biomarker_synergy_analysis(model.predictions = dff,
    models.stable.state = models_ss_2, observed.synergies = c("r-u"), threshold = 0.7))

  # `model.stable.states` between 0 and 1
  models_ss_3 = models_ss
  models_ss_3["test4", "NLK"] = 1.001
  expect_error(biomarker_tp_analysis(model.predictions = dff,
    models.stable.state = models_ss_3, observed.synergies = c("r-u"), threshold = 0.7))
  expect_error(biomarker_mcc_analysis(model.predictions = dff,
    models.stable.state = models_ss_3, observed.synergies = c("r-u"), threshold = 0.7))
  expect_error(biomarker_synergy_analysis(model.predictions = dff,
    models.stable.state = models_ss_3, observed.synergies = c("r-u"), threshold = 0.7))

  # `models.link.operator` can be one of 0, 1 or 0.5
  models_lo_2 = models_lo
  models_lo_2["test3", "MAP2K6"] = 0.2
  models_lo_2["test3", "NLK"] = 0.99
  models_lo_2["test2", "NLK"] = 99
  expect_error(suppressWarnings(biomarker_tp_analysis(model.predictions = dff,
    models.stable.state = models_ss, models.link.operator = models_lo_2,
    observed.synergies = c("r-u"), threshold = 0.7)))
  expect_error(suppressWarnings(biomarker_mcc_analysis(model.predictions = dff,
    models.stable.state = models_ss, models.link.operator = models_lo_2,
    num.of.mcc.classes = 2, observed.synergies = c("r-u"), threshold = 0.7)))
  expect_error(biomarker_synergy_analysis(model.predictions = dff,
    models.stable.state = models_ss, models.link.operator = models_lo_2,
    calculate.subsets.stats = TRUE, observed.synergies = c("r-u", "f-l"), threshold = 0.7))

  # `penalty` between 0 and 1
  expect_error(suppressWarnings(biomarker_tp_analysis(model.predictions = dff,
    models.stable.state = models_ss, observed.synergies = c("r-u"),
    threshold = 0.7, penalty = 100)))
  expect_error(suppressWarnings(biomarker_mcc_analysis(model.predictions = dff,
    models.stable.state = models_ss, observed.synergies = c("r-u"),
    threshold = 0.7, penalty = -0.001)))
  expect_error(biomarker_synergy_analysis(model.predictions = dff,
    models.stable.state = models_ss, observed.synergies = c("r-u"), threshold = 0.7, penalty = 1.1))

  # check results
  res_list = biomarker_synergy_analysis(model.predictions = dff,
    models.stable.state = models_ss, models.link.operator = models_lo,
    calculate.subsets.stats = TRUE, observed.synergies = c("r-u", "f-l"),
    threshold = 0.7, penalty = 0.3)

  expect_equal(length(res_list), 7)
  expect_equal(names(res_list), c("predicted.synergies", "synergy.subset.stats",
    "synergy.comparison.sets", "diff.state.synergies.mat", "activity.biomarkers",
    "diff.link.synergies.mat", "link.operator.biomarkers"))
  expect_equal(res_list$predicted.synergies, c("r-u", "f-l"))
  expect_equal(names(res_list$synergy.subset.stats), c("", "r-u", "f-l", "r-u,f-l"))
  expect_equal(unlist(unname(res_list$synergy.comparison.sets[1, ])),
    c("f-l", "r-u,f-l", "r-u"))
  expect_true(all(res_list$diff.state.synergies.mat >= -1, res_list$diff.state.synergies.mat <= 1))
  expect_true(all(res_list$diff.link.synergies.mat >= -1, res_list$diff.link.synergies.mat <= 1))
  expect_equal(res_list$activity.biomarkers["f-l", "NLK"], 0)
  expect_equal(res_list$link.operator.biomarkers["f-l", "NLK"], -1)
})
