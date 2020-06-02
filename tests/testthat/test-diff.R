# Common data used in all tests

## dff => model.predictions
df = model_predictions_df %>% as.data.frame() # from R/sysdata.rda
dff = df[c(567, 736, 810, 1000, 1000, 567, 736, 736),]
rownames(dff) = c("test", "test2", "test3", "test4", "test5", "test6", "test7", "test8")

models_mcc = c(-0.04, -0.17, 0.15, -0.24, -0.02 , 0.27, -0.42 , 0.38)
names(models_mcc) = rownames(dff)

models_tp = c(1,1,2,1,3,3,2,1)
names(models_tp) = rownames(dff)

## models_ss => models.stable.state, models_lo => models.link.operators
models_dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
models_ss = get_stable_state_from_models_dir(models_dir)
models_lo = get_link_operators_from_models_dir(models_dir)
models_ss = models_ss[, 1:5]
models_lo = models_lo[, 1:5]
df_new = as.data.frame(matrix(c(1,0,1,1,0), ncol = 5, nrow = 5))
colnames(df_new) = colnames(models_ss)
rownames(df_new) = c("test4", "test5", "test6", "test7", "test8")
models_ss = rbind(models_ss, df_new)
models_ss = models_ss[rownames(dff),] # keep the order as in the dff
colnames(df_new) = colnames(models_lo)
models_lo = rbind(models_lo, df_new[-1,])
models_lo = models_lo[rownames(dff),] # keep the order as in the dff
models_ss_cont = models_ss_cont[rownames(dff),] # keep the order as in the dff (from R/sysdata.rda)

context("Testing 'get_avg_{activity/link_operator}_diff_mat_based_on_tp_predictions'")
test_that("it returns proper results", {
  res1 = get_avg_activity_diff_mat_based_on_tp_predictions(models_tp, models_ss)
  res2 = get_avg_activity_diff_mat_based_on_tp_predictions(models_tp, models_ss, penalty = 0.1)
  res3 = get_avg_link_operator_diff_mat_based_on_tp_predictions(models_tp, models_lo)
  res4 = get_avg_link_operator_diff_mat_based_on_tp_predictions(models_tp, models_lo, penalty = 0.1)
  res5 = get_avg_activity_diff_mat_based_on_tp_predictions(models_tp, models_ss_cont, penalty = 0.1)

  expect_equal(dim(res1), c(3,5))
  expect_equal(dim(res2), c(3,5))
  expect_equal(dim(res3), c(3,5))
  expect_equal(dim(res4), c(3,5))
  expect_equal(dim(res5), c(3,5))

  expect_equal(colnames(res1), colnames(models_ss))
  expect_equal(colnames(res2), colnames(models_ss))
  expect_equal(colnames(res3), colnames(models_lo))
  expect_equal(colnames(res4), colnames(models_lo))
  expect_equal(colnames(res5), colnames(models_ss_cont))

  expect_equal(rownames(res1), c("(1,2)", "(1,3)", "(2,3)"))
  expect_equal(rownames(res2), c("(1,2)", "(1,3)", "(2,3)"))
  expect_equal(rownames(res3), c("(1,2)", "(1,3)", "(2,3)"))
  expect_equal(rownames(res4), c("(1,2)", "(1,3)", "(2,3)"))
  expect_equal(rownames(res5), c("(1,2)", "(1,3)", "(2,3)"))

  expect_equal(unname(res1[1,]), c(0.5, 0.5, 0.5, 0.5, 0.25))
  expect_equal(unname(res3[1,]), c(0.5, 0, 0.75, 0.5, 0.25))
  expect_equal(res5["(1,2)", "NLK"], 0.343, tolerance = 0.01)

  # same number of models in each group, same result no matter the penalty
  expect_equal(res1[3,], res2[3,])
  expect_equal(res3[3,], res4[3,])

  # test the penalty effect
  expect_true(all(res2[1,] <= res1[1,]))
  expect_true(all(res4[1,] <= res3[1,]))
})

context("Testing 'get_avg_activity_diff_based_on_tp_predictions'")
test_that("it does proper input check", {
  expect_error(get_avg_activity_diff_based_on_tp_predictions(models_tp, models_ss, 2, 1))
})

context("Testing 'get_avg_{activity/link_operator}_diff_mat_based_on_mcc_clustering'")
test_that("it returns proper results", {
  expect_error(get_avg_activity_diff_mat_based_on_mcc_clustering(models_mcc, models_ss,
    num.of.mcc.classes = 1))

  res1 = get_avg_activity_diff_mat_based_on_mcc_clustering(models_mcc, models_ss,
    num.of.mcc.classes = 3)
  res2 = get_avg_activity_diff_mat_based_on_mcc_clustering(models_mcc, models_ss,
    num.of.mcc.classes = 3, penalty = 0.3)
  res3 = get_avg_activity_diff_mat_based_on_mcc_clustering(models_mcc, models_ss,
    num.of.mcc.classes = 2, penalty = 0.1)
  res4 = get_avg_link_operator_diff_mat_based_on_mcc_clustering(models_mcc, models_lo,
    num.of.mcc.classes = 3, penalty = 0.4)
  res5 = get_avg_link_operator_diff_mat_based_on_mcc_clustering(models_mcc, models_ss_cont,
    num.of.mcc.classes = 3, penalty = 0.1)

  expect_equal(dim(res1), c(3,5))
  expect_equal(dim(res2), c(3,5))
  expect_equal(dim(res3), c(1,5))
  expect_equal(dim(res4), c(3,5))
  expect_equal(dim(res5), c(3,5))

  expect_equal(colnames(res1), colnames(models_ss))
  expect_equal(colnames(res2), colnames(models_ss))
  expect_equal(colnames(res3), colnames(models_ss))
  expect_equal(colnames(res4), colnames(models_lo))
  expect_equal(colnames(res5), colnames(models_ss_cont))

  expect_equal(rownames(res1), c("(1,2)", "(1,3)", "(2,3)"))
  expect_equal(rownames(res2), c("(1,2)", "(1,3)", "(2,3)"))
  expect_equal(rownames(res3), c("(1,2)"))
  expect_equal(rownames(res4), c("(1,2)", "(1,3)", "(2,3)"))
  expect_equal(rownames(res5), c("(1,2)", "(1,3)", "(2,3)"))

  expect_equal(unname(res3[1,]), c(0.0633, 0.0633, 0.0633, 0.0633, -0.0633), tolerance = 0.0001)
  expect_equal(res5["(1,3)", "NLK"], 0.0095, tolerance = 0.0001)

  # same number of models in each group, same result no matter the penalty
  expect_equal(res1[3,], res2[3,])

  # test the penalty effect
  expect_true(all(res1[1,] <= res2[1,]))
  expect_true(all(res1[2,] <= res2[2,]))
})

context("Testing 'get_avg_activity_diff_based_on_mcc_clustering'")
test_that("it does proper input check", {
  expect_error(get_avg_activity_diff_based_on_mcc_clustering(models.mcc = models_mcc,
    models.stable.state = models_ss, mcc.class.ids = c(1,2),
    models.cluster.ids = c(1,1,2,1,2,1,1,2), class.id.low = 2, class.id.high = 1))
})

context("Testing 'get_avg_{activity/link_operator}_diff_mat_based_on_specific_synergy_prediction'")
test_that("it returns proper results", {
  expect_error(get_avg_activity_diff_mat_based_on_specific_synergy_prediction(
    dff, models_ss, predicted.synergies = c("NO!", "YES!")))

  synergies = c("r-u", "f-l")
  res1 = get_avg_activity_diff_mat_based_on_specific_synergy_prediction(model.predictions = dff, models.stable.state = models_ss, predicted.synergies = synergies, penalty = 0.2)
  res2 = get_avg_link_operator_diff_mat_based_on_specific_synergy_prediction(model.predictions = dff, models.link.operator = models_lo, predicted.synergies = synergies)
  res3 = get_avg_activity_diff_mat_based_on_specific_synergy_prediction(model.predictions = dff, models.stable.state = models_ss_cont, predicted.synergies = synergies, penalty = 0.1)

  expect_equal(unname(res1[1,]), c(-0.2904, -0.2904, -0.2904, -0.2904, 0.2904), tolerance = 0.0001)
  expect_equal(unname(res2[2,]), c(0.5, 0.5, 0, 0.5, 0))
  expect_equal(res3["f-l", "NLK"], -0.271, tolerance = 0.01)

  expect_equal(colnames(res1), colnames(models_ss))
  expect_equal(colnames(res2), colnames(models_lo))
  expect_equal(colnames(res3), colnames(models_ss_cont))

  expect_equal(rownames(res1), synergies)
  expect_equal(rownames(res2), synergies)
  expect_equal(rownames(res3), synergies)
})

context("Testing 'get_avg_activity_diff_based_on_specific_synergy_prediction'")
test_that("it does proper input check", {
  # input check
  expect_error(get_avg_activity_diff_based_on_specific_synergy_prediction(
    dff, models_ss, drug.comb = "NO!"))
})

context("Testing 'get_avg_{activity/link_operator}_diff_based_on_synergy_set_cmp'")
test_that("it does proper input checks", {
  # length(synergy.subset) > 0
  expect_error(get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "A-B", synergy.subset.str = "", NULL, NULL)
  )

  # length(synergy.set) > length(synergy.subset)
  expect_error(get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "A-B,C-D", synergy.subset.str = "A-B,C-D,A-C", NULL, NULL)
  )

  # all(synergy.subset %in% synergy.set)
  expect_error(get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "A-B,C-D,E-F", synergy.subset.str = "A-B,E-G", NULL, NULL)
  )
})

test_that("it returns proper results", {
  diff1 = get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "r-u,i-k", synergy.subset.str = "i-k",
    model.predictions = dff, models.stable.state = models_ss)
  diff2 = get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "r-u,i-k,g-o", synergy.subset.str = "r-u,g-o",
    model.predictions = dff, models.stable.state = models_ss)
  diff3 = get_avg_link_operator_diff_based_on_synergy_set_cmp(
    synergy.set.str = "i-k,g-o,w-x,n-s,b-m,c-y", synergy.subset.str = "i-k,c-y",
    model.predictions = dff, models.link.operator = models_lo)
  diff4 = get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "r-u,i-k,g-o", synergy.subset.str = "r-u,g-o",
    model.predictions = dff, models.stable.state = models_ss_cont)

  expect_equal(unname(diff1), c(-0.4, -0.4, -0.4, -0.4, 0.4))
  expect_equal(unname(diff2), c(0.1, 0.1, 0.1, 0.1, -0.1))
  expect_equal(unname(diff3), c(-0.2, -0.6, 0.8, -0.2, 0.4))
  expect_equal(unname(diff4), c(0.007, -0.0234, -0.016, 0.2098, 0.097), tolerance = 0.001)

  expect_equal(names(diff1), colnames(models_ss))
  expect_equal(names(diff2), colnames(models_ss))
  expect_equal(names(diff3), colnames(models_lo))
  expect_equal(names(diff4), colnames(models_ss_cont))
})

context("Testing 'get_vector_diff'")
test_that("it returns proper results", {
  # input check
  expect_error(get_vector_diff(vec1 = c(1,2), vec2 = c(3)))
  expect_error(get_vector_diff(vec1 = c(1,2), vec2 = c(3,4), penalty = 3))

  vec1 = c(1,2,3,2,1)
  vec2 = c(3,2,1,3,3)

  res = c(-2,0,2,-1,-2)
  expect_equal(get_vector_diff(vec1, vec2), res)

  m1 = -3
  m2 = 0
  expect_equal(get_vector_diff(vec1, vec2, m1, m2), res)

  m1 = 1
  m2 = 1
  expect_equal(get_vector_diff(vec1, vec2, m1, m2), res)

  m1 = 100
  m2 = 3
  # default value for penalty is 0
  expect_equal(get_vector_diff(vec1, vec2, m1, m2), res)
  expect_equal(get_vector_diff(vec1, vec2, m1, m2, penalty = 0), res)
  expect_equal(get_vector_diff(vec1, vec2, m1, m2, penalty = 0.1), res * (m2/m1)^0.1)
  expect_equal(get_vector_diff(vec1, vec2, m1, m2, penalty = 0.5), res * (m2/m1)^0.5)
  expect_equal(get_vector_diff(vec1, vec2, m1, m2, penalty = 1), res * (m2/m1))
  expect_equal(get_vector_diff(vec1, vec2, m2, m1, penalty = 1), res * (m2/m1))

  # if m1=m2, penalty does not matter
  expect_equal(get_vector_diff(vec1, vec2, m1 = 1000, m2 = 1000, penalty = 0.5), res)

  names(vec1) = letters[1:5]
  names(vec2) = letters[6:10]
  expect_equal(names(get_vector_diff(vec1, vec2, m1, m2)), names(vec1))
})


