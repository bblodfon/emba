context("Testing 'get_avg_activity_diff_based_on_synergy_set_cmp'")
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
  df = model_predictions_df %>% as.data.frame() # from R/sysdata.rda
  dff = df[c(567, 736, 810, 1000, 1000, 567, 736, 736),]
  rownames(dff) = c("test", "test2", "test3", "test4", "test5", "test6", "test7", "test8")

  models.dir = system.file("extdata", "models", package = "emba", mustWork = TRUE)
  models.ss = get_stable_state_from_models_dir(models.dir)
  models.ss = models.ss[, 1:5]
  models.ss.new = as.data.frame(matrix(c(1,0,1,1,0), ncol = 5, nrow = 5))
  colnames(models.ss.new) = colnames(models.ss)
  rownames(models.ss.new) = c("test4", "test5", "test6", "test7", "test8")
  models.ss = rbind(models.ss, models.ss.new)

  diff1 = get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "r-u,i-k", synergy.subset.str = "i-k",
    model.predictions = dff, models.stable.state = models.ss)
  diff2 = get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "r-u,i-k,g-o", synergy.subset.str = "r-u,g-o",
    model.predictions = dff, models.stable.state = models.ss)
  diff3 = get_avg_activity_diff_based_on_synergy_set_cmp(
    synergy.set.str = "i-k,g-o,w-x,n-s,b-m,c-y", synergy.subset.str = "i-k,c-y",
    model.predictions = dff, models.stable.state = models.ss)

  expect_equal(unname(diff1), c(-0.4, -0.4, -0.4, -0.4, 0.4))
  expect_equal(unname(diff2), c(0.1, 0.1, 0.1, 0.1, -0.1))
  expect_equal(unname(diff3), c(0.4, 0.4, 0.4, 0.4, -0.4))

  expect_equal(names(diff1), colnames(models.ss))
  expect_equal(names(diff2), colnames(models.ss))
  expect_equal(names(diff3), colnames(models.ss))
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


