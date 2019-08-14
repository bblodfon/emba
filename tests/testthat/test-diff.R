context("Testing 'get_avg_activity_diff_based_on_comparing_synergy_sets'")

test_that("it does proper input checks", {
  # length(synergy.subset) > 0
  expect_error(get_avg_activity_diff_based_on_comparing_synergy_sets(
    synergy.set.str = "A-B", synergy.subset.str = "", NULL, NULL)
  )

  # length(synergy.set) > length(synergy.subset)
  expect_error(get_avg_activity_diff_based_on_comparing_synergy_sets(
    synergy.set.str = "A-B,C-D", synergy.subset.str = "A-B,C-D,A-C", NULL, NULL)
  )

  # all(synergy.subset %in% synergy.set)
  expect_error(get_avg_activity_diff_based_on_comparing_synergy_sets(
    synergy.set.str = "A-B,C-D,E-F", synergy.subset.str = "A-B,E-G", NULL, NULL)
  )
})
