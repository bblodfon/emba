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
