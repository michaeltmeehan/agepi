test_that("wpp_age_structure_1year creates default single-year WPP ages", {
  ages <- wpp_age_structure_1year()

  expect_s3_class(ages, "AgeStructure")
  expect_identical(ages$n_age_groups, 101L)
  expect_identical(ages$age_groups[1:5], as.character(0:4))
  expect_identical(ages$age_groups[97:101], c("96", "97", "98", "99", "100+"))
  expect_equal(ages$lower_bounds, 0:100)
  expect_equal(ages$upper_bounds, c(0:99, Inf))
  expect_true(is.infinite(ages$upper_bounds[101]))
})

test_that("wpp_age_structure_1year validates max_age", {
  expect_error(wpp_age_structure_1year(0), "positive")
  expect_error(wpp_age_structure_1year(-1), "positive")
  expect_error(wpp_age_structure_1year(100.5), "whole number")
  expect_error(wpp_age_structure_1year(NA_real_), "finite numeric scalar")
  expect_error(wpp_age_structure_1year(c(99, 100)), "finite numeric scalar")
})

test_that("wpp_age_structure_5year creates default five-year WPP ages", {
  ages <- wpp_age_structure_5year()

  expect_s3_class(ages, "AgeStructure")
  expect_identical(ages$n_age_groups, 21L)
  expect_identical(ages$age_groups[1:4], c("0-4", "5-9", "10-14", "15-19"))
  expect_identical(ages$age_groups[18:21], c("85-89", "90-94", "95-99", "100+"))
  expect_equal(ages$lower_bounds, seq(0, 100, by = 5))
  expect_equal(ages$upper_bounds, c(seq(4, 99, by = 5), Inf))
  expect_true(is.infinite(ages$upper_bounds[21]))
})

test_that("wpp_age_structure_5year validates max_age", {
  expect_error(wpp_age_structure_5year(0), "positive")
  expect_error(wpp_age_structure_5year(98), "multiple of 5")
  expect_error(wpp_age_structure_5year(100.5), "whole number")
  expect_error(wpp_age_structure_5year(NA_real_), "finite numeric scalar")
})
