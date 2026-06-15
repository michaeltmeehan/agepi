test_that("expand_age_rates copies annual rates to nested finer ages", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  rates <- data.frame(
    age_group = coarse_ages$age_groups,
    value = c(0.02, 0.05, 0.1),
    stringsAsFactors = FALSE
  )

  expanded <- expand_age_rates(rates, mapping)

  expect_identical(expanded$age_group, fine_ages$age_groups)
  expect_equal(expanded$value, c(rep(0.02, 5), rep(0.05, 5), 0.1))
})

test_that("expand_age_hazards copies hazards to nested finer ages", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  hazards <- data.frame(
    age_group = coarse_ages$age_groups,
    mortality_rate = c(0.01, 0.03, 0.08),
    stringsAsFactors = FALSE
  )

  expanded <- expand_age_hazards(hazards, mapping, value_col = "mortality_rate")

  expect_identical(expanded$age_group, fine_ages$age_groups)
  expect_equal(expanded$mortality_rate, c(rep(0.01, 5), rep(0.03, 5), 0.08))
})

test_that("expand_age_interval_probabilities composes to the coarse probability", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  probabilities <- data.frame(
    age_group = coarse_ages$age_groups,
    value = c(0.20, 0.10, 0.40),
    stringsAsFactors = FALSE
  )

  expanded <- expand_age_interval_probabilities(probabilities, mapping)

  expect_equal(expanded$value[1:5], rep(1 - (1 - 0.20)^(1 / 5), 5))
  expect_equal(expanded$value[6:10], rep(1 - (1 - 0.10)^(1 / 5), 5))
  expect_equal(expanded$value[11], 0.40)

  composed_first <- 1 - prod(1 - expanded$value[1:5])
  composed_second <- 1 - prod(1 - expanded$value[6:10])
  expect_equal(composed_first, 0.20)
  expect_equal(composed_second, 0.10)
})

test_that("expand_age_interval_probabilities validates probability bounds", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  low <- data.frame(
    age_group = coarse_ages$age_groups,
    value = c(-0.01, 0.1, 0.2),
    stringsAsFactors = FALSE
  )
  high <- data.frame(
    age_group = coarse_ages$age_groups,
    value = c(0.01, 1.1, 0.2),
    stringsAsFactors = FALSE
  )

  expect_error(expand_age_interval_probabilities(low, mapping), "between 0 and 1")
  expect_error(expand_age_interval_probabilities(high, mapping), "between 0 and 1")
})

test_that("rate expansion preserves non-age grouping columns", {
  coarse_ages <- wpp_age_structure_5year(max_age = 5)
  fine_ages <- wpp_age_structure_1year(max_age = 5)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  rates <- data.frame(
    time = rep(c(2020, 2021), each = 4),
    sex = rep(c("female", "male"), each = 2, times = 2),
    compartment = rep("S", times = 8),
    scenario = rep(c("baseline", "high"), times = 4),
    age_group = rep(coarse_ages$age_groups, times = 4),
    rate = seq(0.01, 0.08, by = 0.01),
    stringsAsFactors = FALSE
  )

  expanded <- expand_age_rates(rates, mapping, value_col = "rate")

  expect_identical(names(expanded), c("time", "sex", "compartment", "scenario", "age_group", "rate"))
  expect_equal(nrow(expanded), 24)
  expect_equal(expanded$rate[1:6], c(rep(0.01, 5), 0.02))
  expect_identical(expanded$time[1:6], rep(2020, 6))
  expect_identical(expanded$sex[1:6], rep("female", 6))
})

test_that("rate expansion requires an expansion mapping", {
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  aggregate_mapping <- AgeGridMapping(fine_ages, coarse_ages, open_ended = "include")

  rates <- data.frame(
    age_group = fine_ages$age_groups,
    value = rep(0.01, fine_ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  expect_error(
    expand_age_rates(rates, aggregate_mapping),
    "does not support expansion"
  )
})

test_that("non-nested age grids error clearly before rate expansion", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = c("0-6", "7-9"),
    lower_bounds = c(0, 7),
    upper_bounds = c(6, 9)
  )

  expect_error(
    AgeGridMapping(from_ages, to_ages),
    "exact nested aggregation or expansion|exactly covered"
  )
})
