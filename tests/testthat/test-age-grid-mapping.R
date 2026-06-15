test_that("AgeGridMapping aggregates one-year counts to five-year counts", {
  from_ages <- wpp_age_structure_1year(max_age = 10)
  to_ages <- wpp_age_structure_5year(max_age = 10)
  mapping <- AgeGridMapping(from_ages, to_ages, open_ended = "include")

  counts <- data.frame(
    age_group = from_ages$age_groups,
    value = seq_along(from_ages$age_groups),
    stringsAsFactors = FALSE
  )

  aggregated <- aggregate_age_counts(counts, mapping)

  expect_identical(aggregated$age_group, to_ages$age_groups)
  expect_equal(aggregated$value, c(sum(1:5), sum(6:10), 11))
})

test_that("expand_age_counts uniformly expands five-year counts to one-year counts", {
  from_ages <- wpp_age_structure_5year(max_age = 10)
  to_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(from_ages, to_ages, open_ended = "include")

  counts <- data.frame(
    age_group = from_ages$age_groups,
    value = c(50, 100, 30),
    stringsAsFactors = FALSE
  )

  expanded <- expand_age_counts(counts, mapping)

  expect_identical(expanded$age_group, to_ages$age_groups)
  expect_equal(expanded$value, c(rep(10, 5), rep(20, 5), 30))
  expect_equal(sum(expanded$value), sum(counts$value))
})

test_that("aggregate_age_counts preserves non-age grouping columns", {
  from_ages <- wpp_age_structure_1year(max_age = 5)
  to_ages <- wpp_age_structure_5year(max_age = 5)
  mapping <- AgeGridMapping(from_ages, to_ages, open_ended = "include")

  counts <- data.frame(
    time = rep(c(2020, 2021), each = 12),
    compartment = rep(rep(c("S", "I"), each = 6), times = 2),
    age_group = rep(from_ages$age_groups, times = 4),
    value = 1:24,
    stringsAsFactors = FALSE
  )

  aggregated <- aggregate_age_counts(counts, mapping)

  expect_identical(names(aggregated), c("time", "compartment", "age_group", "value"))
  expect_equal(
    aggregated$value[aggregated$time == 2020 & aggregated$compartment == "S"],
    c(sum(1:5), 6)
  )
  expect_equal(
    aggregated$value[aggregated$time == 2021 & aggregated$compartment == "I"],
    c(sum(19:23), 24)
  )
})

test_that("AgeGridMapping rejects non-nested age grids", {
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
    "exactly covered|maps"
  )
})

test_that("AgeGridMapping requires explicit open-ended policy", {
  from_ages <- wpp_age_structure_1year(max_age = 10)
  to_ages <- wpp_age_structure_5year(max_age = 10)

  expect_error(
    AgeGridMapping(from_ages, to_ages, open_ended = "error"),
    "open_ended = 'include'"
  )

  expect_s3_class(
    AgeGridMapping(from_ages, to_ages, open_ended = "include"),
    "AgeGridMapping"
  )
})

test_that("expansion followed by aggregation returns original nested-bin totals", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  expand_mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")
  aggregate_mapping <- AgeGridMapping(fine_ages, coarse_ages, open_ended = "include")

  counts <- data.frame(
    time = rep(c(0, 1), each = coarse_ages$n_age_groups),
    compartment = rep("S", times = 2 * coarse_ages$n_age_groups),
    age_group = rep(coarse_ages$age_groups, times = 2),
    value = c(50, 100, 30, 60, 110, 40),
    stringsAsFactors = FALSE
  )

  round_trip <- aggregate_age_counts(
    expand_age_counts(counts, expand_mapping),
    aggregate_mapping
  )

  expect_equal(round_trip$value, counts$value)
  expect_identical(round_trip$age_group, counts$age_group)
  expect_identical(round_trip$time, counts$time)
  expect_identical(round_trip$compartment, counts$compartment)
})
