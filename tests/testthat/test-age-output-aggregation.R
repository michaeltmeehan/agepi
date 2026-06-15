output_aggregation_source_ages <- function() {
  wpp_age_structure_1year(max_age = 5)
}

output_aggregation_reporting_ages <- function() {
  wpp_age_structure_5year(max_age = 5)
}

output_aggregation_mapping <- function() {
  AgeGridMapping(
    output_aggregation_source_ages(),
    output_aggregation_reporting_ages(),
    open_ended = "include"
  )
}

test_that("demography trajectory aggregation conserves population totals by time", {
  mapping <- output_aggregation_mapping()
  trajectory <- data.frame(
    time = rep(c(0, 1), each = 6),
    age_group = rep(output_aggregation_source_ages()$age_groups, times = 2),
    population = c(1:6, 11:16),
    stringsAsFactors = FALSE
  )

  aggregated <- aggregate_demography_trajectory_age_grid(trajectory, mapping)

  expect_identical(names(aggregated), c("time", "age_group", "population"))
  expect_identical(aggregated$age_group, rep(output_aggregation_reporting_ages()$age_groups, times = 2))
  expect_equal(aggregated$population, c(sum(1:5), 6, sum(11:15), 16))
  expect_equal(
    tapply(aggregated$population, aggregated$time, sum),
    tapply(trajectory$population, trajectory$time, sum)
  )
})

test_that("epidemic trajectory aggregation conserves compartment totals", {
  mapping <- output_aggregation_mapping()
  trajectory <- data.frame(
    time = rep(c(0, 1), each = 12),
    compartment = rep(rep(c("S", "I"), each = 6), times = 2),
    scenario = rep(c("baseline", "scenario"), each = 12),
    age_group = rep(output_aggregation_source_ages()$age_groups, times = 4),
    value = 1:24,
    stringsAsFactors = FALSE
  )

  aggregated <- aggregate_epidemic_trajectory_age_grid(trajectory, mapping)

  expect_identical(names(aggregated), c("time", "compartment", "scenario", "age_group", "value"))
  expect_equal(
    aggregated$value[aggregated$time == 0 & aggregated$compartment == "S"],
    c(sum(1:5), 6)
  )
  expect_equal(
    aggregate(value ~ time + compartment, aggregated, sum),
    aggregate(value ~ time + compartment, trajectory, sum)
  )
})

test_that("population summary aggregation preserves grouping columns", {
  mapping <- output_aggregation_mapping()
  summary <- data.frame(
    scenario = rep(c("baseline", "shock"), each = 6),
    time = rep(c(0, 1), each = 6),
    age_group = rep(output_aggregation_source_ages()$age_groups, times = 2),
    value = c(1:6, 11:16),
    stringsAsFactors = FALSE
  )

  aggregated <- aggregate_population_summary_age_grid(summary, mapping)

  expect_identical(names(aggregated), c("scenario", "time", "age_group", "value"))
  expect_identical(aggregated$scenario, rep(c("baseline", "shock"), each = 2))
  expect_equal(aggregated$value, c(sum(1:5), 6, sum(11:15), 16))
})

test_that("cumulative-flow aggregation conserves flow values by metadata", {
  mapping <- output_aggregation_mapping()
  cumulative <- data.frame(
    time = rep(c(0, 1), each = 6),
    cumulative_name = "infections",
    transition_id = "S->I",
    from = "S",
    to = "I",
    age_group = rep(output_aggregation_source_ages()$age_groups, times = 2),
    value = c(1:6, 11:16),
    stringsAsFactors = FALSE
  )

  aggregated <- aggregate_cumulative_flows_age_grid(cumulative, mapping)

  expect_identical(
    names(aggregated),
    c("time", "cumulative_name", "transition_id", "from", "to", "age_group", "value")
  )
  expect_equal(aggregated$value, c(sum(1:5), 6, sum(11:15), 16))
  expect_equal(
    aggregate(value ~ time + cumulative_name + transition_id + from + to, aggregated, sum),
    aggregate(value ~ time + cumulative_name + transition_id + from + to, cumulative, sum)
  )
})

test_that("output aggregation wrappers validate required columns", {
  mapping <- output_aggregation_mapping()
  trajectory <- data.frame(
    time = 0,
    age_group = "0",
    value = 1,
    stringsAsFactors = FALSE
  )

  expect_error(
    aggregate_demography_trajectory_age_grid(trajectory, mapping),
    "demography trajectory is missing required column\\(s\\): population"
  )

  expect_error(
    aggregate_cumulative_flows_age_grid(trajectory, mapping),
    "cumulative flows is missing required column"
  )
})

test_that("output aggregation wrappers reject mappings that cannot aggregate", {
  source_ages <- output_aggregation_reporting_ages()
  target_ages <- output_aggregation_source_ages()
  mapping <- AgeGridMapping(source_ages, target_ages, open_ended = "include")
  trajectory <- data.frame(
    time = 0,
    age_group = source_ages$age_groups,
    population = c(10, 5),
    stringsAsFactors = FALSE
  )

  expect_error(
    aggregate_demography_trajectory_age_grid(trajectory, mapping),
    "mapping does not support aggregation"
  )
})
