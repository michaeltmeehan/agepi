test_that("standardise_wpp_mortality creates sorted MortalitySchedule", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality_like <- data.frame(
    year = rep(c(2025, 2020), each = ages$n_age_groups),
    age = rep(rev(ages$age_groups), times = 2),
    mx = rep(0.01, 2 * ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  mortality <- standardise_wpp_mortality(
    mortality_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    mortality_col = "mx"
  )

  expect_s3_class(mortality, "agepi_mortality_schedule")
  expect_identical(mortality$rate_convention, "annual_hazard")
  expect_identical(mortality$data$time[1:ages$n_age_groups], rep(2020, ages$n_age_groups))
  expect_identical(mortality$data$age_group[1:ages$n_age_groups], ages$age_groups)
})

test_that("standardise_wpp_mortality validates input and delegates rate checks", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    mx = rep(0.01, ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  missing_mx <- mortality_like
  missing_mx$mx <- NULL
  expect_error(
    standardise_wpp_mortality(missing_mx, ages, "year", "age", "mx"),
    "missing required column"
  )

  invalid_age <- mortality_like
  invalid_age$age[1] <- "not-an-age"
  expect_error(
    standardise_wpp_mortality(invalid_age, ages, "year", "age", "mx"),
    "Unsupported WPP age label"
  )

  negative_mortality <- mortality_like
  negative_mortality$mx[1] <- -0.01
  expect_error(
    standardise_wpp_mortality(negative_mortality, ages, "year", "age", "mx"),
    "negative"
  )
})

test_that("standardise_wpp_fertility creates FertilitySchedule with partial age coverage", {
  ages <- wpp_age_structure_5year()
  fertility_like <- data.frame(
    year = c(2025, 2020, 2020),
    age = c("25-29", "20-24", "15-19"),
    asfr = c(0.1, 0.08, 0.04),
    stringsAsFactors = FALSE
  )

  fertility <- standardise_wpp_fertility(
    fertility_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    fertility_col = "asfr"
  )

  expect_s3_class(fertility, "agepi_fertility_schedule")
  expect_identical(fertility$rate_convention, "births_per_female_person_year")
  expect_identical(fertility$data$time, c(2020, 2020, 2025))
  expect_identical(fertility$data$age_group, c("15-19", "20-24", "25-29"))
})

test_that("standardise_wpp_fertility validates inputs and non-negative rates", {
  ages <- wpp_age_structure_5year()
  fertility_like <- data.frame(
    year = 2020,
    age = "20-24",
    asfr = 0.08,
    stringsAsFactors = FALSE
  )

  missing_asfr <- fertility_like
  missing_asfr$asfr <- NULL
  expect_error(
    standardise_wpp_fertility(missing_asfr, ages, "year", "age", "asfr"),
    "missing required column"
  )

  negative_fertility <- fertility_like
  negative_fertility$asfr <- -0.01
  expect_error(
    standardise_wpp_fertility(negative_fertility, ages, "year", "age", "asfr"),
    "negative"
  )
})

test_that("standardise_wpp_migration creates count and rate schedules", {
  ages <- wpp_age_structure_5year(max_age = 10)
  migration_like <- data.frame(
    year = rep(c(2025, 2020), each = ages$n_age_groups),
    age = rep(rev(ages$age_groups), times = 2),
    net = c(-10, 0, 3, 5, -2, 1),
    stringsAsFactors = FALSE
  )

  migration_count <- standardise_wpp_migration(
    migration_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    migration_col = "net",
    migration_type = "count"
  )
  expect_s3_class(migration_count, "agepi_migration_schedule")
  expect_identical(migration_count$migration_type, "count")
  expect_true(any(migration_count$data$migration_count < 0))
  expect_identical(migration_count$data$age_group[1:ages$n_age_groups], ages$age_groups)

  migration_rate <- standardise_wpp_migration(
    transform(migration_like, net = net / 1000),
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    migration_col = "net",
    migration_type = "rate"
  )
  expect_s3_class(migration_rate, "agepi_migration_schedule")
  expect_identical(migration_rate$migration_type, "rate")
  expect_true(any(migration_rate$data$migration_rate < 0))
})

test_that("standardise_wpp_migration validates migration_type", {
  ages <- wpp_age_structure_5year(max_age = 10)
  migration_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    net = rep(0, ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  expect_error(
    standardise_wpp_migration(
      migration_like,
      age_structure = ages,
      time_col = "year",
      age_col = "age",
      migration_col = "net",
      migration_type = "invalid"
    ),
    "'arg' should be one of"
  )
})

test_that("WPP age labels map to one-year age structures", {
  ages <- wpp_age_structure_1year(max_age = 3)
  mortality_like <- data.frame(
    year = 2020,
    age = c(0, 1, 2, 3),
    mx = rep(0.01, ages$n_age_groups)
  )

  mortality <- standardise_wpp_mortality(mortality_like, ages, "year", "age", "mx")

  expect_identical(mortality$data$age_group, c("0", "1", "2", "3+"))
})

test_that("WPP age labels map to five-year and lower-bound age structures", {
  ages <- wpp_age_structure_5year(max_age = 10)
  five_year_like <- data.frame(
    year = 2020,
    age = c("0-4", "5-9", "10+"),
    mx = rep(0.01, ages$n_age_groups)
  )
  lower_bound_like <- data.frame(
    year = 2020,
    age = c("0", "5", "10+"),
    mx = rep(0.01, ages$n_age_groups)
  )

  five_year <- standardise_wpp_mortality(five_year_like, ages, "year", "age", "mx")
  lower_bound <- standardise_wpp_mortality(lower_bound_like, ages, "year", "age", "mx")

  expect_identical(five_year$data$age_group, ages$age_groups)
  expect_identical(lower_bound$data$age_group, ages$age_groups)
})

test_that("wpp2024 remains optional for demographic input standardisation", {
  testthat::skip_if_not_installed("wpp2024")
  expect_true(requireNamespace("wpp2024", quietly = TRUE))
})
