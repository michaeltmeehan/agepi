residual_test_age_structure <- function() {
  AgeStructure(
    age_groups = c("0-4", "5+"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, Inf)
  )
}

residual_test_demography <- function(times, population, ages = residual_test_age_structure()) {
  Demography(
    data.frame(
      time = rep(times, each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = length(times)),
      population = population,
      stringsAsFactors = FALSE
    ),
    ages
  )
}

test_that("implied_demographic_residual is zero for a perfect Euler match", {
  ages <- residual_test_age_structure()
  process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = c(0, 1),
    population = c(100, 200, 80, 220),
    ages = ages
  )

  residual <- implied_demographic_residual(observed, process)

  expect_equal(residual$residual_count, c(0, 0))
  expect_equal(residual$residual_rate, c(0, 0))
  expect_equal(residual$predicted_end, c(80, 220))
})

test_that("implied_demographic_residual computes hand-checkable residuals", {
  ages <- residual_test_age_structure()
  process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = c(0, 2),
    population = c(100, 200, 70, 250),
    ages = ages
  )

  residual <- implied_demographic_residual(observed, process)

  expect_equal(residual$model_derivative, c(-20, 20))
  expect_equal(residual$predicted_end, c(60, 240))
  expect_equal(residual$residual_count, c(10, 10))
})

test_that("implied_demographic_residual computes rates and uses NA for zero starts", {
  ages <- residual_test_age_structure()
  process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = c(0, 1),
    population = c(0, 100, 5, 120),
    ages = ages
  )

  residual <- implied_demographic_residual(observed, process, type = "rate")

  expect_false("residual_count" %in% names(residual))
  expect_true("residual_rate" %in% names(residual))
  expect_true(is.na(residual$residual_rate[1]))
  expect_equal(residual$residual_rate[2], 0.2)
})

test_that("implied_demographic_residual handles multiple intervals in order", {
  ages <- residual_test_age_structure()
  process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = c(0, 1, 3),
    population = c(100, 200, 80, 220, 55, 260),
    ages = ages
  )

  residual <- implied_demographic_residual(observed, process, type = "count")

  expect_equal(nrow(residual), 4)
  expect_identical(residual$time_start, c(0, 0, 1, 1))
  expect_identical(residual$time_end, c(1, 1, 3, 3))
  expect_identical(residual$dt, c(1, 1, 2, 2))
  expect_identical(residual$age_group, rep(ages$age_groups, times = 2))
  expect_true("residual_count" %in% names(residual))
  expect_false("residual_rate" %in% names(residual))
})

test_that("implied_demographic_residual errors for age-group mismatches", {
  ages <- residual_test_age_structure()
  other_ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
  process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = c(0, 1),
    population = c(100, 50, 200, 90, 60, 220),
    ages = other_ages
  )

  expect_error(
    implied_demographic_residual(observed, process),
    "observed must use the same age_structure"
  )
})

test_that("implied_demographic_residual requires at least two observed times", {
  ages <- residual_test_age_structure()
  process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = 0,
    population = c(100, 200),
    ages = ages
  )

  expect_error(
    implied_demographic_residual(observed, process),
    "at least two exact time points"
  )
})

test_that("implied_demographic_residual uses observed exact times only", {
  ages <- residual_test_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      mortality_rate = c(0.01, 0.02),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- build_demographic_process(ages, mortality_schedule = mortality)
  observed <- residual_test_demography(
    times = c(1, 2),
    population = c(100, 200, 90, 220),
    ages = ages
  )

  expect_error(
    implied_demographic_residual(observed, process),
    "Exact time 1 is not available.*no interpolation"
  )
})

test_that("implied_demographic_residual is relative to processes that include migration", {
  ages <- residual_test_age_structure()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(5, -10),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- build_demographic_process(
    age_structure = ages,
    migration_schedule = migration,
    mode = "migration"
  )
  observed <- residual_test_demography(
    times = c(0, 1),
    population = c(100, 200, 85, 205),
    ages = ages
  )

  residual <- implied_demographic_residual(observed, process)

  expect_equal(residual$predicted_end, c(85, 210))
  expect_equal(residual$residual_count, c(0, -5))
})

test_that("residual_to_migration_schedule converts interval counts to per-time flows", {
  ages <- residual_test_age_structure()
  residual <- data.frame(
    time_start = c(0, 0, 2, 2),
    time_end = c(2, 2, 4, 4),
    dt = c(2, 2, 2, 2),
    age_group = rep(ages$age_groups, times = 2),
    residual_count = c(10, -4, 6, -2),
    residual_rate = c(0.05, -0.01, 0.03, -0.005),
    stringsAsFactors = FALSE
  )

  migration <- residual_to_migration_schedule(residual, ages, use = "count")

  expect_s3_class(migration, "agepi_migration_schedule")
  expect_identical(migration$migration_type, "count")
  expect_identical(migration$data$time, c(0, 0, 2, 2))
  expect_identical(migration$data$age_group, rep(ages$age_groups, times = 2))
  expect_equal(migration$data$migration_count, c(5, -2, 3, -1))
})

test_that("residual_to_migration_schedule converts residual rates directly", {
  ages <- residual_test_age_structure()
  residual <- data.frame(
    time_start = c(0, 0),
    time_end = c(1, 1),
    dt = c(1, 1),
    age_group = ages$age_groups,
    residual_count = c(10, -4),
    residual_rate = c(0.1, -0.02),
    stringsAsFactors = FALSE
  )

  migration <- residual_to_migration_schedule(residual, ages, use = "rate")

  expect_s3_class(migration, "agepi_migration_schedule")
  expect_identical(migration$migration_type, "rate")
  expect_identical(migration$data$time, c(0, 0))
  expect_equal(migration$data$migration_rate, c(0.1, -0.02))
})

test_that("residual_to_migration_schedule errors for age structure mismatches", {
  ages <- residual_test_age_structure()
  residual <- data.frame(
    time_start = c(0, 0),
    time_end = c(1, 1),
    dt = c(1, 1),
    age_group = c("0-4", "10+"),
    residual_count = c(1, 2),
    residual_rate = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  expect_error(
    residual_to_migration_schedule(residual, ages),
    "not in age_structure"
  )
})

test_that("residual_to_migration_schedule errors for missing columns", {
  ages <- residual_test_age_structure()
  residual <- data.frame(
    time_start = c(0, 0),
    time_end = c(1, 1),
    dt = c(1, 1),
    age_group = ages$age_groups,
    residual_count = c(1, 2),
    stringsAsFactors = FALSE
  )

  expect_error(
    residual_to_migration_schedule(residual, ages),
    "missing required column"
  )
})

test_that("residual_to_migration_schedule rejects NA rates for rate conversion", {
  ages <- residual_test_age_structure()
  residual <- data.frame(
    time_start = c(0, 0),
    time_end = c(1, 1),
    dt = c(1, 1),
    age_group = ages$age_groups,
    residual_count = c(1, 2),
    residual_rate = c(NA_real_, 0.02),
    stringsAsFactors = FALSE
  )

  expect_error(
    residual_to_migration_schedule(residual, ages, use = "rate"),
    "residual_rate cannot contain NA"
  )
  expect_silent(residual_to_migration_schedule(residual, ages, use = "count"))
})

test_that("residual_to_migration_schedule checks dt against interval length", {
  ages <- residual_test_age_structure()
  residual <- data.frame(
    time_start = c(0, 0),
    time_end = c(1, 1),
    dt = c(2, 2),
    age_group = ages$age_groups,
    residual_count = c(1, 2),
    residual_rate = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  expect_error(
    residual_to_migration_schedule(residual, ages, use = "count"),
    "residual dt must equal time_end - time_start"
  )
})

test_that("residual_to_migration_schedule count conversion closes one Euler interval", {
  ages <- residual_test_age_structure()
  base_process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = c(0, 2),
    population = c(100, 200, 70, 250),
    ages = ages
  )

  residual <- implied_demographic_residual(observed, base_process)
  migration <- residual_to_migration_schedule(residual, ages, use = "count")
  process_with_residual_migration <- build_demographic_process(
    age_structure = ages,
    migration_schedule = migration,
    mode = "migration"
  )

  simulated <- simulate_demography(
    process = process_with_residual_migration,
    initial_state = c(100, 200),
    times = c(0, 2)
  )

  expect_equal(simulated$population[simulated$time == 2], c(70, 250))
})

test_that("residual-derived migration closes the intended interval with step lookup", {
  ages <- residual_test_age_structure()
  base_process <- build_demographic_process(ages)
  observed <- residual_test_demography(
    times = c(0, 2),
    population = c(100, 200, 70, 250),
    ages = ages
  )

  residual <- implied_demographic_residual(observed, base_process)
  migration <- residual_to_migration_schedule(residual, ages, use = "count")
  process_with_residual_migration <- build_demographic_process(
    age_structure = ages,
    migration_schedule = migration,
    mode = "migration"
  )

  simulated <- simulate_demography(
    process = process_with_residual_migration,
    initial_state = c(100, 200),
    times = c(0, 2),
    time_policy = "step"
  )

  expect_equal(simulated$population[simulated$time == 2], c(70, 250))
})
