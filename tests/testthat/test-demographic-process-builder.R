builder_age_structure <- function() {
  wpp_age_structure_5year(max_age = 20)
}

builder_mortality <- function(ages, times = 2020) {
  MortalitySchedule(
    data.frame(
      time = rep(times, each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = length(times)),
      mortality_rate = rep(0.01, ages$n_age_groups * length(times)),
      stringsAsFactors = FALSE
    ),
    ages
  )
}

builder_fertility <- function(ages, times = 2020) {
  FertilitySchedule(
    data.frame(
      time = times,
      age_group = rep("15-19", length(times)),
      fertility_rate = rep(0.05, length(times)),
      stringsAsFactors = FALSE
    ),
    ages
  )
}

builder_migration <- function(ages, times = 2020) {
  MigrationSchedule(
    data.frame(
      time = rep(times, each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = length(times)),
      migration_count = rep(0, ages$n_age_groups * length(times)),
      stringsAsFactors = FALSE
    ),
    ages
  )
}

test_that("build_demographic_process creates closed process with ageing only", {
  ages <- builder_age_structure()
  process <- build_demographic_process(ages)

  expect_s3_class(process, "agepi_demographic_process")
  expect_s3_class(process$ageing_operator, "agepi_ageing_operator")
  expect_identical(process$age_structure, ages)
  expect_identical(process$mode, "closed")
  expect_null(process$fertility_schedule)
  expect_null(process$mortality_schedule)
  expect_null(process$migration_schedule)
  expect_silent(validate_demographic_process(process))
})

test_that("build_demographic_process stores fertility and mortality schedules", {
  ages <- builder_age_structure()
  fertility <- builder_fertility(ages, times = c(2020, 2025))
  mortality <- builder_mortality(ages, times = c(2020, 2025))

  process <- build_demographic_process(
    age_structure = ages,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    mode = "closed"
  )

  expect_identical(process$fertility_schedule, fertility)
  expect_identical(process$mortality_schedule, mortality)
  expect_identical(process$times, c(2020, 2025))
  expect_silent(validate_demographic_process(process))
})

test_that("build_demographic_process rejects migration in closed mode", {
  ages <- builder_age_structure()
  migration <- builder_migration(ages)

  expect_error(
    build_demographic_process(
      age_structure = ages,
      migration_schedule = migration,
      mode = "closed"
    ),
    "closed"
  )
})

test_that("build_demographic_process accepts migration in migration mode", {
  ages <- builder_age_structure()
  migration <- builder_migration(ages)

  process <- build_demographic_process(
    age_structure = ages,
    migration_schedule = migration,
    mode = "migration"
  )

  expect_s3_class(process, "agepi_demographic_process")
  expect_identical(process$migration_schedule, migration)
  expect_identical(process$mode, "migration")
  expect_silent(validate_demographic_process(process))
})

test_that("build_demographic_process rejects mismatched schedule age structures", {
  ages <- builder_age_structure()
  other_ages <- wpp_age_structure_5year(max_age = 25)
  mortality <- builder_mortality(other_ages)

  expect_error(
    build_demographic_process(
      age_structure = ages,
      mortality_schedule = mortality
    ),
    "mortality_schedule"
  )
})

test_that("build_demographic_process matches direct construction in key fields", {
  ages <- builder_age_structure()
  fertility <- builder_fertility(ages)
  mortality <- builder_mortality(ages)
  migration <- builder_migration(ages)

  built <- build_demographic_process(
    age_structure = ages,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = "migration"
  )
  direct <- DemographicProcess(
    age_structure = ages,
    ageing_operator = AgeingOperator(ages),
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = "migration"
  )

  expect_identical(built$age_structure$age_groups, direct$age_structure$age_groups)
  expect_identical(built$mode, direct$mode)
  expect_identical(built$fertility_schedule, direct$fertility_schedule)
  expect_identical(built$mortality_schedule, direct$mortality_schedule)
  expect_identical(built$migration_schedule, direct$migration_schedule)
  expect_identical(built$ageing_operator$departure_rate, direct$ageing_operator$departure_rate)
  expect_identical(built$ageing_operator$destination_index, direct$ageing_operator$destination_index)
})
