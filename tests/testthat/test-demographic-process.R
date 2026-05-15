test_process_age_structure <- function() {
  wpp_age_structure_5year(max_age = 20)
}

test_process_mortality <- function(ages, times = 2020) {
  MortalitySchedule(
    data.frame(
      time = rep(times, each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = length(times)),
      mortality_rate = 0.01,
      stringsAsFactors = FALSE
    ),
    ages
  )
}

test_process_fertility <- function(ages, times = 2020) {
  FertilitySchedule(
    data.frame(
      time = rep(times, each = 1),
      age_group = rep("15-19", times = length(times)),
      fertility_rate = 0.05,
      stringsAsFactors = FALSE
    ),
    ages
  )
}

test_process_migration <- function(ages, times = 2020) {
  MigrationSchedule(
    data.frame(
      time = rep(times, each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = length(times)),
      migration_count = 0,
      stringsAsFactors = FALSE
    ),
    ages
  )
}

test_that("DemographicProcess accepts closed process with ageing only", {
  ages <- test_process_age_structure()
  process <- DemographicProcess(age_structure = ages, mode = "closed")

  expect_s3_class(process, "agepi_demographic_process")
  expect_identical(
    names(process),
    c(
      "age_structure",
      "ageing_operator",
      "fertility_schedule",
      "mortality_schedule",
      "migration_schedule",
      "mode",
      "times"
    )
  )
  expect_identical(process$age_structure, ages)
  expect_s3_class(process$ageing_operator, "agepi_ageing_operator")
  expect_null(process$fertility_schedule)
  expect_null(process$mortality_schedule)
  expect_null(process$migration_schedule)
  expect_identical(process$mode, "closed")
  expect_null(process$times)
  expect_silent(validate_demographic_process(process))
})

test_that("DemographicProcess accepts closed process with fertility and mortality", {
  ages <- test_process_age_structure()
  fertility <- test_process_fertility(ages, times = c(2020, 2025))
  mortality <- test_process_mortality(ages, times = c(2020, 2025))

  process <- DemographicProcess(
    age_structure = ages,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    mode = "closed"
  )

  expect_identical(process$fertility_schedule, fertility)
  expect_identical(process$mortality_schedule, mortality)
  expect_identical(process$times, c(2020, 2025))
})

test_that("DemographicProcess records NULL common times for different schedule grids", {
  ages <- test_process_age_structure()
  fertility <- test_process_fertility(ages, times = c(2020, 2025))
  mortality <- test_process_mortality(ages, times = 2020)

  process <- DemographicProcess(
    age_structure = ages,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    mode = "closed"
  )

  expect_null(process$times)
})

test_that("DemographicProcess handles migration mode rules", {
  ages <- test_process_age_structure()
  migration <- test_process_migration(ages)

  expect_error(
    DemographicProcess(
      age_structure = ages,
      migration_schedule = migration,
      mode = "closed"
    ),
    "closed"
  )

  process <- DemographicProcess(
    age_structure = ages,
    migration_schedule = migration,
    mode = "migration"
  )
  expect_identical(process$migration_schedule, migration)
  expect_identical(process$mode, "migration")
})

test_that("DemographicProcess rejects mismatched age structures", {
  ages <- test_process_age_structure()
  other_ages <- wpp_age_structure_5year(max_age = 25)

  expect_error(
    DemographicProcess(
      age_structure = ages,
      ageing_operator = AgeingOperator(other_ages)
    ),
    "ageing_operator"
  )

  expect_error(
    DemographicProcess(
      age_structure = ages,
      mortality_schedule = test_process_mortality(other_ages)
    ),
    "mortality_schedule"
  )
})

test_that("validate_demographic_process rejects corrupted objects", {
  ages <- test_process_age_structure()
  process <- DemographicProcess(age_structure = ages)

  process$mode <- "forced"
  expect_error(validate_demographic_process(process), "mode")

  process <- DemographicProcess(age_structure = ages)
  class(process) <- "list"
  expect_error(validate_demographic_process(process), "agepi_demographic_process")
})
