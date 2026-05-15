test_rate_age_structure <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9", "10-14", "15-19", "20-24"),
    lower_bounds = c(0, 5, 10, 15, 20),
    upper_bounds = c(4, 9, 14, 19, Inf)
  )
}

test_that("FertilitySchedule accepts valid partial age coverage", {
  ages <- test_rate_age_structure()
  fertility <- FertilitySchedule(
    data.frame(
      time = c(2025, 2020, 2020),
      age_group = c("20-24", "20-24", "15-19"),
      fertility_rate = c(0.1, 0.08, 0.04),
      rate_source = "mock",
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_s3_class(fertility, "agepi_fertility_schedule")
  expect_identical(fertility$times, c(2020, 2025))
  expect_identical(fertility$n_times, 2L)
  expect_identical(fertility$n_age_groups, ages$n_age_groups)
  expect_identical(fertility$rate_convention, "births_per_female_person_year")
  expect_identical(fertility$data$time, c(2020, 2020, 2025))
  expect_identical(fertility$data$age_group, c("15-19", "20-24", "20-24"))
  expect_identical(fertility$data$rate_source, rep("mock", 3))
  expect_silent(validate_fertility_schedule(fertility))
})

test_that("FertilitySchedule validates required columns and values", {
  ages <- test_rate_age_structure()
  valid <- data.frame(
    time = 2020,
    age_group = "15-19",
    fertility_rate = 0.04,
    stringsAsFactors = FALSE
  )

  missing_rate <- valid
  missing_rate$fertility_rate <- NULL
  expect_error(FertilitySchedule(missing_rate, ages), "missing required column")

  negative_rate <- valid
  negative_rate$fertility_rate <- -0.1
  expect_error(FertilitySchedule(negative_rate, ages), "negative")

  duplicate_rows <- rbind(valid, valid)
  expect_error(FertilitySchedule(duplicate_rows, ages), "duplicate time-age_group")

  invalid_age <- valid
  invalid_age$age_group <- "25-29"
  expect_error(FertilitySchedule(invalid_age, ages), "not in age_structure")
})

test_that("MortalitySchedule accepts valid full age coverage and sorts rows", {
  ages <- test_rate_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2025, 2020), each = ages$n_age_groups),
      age_group = rep(rev(ages$age_groups), times = 2),
      mortality_rate = rep(0.01, 2 * ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_s3_class(mortality, "agepi_mortality_schedule")
  expect_identical(mortality$times, c(2020, 2025))
  expect_identical(mortality$rate_convention, "annual_hazard")
  expect_identical(mortality$data$time[1:ages$n_age_groups], rep(2020, ages$n_age_groups))
  expect_identical(mortality$data$age_group[1:ages$n_age_groups], ages$age_groups)
  expect_silent(validate_mortality_schedule(mortality))
})

test_that("MortalitySchedule rejects incomplete duplicate and invalid rates", {
  ages <- test_rate_age_structure()
  valid <- data.frame(
    time = 2020,
    age_group = ages$age_groups,
    mortality_rate = rep(0.01, ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  incomplete <- valid[-1, ]
  expect_error(MortalitySchedule(incomplete, ages), "missing age_group")

  duplicate_rows <- rbind(valid, valid[1, ])
  expect_error(MortalitySchedule(duplicate_rows, ages), "duplicate time-age_group")

  negative_rate <- valid
  negative_rate$mortality_rate[1] <- -0.01
  expect_error(MortalitySchedule(negative_rate, ages), "negative")

  infinite_rate <- valid
  infinite_rate$mortality_rate[1] <- Inf
  expect_error(MortalitySchedule(infinite_rate, ages), "finite")
})

test_that("MigrationSchedule accepts rate schedules and allows negative values", {
  ages <- test_rate_age_structure()
  migration <- MigrationSchedule(
    data.frame(
      time = 2020,
      age_group = ages$age_groups,
      migration_rate = c(-0.01, 0, 0.01, -0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_s3_class(migration, "agepi_migration_schedule")
  expect_identical(migration$migration_type, "rate")
  expect_true(any(migration$data$migration_rate < 0))
  expect_silent(validate_migration_schedule(migration))
})

test_that("MigrationSchedule accepts count schedules and sorts rows", {
  ages <- test_rate_age_structure()
  migration <- MigrationSchedule(
    data.frame(
      time = rep(c(2025, 2020), each = ages$n_age_groups),
      age_group = rep(rev(ages$age_groups), times = 2),
      migration_count = seq_len(2 * ages$n_age_groups) - 6,
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_identical(migration$migration_type, "count")
  expect_identical(migration$data$time[1:ages$n_age_groups], rep(2020, ages$n_age_groups))
  expect_identical(migration$data$age_group[1:ages$n_age_groups], ages$age_groups)
})

test_that("MigrationSchedule validates type and coverage", {
  ages <- test_rate_age_structure()
  valid <- data.frame(
    time = 2020,
    age_group = ages$age_groups,
    migration_rate = rep(0, ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  both <- valid
  both$migration_count <- 0
  expect_error(MigrationSchedule(both, ages), "exactly one")

  neither <- valid
  neither$migration_rate <- NULL
  expect_error(MigrationSchedule(neither, ages), "exactly one")

  incomplete <- valid[-1, ]
  expect_error(MigrationSchedule(incomplete, ages), "missing age_group")

  duplicate_rows <- rbind(valid, valid[1, ])
  expect_error(MigrationSchedule(duplicate_rows, ages), "duplicate time-age_group")

  nonfinite <- valid
  nonfinite$migration_rate[1] <- NA_real_
  expect_error(MigrationSchedule(nonfinite, ages), "finite")
})
