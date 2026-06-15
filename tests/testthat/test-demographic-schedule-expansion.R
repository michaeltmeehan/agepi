test_that("fertility schedule expansion copies rates to nested finer ages and preserves time", {
  coarse_ages <- wpp_age_structure_5year(max_age = 25)
  fine_ages <- wpp_age_structure_1year(max_age = 25)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  fertility <- FertilitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = 2),
      age_group = rep(c("15-19", "20-24"), times = 2),
      fertility_rate = c(0.04, 0.08, 0.05, 0.09),
      rate_source = "mock",
      stringsAsFactors = FALSE
    ),
    coarse_ages
  )
  original_data <- fertility$data

  expanded <- expand_fertility_schedule_age_grid(fertility, mapping)

  expect_s3_class(expanded, "agepi_fertility_schedule")
  expect_identical(expanded$age_structure, fine_ages)
  expect_identical(expanded$rate_convention, "births_per_female_person_year")
  expect_identical(expanded$data$time, rep(c(2020, 2025), each = 10))
  expect_identical(expanded$data$age_group, rep(as.character(15:24), times = 2))
  expect_equal(expanded$data$fertility_rate, c(rep(0.04, 5), rep(0.08, 5), rep(0.05, 5), rep(0.09, 5)))
  expect_identical(expanded$data$rate_source, rep("mock", 20))
  expect_identical(fertility$data, original_data)
})

test_that("mortality schedule expansion copies annual hazards and preserves time", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = coarse_ages$n_age_groups),
      age_group = rep(coarse_ages$age_groups, times = 2),
      mortality_rate = c(0.01, 0.02, 0.08, 0.015, 0.025, 0.09),
      source = "mx",
      stringsAsFactors = FALSE
    ),
    coarse_ages
  )

  expanded <- expand_mortality_schedule_age_grid(mortality, mapping)

  expect_s3_class(expanded, "agepi_mortality_schedule")
  expect_identical(expanded$age_structure, fine_ages)
  expect_identical(expanded$rate_convention, "annual_hazard")
  expect_identical(expanded$data$time, rep(c(2020, 2025), each = fine_ages$n_age_groups))
  expect_equal(expanded$data$mortality_rate, c(rep(0.01, 5), rep(0.02, 5), 0.08, rep(0.015, 5), rep(0.025, 5), 0.09))
  expect_identical(expanded$data$source, rep("mx", 22))
})

test_that("migration count schedule expansion uniformly splits net counts and conserves totals", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  migration <- MigrationSchedule(
    data.frame(
      time = rep(c(2020, 2025), each = coarse_ages$n_age_groups),
      age_group = rep(coarse_ages$age_groups, times = 2),
      migration_count = c(-10, 20, 5, -5, 10, -3),
      stringsAsFactors = FALSE
    ),
    coarse_ages
  )

  expanded <- expand_migration_schedule_age_grid(migration, mapping)

  expect_s3_class(expanded, "agepi_migration_schedule")
  expect_identical(expanded$migration_type, "count")
  expect_equal(expanded$data$migration_count[1:11], c(rep(-2, 5), rep(4, 5), 5))
  expect_equal(
    tapply(expanded$data$migration_count, expanded$data$time, sum),
    tapply(migration$data$migration_count, migration$data$time, sum)
  )
})

test_that("migration rate schedule expansion copies rates to nested finer ages", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  mapping <- AgeGridMapping(coarse_ages, fine_ages, open_ended = "include")

  migration <- MigrationSchedule(
    data.frame(
      time = 2020,
      age_group = coarse_ages$age_groups,
      migration_rate = c(-0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    coarse_ages
  )

  expanded <- expand_migration_schedule_age_grid(migration, mapping)

  expect_s3_class(expanded, "agepi_migration_schedule")
  expect_identical(expanded$migration_type, "rate")
  expect_equal(expanded$data$migration_rate, c(rep(-0.01, 5), rep(0.02, 5), 0.03))
})

test_that("schedule expansion validates mapping compatibility", {
  coarse_ages <- wpp_age_structure_5year(max_age = 10)
  fine_ages <- wpp_age_structure_1year(max_age = 10)
  wrong_source_mapping <- AgeGridMapping(fine_ages, fine_ages, open_ended = "include")

  fertility <- FertilitySchedule(
    data.frame(
      time = 2020,
      age_group = "5-9",
      fertility_rate = 0.01,
      stringsAsFactors = FALSE
    ),
    coarse_ages
  )

  expect_error(
    expand_fertility_schedule_age_grid(fertility, wrong_source_mapping),
    "age_structure must match mapping\\$from_age_groups"
  )
})

test_that("non-nested age mappings error clearly before schedule expansion", {
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
