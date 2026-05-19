test_derivative_age_structure <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
}

test_derivative_process <- function(
  ages = test_derivative_age_structure(),
  fertility = NULL,
  mortality = NULL,
  migration = NULL
) {
  DemographicProcess(
    age_structure = ages,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = if (is.null(migration)) "closed" else "migration"
  )
}

test_that("demographic_derivative computes ageing-only flows", {
  ages <- test_derivative_age_structure()
  process <- test_derivative_process(ages)
  state <- c(100, 50, 25)

  derivative <- demographic_derivative(state, 2020, process)

  expect_equal(derivative, c("0-4" = -20, "5-9" = 10, "10+" = 10))
})

test_that("demographic_derivative keeps final open-ended bin terminal", {
  ages <- test_derivative_age_structure()
  process <- test_derivative_process(ages)
  state <- c(0, 50, 100)

  derivative <- demographic_derivative(state, 2020, process)

  expect_equal(unname(derivative["10+"]), 10)
  expect_identical(process$ageing_operator$departure_rate[3], 0)
})

test_that("demographic_derivative adds fertility births to youngest group only", {
  ages <- test_derivative_age_structure()
  fertility <- FertilitySchedule(
    data.frame(
      time = 2020,
      age_group = c("5-9", "10+"),
      fertility_rate = c(0.1, 0.2),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(ages, fertility = fertility)
  state <- c(100, 50, 25)

  derivative <- demographic_derivative(state, 2020, process)

  expect_equal(derivative, c("0-4" = -10, "5-9" = 10, "10+" = 10))
})

test_that("demographic_derivative subtracts mortality rates by age", {
  ages <- test_derivative_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = 2020,
      age_group = ages$age_groups,
      mortality_rate = c(0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(ages, mortality = mortality)
  state <- c(100, 50, 25)

  derivative <- demographic_derivative(state, 2020, process)

  expect_equal(derivative, c("0-4" = -21, "5-9" = 9, "10+" = 9.25))
})

test_that("demographic_derivative applies migration rates including negative rates", {
  ages <- test_derivative_age_structure()
  migration <- MigrationSchedule(
    data.frame(
      time = 2020,
      age_group = ages$age_groups,
      migration_rate = c(0.1, -0.2, 0),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(ages, migration = migration)
  state <- c(100, 50, 25)

  derivative <- demographic_derivative(state, 2020, process)

  expect_equal(derivative, c("0-4" = -10, "5-9" = 0, "10+" = 10))
})

test_that("demographic_derivative applies migration counts including negative counts", {
  ages <- test_derivative_age_structure()
  migration <- MigrationSchedule(
    data.frame(
      time = 2020,
      age_group = ages$age_groups,
      migration_count = c(5, -10, 2),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(ages, migration = migration)
  state <- c(100, 50, 25)

  derivative <- demographic_derivative(state, 2020, process)

  expect_equal(derivative, c("0-4" = -15, "5-9" = 0, "10+" = 12))
})

test_that("demographic_derivative combines fertility ageing mortality and migration", {
  ages <- test_derivative_age_structure()
  fertility <- FertilitySchedule(
    data.frame(
      time = 2020,
      age_group = "5-9",
      fertility_rate = 0.1,
      stringsAsFactors = FALSE
    ),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(
      time = 2020,
      age_group = ages$age_groups,
      mortality_rate = c(0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(
      time = 2020,
      age_group = ages$age_groups,
      migration_count = c(1, -2, 3),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(
    ages,
    fertility = fertility,
    mortality = mortality,
    migration = migration
  )
  state <- c(100, 50, 25)

  derivative <- demographic_derivative(state, 2020, process)

  expect_equal(derivative, c("0-4" = -15, "5-9" = 7, "10+" = 12.25))
})

test_that("demographic_derivative uses exact-time schedule lookup only", {
  ages <- test_derivative_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.01, ages$n_age_groups), rep(0.02, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(ages, mortality = mortality)

  expect_silent(demographic_derivative(c(100, 50, 25), 2020, process))
  expect_error(
    demographic_derivative(c(100, 50, 25), 2022.5, process),
    "Exact time 2022.5 is not available.*no interpolation"
  )
})

test_that("demographic_derivative exact time_policy preserves exact-time errors", {
  ages <- test_derivative_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.01, ages$n_age_groups), rep(0.02, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(ages, mortality = mortality)

  expect_silent(demographic_derivative(c(100, 50, 25), 2020, process, time_policy = "exact"))
  expect_error(
    demographic_derivative(c(100, 50, 25), 2022.5, process, time_policy = "exact"),
    "Exact time 2022.5 is not available.*no interpolation"
  )
})

test_that("demographic_derivative step time_policy is left-continuous", {
  ages <- test_derivative_age_structure()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.01, ages$n_age_groups), rep(0.02, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    mortality_schedule = mortality
  )

  expect_equal(
    demographic_derivative(c(100, 50, 25), 2020, process, time_policy = "step"),
    c("0-4" = -1, "5-9" = -0.5, "10+" = -0.25)
  )
  expect_equal(
    demographic_derivative(c(100, 50, 25), 2022.5, process, time_policy = "step"),
    c("0-4" = -1, "5-9" = -0.5, "10+" = -0.25)
  )
  expect_equal(
    demographic_derivative(c(100, 50, 25), 2025, process, time_policy = "step"),
    c("0-4" = -2, "5-9" = -1, "10+" = -0.5)
  )
})

test_that("demographic_derivative step time_policy rejects times outside schedule coverage", {
  ages <- test_derivative_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(0.01, 2 * ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- test_derivative_process(ages, mortality = mortality)

  expect_error(
    demographic_derivative(c(100, 50, 25), 2019, process, time_policy = "step"),
    "before the first available schedule time 2020"
  )
  expect_error(
    demographic_derivative(c(100, 50, 25), 2026, process, time_policy = "step"),
    "after the final available schedule time 2025"
  )
})

test_that("demographic_derivative step time_policy is consistent across schedules", {
  ages <- test_derivative_age_structure()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  fertility <- FertilitySchedule(
    data.frame(
      time = c(2020, 2025),
      age_group = c("5-9", "5-9"),
      fertility_rate = c(0.1, 0.2),
      stringsAsFactors = FALSE
    ),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.01, ages$n_age_groups), rep(0.02, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      migration_count = c(c(1, 2, 3), c(4, 5, 6)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = "migration"
  )

  expect_equal(
    demographic_derivative(c(100, 50, 25), 2022.5, process, time_policy = "step"),
    c("0-4" = 5, "5-9" = 1.5, "10+" = 2.75)
  )
  expect_equal(
    demographic_derivative(c(100, 50, 25), 2025, process, time_policy = "step"),
    c("0-4" = 12, "5-9" = 4, "10+" = 5.5)
  )
})

test_that("linear time_policy returns exact endpoint schedule values", {
  ages <- test_derivative_age_structure()
  fertility <- FertilitySchedule(
    data.frame(
      time = c(2020, 2025),
      age_group = c("5-9", "5-9"),
      fertility_rate = c(0.1, 0.2),
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_equal(
    fertility_rates_at(fertility, 2020, ages$age_groups, time_policy = "linear"),
    c(0, 0.1, 0)
  )
  expect_equal(
    fertility_rates_at(fertility, 2025, ages$age_groups, time_policy = "linear"),
    c(0, 0.2, 0)
  )
})

test_that("linear time_policy interpolates fertility mortality and migration by age", {
  ages <- test_derivative_age_structure()
  fertility <- FertilitySchedule(
    data.frame(
      time = c(2020, 2025, 2020, 2025),
      age_group = c("5-9", "5-9", "10+", "10+"),
      fertility_rate = c(0.1, 0.3, 0.2, 0.4),
      stringsAsFactors = FALSE
    ),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(0.01, 0.02, 0.03, 0.03, 0.06, 0.09),
      stringsAsFactors = FALSE
    ),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      migration_count = c(1, 2, 3, 5, 6, 7),
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_equal(
    fertility_rates_at(fertility, 2022.5, ages$age_groups, time_policy = "linear"),
    c(0, 0.2, 0.3)
  )
  expect_equal(
    mortality_rates_at(mortality, 2022.5, ages$age_groups, time_policy = "linear"),
    c(0.02, 0.04, 0.06)
  )
  expect_equal(
    migration_values_at(migration, 2022.5, c(100, 50, 25), ages$age_groups, time_policy = "linear"),
    c(3, 4, 5)
  )
})

test_that("linear time_policy rejects extrapolation and inconsistent interpolation inputs", {
  ages <- test_derivative_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(0.01, 2 * ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_error(
    mortality_rates_at(mortality, 2019, ages$age_groups, time_policy = "linear"),
    "before the first available schedule time 2020"
  )
  expect_error(
    mortality_rates_at(mortality, 2026, ages$age_groups, time_policy = "linear"),
    "after the final available schedule time 2025"
  )

  one_time <- MortalitySchedule(
    data.frame(
      time = 2020,
      age_group = ages$age_groups,
      mortality_rate = rep(0.01, ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    ages
  )
  expect_error(
    mortality_rates_at(one_time, 2020.5, ages$age_groups, time_policy = "linear"),
    "requires at least two schedule times"
  )

  fertility <- FertilitySchedule(
    data.frame(
      time = c(2020, 2025),
      age_group = c("5-9", "10+"),
      fertility_rate = c(0.1, 0.2),
      stringsAsFactors = FALSE
    ),
    ages
  )
  expect_error(
    fertility_rates_at(fertility, 2022.5, ages$age_groups, time_policy = "linear"),
    "consistent age-group coverage"
  )
})

test_that("demographic time_policy rejects unknown policies", {
  expect_error(
    validate_demographic_time_policy("nearest"),
    "unsupported time_policy: nearest"
  )
})

test_that("linear time_policy preserves age-group ordering", {
  ages <- test_derivative_age_structure()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(rev(ages$age_groups), times = 2),
      mortality_rate = c(0.03, 0.02, 0.01, 0.09, 0.06, 0.03),
      stringsAsFactors = FALSE
    ),
    ages
  )

  expect_equal(
    mortality_rates_at(mortality, 2022.5, ages$age_groups, time_policy = "linear"),
    c(0.02, 0.04, 0.06)
  )
})

test_that("demographic_derivative uses linear interpolated rates", {
  ages <- test_derivative_age_structure()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(2020, 2025), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.01, ages$n_age_groups), rep(0.03, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    mortality_schedule = mortality
  )

  expect_equal(
    demographic_derivative(c(100, 50, 25), 2022.5, process, time_policy = "linear"),
    c("0-4" = -2, "5-9" = -1, "10+" = -0.5)
  )
})

test_that("demographic_derivative validates state time and process inputs", {
  ages <- test_derivative_age_structure()
  process <- test_derivative_process(ages)

  expect_error(demographic_derivative(c(1, 2), 2020, process), "length")
  expect_error(demographic_derivative(c(1, -2, 3), 2020, process), "non-negative")
  expect_error(demographic_derivative(c(1, Inf, 3), 2020, process), "finite")
  expect_error(demographic_derivative(c(1, 2, 3), c(2020, 2021), process), "single finite numeric")

  invalid_process <- process
  class(invalid_process) <- "list"
  expect_error(
    demographic_derivative(c(1, 2, 3), 2020, invalid_process),
    "agepi_demographic_process"
  )
})

test_that("demographic_derivative returns numeric named age-ordered output", {
  ages <- test_derivative_age_structure()
  process <- test_derivative_process(ages)

  derivative <- demographic_derivative(c(100, 50, 25), 2020, process)

  expect_type(derivative, "double")
  expect_length(derivative, ages$n_age_groups)
  expect_identical(names(derivative), ages$age_groups)
  expect_equal(unname(derivative), c(-20, 10, 10))
})
