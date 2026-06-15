simulate_demography_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
}

simulate_demography_process <- function(
  ages = simulate_demography_test_ages(),
  fertility = NULL,
  mortality = NULL,
  migration = NULL,
  fertility_exposure_fraction = 1
) {
  DemographicProcess(
    age_structure = ages,
    fertility_schedule = fertility,
    fertility_exposure_fraction = fertility_exposure_fraction,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = if (is.null(migration)) "closed" else "migration"
  )
}

simulate_demography_annual_ages <- function() {
  wpp_age_structure_1year(max_age = 2)
}

simulate_demography_zero_annual_mortality <- function(ages = simulate_demography_annual_ages(), time = 0) {
  MortalitySchedule(
    data.frame(
      time = time,
      age_group = ages$age_groups,
      mortality_rate = 0,
      stringsAsFactors = FALSE
    ),
    ages
  )
}

test_that("ageing-only Euler simulation moves people into the open-ended final bin", {
  process <- simulate_demography_process()
  output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1))

  final <- output[output$time == 1, ]

  expect_equal(final$population, c(80, 60, 35))
})

test_that("fertility-only Euler simulation adds births to youngest group", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  fertility <- FertilitySchedule(
    data.frame(
      time = 0,
      age_group = "5-9",
      fertility_rate = 0.1,
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    fertility_schedule = fertility
  )

  output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1))

  expect_equal(output$population[output$time == 1], c(105, 50, 25))
})

test_that("mortality-only Euler simulation subtracts mortality by age", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  mortality <- MortalitySchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      mortality_rate = c(0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    mortality_schedule = mortality
  )

  output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 2))

  expect_equal(output$population[output$time == 2], c(98, 48, 23.5))
})

test_that("migration-count Euler simulation adds counts directly over dt", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(5, -10, 2),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    migration_schedule = migration,
    mode = "migration"
  )

  output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 2))

  expect_equal(output$population[output$time == 2], c(110, 30, 29))
})

test_that("migration-rate Euler simulation multiplies rates by population over dt", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_rate = c(0.1, -0.2, 0),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    migration_schedule = migration,
    mode = "migration"
  )

  output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 2))

  expect_equal(output$population[output$time == 2], c(120, 30, 25))
})

test_that("combined Euler simulation matches hand-computable example", {
  ages <- simulate_demography_test_ages()
  fertility <- FertilitySchedule(
    data.frame(
      time = 0,
      age_group = "5-9",
      fertility_rate = 0.1,
      stringsAsFactors = FALSE
    ),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      mortality_rate = c(0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(1, -2, 3),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- simulate_demography_process(
    ages = ages,
    fertility = fertility,
    mortality = mortality,
    migration = migration
  )

  output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1))

  expect_equal(output$population[output$time == 1], c(85, 57, 37.25))
})

test_that("simulate_demography returns tidy output ordered by time then age group", {
  process <- simulate_demography_process()
  output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1, 2))

  expect_identical(names(output), c("time", "age_group", "population"))
  expect_identical(output$time, rep(c(0, 1, 2), each = 3))
  expect_identical(output$age_group, rep(process$age_structure$age_groups, times = 3))
  expect_equal(output$population[output$time == 0], c(100, 50, 25))
  expect_type(output$population, "double")
})

test_that("simulate_demography default ageing policy matches explicit exponential policy", {
  process <- simulate_demography_process()

  default_output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1))
  exponential_output <- simulate_demography(
    process,
    initial_state = c(100, 50, 25),
    times = c(0, 1),
    ageing_policy = "exponential"
  )

  expect_equal(default_output, exponential_output)
})

test_that("annual cohort simulate_demography shifts a cohort from age 0 to age 1", {
  ages <- simulate_demography_annual_ages()
  process <- simulate_demography_process(
    ages = ages,
    mortality = simulate_demography_zero_annual_mortality(ages)
  )

  output <- simulate_demography(
    process,
    initial_state = c(10, 0, 0),
    times = c(0, 1),
    ageing_policy = "annual_cohort"
  )

  expect_equal(output$population[output$time == 1], c(0, 10, 0))
})

test_that("annual cohort simulate_demography conserves population under ageing only", {
  ages <- simulate_demography_annual_ages()
  process <- simulate_demography_process(
    ages = ages,
    mortality = simulate_demography_zero_annual_mortality(ages)
  )

  output <- simulate_demography(
    process,
    initial_state = c(4, 6, 8),
    times = c(0, 1),
    ageing_policy = "annual_cohort"
  )

  expect_equal(sum(output$population[output$time == 1]), 18)
})

test_that("annual cohort simulate_demography handles the open-ended final age group", {
  ages <- simulate_demography_annual_ages()
  process <- simulate_demography_process(
    ages = ages,
    mortality = simulate_demography_zero_annual_mortality(ages)
  )

  output <- simulate_demography(
    process,
    initial_state = c(0, 7, 11),
    times = c(0, 1),
    ageing_policy = "annual_cohort"
  )

  expect_equal(output$population[output$time == 1], c(0, 0, 18))
})

test_that("annual cohort simulate_demography applies mortality as exp minus hazard survival", {
  ages <- AgeStructure("0+", lower_bounds = 0, upper_bounds = Inf)
  mortality <- MortalitySchedule(
    data.frame(
      time = 0,
      age_group = "0+",
      mortality_rate = log(2),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- simulate_demography_process(ages = ages, mortality = mortality)

  output <- simulate_demography(
    process,
    initial_state = 100,
    times = c(0, 1),
    ageing_policy = "annual_cohort"
  )

  expect_equal(output$population[output$time == 1], 50)
})

test_that("annual cohort simulate_demography applies births to the youngest age group", {
  ages <- simulate_demography_annual_ages()
  fertility <- FertilitySchedule(
    data.frame(
      time = 0,
      age_group = "1",
      fertility_rate = 0.25,
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- simulate_demography_process(
    ages = ages,
    fertility = fertility,
    mortality = simulate_demography_zero_annual_mortality(ages)
  )

  output <- simulate_demography(
    process,
    initial_state = c(0, 20, 0),
    times = c(0, 1),
    ageing_policy = "annual_cohort"
  )

  expect_equal(output$population[output$time == 1], c(5, 0, 20))
})

test_that("annual cohort simulate_demography uses fertility exposure fraction", {
  ages <- simulate_demography_annual_ages()
  fertility <- FertilitySchedule(
    data.frame(
      time = 0,
      age_group = "1",
      fertility_rate = 0.25,
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- simulate_demography_process(
    ages = ages,
    fertility = fertility,
    mortality = simulate_demography_zero_annual_mortality(ages),
    fertility_exposure_fraction = 0.5
  )

  output <- simulate_demography(
    process,
    initial_state = c(0, 20, 0),
    times = c(0, 1),
    ageing_policy = "annual_cohort"
  )

  expect_equal(output$population[output$time == 1], c(2.5, 0, 20))
})

test_that("annual cohort simulate_demography errors for non-annual time steps", {
  ages <- simulate_demography_annual_ages()
  process <- simulate_demography_process(
    ages = ages,
    mortality = simulate_demography_zero_annual_mortality(ages)
  )

  expect_error(
    simulate_demography(
      process,
      initial_state = c(10, 0, 0),
      times = c(0, 0.5, 1),
      ageing_policy = "annual_cohort"
    ),
    "annual time steps"
  )
})

test_that("annual cohort simulate_demography errors for non-1-year age grids", {
  process <- simulate_demography_process()

  expect_error(
    simulate_demography(
      process,
      initial_state = c(100, 50, 25),
      times = c(0, 1),
      ageing_policy = "annual_cohort"
    ),
    "1-year finite age groups"
  )
})

test_that("simulate_demography validates process initial state times schedules and method", {
  process <- simulate_demography_process()

  invalid_process <- process
  class(invalid_process) <- "list"
  expect_error(
    simulate_demography(invalid_process, c(100, 50, 25), c(0, 1)),
    "agepi_demographic_process"
  )
  expect_error(
    simulate_demography(process, c(100, 50), c(0, 1)),
    "initial_state length"
  )
  expect_error(
    simulate_demography(process, c(100, -50, 25), c(0, 1)),
    "non-negative"
  )
  expect_error(
    simulate_demography(process, c(100, Inf, 25), c(0, 1)),
    "finite"
  )
  expect_error(
    simulate_demography(process, c(100, 50, 25), c(0, 0)),
    "strictly increasing"
  )
  expect_error(
    simulate_demography(process, c(100, 50, 25), c(0, NA_real_)),
    "missing"
  )
  expect_error(
    simulate_demography(process, c(100, 50, 25), c(0, 1), method = "rk4"),
    "unsupported simulation method"
  )
  expect_error(
    simulate_demography(process, c(100, 50, 25), c(0, 1), method = c("euler", "deSolve")),
    "method must be a non-missing character scalar"
  )
  ages <- simulate_demography_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      mortality_rate = c(0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    ages
  )
  scheduled_process <- simulate_demography_process(ages = ages, mortality = mortality)
  expect_error(
    simulate_demography(scheduled_process, c(100, 50, 25), c(1, 2)),
    "Exact time 1 is not available.*no interpolation"
  )
})

test_that("simulate_demography step time_policy supports subannual Euler times with annual schedules", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.1, ages$n_age_groups), rep(0.2, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    mortality_schedule = mortality
  )

  output <- simulate_demography(
    process,
    initial_state = c(100, 50, 25),
    times = c(0, 0.5, 1),
    time_policy = "step"
  )

  expect_equal(output$population[output$time == 0.5], c(95, 47.5, 23.75))
  expect_equal(output$population[output$time == 1], c(90.25, 45.125, 22.5625))
})

test_that("simulate_demography linear time_policy supports off-grid Euler left endpoints", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.1, ages$n_age_groups), rep(0.3, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    mortality_schedule = mortality
  )

  output <- simulate_demography(
    process,
    initial_state = c(100, 50, 25),
    times = c(0, 0.5, 1),
    time_policy = "linear"
  )

  expect_equal(output$population[output$time == 0.5], c(95, 47.5, 23.75))
  expect_equal(output$population[output$time == 1], c(85.5, 42.75, 21.375))
})

test_that("simulate_demography exact time_policy still requires exact Euler left endpoints", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(0.1, 2 * ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    mortality_schedule = mortality
  )

  expect_error(
    simulate_demography(
      process,
      initial_state = c(100, 50, 25),
      times = c(0, 0.5, 1),
      time_policy = "exact"
    ),
    "Exact time 0.5 is not available.*no interpolation"
  )
})

test_that("simulate_demography step coverage is validated before Euler stepping", {
  ages <- simulate_demography_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(0.1, 2 * ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- simulate_demography_process(ages = ages, mortality = mortality)

  expect_error(
    simulate_demography(process, c(100, 50, 25), c(-0.5, 0), time_policy = "step"),
    "before the first available schedule time 0"
  )
  expect_silent(
    simulate_demography(process, c(100, 50, 25), c(0, 1, 1.5), time_policy = "step")
  )
  expect_error(
    simulate_demography(process, c(100, 50, 25), c(0, 1, 1.5, 2), time_policy = "step"),
    "after the final available schedule time 1"
  )
})

test_that("simulate_demography linear coverage is validated before Euler stepping", {
  ages <- simulate_demography_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(0.1, 2 * ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- simulate_demography_process(ages = ages, mortality = mortality)

  expect_error(
    simulate_demography(process, c(100, 50, 25), c(-0.5, 0), time_policy = "linear"),
    "before the first available schedule time 0"
  )
  expect_error(
    simulate_demography(process, c(100, 50, 25), c(0, 1, 1.5, 2), time_policy = "linear"),
    "after the final available schedule time 1"
  )
})

test_that("Euler steps that produce negative populations are rejected", {
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(0, -100, 0),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    migration_schedule = migration,
    mode = "migration"
  )

  expect_error(
    simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1)),
    "negative population value"
  )
})

test_that("method deSolve and ode run for demographic-only simulation", {
  skip_if_not_installed("deSolve")
  process <- simulate_demography_process()

  desolve_output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1), method = "deSolve")
  ode_output <- simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1), method = "ode")

  expect_identical(names(desolve_output), c("time", "age_group", "population"))
  expect_equal(ode_output, desolve_output)
  expect_true(all(is.finite(desolve_output$population)))
})

test_that("Euler and deSolve agree for zero demographic derivative", {
  skip_if_not_installed("deSolve")
  ages <- simulate_demography_test_ages()
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  process <- DemographicProcess(age_structure = ages, ageing_operator = ageing)

  euler_output <- simulate_demography(
    process,
    initial_state = c(100, 50, 25),
    times = c(0, 1),
    method = "euler"
  )
  desolve_output <- simulate_demography(
    process,
    initial_state = c(100, 50, 25),
    times = c(0, 1),
    method = "deSolve"
  )

  expect_equal(desolve_output, euler_output)
})

test_that("method deSolve errors clearly when deSolve is unavailable", {
  local_mocked_bindings(desolve_is_available = function() FALSE)
  process <- simulate_demography_process()

  expect_error(
    simulate_demography(process, initial_state = c(100, 50, 25), times = c(0, 1), method = "deSolve", time_policy = "step"),
    "requires the optional deSolve package"
  )
})
