coupled_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
}

coupled_zero_ageing <- function(ages = coupled_test_ages()) {
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  ageing
}

coupled_state <- function(S = c(100, 50, 25), I = c(10, 5, 2), R = c(20, 10, 4)) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 3),
    age_group = rep(c("0-4", "5-9", "10+"), times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

coupled_zero_contacts <- function(ages = coupled_test_ages()) {
  matrix(0, nrow = ages$n_age_groups, ncol = ages$n_age_groups)
}

coupled_run <- function(
  initial_state = coupled_state(),
  process = DemographicProcess(coupled_test_ages()),
  times = c(0, 1),
  model = SIRModel(0),
  age_structure = coupled_test_ages(),
  contact_matrix = coupled_zero_contacts(age_structure),
  ...
) {
  simulate_deterministic(
    initial_state = initial_state,
    times = times,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    demographic_process = process,
    ...
  )
}

test_that("simulate_deterministic output is unchanged when demographic_process is NULL", {
  ages <- coupled_test_ages()
  default_output <- simulate_deterministic(
    initial_state = coupled_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = coupled_zero_contacts(ages)
  )
  explicit_null_output <- simulate_deterministic(
    initial_state = coupled_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = coupled_zero_contacts(ages),
    demographic_process = NULL,
    time_policy = "step"
  )

  expect_equal(explicit_null_output, default_output)
})

test_that("coupled ageing-only simulation moves S I and R independently", {
  output <- coupled_run()

  expect_equal(
    output$value[output$time == 1],
    c(80, 60, 35, 8, 6, 3, 16, 12, 6)
  )
})

test_that("coupled fertility births enter only the youngest susceptible group", {
  ages <- coupled_test_ages()
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
    ageing_operator = coupled_zero_ageing(ages),
    fertility_schedule = fertility
  )

  output <- coupled_run(process = process)

  expect_equal(
    output$value[output$time == 1],
    c(106.5, 50, 25, 10, 5, 2, 20, 10, 4)
  )
})

test_that("coupled fertility exposure uses total S plus I plus R population", {
  ages <- coupled_test_ages()
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
    ageing_operator = coupled_zero_ageing(ages),
    fertility_schedule = fertility
  )
  output <- coupled_run(
    initial_state = coupled_state(S = c(1, 0, 1), I = c(0, 5, 0), R = c(0, 10, 0)),
    process = process
  )

  expect_equal(output$value[output$time == 1], c(2.5, 0, 1, 0, 5, 0, 0, 10, 0))
})

test_that("coupled mortality removes proportionally from S I and R", {
  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      mortality_rate = c(0.1, 0.2, 0.3),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )

  output <- coupled_run(process = process)

  expect_equal(
    output$value[output$time == 1],
    c(90, 40, 17.5, 9, 4, 1.4, 18, 8, 2.8)
  )
})

test_that("coupled migration counts enter susceptible compartments only", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(5, -2, 1),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    migration_schedule = migration,
    mode = "migration"
  )

  output <- coupled_run(process = process)

  expect_equal(
    output$value[output$time == 1],
    c(105, 48, 26, 10, 5, 2, 20, 10, 4)
  )
})

test_that("coupled migration rates use total population then apply to S only", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_rate = c(0.1, 0.2, -0.1),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    migration_schedule = migration,
    mode = "migration"
  )

  output <- coupled_run(process = process)

  expect_equal(
    output$value[output$time == 1],
    c(113, 63, 21.9, 10, 5, 2, 20, 10, 4)
  )
})

test_that("coupled infection and demography combine additively in one Euler update", {
  ages <- AgeStructure("0+", 0, Inf)
  mortality <- MortalitySchedule(
    data.frame(time = 0, age_group = "0+", mortality_rate = 0.1),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    mortality_schedule = mortality
  )
  state <- data.frame(
    compartment = c("S", "I", "R"),
    age_group = "0+",
    value = c(90, 10, 0)
  )

  output <- simulate_deterministic(
    initial_state = state,
    times = c(0, 1),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = matrix(1, nrow = 1, ncol = 1),
    demographic_process = process
  )

  expect_equal(output$value[output$time == 1], c(72, 16, 2))
})

test_that("coupled force of infection denominator uses current S plus I plus R across steps", {
  ages <- AgeStructure("0+", 0, Inf)
  migration <- MigrationSchedule(
    data.frame(time = c(0, 1), age_group = "0+", migration_count = 100),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    migration_schedule = migration,
    mode = "migration"
  )
  state <- data.frame(
    compartment = c("S", "I", "R"),
    age_group = "0+",
    value = c(90, 10, 0)
  )

  output <- simulate_deterministic(
    initial_state = state,
    times = c(0, 1, 2),
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = matrix(1, nrow = 1, ncol = 1),
    demographic_process = process
  )

  expect_equal(output$value[output$time == 1], c(181, 19, 0))
  expect_equal(output$value[output$time == 2], c(263.805, 36.195, 0))
})

test_that("coupled exact time_policy errors for missing exact demographic schedule times", {
  ages <- coupled_test_ages()
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
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )

  expect_error(
    coupled_run(process = process, times = c(0, 0.5, 1), time_policy = "exact"),
    "Exact time 0.5 is not available.*no interpolation"
  )
})

test_that("coupled step time_policy supports subannual Euler times with annual schedules", {
  ages <- coupled_test_ages()
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
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )

  output <- coupled_run(process = process, times = c(0, 0.5, 1), time_policy = "step")

  expect_equal(output$value[output$time == 0.5], c(95, 47.5, 23.75, 9.5, 4.75, 1.9, 19, 9.5, 3.8))
  expect_equal(output$value[output$time == 1], c(90.25, 45.125, 22.5625, 9.025, 4.5125, 1.805, 18.05, 9.025, 3.61))
})

test_that("coupled demographic schedule coverage is validated before stepping", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      migration_count = c(rep(-1000, ages$n_age_groups), rep(0, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    migration_schedule = migration,
    mode = "migration"
  )

  expect_error(
    coupled_run(process = process, times = c(-0.5, 0), time_policy = "step"),
    "before the first available schedule time 0"
  )
})

test_that("coupled simulation rejects deSolve and ode methods", {
  process <- DemographicProcess(coupled_test_ages())

  expect_error(
    coupled_run(process = process, method = "deSolve"),
    "deSolve.*ode.*not supported"
  )
  expect_error(
    coupled_run(process = process, method = "ode"),
    "deSolve.*ode.*not supported"
  )
})

test_that("coupled simulation rejects demographic process age-structure mismatch", {
  infection_ages <- coupled_test_ages()
  demographic_ages <- AgeStructure(
    age_groups = c("0-9", "10+"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, Inf)
  )
  process <- DemographicProcess(demographic_ages)

  expect_error(
    simulate_deterministic(
      initial_state = coupled_state(),
      times = c(0, 1),
      model = SIRModel(0),
      age_structure = infection_ages,
      contact_matrix = coupled_zero_contacts(infection_ages),
      demographic_process = process
    ),
    "demographic_process must use the same age_structure"
  )
})

test_that("coupled Euler updates that produce negative compartments still error", {
  ages <- coupled_test_ages()
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
    ageing_operator = coupled_zero_ageing(ages),
    migration_schedule = migration,
    mode = "migration"
  )

  expect_error(
    coupled_run(process = process),
    "negative compartment value"
  )
})
