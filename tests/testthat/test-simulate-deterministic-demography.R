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

coupled_seir_state <- function(
  S = c(100, 50, 25),
  E = c(3, 2, 1),
  I = c(10, 5, 2),
  R = c(20, 10, 4)
) {
  data.frame(
    compartment = rep(c("S", "E", "I", "R"), each = 3),
    age_group = rep(c("0-4", "5-9", "10+"), times = 4),
    value = c(S, E, I, R),
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

test_that("SEIR demographic coupling runs with Euler", {
  output <- coupled_run(
    initial_state = coupled_seir_state(),
    model = SEIRModel(sigma = 0, gamma = 0),
    time_policy = "linear"
  )

  expect_identical(names(output), c("time", "compartment", "age_group", "value"))
  expect_identical(output$time, rep(c(0, 1), each = 12))
  expect_identical(unique(output$compartment), c("S", "E", "I", "R"))
  expect_true(all(is.finite(output$value)))
  expect_true(all(output$value >= 0))
})

test_that("coupled ageing-only simulation moves S I and R independently", {
  output <- coupled_run()

  expect_equal(
    output$value[output$time == 1],
    c(80, 60, 35, 8, 6, 3, 16, 12, 6)
  )
})

test_that("coupled ageing remains compartment-wise under proportional migration", {
  output <- coupled_run(migration_policy = "proportional")

  expect_equal(
    output$value[output$time == 1],
    c(80, 60, 35, 8, 6, 3, 16, 12, 6)
  )
})

test_that("SEIR coupled ageing-only simulation moves S E I and R independently", {
  output <- coupled_run(
    initial_state = coupled_seir_state(),
    model = SEIRModel(sigma = 0, gamma = 0)
  )

  expect_equal(
    output$value[output$time == 1],
    c(80, 60, 35, 2.4, 2.2, 1.4, 8, 6, 3, 16, 12, 6)
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

test_that("coupled fertility births remain S-only under proportional migration", {
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

  output <- coupled_run(process = process, migration_policy = "proportional")

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

test_that("SEIR coupled fertility births enter only youngest S using total S plus E plus I plus R exposure", {
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
    initial_state = coupled_seir_state(),
    process = process,
    model = SEIRModel(sigma = 0, gamma = 0)
  )
  expect_equal(
    output$value[output$time == 1],
    c(106.7, 50, 25, 3, 2, 1, 10, 5, 2, 20, 10, 4)
  )

  exposure_output <- coupled_run(
    initial_state = coupled_seir_state(
      S = c(1, 0, 1),
      E = c(0, 2, 0),
      I = c(0, 5, 0),
      R = c(0, 10, 0)
    ),
    process = process,
    model = SEIRModel(sigma = 0, gamma = 0)
  )
  expect_equal(
    exposure_output$value[exposure_output$time == 1],
    c(2.7, 0, 1, 0, 2, 0, 0, 5, 0, 0, 10, 0)
  )
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

test_that("coupled mortality remains compartment-wise under proportional migration", {
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

  output <- coupled_run(process = process, migration_policy = "proportional")

  expect_equal(
    output$value[output$time == 1],
    c(90, 40, 17.5, 9, 4, 1.4, 18, 8, 2.8)
  )
})

test_that("SEIR coupled mortality removes proportionally from S E I and R", {
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

  output <- coupled_run(
    initial_state = coupled_seir_state(),
    process = process,
    model = SEIRModel(sigma = 0, gamma = 0)
  )

  expect_equal(
    output$value[output$time == 1],
    c(90, 40, 17.5, 2.7, 1.6, 0.7, 9, 4, 1.4, 18, 8, 2.8)
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

test_that("omitted migration_policy matches susceptible migration_policy", {
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

  default_output <- coupled_run(process = process)
  susceptible_output <- coupled_run(process = process, migration_policy = "susceptible")

  expect_equal(default_output, susceptible_output)
})

test_that("proportional migration distributes positive counts across SIR compartment shares", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(13, 13, 31),
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

  output <- coupled_run(process = process, migration_policy = "proportional")

  expect_equal(
    output$value[output$time == 1],
    c(110, 60, 50, 11, 6, 4, 22, 12, 8)
  )
})

test_that("proportional migration distributes positive counts across SEIR compartment shares", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(13.3, 13.4, 32),
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

  output <- coupled_run(
    initial_state = coupled_seir_state(),
    process = process,
    model = SEIRModel(sigma = 0, gamma = 0),
    migration_policy = "proportional"
  )

  expect_equal(
    output$value[output$time == 1],
    c(110, 60, 50, 3.3, 2.4, 2, 11, 6, 4, 22, 12, 8)
  )
})

test_that("proportional migration distributes negative counts across compartment shares", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(-13, -13, -31),
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

  output <- coupled_run(process = process, migration_policy = "proportional")

  expect_equal(
    output$value[output$time == 1],
    c(90, 40, 0, 9, 4, 0, 18, 8, 0)
  )
})

test_that("proportional migration allows zero population ages with zero migration", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(0, 13, 31),
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

  derivative <- compartment_demographic_derivative(
    state_vector = state_long_to_vector(
      coupled_state(S = c(0, 50, 25), I = c(0, 5, 2), R = c(0, 10, 4)),
      ages,
      c("S", "I", "R")
    ),
    time = 0,
    model = SIRModel(gamma = 0),
    age_structure = ages,
    demographic_process = process,
    migration_policy = "proportional"
  )

  expect_equal(
    derivative,
    c(0, 10, 25, 0, 1, 2, 0, 2, 4)
  )
})

test_that("proportional migration errors for zero population ages with non-zero migration", {
  ages <- coupled_test_ages()
  migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(1, 0, 0),
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
    compartment_demographic_derivative(
      state_vector = state_long_to_vector(
        coupled_state(S = c(0, 50, 25), I = c(0, 5, 2), R = c(0, 10, 4)),
        ages,
        c("S", "I", "R")
      ),
      time = 0,
      model = SIRModel(gamma = 0),
      age_structure = ages,
      demographic_process = process,
      migration_policy = "proportional"
    ),
    "cannot allocate non-zero net migration.*total population is zero"
  )
})

test_that("migration_policy error rejects non-zero migration and allows zero migration", {
  ages <- coupled_test_ages()
  zero_migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(0, 0, 0),
      stringsAsFactors = FALSE
    ),
    ages
  )
  non_zero_migration <- MigrationSchedule(
    data.frame(
      time = 0,
      age_group = ages$age_groups,
      migration_count = c(0, 1, 0),
      stringsAsFactors = FALSE
    ),
    ages
  )

  zero_process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    migration_schedule = zero_migration,
    mode = "migration"
  )
  non_zero_process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    migration_schedule = non_zero_migration,
    mode = "migration"
  )

  expect_silent(coupled_run(process = zero_process, migration_policy = "error"))
  expect_error(
    coupled_run(process = non_zero_process, migration_policy = "error"),
    "does not allow non-zero net migration.*ambiguous"
  )
})

test_that("SEIR coupled migration counts enter susceptible compartments only", {
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

  output <- coupled_run(
    initial_state = coupled_seir_state(),
    process = process,
    model = SEIRModel(sigma = 0, gamma = 0)
  )

  expect_equal(
    output$value[output$time == 1],
    c(105, 48, 26, 3, 2, 1, 10, 5, 2, 20, 10, 4)
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

test_that("SEIR coupled migration rates use total population then apply to S only", {
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

  output <- coupled_run(
    initial_state = coupled_seir_state(),
    process = process,
    model = SEIRModel(sigma = 0, gamma = 0)
  )

  expect_equal(
    output$value[output$time == 1],
    c(113.3, 63.4, 21.8, 3, 2, 1, 10, 5, 2, 20, 10, 4)
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

test_that("SEIR coupled infection and demography combine additively without changing total through disease transitions", {
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
    compartment = c("S", "E", "I", "R"),
    age_group = "0+",
    value = c(90, 5, 10, 0)
  )

  output <- simulate_deterministic(
    initial_state = state,
    times = c(0, 1),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = ages,
    contact_matrix = matrix(1, nrow = 1, ncol = 1),
    demographic_process = process
  )
  final_values <- output$value[output$time == 1]

  expect_equal(
    final_values,
    c(72.4285714285714, 11.5714285714286, 8.5, 2)
  )
  expect_equal(sum(final_values), 94.5)
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

test_that("SEIR coupled force of infection depends on I rather than E", {
  ages <- AgeStructure("0+", 0, Inf)
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages)
  )
  state_with_exposed <- data.frame(
    compartment = c("S", "E", "I", "R"),
    age_group = "0+",
    value = c(90, 500, 10, 0)
  )
  state_without_exposed <- data.frame(
    compartment = c("S", "E", "I", "R"),
    age_group = "0+",
    value = c(90, 0, 10, 500)
  )
  state_without_infectious <- data.frame(
    compartment = c("S", "E", "I", "R"),
    age_group = "0+",
    value = c(90, 500, 0, 10)
  )

  run <- function(state) {
    simulate_deterministic(
      initial_state = state,
      times = c(0, 1),
      model = SEIRModel(sigma = 0, gamma = 0),
      age_structure = ages,
      contact_matrix = matrix(1, nrow = 1, ncol = 1),
      demographic_process = process
    )
  }

  final_with_exposed <- run(state_with_exposed)$value[5:8]
  final_without_exposed <- run(state_without_exposed)$value[5:8]
  final_without_infectious <- run(state_without_infectious)$value[5:8]

  expect_equal(final_with_exposed[1], final_without_exposed[1])
  expect_equal(
    final_with_exposed[2] - state_with_exposed$value[2],
    final_without_exposed[2] - state_without_exposed$value[2]
  )
  expect_true(final_with_exposed[2] > state_with_exposed$value[2])
  expect_equal(final_without_infectious[c(1, 2)], c(90, 500))
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

test_that("coupled linear time_policy follows demographic derivative interpolation", {
  ages <- coupled_test_ages()
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
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )

  output <- coupled_run(process = process, times = c(0, 0.5, 1), time_policy = "linear")

  expect_equal(output$value[output$time == 0.5], c(95, 47.5, 23.75, 9.5, 4.75, 1.9, 19, 9.5, 3.8))
  expect_equal(output$value[output$time == 1], c(85.5, 42.75, 21.375, 8.55, 4.275, 1.71, 17.1, 8.55, 3.42))
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

test_that("SIR demography with deSolve and linear time_policy runs", {
  skip_if_not_installed("deSolve")

  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1, 2), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 3),
      mortality_rate = c(
        rep(0.01, ages$n_age_groups),
        rep(0.02, ages$n_age_groups),
        rep(0.03, ages$n_age_groups)
      ),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )
  times <- c(0, 0.25, 0.5, 1)

  euler_output <- coupled_run(process = process, times = times, time_policy = "linear")
  desolve_output <- coupled_run(
    process = process,
    times = times,
    method = "deSolve",
    time_policy = "linear"
  )

  expect_identical(names(desolve_output), names(euler_output))
  expect_identical(desolve_output$time, euler_output$time)
  expect_identical(desolve_output$compartment, euler_output$compartment)
  expect_identical(desolve_output$age_group, euler_output$age_group)
  expect_true(all(is.finite(desolve_output$value)))
  expect_true(all(desolve_output$value >= 0))
})

test_that("Euler and deSolve both honour proportional migration_policy", {
  skip_if_not_installed("deSolve")

  ages <- AgeStructure("0+", 0, Inf)
  migration <- MigrationSchedule(
    data.frame(time = c(0, 1), age_group = "0+", migration_count = 100),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    migration_schedule = migration,
    mode = "migration"
  )
  state <- data.frame(
    compartment = c("S", "I", "R"),
    age_group = "0+",
    value = c(90, 10, 0)
  )
  run <- function(method) {
    simulate_deterministic(
      initial_state = state,
      times = c(0, 1),
      model = SIRModel(gamma = 0),
      age_structure = ages,
      contact_matrix = matrix(0, nrow = 1, ncol = 1),
      demographic_process = process,
      method = method,
      time_policy = "step",
      migration_policy = "proportional"
    )
  }

  euler_output <- run("euler")
  desolve_output <- run("deSolve")

  expect_equal(euler_output$value[euler_output$time == 1], c(180, 20, 0))
  expect_equal(desolve_output$value[desolve_output$time == 1], c(180, 20, 0), tolerance = 1e-6)
})

test_that("invalid migration_policy errors clearly", {
  expect_error(
    coupled_run(migration_policy = "banana"),
    "unsupported migration_policy: banana"
  )
})

test_that("SIR demography with deSolve linear policy respects final schedule boundary", {
  skip_if_not_installed("deSolve")

  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.01, ages$n_age_groups), rep(0.02, ages$n_age_groups)),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )
  times <- c(0, 0.25, 0.5, 1)

  output <- coupled_run(
    process = process,
    times = times,
    method = "deSolve",
    time_policy = "linear"
  )

  expect_identical(output$time, rep(times, each = 9))
  expect_true(all(is.finite(output$value)))
  expect_true(all(output$value >= 0))
})

test_that("SIR demography with deSolve preserves off-grid requested output times under linear policy", {
  skip_if_not_installed("deSolve")

  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1, 2), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 3),
      mortality_rate = c(
        rep(0.01, ages$n_age_groups),
        rep(0.02, ages$n_age_groups),
        rep(0.03, ages$n_age_groups)
      ),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )
  times <- c(0, 0.125, 0.5, 0.875, 1)

  output <- coupled_run(
    process = process,
    times = times,
    method = "deSolve",
    time_policy = "linear"
  )

  expect_identical(output$time, rep(times, each = 9))
  expect_true(all(is.finite(output$value)))
  expect_true(all(output$value >= 0))
})

test_that("SIR demography with deSolve supports step time_policy", {
  skip_if_not_installed("deSolve")

  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1, 2), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 3),
      mortality_rate = c(
        rep(0.01, ages$n_age_groups),
        rep(0.02, ages$n_age_groups),
        rep(0.03, ages$n_age_groups)
      ),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )

  output <- coupled_run(
    process = process,
    times = c(0, 0.5, 1),
    method = "deSolve",
    time_policy = "step"
  )

  expect_identical(output$time, rep(c(0, 0.5, 1), each = 9))
  expect_true(all(is.finite(output$value)))
  expect_true(all(output$value >= 0))
})

test_that("SIR demography with deSolve exact time_policy errors on off-grid times", {
  skip_if_not_installed("deSolve")

  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(0.01, 2 * ages$n_age_groups),
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
    coupled_run(process = process, times = c(0, 0.5, 1), method = "deSolve", time_policy = "exact"),
    "Exact time 0.5 is not available.*no interpolation"
  )
})

test_that("SEIR demography with deSolve runs through the coupled derivative path", {
  skip_if_not_installed("deSolve")

  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1, 2), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 3),
      mortality_rate = c(
        rep(0.01, ages$n_age_groups),
        rep(0.02, ages$n_age_groups),
        rep(0.03, ages$n_age_groups)
      ),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    mortality_schedule = mortality
  )
  times <- c(0, 0.25, 0.5, 1)

  output <- coupled_run(
    initial_state = coupled_seir_state(),
    process = process,
    times = times,
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    method = "deSolve",
    time_policy = "linear"
  )

  expect_identical(names(output), c("time", "compartment", "age_group", "value"))
  expect_identical(output$time, rep(times, each = 12))
  expect_identical(output$compartment, rep(c("S", "S", "S", "E", "E", "E", "I", "I", "I", "R", "R", "R"), times = length(times)))
  expect_true(all(is.finite(output$value)))
  expect_true(all(output$value >= 0))
})

test_that("SEIR demography with deSolve linear policy respects final schedule boundary", {
  skip_if_not_installed("deSolve")

  ages <- coupled_test_ages()
  fertility <- FertilitySchedule(
    data.frame(
      time = c(0, 1),
      age_group = c("10+", "10+"),
      fertility_rate = c(0.025, 0.024),
      stringsAsFactors = FALSE
    ),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(0.004, 0.002, 0.015, 0.004, 0.002, 0.016),
      stringsAsFactors = FALSE
    ),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      migration_count = c(2, -1, 1, 2, -1, 1),
      stringsAsFactors = FALSE
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = coupled_zero_ageing(ages),
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = "migration"
  )
  times <- c(0, 0.25, 0.5, 1)

  output <- coupled_run(
    initial_state = coupled_seir_state(
      S = c(495, 597, 898),
      E = c(3, 2, 1),
      I = c(2, 1, 1),
      R = c(0, 0, 0)
    ),
    process = process,
    times = times,
    model = SEIRModel(sigma = 0.4, gamma = 0.25),
    contact_matrix = matrix(
      c(4, 2, 1, 2, 5, 2, 1, 2, 4),
      nrow = ages$n_age_groups,
      byrow = TRUE
    ),
    beta = 0.08,
    method = "deSolve",
    time_policy = "linear"
  )

  expect_identical(output$time, rep(times, each = 12))
  expect_true(all(is.finite(output$value)))
  expect_true(all(output$value >= 0))
})

test_that("coupled linear policy still errors when requested outputs exceed schedule support", {
  ages <- coupled_test_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = c(rep(0.01, ages$n_age_groups), rep(0.02, ages$n_age_groups)),
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
    coupled_run(process = process, times = c(0, 1, 1.25), time_policy = "linear"),
    "after the final available schedule time 1"
  )
  skip_if_not_installed("deSolve")
  expect_error(
    coupled_run(process = process, times = c(0, 1, 1.25), method = "deSolve", time_policy = "linear"),
    "after the final available schedule time 1"
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
