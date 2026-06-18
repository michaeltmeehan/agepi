annual_split_ages <- function(max_age = 3) {
  wpp_age_structure_1year(max_age = max_age)
}

annual_split_state <- function(ages = annual_split_ages(), S = NULL, I = NULL, R = NULL) {
  if (is.null(S)) {
    S <- c(100, 80, 60, 40)
  }
  if (is.null(I)) {
    I <- c(10, 8, 6, 4)
  }
  if (is.null(R)) {
    R <- c(5, 4, 3, 2)
  }
  data.frame(
    compartment = rep(c("S", "I", "R"), each = ages$n_age_groups),
    age_group = rep(ages$age_groups, times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

annual_split_zero_contacts <- function(ages) {
  matrix(0, nrow = ages$n_age_groups, ncol = ages$n_age_groups)
}

annual_split_totals <- function(output) {
  rows <- aggregate(value ~ time + age_group, output, sum)
  rows[order(rows$time, rows$age_group), ]
}

annual_split_process <- function(ages,
                                 fertility = NULL,
                                 mortality = NULL,
                                 migration = NULL,
                                 fertility_exposure_fraction = 1) {
  DemographicProcess(
    age_structure = ages,
    fertility_schedule = fertility,
    fertility_exposure_fraction = fertility_exposure_fraction,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = if (is.null(migration)) "closed" else "migration"
  )
}

annual_split_run_zero_epidemic <- function(process, initial_state, ages, ...) {
  simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 1, 2),
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = annual_split_zero_contacts(ages),
    beta = 0,
    method = "euler",
    demographic_process = process,
    ageing_policy = "annual_cohort",
    time_policy = "step",
    ...
  )
}

test_that("annual cohort compartment step applies births to S and ageing/mortality to every compartment", {
  ages <- annual_split_ages()
  fertility <- FertilitySchedule(
    data.frame(time = 0, age_group = "1", fertility_rate = 0.5),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(time = 0, age_group = ages$age_groups, mortality_rate = log(2)),
    ages
  )
  process <- annual_split_process(ages, fertility = fertility, mortality = mortality)
  state <- state_long_to_vector(annual_split_state(ages), ages, c("S", "I", "R"))

  stepped <- compartment_annual_cohort_demographic_step(
    state_vector = state,
    time = 0,
    model = SIRModel(gamma = 0),
    age_structure = ages,
    demographic_process = process
  )

  expect_equal(
    stepped,
    c(46, 50, 40, 50, 0, 5, 4, 5, 0, 2.5, 2, 2.5)
  )
})

test_that("annual cohort coupling accounts for births deaths ageing and susceptible migration", {
  ages <- wpp_age_structure_1year(max_age = 2)
  fertility <- FertilitySchedule(
    data.frame(time = 0, age_group = "1", fertility_rate = 0.2),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(time = 0, age_group = ages$age_groups, mortality_rate = log(2)),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(time = 0, age_group = ages$age_groups, migration_count = c(2, -1, 3)),
    ages
  )
  process <- annual_split_process(
    ages,
    fertility = fertility,
    mortality = mortality,
    migration = migration
  )
  initial <- annual_split_state(
    ages,
    S = c(10, 20, 30),
    I = c(1, 2, 3),
    R = c(4, 5, 6)
  )

  output <- simulate_deterministic(
    initial_state = initial,
    times = c(0, 1),
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = annual_split_zero_contacts(ages),
    beta = 0,
    method = "euler",
    demographic_process = process,
    ageing_policy = "annual_cohort",
    time_policy = "step"
  )

  final <- output[output$time == 1, ]
  expect_equal(
    final$value,
    c(7.4, 4, 28, 0, 0.5, 2.5, 0, 2, 5.5)
  )
  expect_equal(sum(final$value), sum(initial$value) + 5.4 - 40.5 + 4)
  expect_true(all(final$value >= 0))
})

test_that("annual cohort zero-epidemic births-only coupling matches simulate_demography", {
  ages <- annual_split_ages()
  fertility <- FertilitySchedule(
    data.frame(time = c(0, 1), age_group = "1", fertility_rate = 0.25),
    ages
  )
  process <- annual_split_process(ages, fertility = fertility)
  initial <- annual_split_state(ages, I = rep(0, 4), R = rep(0, 4))

  coupled <- annual_split_run_zero_epidemic(process, initial, ages)
  demography <- simulate_demography(
    process,
    initial_state = initial$value[initial$compartment == "S"],
    times = c(0, 1, 2),
    ageing_policy = "annual_cohort",
    time_policy = "step"
  )

  expect_equal(annual_split_totals(coupled)$value, demography$population)
})

test_that("annual cohort zero-epidemic mortality-only coupling matches simulate_demography", {
  ages <- annual_split_ages()
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(c(0.1, 0.2, 0.3, 0.4), times = 2)
    ),
    ages
  )
  process <- annual_split_process(ages, mortality = mortality)
  initial <- annual_split_state(ages, I = rep(0, 4), R = rep(0, 4))

  coupled <- annual_split_run_zero_epidemic(process, initial, ages)
  demography <- simulate_demography(
    process,
    initial_state = initial$value[initial$compartment == "S"],
    times = c(0, 1, 2),
    ageing_policy = "annual_cohort",
    time_policy = "step"
  )

  expect_equal(annual_split_totals(coupled)$value, demography$population)
})

test_that("annual cohort zero-epidemic ageing-only coupling matches simulate_demography", {
  ages <- annual_split_ages()
  process <- annual_split_process(ages)
  initial <- annual_split_state(ages, I = rep(0, 4), R = rep(0, 4))

  coupled <- annual_split_run_zero_epidemic(process, initial, ages)
  demography <- simulate_demography(
    process,
    initial_state = initial$value[initial$compartment == "S"],
    times = c(0, 1, 2),
    ageing_policy = "annual_cohort"
  )

  expect_equal(annual_split_totals(coupled)$value, demography$population)
})

test_that("annual cohort coupling applies WPP-style age-specific migration after ageing", {
  ages <- annual_split_ages()
  migration_like <- data.frame(
    year = rep(0, ages$n_age_groups),
    age = ages$age_groups,
    net_migration = c(2, 3, -1, 4),
    stringsAsFactors = FALSE
  )
  migration <- standardise_wpp_migration(
    migration_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    migration_col = "net_migration",
    migration_type = "count"
  )
  process <- annual_split_process(ages, migration = migration)
  initial <- annual_split_state(
    ages,
    S = c(100, 80, 60, 40),
    I = c(10, 8, 6, 4),
    R = c(5, 4, 3, 2)
  )

  output <- simulate_deterministic(
    initial_state = initial,
    times = c(0, 1),
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = annual_split_zero_contacts(ages),
    beta = 0,
    method = "euler",
    demographic_process = process,
    ageing_policy = "annual_cohort",
    time_policy = "step"
  )
  final <- output[output$time == 1, ]
  final_by_compartment <- split(final$value, final$compartment)
  initial_total <- sum(initial$value)
  final_total <- sum(final$value)

  expect_equal(final_total - initial_total, sum(migration_like$net_migration))
  expect_equal(final_by_compartment$S, c(2, 103, 79, 104))
  expect_equal(final_by_compartment$I, c(0, 10, 8, 10))
  expect_equal(final_by_compartment$R, c(0, 5, 4, 5))
  expect_identical(final$age_group, rep(ages$age_groups, times = 3))
  expect_true(all(final$value >= 0))
})

test_that("annual cohort zero-epidemic fertility and mortality coupling matches simulate_demography", {
  ages <- annual_split_ages()
  fertility <- FertilitySchedule(
    data.frame(time = c(0, 1), age_group = "1", fertility_rate = 0.25),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(c(0, 1), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 2),
      mortality_rate = rep(c(0.1, 0.2, 0.3, 0.4), times = 2)
    ),
    ages
  )
  process <- annual_split_process(ages, fertility = fertility, mortality = mortality)
  initial <- annual_split_state(ages, I = rep(0, 4), R = rep(0, 4))

  coupled <- annual_split_run_zero_epidemic(process, initial, ages)
  demography <- simulate_demography(
    process,
    initial_state = initial$value[initial$compartment == "S"],
    times = c(0, 1, 2),
    ageing_policy = "annual_cohort",
    time_policy = "step"
  )

  expect_equal(annual_split_totals(coupled)$value, demography$population)
})

test_that("annual cohort multi-compartment no-transition totals match simulate_demography", {
  ages <- annual_split_ages()
  process <- annual_split_process(ages)
  model <- CompartmentModel(compartments = c("S", "A", "B"), infectious_compartments = character())
  initial <- data.frame(
    compartment = rep(c("S", "A", "B"), each = ages$n_age_groups),
    age_group = rep(ages$age_groups, times = 3),
    value = c(100, 80, 60, 40, 10, 8, 6, 4, 5, 4, 3, 2),
    stringsAsFactors = FALSE
  )

  coupled <- simulate_deterministic(
    initial_state = initial,
    times = c(0, 1, 2),
    model = model,
    age_structure = ages,
    contact_matrix = annual_split_zero_contacts(ages),
    beta = 0,
    method = "euler",
    demographic_process = process,
    ageing_policy = "annual_cohort"
  )
  initial_population <- tapply(initial$value, initial$age_group, sum)
  demography <- simulate_demography(
    process,
    initial_state = initial_population[ages$age_groups],
    times = c(0, 1, 2),
    ageing_policy = "annual_cohort"
  )

  expect_equal(annual_split_totals(coupled)$value, demography$population)
})

test_that("annual split with no demography matches existing deSolve SIR and SEIR trajectories", {
  skip_if_not_installed("deSolve")
  ages <- annual_split_ages()
  run <- function(model, state) {
    args <- list(
      initial_state = state,
      times = c(0, 1, 2),
      model = model,
      age_structure = ages,
      contact_matrix = matrix(1, nrow = ages$n_age_groups, ncol = ages$n_age_groups),
      beta = 0.2,
      method = "deSolve"
    )
    list(
      baseline = do.call(simulate_deterministic, args),
      split = do.call(simulate_deterministic, c(args, list(ageing_policy = "annual_cohort")))
    )
  }

  sir <- run(SIRModel(gamma = 0.3), annual_split_state(ages))
  seir <- run(
    SEIRModel(sigma = 0.4, gamma = 0.3),
    data.frame(
      compartment = rep(c("S", "E", "I", "R"), each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = 4),
      value = c(100, 80, 60, 40, 3, 2, 1, 1, 10, 8, 6, 4, 5, 4, 3, 2)
    )
  )

  expect_equal(sir$split, sir$baseline, tolerance = 1e-6)
  expect_equal(seir$split, seir$baseline, tolerance = 1e-6)
})

test_that("annual split with no demography matches existing age-structured cumulative flows", {
  skip_if_not_installed("deSolve")
  ages <- annual_split_ages()
  args <- list(
    initial_state = annual_split_state(ages),
    times = c(0, 1, 2),
    model = SIRModel(gamma = 0.3),
    age_structure = ages,
    contact_matrix = matrix(c(
      2, 1, 0, 0,
      1, 2, 1, 0,
      0, 1, 2, 1,
      0, 0, 1, 2
    ), nrow = ages$n_age_groups, byrow = TRUE),
    beta = 0.2,
    method = "deSolve",
    cumulative_flows = data.frame(name = c("infections", "recoveries"), from = c("S", "I"), to = c("I", "R"))
  )

  baseline <- do.call(simulate_deterministic, args)
  split <- do.call(simulate_deterministic, c(args, list(ageing_policy = "annual_cohort")))

  expect_equal(split$trajectory, baseline$trajectory, tolerance = 1e-6)
  expect_equal(split$cumulative, baseline$cumulative, tolerance = 1e-6)
})

test_that("annual cohort SIR demography sanity run has expected demographic conventions", {
  ages <- annual_split_ages()
  fertility <- FertilitySchedule(data.frame(time = 0, age_group = "1", fertility_rate = 0.1), ages)
  mortality <- MortalitySchedule(
    data.frame(time = 0, age_group = ages$age_groups, mortality_rate = 0.05),
    ages
  )
  process <- annual_split_process(ages, fertility = fertility, mortality = mortality)

  output <- simulate_deterministic(
    initial_state = annual_split_state(ages),
    times = c(0, 1),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = matrix(1, nrow = ages$n_age_groups, ncol = ages$n_age_groups),
    beta = 0.1,
    method = "deSolve",
    demographic_process = process,
    ageing_policy = "annual_cohort",
    time_policy = "step"
  )

  initial_total <- sum(output$value[output$time == 0])
  final_total <- sum(output$value[output$time == 1])
  expect_true(final_total > 0)
  expect_true(final_total < initial_total + 100)
  expect_true(all(output$value >= 0))
  expect_gt(output$value[output$time == 1 & output$compartment == "S" & output$age_group == "0"], 0)
})

test_that("annual cohort output_age_structure aggregates trajectories and cumulative flows", {
  skip_if_not_installed("deSolve")
  ages <- wpp_age_structure_1year(max_age = 5)
  reporting_ages <- wpp_age_structure_5year(max_age = 5)
  process <- annual_split_process(ages)
  initial <- annual_split_state(
    ages,
    S = c(100, 90, 80, 70, 60, 50),
    I = c(10, 9, 8, 7, 6, 5),
    R = c(1, 2, 3, 4, 5, 6)
  )

  internal <- simulate_deterministic(
    initial_state = initial,
    times = c(0, 1),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = matrix(1, nrow = ages$n_age_groups, ncol = ages$n_age_groups),
    beta = 0.05,
    method = "deSolve",
    demographic_process = process,
    ageing_policy = "annual_cohort",
    cumulative_flows = data.frame(name = "infections", from = "S", to = "I")
  )
  aggregated <- simulate_deterministic(
    initial_state = initial,
    times = c(0, 1),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = matrix(1, nrow = ages$n_age_groups, ncol = ages$n_age_groups),
    beta = 0.05,
    method = "deSolve",
    demographic_process = process,
    ageing_policy = "annual_cohort",
    cumulative_flows = data.frame(name = "infections", from = "S", to = "I"),
    output_age_structure = reporting_ages
  )

  mapping <- AgeGridMapping(ages, reporting_ages, open_ended = "include")
  expect_equal(aggregated$trajectory, aggregate_epidemic_trajectory_age_grid(internal$trajectory, mapping))
  expect_equal(aggregated$cumulative, aggregate_cumulative_flows_age_grid(internal$cumulative, mapping))
  expect_equal(
    aggregate(value ~ time + compartment, aggregated$trajectory, sum),
    aggregate(value ~ time + compartment, internal$trajectory, sum),
    tolerance = 1e-6
  )
  expect_equal(
    aggregate(value ~ time + cumulative_name + transition_id + from + to, aggregated$cumulative, sum),
    aggregate(value ~ time + cumulative_name + transition_id + from + to, internal$cumulative, sum),
    tolerance = 1e-6
  )
})

test_that("annual cohort coupling is opt-in and enforces annual complete-grid rules", {
  ages <- AgeStructure(c("0-4", "5+"), c(0, 5), c(4, Inf))
  process <- DemographicProcess(ages)
  state <- data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 3),
    value = c(100, 50, 10, 5, 0, 0),
    stringsAsFactors = FALSE
  )

  default_output <- simulate_deterministic(
    initial_state = state,
    times = c(0, 1),
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = matrix(0, nrow = 2, ncol = 2),
    demographic_process = process,
    method = "euler"
  )
  explicit_output <- simulate_deterministic(
    initial_state = state,
    times = c(0, 1),
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = matrix(0, nrow = 2, ncol = 2),
    demographic_process = process,
    method = "euler",
    ageing_policy = "exponential"
  )

  expect_equal(default_output, explicit_output)
  expect_error(
    simulate_deterministic(
      initial_state = state,
      times = c(0, 1),
      model = SIRModel(gamma = 0),
      age_structure = ages,
      contact_matrix = matrix(0, nrow = 2, ncol = 2),
      demographic_process = process,
      method = "euler",
      ageing_policy = "annual_cohort"
    ),
    "requires 1-year finite age groups"
  )
})
