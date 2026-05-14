simulate_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

simulate_test_state <- function(S = c(90, 180), I = c(10, 20), R = c(0, 0)) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

simulate_test_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

simulate_test_run <- function(initial_state = simulate_test_state(), times = c(0, 0.1)) {
  simulate_deterministic(
    initial_state = initial_state,
    times = times,
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts()
  )
}

test_that("simulate_deterministic returns tidy long-form output columns", {
  output <- simulate_test_run()

  expect_identical(names(output), c("time", "compartment", "age_group", "value"))
  expect_type(output$value, "double")
})

test_that("simulate_deterministic includes the initial state at the first time point", {
  output <- simulate_test_run()
  initial_rows <- output[output$time == 0, ]

  expect_equal(initial_rows$value, c(90, 180, 10, 20, 0, 0))
  expect_identical(initial_rows$compartment, c("S", "S", "I", "I", "R", "R"))
  expect_identical(initial_rows$age_group, c("0-4", "5-9", "0-4", "5-9", "0-4", "5-9"))
})

test_that("simulate_deterministic orders time outermost, compartment next, age group innermost", {
  times <- c(-1, 0.1, 0.25)
  output <- simulate_test_run(times = times)

  expect_identical(output$time, rep(times, each = 6))
  expect_identical(output$compartment, rep(c("S", "S", "I", "I", "R", "R"), times = 3))
  expect_identical(output$age_group, rep(c("0-4", "5-9", "0-4", "5-9", "0-4", "5-9"), times = 3))
})

test_that("one-step Euler update matches hand-calculated SIR example", {
  output <- simulate_test_run(times = c(0, 0.1))
  final_rows <- output[output$time == 0.1, ]

  expect_equal(final_rows$value, c(87.3, 167.4, 12.5, 32.2, 0.2, 0.4))
})

test_that("method euler remains the default", {
  default_output <- simulate_test_run(times = c(0, 0.1, 0.2))
  explicit_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "euler"
  )

  expect_equal(default_output, explicit_output)
})

test_that("multi-step simulation runs end-to-end", {
  output <- simulate_test_run(times = seq(0, 0.3, by = 0.1))

  expect_equal(nrow(output), 24)
  expect_equal(unique(output$time), c(0, 0.1, 0.2, 0.3))
  expect_true(all(is.finite(output$value)))
})

test_that("non-uniform time steps are handled correctly", {
  output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.3),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts()
  )

  step_two <- output[output$time == 0.3, ]
  expect_equal(
    step_two$value,
    c(80.12394, 133.28388, 19.17606, 65.02812, 0.7, 1.688),
    tolerance = 1e-10
  )
})

test_that("total population is conserved within each age group for closed SIR dynamics", {
  output <- simulate_test_run(times = seq(0, 0.3, by = 0.1))

  totals <- aggregate(value ~ time + age_group, output, sum)
  expect_equal(totals$value[totals$age_group == "0-4"], rep(100, 4))
  expect_equal(totals$value[totals$age_group == "5-9"], rep(200, 4))

  overall_totals <- aggregate(value ~ time, output, sum)
  expect_equal(overall_totals$value, rep(300, 4))
})

test_that("vector and long-form initial states produce equivalent output", {
  ages <- simulate_test_ages()
  state_long <- simulate_test_state()
  state_vector <- state_long_to_vector(state_long, ages, c("S", "I", "R"))

  long_output <- simulate_test_run(initial_state = state_long, times = c(0, 0.1, 0.2))
  vector_output <- simulate_deterministic(
    initial_state = state_vector,
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = simulate_test_contacts()
  )

  expect_equal(vector_output, long_output)
})

test_that("simulate_deterministic ignores state vector names intentionally", {
  named_state <- c(
    not_s_0_4 = 90, not_s_5_9 = 180,
    not_i_0_4 = 10, not_i_5_9 = 20,
    not_r_0_4 = 0, not_r_5_9 = 0
  )

  output <- simulate_deterministic(
    initial_state = named_state,
    times = c(0, 0.1),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts()
  )

  expect_equal(output$value[output$time == 0.1], c(87.3, 167.4, 12.5, 32.2, 0.2, 0.4))
})

test_that("simulate_deterministic ignores state vector names with deSolve", {
  skip_if_not_installed("deSolve")

  named_state <- c(
    not_s_0_4 = 90, not_s_5_9 = 180,
    not_i_0_4 = 10, not_i_5_9 = 20,
    not_r_0_4 = 0, not_r_5_9 = 0
  )
  unnamed_state <- as.numeric(named_state)

  named_output <- simulate_deterministic(
    initial_state = named_state,
    times = c(0, 0.1),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve"
  )
  unnamed_output <- simulate_deterministic(
    initial_state = unnamed_state,
    times = c(0, 0.1),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve"
  )

  expect_equal(named_output, unnamed_output)
})

test_that("method deSolve works when deSolve is installed", {
  skip_if_not_installed("deSolve")

  output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve"
  )

  expect_identical(names(output), c("time", "compartment", "age_group", "value"))
  expect_equal(unique(output$time), c(0, 0.1, 0.2))
  expect_true(all(is.finite(output$value)))
})

test_that("method ode aliases the optional deSolve backend", {
  skip_if_not_installed("deSolve")

  desolve_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve"
  )
  ode_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "ode"
  )

  expect_equal(ode_output, desolve_output)
})

test_that("method deSolve errors clearly when deSolve is unavailable", {
  local_mocked_bindings(desolve_is_available = function() FALSE)

  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 0.1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = simulate_test_contacts(),
      method = "deSolve"
    ),
    "requires the optional deSolve package"
  )
})

test_that("method ode errors clearly when deSolve is unavailable", {
  local_mocked_bindings(desolve_is_available = function() FALSE)

  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 0.1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = simulate_test_contacts(),
      method = "ode"
    ),
    "requires the optional deSolve package"
  )
})

test_that("deSolve output has the same columns and ordering as Euler output", {
  skip_if_not_installed("deSolve")

  times <- c(-1, 0.1, 0.25)
  euler_output <- simulate_test_run(times = times)
  desolve_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = times,
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve"
  )

  expect_identical(names(desolve_output), names(euler_output))
  expect_identical(desolve_output$time, euler_output$time)
  expect_identical(desolve_output$compartment, euler_output$compartment)
  expect_identical(desolve_output$age_group, euler_output$age_group)
})

test_that("deSolve approximately agrees with Euler for small time steps", {
  skip_if_not_installed("deSolve")

  times <- seq(0, 0.01, by = 0.001)
  euler_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = times,
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.001,
    method = "euler"
  )
  desolve_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = times,
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.001,
    method = "deSolve"
  )

  expect_equal(desolve_output$value, euler_output$value, tolerance = 1e-4)
})

test_that("invalid times are rejected", {
  expect_error(simulate_test_run(times = 0), "at least two")
  expect_error(simulate_test_run(times = c("0", "1")), "numeric vector")
  expect_error(simulate_test_run(times = c(0, NA_real_)), "missing")
  expect_error(simulate_test_run(times = c(0, Inf)), "non-finite")
  expect_error(simulate_test_run(times = c(0, 0)), "strictly increasing")
  expect_error(simulate_test_run(times = c(0, 2, 1)), "strictly increasing")
})

test_that("unsupported method is rejected", {
  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = simulate_test_contacts(),
      method = "rk4"
    ),
    "unsupported simulation method"
  )
  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = simulate_test_contacts(),
      method = NA_character_
    ),
    "method must be a non-missing character scalar"
  )
})

test_that("negative initial state values are rejected", {
  state <- simulate_test_state()
  state$value[1] <- -1

  expect_error(
    simulate_test_run(initial_state = state),
    "initial_state values cannot be negative|state values cannot be negative"
  )
})

test_that("non-finite initial state values are rejected", {
  state <- simulate_test_state()
  state$value[1] <- Inf

  expect_error(
    simulate_test_run(initial_state = state),
    "initial_state values cannot contain non-finite values"
  )
})

test_that("invalid initial state shape is rejected", {
  expect_error(
    simulate_test_run(initial_state = list(S = 1)),
    "initial_state must be a long-form data frame or a numeric vector"
  )
  expect_error(
    simulate_test_run(initial_state = c(1, 2, 3)),
    "state_vector length"
  )
})

test_that("missing and duplicate initial state rows are rejected", {
  missing_row <- simulate_test_state()
  missing_row <- missing_row[-1, ]
  duplicate_row <- rbind(simulate_test_state(), simulate_test_state()[1, ])

  expect_error(
    simulate_test_run(initial_state = missing_row),
    "missing compartment-age row"
  )
  expect_error(
    simulate_test_run(initial_state = duplicate_row),
    "duplicate compartment-age rows"
  )
})

test_that("zero total population in an age group is rejected", {
  expect_error(
    simulate_test_run(
      initial_state = simulate_test_state(S = c(0, 180), I = c(0, 20), R = c(0, 0))
    ),
    "state population must be positive"
  )
})

test_that("Euler steps that produce negative values are rejected", {
  expect_error(
    simulate_test_run(times = c(0, 2)),
    "Euler step produced negative compartment value"
  )
})

test_that("simulate_deterministic does not accept contact schedules as contact matrices", {
  contact_schedule <- ContactSchedule(
    list("0" = simulate_test_contacts()),
    simulate_test_ages()
  )

  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 0.1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = contact_schedule
    ),
    "contact_matrix must be a numeric matrix"
  )
})

test_that("simulate_deterministic does not accept demography objects as initial state", {
  ages <- simulate_test_ages()
  demography <- Demography(
    data.frame(
      time = c(0, 0),
      age_group = ages$age_groups,
      population = c(100, 200)
    ),
    ages
  )

  expect_error(
    simulate_deterministic(
      initial_state = demography,
      times = c(0, 0.1),
      model = SIRModel(gamma = 0.2),
      age_structure = ages,
      contact_matrix = simulate_test_contacts()
    ),
    "initial_state must be a long-form data frame or a numeric vector"
  )
})

test_that("invalid contact matrix dimensions are rejected via transition validation", {
  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 0.1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = diag(3)
    ),
    "contact_matrix dimensions"
  )
})

test_that("invalid susceptibility and infectiousness are rejected via transition validation", {
  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 0.1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = simulate_test_contacts(),
      susceptibility = c(1, -1)
    ),
    "susceptibility cannot contain negative"
  )
  expect_error(
    simulate_deterministic(
      initial_state = simulate_test_state(),
      times = c(0, 0.1),
      model = SIRModel(gamma = 0.2),
      age_structure = simulate_test_ages(),
      contact_matrix = simulate_test_contacts(),
      infectiousness = c(1, 1, 1)
    ),
    "infectiousness length"
  )
})
