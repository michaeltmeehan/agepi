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

simulate_test_seir_state <- function(
  S = c(90, 180),
  E = c(5, 15),
  I = c(10, 20),
  R = c(0, 0)
) {
  data.frame(
    compartment = rep(c("S", "E", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 4),
    value = c(S, E, I, R),
    stringsAsFactors = FALSE
  )
}

simulate_test_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

simulate_onset_model <- function() {
  transitions <- data.frame(
    from = c("Lr", "Ld", "I"),
    to = c("I", "I", "T"),
    rate = c(0.2, 0.1, 0.3),
    stringsAsFactors = FALSE
  )

  CompartmentModel(
    compartments = c("S", "Lr", "Ld", "I", "T"),
    infection_transitions = data.frame(from = "S", to = "Lr", stringsAsFactors = FALSE),
    transitions = transitions,
    infectious_compartments = "I"
  )
}

simulate_onset_state <- function(
  S = c(100, 200),
  Lr = c(20, 30),
  Ld = c(40, 50),
  I = c(10, 15),
  T = c(0, 0)
) {
  data.frame(
    compartment = rep(c("S", "Lr", "Ld", "I", "T"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 5),
    value = c(S, Lr, Ld, I, T),
    stringsAsFactors = FALSE
  )
}

cumulative_values_by_time_age <- function(cumulative, cumulative_name) {
  rows <- cumulative[cumulative$cumulative_name == cumulative_name, ]
  rows <- rows[order(rows$time, rows$age_group), ]
  rows$value
}

simulate_test_run <- function(initial_state = simulate_test_state(), times = c(0, 0.1), ...) {
  simulate_deterministic(
    initial_state = initial_state,
    times = times,
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "euler",
    ...
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

test_that("simulate_deterministic works for SEIR with Euler", {
  output <- simulate_deterministic(
    initial_state = simulate_test_seir_state(),
    times = c(0, 0.1),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "euler"
  )

  final_rows <- output[output$time == 0.1, ]
  expect_identical(final_rows$compartment, c("S", "S", "E", "E", "I", "I", "R", "R"))
  expect_equal(
    final_rows$value,
    c(
      87.448504983389,
      168.159468438538,
      7.401495016611,
      26.390531561462,
      9.95,
      20.05,
      0.2,
      0.4
    )
  )
})

test_that("simulate_deterministic supports generic age-specific downstream transitions", {
  ages <- simulate_test_ages()
  transitions <- data.frame(
    from = c("E", "E", "IP", "IC", "IS"),
    to = c("IP", "IS", "IC", "R", "R"),
    stringsAsFactors = FALSE
  )
  transitions$rate <- I(list(
    c("0-4" = 0.2, "5-9" = 0.35),
    c("0-4" = 0.3, "5-9" = 0.15),
    c("0-4" = 0.25, "5-9" = 0.2),
    0.1,
    c("5-9" = 0.4, "0-4" = 0.3)
  ))
  model <- CompartmentModel(
    compartments = c("S", "E", "IP", "IC", "IS", "R"),
    infection_transitions = data.frame(from = "S", to = "E", stringsAsFactors = FALSE),
    transitions = transitions,
    infectious_compartments = c("IP", "IC", "IS"),
    infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "E", "IP", "IC", "IS", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 6),
    value = c(90, 180, 10, 20, 5, 10, 2, 4, 3, 6, 0, 0),
    stringsAsFactors = FALSE
  )

  output <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 0.1),
    model = model,
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    beta = 0,
    method = "euler"
  )

  final_rows <- output[output$time == 0.1, ]
  expect_equal(
    final_rows$value,
    c(90, 180, 9.5, 19, 5.075, 10.5, 2.105, 4.16, 3.21, 6.06, 0.11, 0.28)
  )
  expect_equal(aggregate(value ~ age_group, final_rows, sum)$value, c(110, 220))
})


test_that("method defaults to deSolve when available and Euler otherwise", {
  default_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts()
  )
  explicit_output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = if (desolve_is_available()) "deSolve" else "euler"
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
    contact_matrix = simulate_test_contacts(),
    method = "euler"
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
    contact_matrix = simulate_test_contacts(),
    method = "euler"
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
    contact_matrix = simulate_test_contacts(),
    method = "euler"
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

test_that("SIR infection-only method deSolve runs", {
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

test_that("SEIR infection-only method deSolve runs", {
  skip_if_not_installed("deSolve")

  output <- simulate_deterministic(
    initial_state = simulate_test_seir_state(),
    times = c(0, 0.1, 0.2),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve"
  )

  expect_identical(names(output), c("time", "compartment", "age_group", "value"))
  expect_equal(unique(output$time), c(0, 0.1, 0.2))
  expect_identical(unique(output$compartment), c("S", "E", "I", "R"))
  expect_true(all(is.finite(output$value)))
})

test_that("method ode aliases the deSolve backend", {
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
    "requires the deSolve package"
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
    "requires the deSolve package"
  )
})

test_that("deSolve output has the same tidy structure as Euler output", {
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

test_that("deSolve preserves total population for infection-only SIR", {
  skip_if_not_installed("deSolve")

  output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = seq(0, 1, by = 0.1),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "deSolve"
  )

  totals <- aggregate(value ~ time, output, sum)
  expect_equal(totals$value, rep(300, length(unique(output$time))), tolerance = 1e-8)
})

test_that("deSolve preserves total population for infection-only SEIR", {
  skip_if_not_installed("deSolve")

  output <- simulate_deterministic(
    initial_state = simulate_test_seir_state(),
    times = seq(0, 1, by = 0.1),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "deSolve"
  )

  totals <- aggregate(value ~ time, output, sum)
  expect_equal(totals$value, rep(320, length(unique(output$time))), tolerance = 1e-8)
})

test_that("deSolve preserves requested time vector ordering", {
  skip_if_not_installed("deSolve")

  times <- c(0, 0.05, 0.2, 0.35)
  output <- simulate_deterministic(
    initial_state = simulate_test_seir_state(),
    times = times,
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve"
  )

  expect_identical(output$time, rep(times, each = 8))
})

test_that("deSolve state values are finite and non-negative for well-behaved examples", {
  skip_if_not_installed("deSolve")

  output <- simulate_deterministic(
    initial_state = simulate_test_seir_state(),
    times = seq(0, 0.5, by = 0.1),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.001,
    method = "deSolve"
  )

  expect_true(all(is.finite(output$value)))
  expect_true(all(output$value >= 0))
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

test_that("simulate_deterministic output is unchanged when cumulative_flows is NULL", {
  without_argument <- simulate_test_run(times = c(0, 0.1, 0.2))
  explicit_null <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "euler",
    cumulative_flows = NULL
  )

  expect_equal(explicit_null, without_argument)
})

test_that("simulate_deterministic tracks basic cumulative infection flows", {
  ages <- simulate_test_ages()
  output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = seq(0, 0.2, by = 0.1),
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    method = "euler",
    cumulative_flows = list(infections = list(from = "S", to = "I"))
  )

  expect_named(output, c("trajectory", "cumulative"))
  expect_identical(names(output$trajectory), c("time", "compartment", "age_group", "value"))
  expect_identical(
    names(output$cumulative),
    c("time", "cumulative_name", "transition_id", "from", "to", "age_group", "value")
  )
  expect_true(all(output$cumulative$cumulative_name == "infections"))
  expect_true(all(output$cumulative$value >= 0))

  for (age_group in ages$age_groups) {
    values <- output$cumulative$value[output$cumulative$age_group == age_group]
    expect_true(all(diff(values) >= 0))
    initial_s <- simulate_test_state()$value[
      simulate_test_state()$compartment == "S" & simulate_test_state()$age_group == age_group
    ]
    final_s <- output$trajectory$value[
      output$trajectory$time == 0.2 &
        output$trajectory$compartment == "S" &
        output$trajectory$age_group == age_group
    ]
    expect_equal(tail(values, 1), initial_s - final_s)
  }
})

test_that("simulate_deterministic cumulative flows support generic clinical and subclinical transitions", {
  ages <- simulate_test_ages()
  transitions <- data.frame(
    from = c("E", "E", "IP", "IC", "IS"),
    to = c("IP", "IS", "IC", "R", "R"),
    stringsAsFactors = FALSE
  )
  transitions$rate <- I(list(
    c("0-4" = 0.2, "5-9" = 0.35),
    c("0-4" = 0.3, "5-9" = 0.15),
    c("0-4" = 0.25, "5-9" = 0.2),
    0.1,
    c("5-9" = 0.4, "0-4" = 0.3)
  ))
  model <- CompartmentModel(
    compartments = c("S", "E", "IP", "IC", "IS", "R"),
    infection_transitions = data.frame(from = "S", to = "E", stringsAsFactors = FALSE),
    transitions = transitions,
    infectious_compartments = c("IP", "IC", "IS"),
    infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "E", "IP", "IC", "IS", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 6),
    value = c(90, 180, 10, 20, 5, 10, 2, 4, 3, 6, 0, 0),
    stringsAsFactors = FALSE
  )

  output <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 0.1),
    model = model,
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "euler",
    cumulative_flows = data.frame(
      name = c("exposures", "clinical", "subclinical"),
      from = c("S", "E", "E"),
      to = c("E", "IP", "IS"),
      stringsAsFactors = FALSE
    )
  )

  cumulative <- output$cumulative
  expect_identical(
    names(cumulative),
    c("time", "cumulative_name", "transition_id", "from", "to", "age_group", "value")
  )
  expect_equal(unique(cumulative$cumulative_name), c("exposures", "clinical", "subclinical"))
  expect_equal(unique(cumulative$age_group), ages$age_groups)
  expect_true(all(c("infection:S->E", "transition:E->IP", "transition:E->IS") %in% cumulative$transition_id))
})

test_that("simulate_deterministic aggregates multiple transitions into one cumulative flow", {
  output <- simulate_deterministic(
    initial_state = simulate_onset_state(),
    times = seq(0, 0.3, by = 0.1),
    model = simulate_onset_model(),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "euler",
    cumulative_flows = list(
      infections = list(from = "S", to = "Lr"),
      disease_onset = list(from = c("Lr", "Ld"), to = c("I", "I")),
      treatment = list(from = "I", to = "T")
    )
  )

  cumulative <- output$cumulative
  expect_equal(unique(cumulative$cumulative_name), c("infections", "disease_onset", "treatment"))
  expect_equal(
    unique(cumulative$transition_id[cumulative$cumulative_name == "disease_onset"]),
    "transition:Lr->I,transition:Ld->I"
  )
  expect_equal(nrow(cumulative), length(seq(0, 0.3, by = 0.1)) * 3 * 2)
  expect_true(all(cumulative_values_by_time_age(cumulative, "disease_onset") >= 0))
})

test_that("multi-transition cumulative flow equals the corresponding separate counters", {
  times <- seq(0, 0.3, by = 0.1)
  aggregate_output <- simulate_deterministic(
    initial_state = simulate_onset_state(),
    times = times,
    model = simulate_onset_model(),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "euler",
    cumulative_flows = list(
      disease_onset = list(from = c("Lr", "Ld"), to = c("I", "I"))
    )
  )
  separate_output <- simulate_deterministic(
    initial_state = simulate_onset_state(),
    times = times,
    model = simulate_onset_model(),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "euler",
    cumulative_flows = list(
      recent_onset = list(from = "Lr", to = "I"),
      remote_onset = list(from = "Ld", to = "I")
    )
  )

  aggregate_values <- cumulative_values_by_time_age(aggregate_output$cumulative, "disease_onset")
  separate_values <- cumulative_values_by_time_age(separate_output$cumulative, "recent_onset") +
    cumulative_values_by_time_age(separate_output$cumulative, "remote_onset")

  expect_equal(aggregate_values, separate_values)
})

test_that("cumulative derivative equals the transition flow used by compartment derivative", {
  ages <- simulate_test_ages()
  state <- state_long_to_vector(simulate_test_state(), ages, c("S", "I", "R"))
  spec <- prepare_deterministic_cumulative_flows(
    cumulative_flows = list(infections = list(from = "S", to = "I")),
    state_vector = state,
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    beta = 1,
    susceptibility = NULL,
    infectiousness = NULL
  )

  derivative <- deterministic_derivative_augmented(
    state_vector = c(state, 0, 0),
    ordinary_state_length = length(state),
    time = 0,
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    beta = 1,
    susceptibility = NULL,
    infectiousness = NULL,
    cumulative_spec = spec
  )
  rates <- transition_rates(state, SIRModel(gamma = 0.2), ages, simulate_test_contacts())

  expect_equal(tail(derivative, 2), rates$rate[rates$transition_id == "infection:S->I"])
  expect_equal(head(derivative, length(state)), rates_to_derivative(rates, c("S", "I", "R"), ages)$derivative)
})

test_that("cumulative_flows work with deSolve when available", {
  skip_if_not_installed("deSolve")

  output <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = simulate_test_ages(),
    contact_matrix = simulate_test_contacts(),
    method = "deSolve",
    cumulative_flows = list(infections = list(from = "S", to = "I"))
  )

  expect_named(output, c("trajectory", "cumulative"))
  expect_true(all(is.finite(output$cumulative$value)))
  expect_true(all(output$cumulative$value >= -1e-8))
})

test_that("cumulative states do not alter ordinary trajectories or population totals", {
  ages <- simulate_test_ages()
  baseline <- simulate_test_run(times = seq(0, 0.2, by = 0.1))
  with_cumulative <- simulate_deterministic(
    initial_state = simulate_test_state(),
    times = seq(0, 0.2, by = 0.1),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    method = "euler",
    cumulative_flows = list(infections = list(from = "S", to = "I"))
  )

  expect_equal(with_cumulative$trajectory, baseline)
  expect_equal(total_population(with_cumulative$trajectory), total_population(baseline))
})

test_that("cumulative_flows do not contaminate force-of-infection denominators", {
  ages <- simulate_test_ages()
  state <- state_long_to_vector(simulate_test_state(), ages, c("S", "I", "R"))
  spec <- prepare_deterministic_cumulative_flows(
    cumulative_flows = list(infections = list(from = "S", to = "I")),
    state_vector = state,
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    beta = 1,
    susceptibility = NULL,
    infectiousness = NULL
  )

  low_counter <- deterministic_derivative_augmented(
    c(state, 0, 0), length(state), 0, SIRModel(gamma = 0.2), ages,
    simulate_test_contacts(), 1, NULL, NULL, cumulative_spec = spec
  )
  high_counter <- deterministic_derivative_augmented(
    c(state, 1e9, 1e9), length(state), 0, SIRModel(gamma = 0.2), ages,
    simulate_test_contacts(), 1, NULL, NULL, cumulative_spec = spec
  )

  expect_equal(high_counter, low_counter)
})

test_that("cumulative_flows with demography preserves the ordinary trajectory", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5+"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, Inf)
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 3),
    value = c(90, 180, 10, 20, 0, 0),
    stringsAsFactors = FALSE
  )
  process <- DemographicProcess(ages)
  baseline <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    method = "euler",
    demographic_process = process
  )
  output <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    method = "euler",
    demographic_process = process,
    cumulative_flows = list(
      infections = list(from = "S", to = "I"),
      removals = list(from = "I", to = "R")
    )
  )

  expect_named(output, c("trajectory", "cumulative"))
  expect_equal(output$trajectory, baseline)
  expect_true(all(output$cumulative$value >= 0))
})

test_that("multi-transition cumulative_flows work with demography", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, Inf)
  )
  times <- seq(0, 0.3, by = 0.1)
  aggregate_output <- simulate_deterministic(
    initial_state = simulate_onset_state(),
    times = times,
    model = simulate_onset_model(),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "euler",
    demographic_process = DemographicProcess(ages),
    cumulative_flows = list(
      disease_onset = list(from = c("Lr", "Ld"), to = c("I", "I"))
    )
  )
  separate_output <- simulate_deterministic(
    initial_state = simulate_onset_state(),
    times = times,
    model = simulate_onset_model(),
    age_structure = ages,
    contact_matrix = simulate_test_contacts(),
    beta = 0.01,
    method = "euler",
    demographic_process = DemographicProcess(ages),
    cumulative_flows = list(
      recent_onset = list(from = "Lr", to = "I"),
      remote_onset = list(from = "Ld", to = "I")
    )
  )

  aggregate_values <- cumulative_values_by_time_age(aggregate_output$cumulative, "disease_onset")
  separate_values <- cumulative_values_by_time_age(separate_output$cumulative, "recent_onset") +
    cumulative_values_by_time_age(separate_output$cumulative, "remote_onset")

  expect_equal(aggregate_values, separate_values)
  expect_true(all(aggregate_values >= 0))
})

test_that("invalid cumulative_flows inputs error clearly", {
  expect_error(
    simulate_test_run(cumulative_flows = list(list(from = "S", to = "I"))),
    "list entries must be named"
  )
  expect_error(
    simulate_test_run(cumulative_flows = list(bad = list(from = "X", to = "I"))),
    "unknown source compartment"
  )
  expect_error(
    simulate_test_run(cumulative_flows = data.frame(name = "x", from = "S")),
    "missing required column"
  )
  expect_error(
    simulate_deterministic(
      initial_state = simulate_onset_state(),
      times = c(0, 0.1),
      model = simulate_onset_model(),
      age_structure = simulate_test_ages(),
      contact_matrix = simulate_test_contacts(),
      method = "euler",
      cumulative_flows = list(disease_onset = list(from = c("Lr", "Ld"), to = "I"))
    ),
    "same length"
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
