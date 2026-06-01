stochastic_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

stochastic_test_state <- function(S = c(90, 180), I = c(10, 20), R = c(0, 0)) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

stochastic_test_seir_state <- function(
  S = c(90, 180),
  E = c(5, 10),
  I = c(5, 10),
  R = c(0, 0)
) {
  data.frame(
    compartment = rep(c("S", "E", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 4),
    value = c(S, E, I, R),
    stringsAsFactors = FALSE
  )
}

stochastic_test_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

stochastic_test_run <- function(
  initial_state = stochastic_test_state(),
  times = seq(0, 1, by = 0.25),
  seed = 42,
  return_events = FALSE,
  method = "gillespie",
  beta = 0.03,
  gamma = 0.2,
  demographic_process = NULL
) {
  simulate_stochastic(
    initial_state = initial_state,
    times = times,
    model = SIRModel(gamma = gamma),
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = beta,
    method = method,
    seed = seed,
    return_events = return_events,
    demographic_process = demographic_process
  )
}

stochastic_test_generic_sir_model <- function(gamma = 0.2) {
  CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I"),
    transitions = data.frame(from = "I", to = "R", rate = gamma),
    infectious_compartments = "I"
  )
}

stochastic_test_generic_seir_model <- function(sigma = 0.4, gamma = 0.2) {
  CompartmentModel(
    compartments = c("S", "E", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "E"),
    transitions = data.frame(
      from = c("E", "I"),
      to = c("I", "R"),
      rate = c(sigma, gamma)
    ),
    infectious_compartments = "I"
  )
}

stochastic_test_seir_run <- function(
  initial_state = stochastic_test_seir_state(),
  times = seq(0, 1, by = 0.25),
  seed = 42,
  return_events = FALSE,
  beta = 0.03,
  sigma = 0.4,
  gamma = 0.2,
  demographic_process = NULL
) {
  simulate_stochastic(
    initial_state = initial_state,
    times = times,
    model = SEIRModel(sigma = sigma, gamma = gamma),
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = beta,
    seed = seed,
    return_events = return_events,
    demographic_process = demographic_process
  )
}

test_that("simulate_stochastic is reproducible with a fixed seed", {
  first <- stochastic_test_run(seed = 12, return_events = TRUE)
  second <- stochastic_test_run(seed = 12, return_events = TRUE)

  expect_equal(first, second)
})

test_that("simulate_stochastic includes requested output times", {
  times <- c(0, 0.2, 0.7, 1.3)
  output <- stochastic_test_run(times = times)

  expect_identical(unique(output$time), times)
  expect_identical(names(output), c("time", "compartment", "age_group", "value"))
})

test_that("infection-only stochastic SIR conserves total population", {
  output <- stochastic_test_run(times = seq(0, 2, by = 0.5), seed = 99)

  totals <- aggregate(value ~ time + age_group, output, sum)
  expect_equal(totals$value[totals$age_group == "0-4"], rep(100, 5))
  expect_equal(totals$value[totals$age_group == "5-9"], rep(200, 5))

  overall_totals <- aggregate(value ~ time, output, sum)
  expect_equal(overall_totals$value, rep(300, 5))
})

test_that("zero initial infections produces no stochastic events", {
  result <- stochastic_test_run(
    initial_state = stochastic_test_state(I = c(0, 0)),
    times = seq(0, 1, by = 0.25),
    return_events = TRUE
  )

  expect_equal(nrow(result$events), 0)
  expect_equal(result$trajectory$value[result$trajectory$compartment == "I"], rep(0, 10))
})

test_that("zero total event rate carries state forward to requested times", {
  initial_state <- stochastic_test_state(I = c(0, 0))
  result <- stochastic_test_run(
    initial_state = initial_state,
    times = c(0, 0.5, 2),
    beta = 0,
    gamma = 0,
    return_events = TRUE
  )

  expected_values <- rep(initial_state$value, times = 3)
  expect_equal(result$trajectory$value, expected_values)
  expect_equal(result$trajectory$time, rep(c(0, 0.5, 2), each = 6))
  expect_identical(names(result$events), c("time", "event", "age_group", "from", "to", "rate"))
  expect_equal(nrow(result$events), 0)
})

test_that("event log records recovery-only semantics when beta is zero", {
  result <- stochastic_test_run(
    initial_state = stochastic_test_state(S = c(90, 180), I = c(1, 1), R = c(0, 0)),
    times = c(0, 100),
    beta = 0,
    gamma = 10,
    seed = 7,
    return_events = TRUE
  )

  expect_identical(names(result$events), c("time", "event", "age_group", "from", "to", "rate"))
  expect_true(nrow(result$events) > 0)
  expect_true(all(result$events$event == "recovery"))
  expect_true(all(result$events$from == "I"))
  expect_true(all(result$events$to == "R"))
  expect_true(all(diff(result$events$time) >= 0))
})

test_that("event log records infection-only semantics when gamma is zero", {
  result <- stochastic_test_run(
    initial_state = stochastic_test_state(S = c(1, 1), I = c(10, 10), R = c(0, 0)),
    times = c(0, 100),
    beta = 10,
    gamma = 0,
    seed = 8,
    return_events = TRUE
  )

  expect_true(nrow(result$events) > 0)
  expect_true(all(result$events$event == "infection"))
  expect_true(all(result$events$from == "S"))
  expect_true(all(result$events$to == "I"))
  expect_true(all(diff(result$events$time) >= 0))
})

test_that("stochastic SEIR is reproducible with a fixed seed", {
  first <- stochastic_test_seir_run(seed = 21, return_events = TRUE)
  second <- stochastic_test_seir_run(seed = 21, return_events = TRUE)

  expect_equal(first, second)
})

test_that("stochastic SEIR includes requested output times exactly", {
  times <- c(0, 0.2, 0.7, 1.3)
  output <- stochastic_test_seir_run(times = times)

  expect_identical(unique(output$time), times)
  expect_identical(names(output), c("time", "compartment", "age_group", "value"))
})

test_that("stochastic SEIR conserves total population exactly", {
  output <- stochastic_test_seir_run(times = seq(0, 2, by = 0.5), seed = 31)

  totals <- aggregate(value ~ time + age_group, output, sum)
  expect_equal(totals$value[totals$age_group == "0-4"], rep(100, 5))
  expect_equal(totals$value[totals$age_group == "5-9"], rep(200, 5))

  overall_totals <- aggregate(value ~ time, output, sum)
  expect_equal(overall_totals$value, rep(300, 5))
})

test_that("zero exposed and infected SEIR state produces no events and carries state forward", {
  initial_state <- stochastic_test_seir_state(E = c(0, 0), I = c(0, 0))
  result <- stochastic_test_seir_run(
    initial_state = initial_state,
    times = c(0, 0.5, 2),
    return_events = TRUE
  )

  expect_equal(nrow(result$events), 0)
  expect_equal(result$trajectory$value, rep(initial_state$value, times = 3))
})

test_that("SEIR beta zero allows progression and recovery but no new infections", {
  result <- stochastic_test_seir_run(
    initial_state = stochastic_test_seir_state(S = c(90, 180), E = c(1, 1), I = c(1, 1)),
    times = c(0, 100),
    beta = 0,
    sigma = 10,
    gamma = 10,
    seed = 32,
    return_events = TRUE
  )

  expect_true(nrow(result$events) > 0)
  expect_false(any(result$events$event == "infection"))
  expect_true(all(result$events$event %in% c("progression", "recovery")))
})

test_that("SEIR sigma zero prevents exposed to infected progression", {
  result <- stochastic_test_seir_run(
    initial_state = stochastic_test_seir_state(S = c(0, 0), E = c(2, 2), I = c(0, 0)),
    times = c(0, 100),
    beta = 10,
    sigma = 0,
    gamma = 10,
    seed = 33,
    return_events = TRUE
  )

  expect_equal(nrow(result$events), 0)
  exposed <- result$trajectory$value[result$trajectory$compartment == "E"]
  expect_equal(exposed, rep(c(2, 2), times = 2))
})

test_that("SEIR gamma zero allows infection and progression but no recovery", {
  result <- stochastic_test_seir_run(
    initial_state = stochastic_test_seir_state(S = c(2, 2), E = c(2, 2), I = c(4, 4), R = c(0, 0)),
    times = c(0, 100),
    beta = 10,
    sigma = 10,
    gamma = 0,
    seed = 34,
    return_events = TRUE
  )

  expect_true(nrow(result$events) > 0)
  expect_false(any(result$events$event == "recovery"))
  expect_true(all(result$events$event %in% c("infection", "progression")))
})

test_that("SEIR event times are nondecreasing", {
  result <- stochastic_test_seir_run(times = c(0, 20), seed = 35, return_events = TRUE)

  expect_true(nrow(result$events) > 0)
  expect_true(all(diff(result$events$time) >= 0))
})

test_that("SEIR event log uses stable column order and event labels", {
  result <- stochastic_test_seir_run(times = c(0, 20), seed = 36, return_events = TRUE)

  expect_identical(names(result$events), c("time", "event", "age_group", "from", "to", "rate"))
  expect_true(all(result$events$event %in% c("infection", "progression", "recovery")))
})

test_that("generic CompartmentModel SIR path matches specialised SIR when seeded", {
  specialised <- stochastic_test_run(seed = 41, return_events = TRUE)
  generic <- simulate_stochastic(
    initial_state = stochastic_test_state(),
    times = seq(0, 1, by = 0.25),
    model = stochastic_test_generic_sir_model(gamma = 0.2),
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = 0.03,
    seed = 41,
    return_events = TRUE
  )

  expect_equal(generic, specialised)
})

test_that("generic CompartmentModel SEIR path matches specialised SEIR when seeded", {
  specialised <- stochastic_test_seir_run(seed = 42, return_events = TRUE)
  generic <- simulate_stochastic(
    initial_state = stochastic_test_seir_state(),
    times = seq(0, 1, by = 0.25),
    model = stochastic_test_generic_seir_model(sigma = 0.4, gamma = 0.2),
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = 0.03,
    seed = 42,
    return_events = TRUE
  )

  expect_equal(generic, specialised)
})

test_that("generic stochastic trajectories stay nonnegative and fixed population", {
  result <- simulate_stochastic(
    initial_state = stochastic_test_seir_state(),
    times = seq(0, 2, by = 0.5),
    model = stochastic_test_generic_seir_model(sigma = 0.4, gamma = 0.2),
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = 0.03,
    seed = 43
  )

  expect_true(all(result$value >= 0))

  totals <- aggregate(value ~ time + age_group, result, sum)
  expect_equal(totals$value[totals$age_group == "0-4"], rep(100, 5))
  expect_equal(totals$value[totals$age_group == "5-9"], rep(200, 5))
})

test_that("generic event log keeps the stochastic event structure", {
  result <- simulate_stochastic(
    initial_state = stochastic_test_state(S = c(1, 1), I = c(4, 4), R = c(0, 0)),
    times = c(0, 20),
    model = stochastic_test_generic_sir_model(gamma = 1),
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = 1,
    seed = 44,
    return_events = TRUE
  )

  expect_identical(names(result$events), c("time", "event", "age_group", "from", "to", "rate"))
  expect_true(nrow(result$events) > 0)
  expect_true(all(result$events$from %in% c("S", "I")))
  expect_true(all(result$events$to %in% c("I", "R")))
  expect_true(all(diff(result$events$time) >= 0))
})

test_that("simulate_stochastic rejects unsupported model structures and method", {
  unsupported_model <- list(
    model_type = "Unsupported",
    compartments = c("S", "I", "R"),
    transitions = data.frame(from = "S", to = "I")
  )
  class(unsupported_model) <- "DiseaseModel"

  expect_error(
    simulate_stochastic(
      initial_state = stochastic_test_state(),
      times = c(0, 1),
      model = unsupported_model,
      age_structure = stochastic_test_ages(),
      contact_matrix = stochastic_test_contacts()
    ),
    "unsupported disease model type"
  )

  bad_rate_model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    transitions = data.frame(from = "I", to = "R", rate = I(list(c(0.1, 0.2, 0.3)))),
    infectious_compartments = character()
  )

  expect_error(
    simulate_stochastic(
      initial_state = stochastic_test_state(),
      times = c(0, 1),
      model = bad_rate_model,
      age_structure = stochastic_test_ages(),
      contact_matrix = stochastic_test_contacts()
    ),
    "transition rate for I->R length must be 1 or match the number of age groups"
  )

  expect_error(
    stochastic_test_run(seed = 1, return_events = FALSE, times = c(0, 1), method = "tau"),
    "unsupported stochastic simulation method"
  )
})

test_that("simulate_stochastic rejects demographic processes", {
  demographic_process <- list(not = "a supported stochastic input")

  expect_error(
    stochastic_test_run(
      seed = 1,
      times = c(0, 1),
      demographic_process = demographic_process
    ),
    "demographic_process must be NULL"
  )
})
