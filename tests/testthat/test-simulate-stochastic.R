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
  demographic_process = NULL,
  cumulative_flows = NULL
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
    demographic_process = demographic_process,
    cumulative_flows = cumulative_flows
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

stochastic_vaccine_ages <- function() {
  AgeStructure(
    age_groups = c("0-4"),
    lower_bounds = c(0),
    upper_bounds = c(4)
  )
}

stochastic_vaccine_state <- function(S = 20, V = 20, I = 10) {
  data.frame(
    compartment = c("S", "V", "I"),
    age_group = c("0-4", "0-4", "0-4"),
    value = c(S, V, I),
    stringsAsFactors = FALSE
  )
}

stochastic_vaccine_model <- function(susceptibility = list(1, 0.5)) {
  infection_transitions <- data.frame(
    from = c("S", "V"),
    to = c("I", "I"),
    stringsAsFactors = FALSE
  )
  infection_transitions$susceptibility <- I(if (is.null(susceptibility)) list(1, 1) else susceptibility)

  CompartmentModel(
    compartments = c("S", "V", "I"),
    infection_transitions = infection_transitions,
    transitions = data.frame(
      from = character(),
      to = character(),
      rate = numeric(),
      stringsAsFactors = FALSE
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
  expect_identical(names(result$events), c("time", "event", "transition_type", "transition_id", "age_group", "from", "to", "rate"))
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

  expect_identical(names(result$events), c("time", "event", "transition_type", "transition_id", "age_group", "from", "to", "rate"))
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

  expect_identical(names(result$events), c("time", "event", "transition_type", "transition_id", "age_group", "from", "to", "rate"))
  expect_true(all(result$events$event %in% c("infection", "progression", "recovery")))
})

test_that("stochastic CompartmentModel supports outflow transitions", {
  ages <- stochastic_test_ages()
  model <- CompartmentModel(
    compartments = c("S", "I"),
    outflows = data.frame(from = "I", rate = 5, stringsAsFactors = FALSE),
    infectious_compartments = character()
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "I"), each = 2),
    age_group = rep(ages$age_groups, times = 2),
    value = c(0, 0, 10, 12),
    stringsAsFactors = FALSE
  )

  result <- simulate_stochastic(
    initial_state = initial_state,
    times = c(0, 1),
    model = model,
    age_structure = ages,
    contact_matrix = stochastic_test_contacts(),
    beta = 0,
    seed = 11,
    return_events = TRUE,
    cumulative_flows = list(removals = list(from = "I", to = NA_character_))
  )

  expect_true(nrow(result$events) > 0)
  expect_true(all(result$events$event == "outflow"))
  expect_true(all(is.na(result$events$to)))
  expect_equal(unique(result$cumulative$transition_id), "outflow:I")
  expect_true(all(result$cumulative$value >= 0))
  totals <- aggregate(value ~ time, result$trajectory, sum)
  expect_true(all(diff(totals$value) <= 0))
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

test_that("generic stochastic simulation inherits age-specific downstream transitions", {
  ages <- stochastic_test_ages()
  transitions <- data.frame(
    from = c("E", "E", "IP", "IC", "IS"),
    to = c("IP", "IS", "IC", "R", "R")
  )
  transitions$rate <- I(list(
    c("0-4" = 2, "5-9" = 3),
    c("0-4" = 4, "5-9" = 1),
    c("5-9" = 2, "0-4" = 1),
    1,
    c("0-4" = 2, "5-9" = 2)
  ))
  model <- CompartmentModel(
    compartments = c("S", "E", "IP", "IC", "IS", "R"),
    infection_transitions = data.frame(from = "S", to = "E"),
    transitions = transitions,
    infectious_compartments = c("IP", "IC", "IS"),
    infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "E", "IP", "IC", "IS", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 6),
    value = c(90, 180, 4, 4, 1, 1, 0, 0, 1, 1, 0, 0),
    stringsAsFactors = FALSE
  )

  result <- simulate_stochastic(
    initial_state = initial_state,
    times = c(0, 10),
    model = model,
    age_structure = ages,
    contact_matrix = stochastic_test_contacts(),
    beta = 0,
    seed = 45,
    return_events = TRUE
  )

  expect_true(nrow(result$events) > 0)
  expect_false(any(result$events$event == "infection"))
  expect_true(all(result$events$from %in% c("E", "IP", "IC", "IS")))
  expect_true(all(result$trajectory$value >= 0))

  final_rows <- result$trajectory[result$trajectory$time == 10, ]
  totals <- aggregate(value ~ age_group, final_rows, sum)
  expect_equal(totals$value, c(96, 186))
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

  expect_identical(names(result$events), c("time", "event", "transition_type", "transition_id", "age_group", "from", "to", "rate"))
  expect_true(nrow(result$events) > 0)
  expect_true(all(result$events$transition_id %in% c("infection:S->I", "transition:I->R")))
  expect_true(all(result$events$from %in% c("S", "I")))
  expect_true(all(result$events$to %in% c("I", "R")))
  expect_true(all(diff(result$events$time) >= 0))
})

test_that("stochastic simulation preserves legacy infection-transition semantics when susceptibility is omitted", {
  without_susceptibility <- simulate_stochastic(
    initial_state = stochastic_vaccine_state(),
    times = c(0, 5),
    model = stochastic_vaccine_model(NULL),
    age_structure = stochastic_vaccine_ages(),
    contact_matrix = matrix(1, nrow = 1, ncol = 1),
    beta = 0.6,
    seed = 61,
    return_events = TRUE
  )
  explicit_one <- simulate_stochastic(
    initial_state = stochastic_vaccine_state(),
    times = c(0, 5),
    model = stochastic_vaccine_model(list(1, 1)),
    age_structure = stochastic_vaccine_ages(),
    contact_matrix = matrix(1, nrow = 1, ncol = 1),
    beta = 0.6,
    seed = 61,
    return_events = TRUE
  )

  expect_equal(without_susceptibility, explicit_one)
})

test_that("stochastic simulation treats scalar and list infectiousness weights equivalently", {
  ages <- stochastic_test_ages()
  state <- stochastic_test_state()
  scalar_model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    infectious_compartments = "I",
    infectiousness_weights = 0.5
  )
  list_model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    infectious_compartments = "I",
    infectiousness_weights = list(0.5)
  )
  expanded_model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    infectious_compartments = "I",
    infectiousness_weights = list(c("0-4" = 0.5, "5-9" = 0.5))
  )

  scalar_output <- simulate_stochastic(
    initial_state = state,
    times = c(0, 5),
    model = scalar_model,
    age_structure = ages,
    contact_matrix = diag(2),
    beta = 0.03,
    seed = 61,
    return_events = TRUE
  )
  list_output <- simulate_stochastic(
    initial_state = state,
    times = c(0, 5),
    model = list_model,
    age_structure = ages,
    contact_matrix = diag(2),
    beta = 0.03,
    seed = 61,
    return_events = TRUE
  )
  expanded_output <- simulate_stochastic(
    initial_state = state,
    times = c(0, 5),
    model = expanded_model,
    age_structure = ages,
    contact_matrix = diag(2),
    beta = 0.03,
    seed = 61,
    return_events = TRUE
  )

  expect_equal(list_output, scalar_output)
  expect_equal(expanded_output, scalar_output)
})

test_that("stochastic simulation never generates infections from zero-susceptibility compartments", {
  result <- simulate_stochastic(
    initial_state = stochastic_vaccine_state(S = 30, V = 30, I = 10),
    times = c(0, 5),
    model = stochastic_vaccine_model(list(1, 0)),
    age_structure = stochastic_vaccine_ages(),
    contact_matrix = matrix(1, nrow = 1, ncol = 1),
    beta = 1,
    seed = 62,
    return_events = TRUE
  )

  infection_events <- result$events[result$events$event == "infection", ]
  expect_true(nrow(infection_events) > 0)
  expect_true(all(infection_events$from == "S"))
  expect_false(any(infection_events$from == "V"))
})

test_that("simulate_stochastic respects model beta, explicit zero overrides, and non-infectious models", {
  model <- SIRModel(gamma = 0.2, beta = 0.4)

  default_output <- simulate_stochastic(
    initial_state = stochastic_test_state(),
    times = c(0, 5),
    model = model,
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    seed = 61,
    return_events = TRUE
  )
  explicit_model_output <- simulate_stochastic(
    initial_state = stochastic_test_state(),
    times = c(0, 5),
    model = model,
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = 0.4,
    seed = 61,
    return_events = TRUE
  )
  zero_override <- simulate_stochastic(
    initial_state = stochastic_test_state(),
    times = c(0, 5),
    model = model,
    age_structure = stochastic_test_ages(),
    contact_matrix = stochastic_test_contacts(),
    beta = 0,
    seed = 61,
    return_events = TRUE
  )

  expect_identical(model$beta, 0.4)
  expect_equal(default_output, explicit_model_output)
  expect_equal(sum(zero_override$events$event == "infection"), 0)

  non_infectious_model <- CompartmentModel(
    compartments = c("S", "R"),
    transitions = data.frame(from = "S", to = "R", rate = 0.1, stringsAsFactors = FALSE),
    infectious_compartments = character()
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 2),
    value = c(100, 80, 0, 0),
    stringsAsFactors = FALSE
  )

  baseline <- simulate_stochastic(
    initial_state = initial_state,
    times = c(0, 2),
    model = non_infectious_model,
    age_structure = stochastic_test_ages(),
    contact_matrix = matrix(0, nrow = 2, ncol = 2),
    seed = 91,
    return_events = TRUE
  )
  explicit_beta <- simulate_stochastic(
    initial_state = initial_state,
    times = c(0, 2),
    model = non_infectious_model,
    age_structure = stochastic_test_ages(),
    contact_matrix = matrix(0, nrow = 2, ncol = 2),
    beta = 0.7,
    seed = 91,
    return_events = TRUE
  )

  expect_equal(baseline, explicit_beta)
})

test_that("simulate_stochastic output is unchanged when cumulative_flows is NULL", {
  baseline <- stochastic_test_run(seed = 51)
  explicit_null <- stochastic_test_run(seed = 51, cumulative_flows = NULL)

  expect_equal(explicit_null, baseline)
  expect_identical(names(explicit_null), c("time", "compartment", "age_group", "value"))
})

test_that("stochastic cumulative infection counts match realised events", {
  result <- stochastic_test_run(
    initial_state = stochastic_test_state(S = c(3, 3), I = c(5, 5), R = c(0, 0)),
    times = c(0, 1, 5, 20),
    beta = 2,
    gamma = 0,
    seed = 52,
    return_events = TRUE,
    cumulative_flows = list(infections = list(from = "S", to = "I"))
  )

  expect_named(result, c("trajectory", "events", "cumulative"))
  expect_identical(
    names(result$cumulative),
    c("time", "cumulative_name", "transition_id", "from", "to", "age_group", "value")
  )
  expect_true(all(result$cumulative$cumulative_name == "infections"))
  expect_true(all(result$cumulative$transition_id == "infection:S->I"))

  for (age_group in stochastic_test_ages()$age_groups) {
    values <- result$cumulative$value[result$cumulative$age_group == age_group]
    expect_true(all(diff(values) >= 0))
    expect_equal(
      tail(values, 1),
      sum(result$events$transition_id == "infection:S->I" & result$events$age_group == age_group)
    )
  }
})

test_that("stochastic cumulative outputs support generic clinical and subclinical flows", {
  ages <- stochastic_test_ages()
  transitions <- data.frame(
    from = c("E", "E", "IP", "IC", "IS"),
    to = c("IP", "IS", "IC", "R", "R"),
    stringsAsFactors = FALSE
  )
  transitions$rate <- I(list(
    c("0-4" = 2, "5-9" = 3),
    c("0-4" = 4, "5-9" = 1),
    1,
    1,
    1
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
    value = c(10, 10, 2, 2, 3, 3, 0, 0, 1, 1, 0, 0),
    stringsAsFactors = FALSE
  )

  result <- simulate_stochastic(
    initial_state = initial_state,
    times = c(0, 0.5, 2),
    model = model,
    age_structure = ages,
    contact_matrix = stochastic_test_contacts(),
    beta = 0.2,
    seed = 53,
    return_events = TRUE,
    cumulative_flows = data.frame(
      name = c("exposures", "clinical", "subclinical"),
      from = c("S", "E", "E"),
      to = c("E", "IP", "IS"),
      stringsAsFactors = FALSE
    )
  )

  expect_identical(
    names(result$cumulative),
    c("time", "cumulative_name", "transition_id", "from", "to", "age_group", "value")
  )
  expect_equal(unique(result$cumulative$cumulative_name), c("exposures", "clinical", "subclinical"))
  expect_equal(unique(result$cumulative$age_group), ages$age_groups)
  expect_true(all(c("infection:S->E", "transition:E->IP", "transition:E->IS") %in% result$cumulative$transition_id))
})

test_that("stochastic cumulative output remains zero for selected flows with no events", {
  result <- stochastic_test_run(
    initial_state = stochastic_test_state(S = c(10, 10), I = c(0, 0), R = c(0, 0)),
    times = c(0, 1, 2),
    beta = 0,
    gamma = 0,
    seed = 54,
    cumulative_flows = list(infections = list(from = "S", to = "I"))
  )

  expect_named(result, c("trajectory", "cumulative"))
  expect_equal(result$cumulative$value, rep(0, 6))
  expect_equal(result$cumulative$age_group, rep(stochastic_test_ages()$age_groups, times = 3))
})

test_that("stochastic cumulative counts use inclusive output-time convention", {
  events <- data.frame(
    time = c(0.5, 1),
    event = c("infection", "infection"),
    transition_id = c("infection:S->I", "infection:S->I"),
    age_group = c("0-4", "0-4"),
    from = c("S", "S"),
    to = c("I", "I"),
    rate = c(1, 1),
    stringsAsFactors = FALSE
  )
  spec <- list(
    output_order = data.frame(
      cumulative_name = "infections",
      transition_id = "infection:S->I",
      from = "S",
      to = "I",
      age_group = c("0-4", "5-9"),
      stringsAsFactors = FALSE
    )
  )

  cumulative <- stochastic_cumulative_output(events, times = c(0, 0.5, 1), cumulative_spec = spec)

  expect_equal(cumulative$value[cumulative$age_group == "0-4"], c(0, 1, 2))
  expect_equal(cumulative$value[cumulative$age_group == "5-9"], c(0, 0, 0))
})

test_that("cumulative_flows do not alter stochastic trajectories or event sequence", {
  baseline <- stochastic_test_run(
    initial_state = stochastic_test_state(S = c(3, 3), I = c(5, 5), R = c(0, 0)),
    times = c(0, 1, 5, 20),
    beta = 2,
    gamma = 0.5,
    seed = 55,
    return_events = TRUE
  )
  with_cumulative <- stochastic_test_run(
    initial_state = stochastic_test_state(S = c(3, 3), I = c(5, 5), R = c(0, 0)),
    times = c(0, 1, 5, 20),
    beta = 2,
    gamma = 0.5,
    seed = 55,
    return_events = TRUE,
    cumulative_flows = list(infections = list(from = "S", to = "I"))
  )

  expect_equal(with_cumulative$trajectory, baseline$trajectory)
  expect_equal(with_cumulative$events, baseline$events)
  expect_equal(
    aggregate(value ~ time + age_group, with_cumulative$trajectory, sum),
    aggregate(value ~ time + age_group, baseline$trajectory, sum)
  )
})

test_that("invalid stochastic cumulative_flows inputs error clearly", {
  expect_error(
    stochastic_test_run(cumulative_flows = list(list(from = "S", to = "I"))),
    "cumulative_flows list entries must be named"
  )
  expect_error(
    stochastic_test_run(cumulative_flows = list(bad = list(from = "X", to = "I"))),
    "matched no transitions"
  )
  expect_silent(
    stochastic_test_run(cumulative_flows = data.frame(name = "x", from = "S"))
  )
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
