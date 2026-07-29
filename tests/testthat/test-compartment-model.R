generic_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

generic_test_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

generic_sir_state <- function(S = c(90, 180), I = c(10, 20), R = c(0, 0)) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

generic_seir_state <- function(
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

generic_sir_model <- function(gamma = 0.2) {
  CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = gamma, stringsAsFactors = FALSE),
    infectious_compartments = "I"
  )
}

generic_seir_model <- function(sigma = 0.3, gamma = 0.2) {
  CompartmentModel(
    compartments = c("S", "E", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "E", stringsAsFactors = FALSE),
    transitions = data.frame(
      from = c("E", "I"),
      to = c("I", "R"),
      rate = c(sigma, gamma),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = "I"
  )
}

generic_clinical_subclinical_model <- function(ages = generic_test_ages()) {
  sigma <- 0.5
  y_age <- c("0-4" = 0.4, "5-9" = 0.7)
  sigma_ip <- sigma * y_age
  sigma_is <- sigma * (1 - y_age)
  gamma_p <- c("5-9" = 0.25, "0-4" = 0.2)
  gamma_c <- 0.1
  gamma_s <- c("0-4" = 0.3, "5-9" = 0.4)
  transitions <- data.frame(
    from = c("E", "E", "IP", "IC", "IS"),
    to = c("IP", "IS", "IC", "R", "R"),
    stringsAsFactors = FALSE
  )
  transitions$rate <- I(list(sigma_ip, sigma_is, gamma_p, gamma_c, gamma_s))

  CompartmentModel(
    compartments = c("S", "E", "IP", "IC", "IS", "R"),
    infection_transitions = data.frame(from = "S", to = "E", stringsAsFactors = FALSE),
    transitions = transitions,
    infectious_compartments = c("IP", "IC", "IS"),
    infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
  )
}

test_that("CompartmentModel constructs a valid generic disease model", {
  model <- generic_seir_model()

  expect_s3_class(model, "DiseaseModel")
  expect_identical(model$model_type, "CompartmentModel")
  expect_identical(model$compartments, c("S", "E", "I", "R"))
  expect_identical(model$birth_compartment, "S")
  expect_identical(model$migration_compartment, "S")
  expect_silent(validate_disease_model(model))
})

test_that("CompartmentModel validates compartment names and transition definitions", {
  expect_error(
    CompartmentModel(
      compartments = c("S", "I", "I"),
      infection_transitions = data.frame(from = "S", to = "I")
    ),
    "compartments must be unique"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "I"),
      infection_transitions = data.frame(from = "S", stringsAsFactors = FALSE)
    ),
    "infection_transitions is missing required column"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "I"),
      transitions = data.frame(from = "I", to = "R", rate = 0.2)
    ),
    "unknown destination compartment"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "I"),
      transitions = data.frame(from = "I", to = "S", rate = -0.2)
    ),
    "cannot contain negative"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "I"),
      infection_transitions = data.frame(from = "S", to = "I"),
      infectious_compartments = "E"
    ),
    "infectious_compartments contains unknown compartment"
  )
})

test_that("CompartmentModel supports explicit outflows", {
  model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    outflows = data.frame(from = "I", rate = 0.05, stringsAsFactors = FALSE),
    infectious_compartments = "I"
  )

  expect_identical(model$transitions$from, c("I", "I"))
  expect_identical(model$transitions$to, c("R", NA_character_))
  expect_identical(model$transitions$transition_type, c("internal", "outflow"))
  expect_identical(model$transitions$transition_id, c("transition:I->R", "outflow:I"))
  expect_silent(validate_disease_model(model))
})

test_that("generic SIR transition rates match SIRModel transition rates", {
  ages <- generic_test_ages()
  expected <- transition_rates(
    state = generic_sir_state(),
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )
  observed <- transition_rates(
    state = generic_sir_state(),
    model = generic_sir_model(gamma = 0.2),
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )

  expect_equal(observed, expected)
})

test_that("generic SEIR transition rates and derivatives match SEIRModel", {
  ages <- generic_test_ages()
  expected_rates <- transition_rates(
    state = generic_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )
  observed_rates <- transition_rates(
    state = generic_seir_state(),
    model = generic_seir_model(sigma = 0.3, gamma = 0.2),
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )

  expect_equal(observed_rates, expected_rates)
  expect_equal(
    rates_to_derivative(observed_rates, c("S", "E", "I", "R"), ages),
    rates_to_derivative(expected_rates, c("S", "E", "I", "R"), ages)
  )
})

test_that("generic SIR trajectory matches existing one-age-group SIR", {
  ages <- AgeStructure("all", 0, Inf)
  contacts <- matrix(2, nrow = 1, ncol = 1)
  initial_state <- data.frame(
    compartment = c("S", "I", "R"),
    age_group = "all",
    value = c(99, 1, 0),
    stringsAsFactors = FALSE
  )

  expected <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 0.25, 0.5),
    model = SIRModel(gamma = 0.1),
    age_structure = ages,
    contact_matrix = contacts,
    beta = 0.2,
    method = "euler"
  )
  observed <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 0.25, 0.5),
    model = generic_sir_model(gamma = 0.1),
    age_structure = ages,
    contact_matrix = contacts,
    beta = 0.2,
    method = "euler"
  )

  expect_equal(observed, expected)
})

test_that("generic model supports age-specific per-capita rates", {
  ages <- generic_test_ages()
  transitions <- data.frame(from = "I", to = "R", stringsAsFactors = FALSE)
  transitions$rate <- I(list(c("5-9" = 0.3, "0-4" = 0.1)))
  model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    transitions = transitions
  )

  rates <- transition_rates(
    state = generic_sir_state(),
    model = model,
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )

  expect_equal(rates$rate, c(1, 6))
})

test_that("generic model represents competing age-specific transitions from one source", {
  ages <- generic_test_ages()
  transitions <- data.frame(
    from = c("E", "E"),
    to = c("IP", "IS"),
    stringsAsFactors = FALSE
  )
  transitions$rate <- I(list(
    c("0-4" = 0.2, "5-9" = 0.4),
    c("5-9" = 0.1, "0-4" = 0.3)
  ))
  model <- CompartmentModel(
    compartments = c("S", "E", "IP", "IS", "R"),
    transitions = transitions,
    infectious_compartments = character()
  )
  state <- data.frame(
    compartment = rep(c("S", "E", "IP", "IS", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 5),
    value = c(90, 180, 5, 15, 0, 0, 0, 0, 0, 0),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )

  expect_equal(rates$from, c("E", "E", "E", "E"))
  expect_equal(rates$to, c("IP", "IS", "IP", "IS"))
  expect_equal(rates$age_group, c("0-4", "0-4", "5-9", "5-9"))
  expect_equal(rates$rate, c(1, 1.5, 6, 1.5))

  derivative <- rates_to_derivative(rates, model$compartments, ages)
  expect_equal(
    derivative$derivative,
    c(0, 0, -2.5, -7.5, 1, 6, 1.5, 1.5, 0, 0)
  )
})

test_that("generic clinical/subclinical fixture is expressible with age-specific rates", {
  ages <- generic_test_ages()
  model <- generic_clinical_subclinical_model(ages)
  state <- data.frame(
    compartment = rep(c("S", "E", "IP", "IC", "IS", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 6),
    value = c(90, 180, 10, 20, 5, 10, 2, 4, 3, 6, 0, 0),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = generic_test_contacts(),
    beta = 0
  )

  expect_equal(model$infectiousness_weights, c(IP = 1, IC = 1, IS = 0.5))
  e_rates <- rates[rates$from == "E", c("to", "age_group", "rate")]
  row.names(e_rates) <- NULL
  expect_equal(
    e_rates,
    data.frame(
      to = c("IP", "IS", "IP", "IS"),
      age_group = c("0-4", "0-4", "5-9", "5-9"),
      rate = c(2, 3, 7, 3),
      stringsAsFactors = FALSE
    )
  )
})

test_that("generic age-specific transition-rate vectors are validated at expansion time", {
  ages <- generic_test_ages()
  make_model <- function(rate) {
    transitions <- data.frame(from = "I", to = "R", stringsAsFactors = FALSE)
    transitions$rate <- I(list(rate))
    CompartmentModel(
      compartments = c("S", "I", "R"),
      transitions = transitions,
      infectious_compartments = character()
    )
  }
  expect_bad_rate <- function(rate, pattern) {
    expect_error(
      transition_rates(
        state = generic_sir_state(),
        model = make_model(rate),
        age_structure = ages,
        contact_matrix = generic_test_contacts()
      ),
      pattern
    )
  }

  expect_bad_rate(c(0.1, 0.3), "must be named by age group")
  expect_bad_rate(c("0-4" = 0.1, older = 0.3), "unknown age_group")
  expect_bad_rate(c("0-4" = 0.1, "10-14" = 0.3), "unknown age_group")
  expect_bad_rate(c("0-4" = 0.1, "5-9" = -0.3), "cannot contain negative")
  expect_bad_rate(c("0-4" = 0.1, "5-9" = Inf), "finite non-missing")
  expect_bad_rate(c("0-4" = 0.1, "5-9" = 0.2, "10-14" = 0.3), "length must be 1 or match")
})

test_that("generic model can use arbitrary infectious compartments", {
  ages <- generic_test_ages()
  model <- CompartmentModel(
    compartments = c("S", "A", "R"),
    infection_transitions = data.frame(from = "S", to = "A", stringsAsFactors = FALSE),
    transitions = data.frame(from = "A", to = "R", rate = 0, stringsAsFactors = FALSE),
    infectious_compartments = "A"
  )
  state <- data.frame(
    compartment = rep(c("S", "A", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 3),
    value = c(90, 180, 10, 20, 0, 0),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )

  expect_equal(rates$rate, c(27, 0, 126, 0))
})

test_that("generic model weights multiple infectious compartments", {
  ages <- generic_test_ages()
  model <- CompartmentModel(
    compartments = c("S", "IP", "IC", "R"),
    infection_transitions = data.frame(from = "S", to = "IP", stringsAsFactors = FALSE),
    transitions = data.frame(
      from = c("IP", "IC"),
      to = c("IC", "R"),
      rate = c(0, 0),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = c("IP", "IC"),
    infectiousness_weights = c(IC = 0.5, IP = 1)
  )
  state <- data.frame(
    compartment = rep(c("S", "IP", "IC", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 4),
    value = c(70, 140, 10, 20, 20, 40, 0, 0),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )

  expect_equal(model$infectiousness_weights, c(IP = 1, IC = 0.5))
  expect_equal(rates$rate, c(42, 0, 0, 196, 0, 0))
})

test_that("generic model supports COVID-like relative infectiousness weights", {
  ages <- generic_test_ages()
  model <- CompartmentModel(
    compartments = c("S", "IP", "IC", "IS", "R"),
    infection_transitions = data.frame(from = "S", to = "IP", stringsAsFactors = FALSE),
    transitions = data.frame(
      from = c("IP", "IC", "IS"),
      to = c("IC", "R", "R"),
      rate = c(0, 0, 0),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = c("IP", "IC", "IS"),
    infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
  )
  state <- data.frame(
    compartment = rep(c("S", "IP", "IC", "IS", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 5),
    value = c(90, 180, 10, 20, 5, 10, 8, 16, 0, 0),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = generic_test_contacts()
  )

  expect_equal(rates$rate, c(45.398230088496, 0, 0, 0, 211.858407079646, 0, 0, 0))
})

test_that("generic model rejects malformed infectiousness weights", {
  expect_error(
    CompartmentModel(
      compartments = c("S", "IP", "IC"),
      infection_transitions = data.frame(from = "S", to = "IP"),
      infectious_compartments = c("IP", "IC"),
      infectiousness_weights = c(1, 0.5)
    ),
    "must be named"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "IP", "IC"),
      infection_transitions = data.frame(from = "S", to = "IP"),
      infectious_compartments = c("IP", "IC"),
      infectiousness_weights = c(IP = 1, IS = 0.5)
    ),
    "unknown infectious compartment"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "IP", "IC"),
      infection_transitions = data.frame(from = "S", to = "IP"),
      infectious_compartments = c("IP", "IC"),
      infectiousness_weights = c(IP = -1, IC = 0.5)
    ),
    "cannot contain negative"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "IP", "IC"),
      infection_transitions = data.frame(from = "S", to = "IP"),
      infectious_compartments = c("IP", "IC"),
      infectiousness_weights = c(IP = Inf, IC = 0.5)
    ),
    "finite numeric vector"
  )
  expect_error(
    CompartmentModel(
      compartments = c("S", "IP", "IC"),
      infection_transitions = data.frame(from = "S", to = "IP"),
      infectious_compartments = c("IP", "IC"),
      infectiousness_weights = c(IP = 0, IC = 0)
    ),
    "at least one positive"
  )
})

test_that("generic deterministic simulation supports multiple infectious compartments", {
  ages <- AgeStructure("all", 0, Inf)
  model <- CompartmentModel(
    compartments = c("S", "IP", "IS", "R"),
    infection_transitions = data.frame(from = "S", to = "IP", stringsAsFactors = FALSE),
    transitions = data.frame(
      from = c("IP", "IS"),
      to = c("IS", "R"),
      rate = c(0.2, 0.1),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = c("IP", "IS"),
    infectiousness_weights = c(IP = 1, IS = 0.5)
  )
  state <- data.frame(
    compartment = c("S", "IP", "IS", "R"),
    age_group = "all",
    value = c(90, 5, 10, 0),
    stringsAsFactors = FALSE
  )

  output <- simulate_deterministic(
    initial_state = state,
    times = c(0, 0.1),
    model = model,
    age_structure = ages,
    contact_matrix = matrix(2, nrow = 1),
    beta = 0.2,
    method = "euler"
  )

  expect_equal(output$value[output$time == 0.1], c(89.657142857143, 5.242857142857, 10, 0.1))
})

test_that("generic contact orientation remains recipient-row source-column", {
  ages <- generic_test_ages()
  state <- generic_sir_state(S = c(100, 100), I = c(0, 10), R = c(0, 0))
  contacts <- matrix(c(
    0, 5,
    0, 0
  ), nrow = 2, byrow = TRUE)

  rates <- transition_rates(
    state = state,
    model = generic_sir_model(gamma = 0),
    age_structure = ages,
    contact_matrix = contacts
  )

  expect_equal(rates$rate, c(500 / 11, 0, 0, 0))
})

test_that("generic demographic coupling births and migration enter configured compartments", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5+"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, Inf)
  )
  ageing <- AgeingOperator(ages)
  ageing$departure_rate[] <- 0
  fertility <- FertilitySchedule(
    data.frame(time = 0, age_group = "5+", fertility_rate = 0.1),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(time = 0, age_group = ages$age_groups, migration_count = c(2, 3)),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    fertility_schedule = fertility,
    migration_schedule = migration,
    mode = "migration"
  )
  model <- CompartmentModel(
    compartments = c("U", "I"),
    transitions = data.frame(from = "I", to = "U", rate = 0),
    birth_compartment = "U",
    migration_compartment = "U"
  )
  initial_state <- data.frame(
    compartment = rep(c("U", "I"), each = 2),
    age_group = rep(ages$age_groups, times = 2),
    value = c(10, 20, 1, 2),
    stringsAsFactors = FALSE
  )

  output <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 1),
    model = model,
    age_structure = ages,
    contact_matrix = diag(2),
    demographic_process = process,
    time_policy = "exact",
    method = "euler"
  )

  expect_equal(output$value[output$time == 1], c(14.2, 23, 1, 2))
})
