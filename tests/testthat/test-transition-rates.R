test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

test_state <- function(S = c(90, 180), I = c(10, 20), R = c(0, 0)) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

test_seir_state <- function(
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

test_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

test_generic_vaccine_model <- function(susceptibility = list(1, 0.5)) {
  infection_transitions <- data.frame(
    from = c("S", "V"),
    to = c("I", "I"),
    susceptibility = if (is.list(susceptibility)) I(susceptibility) else susceptibility,
    stringsAsFactors = FALSE
  )

  CompartmentModel(
    compartments = c("S", "V", "I", "E", "R"),
    infection_transitions = infection_transitions,
    transitions = data.frame(
      from = c("I", "E"),
      to = c("R", "R"),
      rate = c(0.2, 0.1),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = "I"
  )
}

test_generic_vaccine_state <- function(
  S = c(100, 100),
  V = c(100, 100),
  I = c(10, 10),
  E = c(0, 0),
  R = c(0, 0)
) {
  data.frame(
    compartment = rep(c("S", "V", "I", "E", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 5),
    value = c(S, V, I, E, R),
    stringsAsFactors = FALSE
  )
}

test_that("transition_rates computes a manually checkable SIR example", {
  rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expected <- data.frame(
    from = c("S", "I", "S", "I"),
    to = c("I", "R", "I", "R"),
    age_group = c("0-4", "0-4", "5-9", "5-9"),
    rate = c(27, 2, 126, 4),
    transition_id = c(
      "infection:S->I",
      "transition:I->R",
      "infection:S->I",
      "transition:I->R"
    ),
    transition_label = rep(NA_character_, 4),
    transition_type = c("infection", "transition", "infection", "transition"),
    stringsAsFactors = FALSE
  )

  expect_equal(rates, expected)
})

test_that("transition_rates computes SEIR infection, progression, and recovery rates", {
  rates <- transition_rates(
    state = test_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expected <- data.frame(
    from = c("S", "E", "I", "S", "E", "I"),
    to = c("E", "I", "R", "E", "I", "R"),
    age_group = c("0-4", "0-4", "0-4", "5-9", "5-9", "5-9"),
    rate = c(25.514950166113, 1.5, 2, 118.405315614618, 4.5, 4),
    transition_id = c(
      "infection:S->E",
      "transition:E->I",
      "transition:I->R",
      "infection:S->E",
      "transition:E->I",
      "transition:I->R"
    ),
    transition_label = rep(NA_character_, 6),
    transition_type = c("infection", "transition", "transition", "infection", "transition", "transition"),
    stringsAsFactors = FALSE
  )

  expect_equal(rates, expected)
})

test_that("SIR and SEIR constructors match explicit unit susceptibility and infectiousness", {
  ages <- test_ages()
  contacts <- test_contacts()

  sir_rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.2, beta = 0.4),
    age_structure = ages,
    contact_matrix = contacts
  )
  generic_sir_rates <- transition_rates(
    state = test_state(),
    model = CompartmentModel(
      compartments = c("S", "I", "R"),
      infection_transitions = data.frame(
        from = "S",
        to = "I",
        susceptibility = 1,
        stringsAsFactors = FALSE
      ),
      transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
      infectious_compartments = "I",
      infectiousness_weights = 1,
      beta = 0.4
    ),
    age_structure = ages,
    contact_matrix = contacts
  )

  seir_rates <- transition_rates(
    state = test_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2, beta = 0.4),
    age_structure = ages,
    contact_matrix = contacts
  )
  generic_seir_rates <- transition_rates(
    state = test_seir_state(),
    model = CompartmentModel(
      compartments = c("S", "E", "I", "R"),
      infection_transitions = data.frame(
        from = "S",
        to = "E",
        susceptibility = 1,
        stringsAsFactors = FALSE
      ),
      transitions = data.frame(
        from = c("E", "I"),
        to = c("I", "R"),
        rate = c(0.3, 0.2),
        stringsAsFactors = FALSE
      ),
      infectious_compartments = "I",
      infectiousness_weights = 1,
      beta = 0.4
    ),
    age_structure = ages,
    contact_matrix = contacts
  )

  expect_equal(sir_rates, generic_sir_rates)
  expect_equal(seir_rates, generic_seir_rates)
})

test_that("transition_rates honors explicit beta overrides and ignores beta for non-infectious models", {
  ages <- test_ages()
  contacts <- test_contacts()
  infectious_model <- SIRModel(gamma = 0.2, beta = 0.4)

  default_rates <- transition_rates(
    state = test_state(),
    model = infectious_model,
    age_structure = ages,
    contact_matrix = contacts
  )
  zero_override <- transition_rates(
    state = test_state(),
    model = infectious_model,
    age_structure = ages,
    contact_matrix = contacts,
    beta = 0
  )
  other_override <- transition_rates(
    state = test_state(),
    model = infectious_model,
    age_structure = ages,
    contact_matrix = contacts,
    beta = 0.1
  )
  no_model_beta <- infectious_model
  no_model_beta$beta <- NULL

  expect_true(all(zero_override$rate[zero_override$transition_id == "infection:S->I"] == 0))
  expect_equal(
    zero_override$rate[zero_override$transition_id == "transition:I->R"],
    default_rates$rate[default_rates$transition_id == "transition:I->R"]
  )
  expect_equal(
    other_override$rate[other_override$transition_id == "infection:S->I"],
    default_rates$rate[default_rates$transition_id == "infection:S->I"] * 0.25
  )
  expect_error(
    transition_rates(
      state = test_state(),
      model = no_model_beta,
      age_structure = ages,
      contact_matrix = contacts
    ),
    "beta is required"
  )

  non_infectious_model <- CompartmentModel(
    compartments = c("S", "R"),
    transitions = data.frame(from = "S", to = "R", rate = 0.1, stringsAsFactors = FALSE),
    infectious_compartments = character()
  )
  non_infectious_state <- data.frame(
    compartment = rep(c("S", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 2),
    value = c(100, 80, 0, 0),
    stringsAsFactors = FALSE
  )
  non_infectious_default <- transition_rates(
    state = non_infectious_state,
    model = non_infectious_model,
    age_structure = ages,
    contact_matrix = contacts
  )
  non_infectious_override <- transition_rates(
    state = non_infectious_state,
    model = non_infectious_model,
    age_structure = ages,
    contact_matrix = contacts,
    beta = 0.7
  )

  expect_equal(non_infectious_default, non_infectious_override)
  expect_identical(infectious_model$beta, 0.4)
})

test_that("SEIR force of infection depends on I rather than E", {
  rates_without_exposed <- transition_rates(
    state = test_seir_state(E = c(0, 0), I = c(10, 20), R = c(5, 15)),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )
  rates_with_exposed <- transition_rates(
    state = test_seir_state(E = c(5, 15), I = c(10, 20), R = c(0, 0)),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )
  rates_without_infectious <- transition_rates(
    state = test_seir_state(E = c(500, 1000), I = c(0, 0)),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_true(all(rates_without_exposed$rate[rates_without_exposed$from == "S"] > 0))
  expect_equal(
    rates_with_exposed$rate[rates_with_exposed$from == "S"],
    rates_without_exposed$rate[rates_without_exposed$from == "S"]
  )
  expect_equal(rates_without_infectious$rate[rates_without_infectious$from == "S"], c(0, 0))
})

test_that("transition_rates accepts compartment-major state vectors", {
  ages <- test_ages()
  state_vector <- state_long_to_vector(test_state(), ages, c("S", "I", "R"))

  rates <- transition_rates(
    state = state_vector,
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = test_contacts()
  )

  expect_equal(rates$rate, c(27, 2, 126, 4))
})

test_that("transition_rates ignores state vector names intentionally", {
  ages <- test_ages()
  state_vector <- c(
    wrong_1 = 90, wrong_2 = 180,
    wrong_3 = 10, wrong_4 = 20,
    wrong_5 = 0, wrong_6 = 0
  )

  rates <- transition_rates(
    state = state_vector,
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = test_contacts()
  )

  expect_equal(rates$rate, c(27, 2, 126, 4))
})

test_that("transition_rates rejects state vectors with wrong length", {
  expect_error(
    transition_rates(c(90, 180, 10), SIRModel(0.2), test_ages(), test_contacts()),
    "state_vector length"
  )
})

test_that("beta scaling affects only infection rates", {
  rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts(),
    beta = 2
  )

  expect_equal(rates$rate, c(54, 2, 252, 4))
})

test_that("gamma scaling affects only recovery rates", {
  rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.4),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_equal(rates$rate, c(27, 4, 126, 8))
})

test_that("susceptibility modifies infection rates by recipient age group", {
  rates <- transition_rates(
    state = test_state(),
    model = CompartmentModel(
      compartments = c("S", "I", "R"),
      infection_transitions = data.frame(
        from = "S",
        to = "I",
        susceptibility = I(list(c(0.5, 2))),
        stringsAsFactors = FALSE
      ),
      transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
      infectious_compartments = "I"
    ),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_equal(rates$rate, c(13.5, 2, 252, 4))
})

test_that("infectiousness modifies infection pressure from source age groups", {
  rates <- transition_rates(
    state = test_state(),
    model = CompartmentModel(
      compartments = c("S", "I", "R"),
      infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
      transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
      infectious_compartments = "I",
      infectiousness_weights = list(c(2, 0.5))
    ),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_equal(rates$rate, c(40.5, 2, 144, 4))
})

test_that("scalar infectiousness weights match scalar list weights and explicit age replication", {
  ages <- test_ages()
  state <- test_state()
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

  scalar_rates <- transition_rates(state, scalar_model, ages, diag(2))
  list_rates <- transition_rates(state, list_model, ages, diag(2))
  expanded_rates <- transition_rates(state, expanded_model, ages, diag(2))

  expect_equal(list_rates, scalar_rates)
  expect_equal(expanded_rates, scalar_rates)
})

test_that("named infectiousness weights are reordered by infectious compartment and age group", {
  ages <- test_ages()
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
    infectiousness_weights = list(
      IC = c("5-9" = 0.5, "0-4" = 0.25),
      IP = c("5-9" = 1, "0-4" = 2)
    )
  )
  state <- data.frame(
    compartment = rep(c("S", "IP", "IC", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 4),
    value = c(90, 180, 10, 20, 30, 40, 0, 0),
    stringsAsFactors = FALSE
  )

  expect_identical(names(model$infectiousness_weights), c("IP", "IC"))

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = diag(2)
  )

  effective_infectious <- c(
    "0-4" = 2 * 10 + 0.25 * 30,
    "5-9" = 1 * 20 + 0.5 * 40
  )
  population <- c("0-4" = 90 + 10 + 30, "5-9" = 180 + 20 + 40)
  expected_infection <- c(
    90 * effective_infectious["0-4"] / population["0-4"],
    180 * effective_infectious["5-9"] / population["5-9"]
  )

  expect_equal(rates$rate[rates$from == "S"], as.numeric(expected_infection))
})

test_that("named infectiousness weights are reordered across infectious compartments and age groups", {
  ages <- AgeStructure(
    age_groups = c("child", "adult"),
    lower_bounds = c(0, 18),
    upper_bounds = c(17, 99)
  )

  # Both the infectious-compartment order and the age-group order are
  # intentionally scrambled here to exercise simultaneous reordering.
  model <- CompartmentModel(
    compartments = c("S", "E", "I1", "I2"),
    infection_transitions = data.frame(from = "S", to = "E", stringsAsFactors = FALSE),
    transitions = data.frame(
      from = character(),
      to = character(),
      rate = numeric(),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = c("I1", "I2"),
    infectiousness_weights = list(
      I2 = c(adult = 1.0, child = 0.5),
      I1 = c(adult = 0.4, child = 0.2)
    )
  )
  state <- data.frame(
    compartment = rep(c("S", "E", "I1", "I2"), each = 2),
    age_group = rep(ages$age_groups, times = 4),
    value = c(140, 160, 0, 0, 10, 20, 30, 40),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = diag(2)
  )

  expected_effective_infectious <- c(child = 17, adult = 48)
  expected_rates <- c(
    child = 140 * expected_effective_infectious["child"] / 180,
    adult = 160 * expected_effective_infectious["adult"] / 220
  )

  expect_identical(names(model$infectiousness_weights), c("I1", "I2"))
  expect_equal(rates$rate[rates$from == "S"], as.numeric(expected_rates))
})

test_that("zero infectiousness removes transmission and weights above one scale rates", {
  ages <- test_ages()
  state <- test_state()
  zero_model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    infectious_compartments = "I",
    infectiousness_weights = 0
  )
  scaled_model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    infectious_compartments = "I",
    infectiousness_weights = 2.5
  )

  zero_rates <- transition_rates(state, zero_model, ages, diag(2))
  scaled_rates <- transition_rates(state, scaled_model, ages, diag(2))

  expect_equal(zero_rates$rate[zero_rates$from == "S"], c(0, 0))
  expect_equal(scaled_rates$rate[scaled_rates$from == "S"], c(22.5, 45))
})

test_that("age-specific susceptibility and infectiousness combine in the infection-flow formula", {
  ages <- test_ages()
  model <- CompartmentModel(
    compartments = c("S", "V", "IP", "IC", "R"),
    infection_transitions = data.frame(
      from = c("S", "V"),
      to = c("IP", "IP"),
      susceptibility = I(list(
        c("0-4" = 0.5, "5-9" = 2),
        c("0-4" = 1, "5-9" = 0.25)
      )),
      stringsAsFactors = FALSE
    ),
    transitions = data.frame(
      from = c("IP", "IC"),
      to = c("R", "R"),
      rate = c(0, 0),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = c("IP", "IC"),
    infectiousness_weights = list(
      c("0-4" = 2, "5-9" = 1),
      c("0-4" = 0.25, "5-9" = 0.5)
    )
  )
  state <- data.frame(
    compartment = rep(c("S", "V", "IP", "IC", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 5),
    value = c(50, 100, 40, 80, 10, 20, 30, 40, 0, 0),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = matrix(1, nrow = 2, ncol = 2)
  )

  effective_infectious <- c(
    "0-4" = 2 * 10 + 0.25 * 30,
    "5-9" = 1 * 20 + 0.5 * 40
  )
  population <- c("0-4" = 50 + 40 + 10 + 30, "5-9" = 100 + 80 + 20 + 40)
  lambda <- sum(effective_infectious / population)
  expected <- c(
    50 * 0.5 * lambda,
    40 * 1 * lambda,
    100 * 2 * lambda,
    80 * 0.25 * lambda
  )

  expect_equal(rates$rate[rates$from %in% c("S", "V")], as.numeric(expected))
})

test_that("duplicate age-group names in infectiousness weights are rejected at expansion time", {
  ages <- test_ages()
  model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0, stringsAsFactors = FALSE),
    infectious_compartments = "I",
    infectiousness_weights = list(c("0-4" = 0.2, "0-4" = 0.4))
  )

  expect_error(
    transition_rates(
      state = test_state(),
      model = model,
      age_structure = ages,
      contact_matrix = test_contacts()
    ),
    "duplicate age group"
  )
})

test_that("legacy infection transitions default susceptibility to one", {
  model_without_susceptibility <- CompartmentModel(
    compartments = c("S", "V", "I", "E", "R"),
    infection_transitions = data.frame(
      from = c("S", "V"),
      to = c("I", "I"),
      stringsAsFactors = FALSE
    ),
    transitions = data.frame(
      from = c("I", "E"),
      to = c("R", "R"),
      rate = c(0.2, 0.1),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = "I"
  )

  model_with_scalar_one <- test_generic_vaccine_model(list(1, 1))

  without_susceptibility <- transition_rates(
    state = test_generic_vaccine_state(),
    model = model_without_susceptibility,
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )
  with_scalar_one <- transition_rates(
    state = test_generic_vaccine_state(),
    model = model_with_scalar_one,
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )

  expect_equal(without_susceptibility, with_scalar_one)
})

test_that("scalar susceptibility is expanded across all age groups", {
  rates <- transition_rates(
    state = test_generic_vaccine_state(),
    model = test_generic_vaccine_model(0.5),
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )

  expect_equal(
    rates$rate[rates$from == "S"],
    c(2.380952380952, 2.380952380952)
  )
  expect_equal(
    rates$rate[rates$from == "V"],
    c(2.380952380952, 2.380952380952)
  )
})

test_that("named age-specific susceptibility vectors are reordered by age group", {
  rates <- transition_rates(
    state = test_generic_vaccine_state(),
    model = test_generic_vaccine_model(list(c("5-9" = 0.2, "0-4" = 0.5), c("0-4" = 1, "5-9" = 1))),
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )

  expect_equal(
    rates$rate[rates$from == "S"],
    c(2.380952380952, 0.952380952381)
  )
  expect_equal(
    rates$rate[rates$from == "V"],
    c(4.761904761905, 4.761904761905)
  )
})

test_that("multiple susceptible compartments can target different infection destinations", {
  model <- CompartmentModel(
    compartments = c("S", "V", "I", "E", "R"),
    infection_transitions = {
      data.frame(
        from = c("S", "V"),
        to = c("I", "E"),
        susceptibility = I(list(1, 0.5)),
        stringsAsFactors = FALSE
      )
    },
    transitions = data.frame(
      from = c("I", "E"),
      to = c("R", "R"),
      rate = c(0.2, 0.1),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = "I"
  )

  rates <- transition_rates(
    state = test_generic_vaccine_state(),
    model = model,
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )

  expect_equal(
    paste(rates$from, rates$to, sep = "->"),
    c("S->I", "V->E", "I->R", "E->R", "S->I", "V->E", "I->R", "E->R")
  )
  expect_equal(rates$rate[rates$from == "S"], c(4.761904761905, 4.761904761905))
  expect_equal(rates$rate[rates$from == "V"], c(2.380952380952, 2.380952380952))
})

test_that("zero susceptibility produces no infections", {
  rates <- transition_rates(
    state = test_generic_vaccine_state(),
    model = test_generic_vaccine_model(list(1, 0)),
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )

  expect_equal(rates$rate[rates$from == "V"], c(0, 0))
})

test_that("susceptibility values greater than one scale incidence", {
  base_rates <- transition_rates(
    state = test_generic_vaccine_state(),
    model = test_generic_vaccine_model(list(1, 1)),
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )
  scaled_rates <- transition_rates(
    state = test_generic_vaccine_state(),
    model = test_generic_vaccine_model(list(2.5, 1)),
    age_structure = test_ages(),
    contact_matrix = diag(2),
    beta = 1
  )

  expect_equal(
    scaled_rates$rate[scaled_rates$from == "S"],
    base_rates$rate[base_rates$from == "S"] * 2.5
  )
})

test_that("infection transition susceptibility validation is strict", {
  expect_error(
    CompartmentModel(
      compartments = c("S", "V", "I", "E", "R"),
      infection_transitions = {
        data.frame(
          from = c("S", "V"),
          to = c("I", "I"),
          susceptibility = I(list(1, c("0-4" = -0.1, "5-9" = 0.2))),
          stringsAsFactors = FALSE
        )
      },
      transitions = data.frame(
        from = c("I", "E"),
        to = c("R", "R"),
        rate = c(0.2, 0.1),
        stringsAsFactors = FALSE
      ),
      infectious_compartments = "I"
    ),
    "cannot contain negative values"
  )

  expect_error(
    transition_rates(
      state = test_generic_vaccine_state(),
      model = test_generic_vaccine_model(list(1, c(0.1, 0.2, 0.3))),
      age_structure = test_ages(),
      contact_matrix = diag(2),
      beta = 1
    ),
    "length must be 1 or match the number of age groups"
  )

  expect_error(
    transition_rates(
      state = test_generic_vaccine_state(),
      model = test_generic_vaccine_model(list(c("0-4" = 1, "unknown" = 2), 1)),
      age_structure = test_ages(),
      contact_matrix = diag(2),
      beta = 1
    ),
    "unknown age_group value"
  )

  expect_error(
    transition_rates(
      state = test_generic_vaccine_state(),
      model = test_generic_vaccine_model(list(c("0-4" = 1, 2), 1)),
      age_structure = test_ages(),
      contact_matrix = diag(2),
      beta = 1
    ),
    "must be non-empty"
  )

  expect_error(
    CompartmentModel(
      compartments = c("S", "V", "I", "E", "R"),
      infection_transitions = {
        data.frame(
          from = c("S", "V"),
          to = c("I", "I"),
          susceptibility = I(list(1, "bad")),
          stringsAsFactors = FALSE
        )
      },
      transitions = data.frame(
        from = c("I", "E"),
        to = c("R", "R"),
        rate = c(0.2, 0.1),
        stringsAsFactors = FALSE
      ),
      infectious_compartments = "I"
    ),
    "must be finite numeric value"
  )

  expect_error(
    CompartmentModel(
      compartments = c("S", "V", "I", "E", "R"),
      infection_transitions = data.frame(
        from = c("S", "S"),
        to = c("I", "I"),
        stringsAsFactors = FALSE
      ),
      transitions = data.frame(
        from = c("I", "E"),
        to = c("R", "R"),
        rate = c(0.2, 0.1),
        stringsAsFactors = FALSE
      ),
      infectious_compartments = "I"
    ),
    "duplicate transition"
  )
})

test_that("zero infectious counts produce zero infection and zero recovery rates", {
  rates <- transition_rates(
    state = test_state(I = c(0, 0)),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_equal(rates$rate, c(0, 0, 0, 0))
})

test_that("recovery rates still depend on infectious counts when some groups have no infection pressure", {
  rates <- transition_rates(
    state = test_state(I = c(0, 20)),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = matrix(c(
      1, 0,
      0, 0
    ), nrow = 2, byrow = TRUE)
  )

  expect_equal(rates$rate, c(0, 0, 0, 4))
})

test_that("zero susceptible counts produce zero infection rates", {
  rates <- transition_rates(
    state = test_state(S = c(0, 0)),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_equal(rates$rate, c(0, 2, 0, 4))
})

test_that("transition_rates returns the expected output columns", {
  rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_identical(
    names(rates),
    c("from", "to", "age_group", "rate", "transition_id", "transition_label", "transition_type")
  )
})

test_that("generic infection transitions expand one row per age group", {
  ages <- test_ages()
  model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(
      from = "S",
      to = "I",
      susceptibility = I(list(c(0.5, 2))),
      stringsAsFactors = FALSE
    ),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    infectious_compartments = "I"
  )

  rates <- transition_rates(
    state = test_state(S = c(90, 180), I = c(10, 20), R = c(0, 0)),
    model = model,
    age_structure = ages,
    contact_matrix = test_contacts()
  )

  expect_equal(
    rates[, c("from", "to", "age_group", "transition_id", "transition_type")],
    data.frame(
      from = c("S", "I", "S", "I"),
      to = c("I", "R", "I", "R"),
      age_group = c("0-4", "0-4", "5-9", "5-9"),
      transition_id = c("infection:S->I", "transition:I->R", "infection:S->I", "transition:I->R"),
      transition_type = c("infection", "transition", "infection", "transition"),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(rates$rate[rates$transition_type == "transition"], c(2, 4))
  expect_true(all(rates$rate[rates$transition_type == "infection"] >= 0))
})

test_that("transition_rates includes stable logical transition identifiers", {
  first <- transition_rates(
    state = test_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )
  second <- transition_rates(
    state = test_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_identical(first$transition_id, second$transition_id)
  expect_identical(
    first$transition_id,
    c(
      "infection:S->E",
      "transition:E->I",
      "transition:I->R",
      "infection:S->E",
      "transition:E->I",
      "transition:I->R"
    )
  )
  expect_identical(
    unique(first$transition_id[first$from == "S" & first$to == "E"]),
    "infection:S->E"
  )
})

test_that("generic transition_ids distinguish infection and ordinary transitions", {
  model <- CompartmentModel(
    compartments = c("S", "E", "IP", "IS", "R"),
    infection_transitions = data.frame(from = "S", to = "E", stringsAsFactors = FALSE),
    transitions = data.frame(
      from = c("E", "E", "IP", "IS"),
      to = c("IP", "IS", "R", "R"),
      rate = c(0.2, 0.3, 0.1, 0.1),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = "IP"
  )
  state <- data.frame(
    compartment = rep(c("S", "E", "IP", "IS", "R"), each = 2),
    age_group = rep(test_ages()$age_groups, times = 5),
    value = c(90, 180, 5, 10, 2, 4, 1, 2, 0, 0),
    stringsAsFactors = FALSE
  )

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = test_ages(),
    contact_matrix = test_contacts(),
    beta = 0
  )

  expect_identical(
    unique(rates$transition_id),
    c(
      "infection:S->E",
      "transition:E->IP",
      "transition:E->IS",
      "transition:IP->R",
      "transition:IS->R"
    )
  )
  expect_identical(
    unique(rates$transition_id[rates$from == "E" & rates$to == "IP"]),
    "transition:E->IP"
  )
})

test_that("transition_rates orders age groups outermost and transitions innermost", {
  rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_identical(rates$age_group, c("0-4", "0-4", "5-9", "5-9"))
  expect_identical(paste(rates$from, rates$to, sep = "->"), c("S->I", "I->R", "S->I", "I->R"))
})

test_that("cumulative-flow validation accepts named list specifications", {
  rates <- transition_rates(
    state = test_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  flows <- validate_cumulative_flows(
    cumulative_flows = list(
      infections = list(from = "S", to = "E"),
      progressions = list(from = "E", to = "I")
    ),
    transition_rate_table = rates
  )

  expect_equal(
    flows[, c("cumulative_name", "transition_id", "from", "to")],
    data.frame(
      cumulative_name = c("infections", "progressions"),
      transition_id = c("infection:S->E", "transition:E->I"),
      from = c("S", "E"),
      to = c("E", "I"),
      stringsAsFactors = FALSE
    )
  )
  expect_true(all(is.na(flows$transition_label)))
  expect_equal(flows$transition_type, c("infection", "transition"))
})

test_that("cumulative-flow validation accepts mixed named list specifications", {
  rates <- data.frame(
    transition_id = c("infection:Mtb.naive->Incipient", "transition:Contained->Incipient"),
    from = c("Mtb.naive", "Contained"),
    to = c("Incipient", "Incipient"),
    age_group = c("0-4", "0-4"),
    rate = c(1, 2),
    stringsAsFactors = FALSE
  )

  flows <- validate_cumulative_flows(
    cumulative_flows = list(
      infections_other = list(from = c("Mtb.naive", "Contained"), to = c("Incipient", "Incipient")),
      infections_contained = list(transition_id = "infection:Mtb.naive->Incipient")
    ),
    transition_rate_table = rates
  )

  expect_equal(
    flows[, c("cumulative_name", "transition_id", "from", "to")],
    data.frame(
      cumulative_name = c("infections_other", "infections_other", "infections_contained"),
      transition_id = c(NA_character_, NA_character_, "infection:Mtb.naive->Incipient"),
      from = c("Mtb.naive", "Contained", "Mtb.naive"),
      to = c("Incipient", "Incipient", "Incipient"),
      stringsAsFactors = FALSE
    )
  )
})

test_that("cumulative-flow validation accepts explicit outflow selectors", {
  ages <- test_ages()
  model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    outflows = data.frame(from = "I", rate = 0.1, stringsAsFactors = FALSE),
    infectious_compartments = character()
  )
  rates <- transition_rates(
    state = data.frame(
      compartment = rep(c("S", "I", "R"), each = 2),
      age_group = rep(ages$age_groups, times = 3),
      value = c(0, 0, 10, 20, 0, 0),
      stringsAsFactors = FALSE
    ),
    model = model,
    age_structure = ages,
    contact_matrix = test_contacts()
  )

  flows <- validate_cumulative_flows(
    cumulative_flows = list(
      removals = list(from = "I", to = NA_character_)
    ),
    transition_rate_table = rates
  )

  expect_equal(
    flows[, c("cumulative_name", "transition_id", "from", "to")],
    data.frame(
      cumulative_name = "removals",
      transition_id = "outflow:I",
      from = "I",
      to = NA_character_,
      stringsAsFactors = FALSE
    )
  )
  expect_equal(flows$transition_type, "outflow")
})

test_that("cumulative-flow validation accepts multi-transition named list specifications", {
  rates <- data.frame(
    transition_id = c("infection:S->Lr", "transition:Lr->I", "transition:Ld->I", "transition:I->T"),
    from = c("S", "Lr", "Ld", "I"),
    to = c("Lr", "I", "I", "T"),
    age_group = "0-4",
    rate = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  flows <- validate_cumulative_flows(
    cumulative_flows = list(
      disease_onset = list(from = c("Lr", "Ld"), to = c("I", "I"))
    ),
    transition_rate_table = rates
  )

  expect_equal(
    flows[, c("cumulative_name", "transition_id", "from", "to")],
    data.frame(
      cumulative_name = c("disease_onset", "disease_onset"),
      transition_id = c("transition:Lr->I", "transition:Ld->I"),
      from = c("Lr", "Ld"),
      to = c("I", "I"),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(flows$selector_index, c(1L, 2L))
})

test_that("cumulative-flow validation accepts data-frame specifications", {
  rates <- transition_rates(
    state = test_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  flows <- validate_cumulative_flows(
    cumulative_flows = data.frame(
      name = c("infections", "recoveries"),
      from = c("S", "I"),
      to = c("E", "R"),
      stringsAsFactors = FALSE
    ),
    transition_rate_table = rates
  )

  expect_equal(flows$cumulative_name, c("infections", "recoveries"))
  expect_equal(flows$transition_id, c("infection:S->E", "transition:I->R"))
})

test_that("cumulative-flow validation rejects malformed specifications", {
  rates <- transition_rates(
    state = test_seir_state(),
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_error(
    validate_cumulative_flows(list(list(from = "S", to = "E")), rates),
    "list entries must be named"
  )
  expect_error(
    validate_cumulative_flows(
      list(infections = list(from = "S", to = "E"), infections = list(from = "E", to = "I")),
      rates
    ),
    "must be unique"
  )
  expect_equal(
    validate_cumulative_flows(list(infections = list(from = "S")), rates)$transition_id,
    "infection:S->E"
  )
  expect_error(
    validate_cumulative_flows(list(infections = list(from = "X", to = "E")), rates),
    "matched no transitions"
  )
  expect_error(
    validate_cumulative_flows(list(infections = list(from = "S", to = "X")), rates),
    "matched no transitions"
  )
  expect_error(
    validate_cumulative_flows(list(infections = list(from = "E", to = "R")), rates),
    "matched no transitions"
  )
  expect_error(
    validate_cumulative_flows(list(onsets = list(from = c("E", "I"), to = "R")), rates),
    "matched no transitions"
  )
  expect_error(
    validate_cumulative_flows(list(onsets = list(from = character(), to = character())), rates),
    "length 0"
  )
})

test_that("cumulative-flow validation rejects ambiguous from-to matches", {
  rates <- data.frame(
    transition_id = c("transition:E->I:a", "transition:E->I:b"),
    from = c("E", "E"),
    to = c("I", "I"),
    age_group = c("0-4", "0-4"),
    rate = c(1, 2),
    stringsAsFactors = FALSE
  )

  expect_error(
    validate_cumulative_flows(list(progressions = list(from = "E", to = "I")), rates),
    "matched multiple transitions"
  )
})

test_that("transition_rates rejects missing compartment-age combinations", {
  state <- test_state()
  state <- state[-1, ]

  expect_error(
    transition_rates(state, SIRModel(0.2), test_ages(), test_contacts()),
    "missing compartment-age row"
  )
})

test_that("transition_rates rejects duplicate compartment-age rows", {
  state <- rbind(test_state(), test_state()[1, ])

  expect_error(
    transition_rates(state, SIRModel(0.2), test_ages(), test_contacts()),
    "duplicate compartment-age rows"
  )
})

test_that("transition_rates rejects negative state values", {
  state <- test_state()
  state$value[1] <- -1

  expect_error(
    transition_rates(state, SIRModel(0.2), test_ages(), test_contacts()),
    "state values cannot be negative"
  )
})

test_that("transition_rates rejects missing and non-numeric long-form state values", {
  missing_state <- test_state()
  missing_state$value[1] <- NA_real_
  expect_error(
    transition_rates(missing_state, SIRModel(0.2), test_ages(), test_contacts()),
    "value must be numeric and cannot contain missing"
  )

  non_numeric_state <- test_state()
  non_numeric_state$value <- as.character(non_numeric_state$value)
  expect_error(
    transition_rates(non_numeric_state, SIRModel(0.2), test_ages(), test_contacts()),
    "value must be numeric"
  )
})

test_that("transition_rates rejects long-form state rows outside the model", {
  extra_compartment <- rbind(
    test_state(),
    data.frame(compartment = "E", age_group = "0-4", value = 1, stringsAsFactors = FALSE)
  )
  expect_error(
    transition_rates(extra_compartment, SIRModel(0.2), test_ages(), test_contacts()),
    "unknown compartment"
  )

  extra_age <- rbind(
    test_state(),
    data.frame(compartment = "S", age_group = "10-14", value = 1, stringsAsFactors = FALSE)
  )
  expect_error(
    transition_rates(extra_age, SIRModel(0.2), test_ages(), test_contacts()),
    "unknown age_group"
  )
})

test_that("transition_rates returns numeric non-negative rates", {
  rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts()
  )

  expect_type(rates$rate, "double")
  expect_true(all(rates$rate >= 0))
})

test_that("transition_rates rejects zero total population by age group", {
  expect_error(
    transition_rates(
      state = test_state(S = c(0, 180), I = c(0, 20), R = c(0, 0)),
      model = SIRModel(gamma = 0.2),
      age_structure = test_ages(),
      contact_matrix = test_contacts()
    ),
    "state population must be positive"
  )
})

test_that("transition_rates rejects unsupported model types", {
  model <- SIRModel(0.2)
  model$model_type <- "SIRS"

  expect_error(
    transition_rates(test_state(), model, test_ages(), test_contacts()),
    "unsupported disease model type"
  )
})

test_that("transition_rates rejects contact matrix dimension mismatch", {
  expect_error(
    transition_rates(test_state(), SIRModel(0.2), test_ages(), diag(3)),
    "contact_matrix dimensions"
  )
})
