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
    stringsAsFactors = FALSE
  )

  expect_equal(rates, expected)
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
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts(),
    susceptibility = c(0.5, 2)
  )

  expect_equal(rates$rate, c(13.5, 2, 252, 4))
})

test_that("infectiousness modifies infection pressure from source age groups", {
  rates <- transition_rates(
    state = test_state(),
    model = SIRModel(gamma = 0.2),
    age_structure = test_ages(),
    contact_matrix = test_contacts(),
    infectiousness = c(2, 0.5)
  )

  expect_equal(rates$rate, c(40.5, 2, 144, 4))
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

  expect_identical(names(rates), c("from", "to", "age_group", "rate", "transition_id"))
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
    flows,
    data.frame(
      cumulative_name = c("infections", "progressions"),
      transition_id = c("infection:S->E", "transition:E->I"),
      from = c("S", "E"),
      to = c("E", "I"),
      stringsAsFactors = FALSE
    )
  )
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
    flows,
    data.frame(
      cumulative_name = c("disease_onset", "disease_onset"),
      transition_id = c("transition:Lr->I", "transition:Ld->I"),
      from = c("Lr", "Ld"),
      to = c("I", "I"),
      stringsAsFactors = FALSE
    )
  )
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
  expect_error(
    validate_cumulative_flows(list(infections = list(from = "S")), rates),
    "must include from and to"
  )
  expect_error(
    validate_cumulative_flows(list(infections = list(from = "X", to = "E")), rates),
    "unknown source compartment"
  )
  expect_error(
    validate_cumulative_flows(list(infections = list(from = "S", to = "X")), rates),
    "unknown destination compartment"
  )
  expect_error(
    validate_cumulative_flows(list(infections = list(from = "E", to = "R")), rates),
    "does not match a declared transition"
  )
  expect_error(
    validate_cumulative_flows(list(onsets = list(from = c("E", "I"), to = "R")), rates),
    "same length"
  )
  expect_error(
    validate_cumulative_flows(list(onsets = list(from = character(), to = character())), rates),
    "non-empty character value"
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
    "ambiguous"
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
