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
    stringsAsFactors = FALSE
  )

  expect_equal(rates, expected)
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

  expect_identical(names(rates), c("from", "to", "age_group", "rate"))
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
  model$model_type <- "SEIR"

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
