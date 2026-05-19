derivative_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

derivative_test_rates <- function() {
  data.frame(
    from = c("S", "I", "S", "I"),
    to = c("I", "R", "I", "R"),
    age_group = c("0-4", "0-4", "5-9", "5-9"),
    rate = c(27, 2, 126, 4),
    stringsAsFactors = FALSE
  )
}

test_that("rates_to_derivative computes a simple two-age-group SIR derivative", {
  derivative <- rates_to_derivative(
    derivative_test_rates(),
    compartments = c("S", "I", "R"),
    age_structure = derivative_test_ages()
  )

  expected <- data.frame(
    compartment = c("S", "S", "I", "I", "R", "R"),
    age_group = c("0-4", "5-9", "0-4", "5-9", "0-4", "5-9"),
    derivative = c(-27, -126, 25, 122, 2, 4),
    stringsAsFactors = FALSE
  )

  expect_equal(derivative, expected)
})

test_that("transition_rates output can be passed directly to rates_to_derivative", {
  ages <- derivative_test_ages()
  state <- data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(90, 180, 10, 20, 0, 0),
    stringsAsFactors = FALSE
  )
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  rates <- transition_rates(
    state = state,
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = contact_matrix
  )
  derivative <- rates_to_derivative(rates, c("S", "I", "R"), ages)

  expect_equal(derivative$derivative, c(-27, -126, 25, 122, 2, 4))
})

test_that("rates_to_derivative returns exactly the expected output columns", {
  derivative <- rates_to_derivative(
    derivative_test_rates(),
    c("S", "I", "R"),
    derivative_test_ages()
  )

  expect_identical(names(derivative), c("compartment", "age_group", "derivative"))
})

test_that("rates_to_derivative orders compartments outermost and age groups innermost", {
  derivative <- rates_to_derivative(
    derivative_test_rates(),
    c("S", "I", "R"),
    derivative_test_ages()
  )

  expect_identical(derivative$compartment, c("S", "S", "I", "I", "R", "R"))
  expect_identical(derivative$age_group, c("0-4", "5-9", "0-4", "5-9", "0-4", "5-9"))
})

test_that("rates_to_derivative conserves mass within each age group for closed SIR rates", {
  derivative <- rates_to_derivative(
    derivative_test_rates(),
    c("S", "I", "R"),
    derivative_test_ages()
  )

  age_sums <- tapply(derivative$derivative, derivative$age_group, sum)
  expect_equal(as.numeric(age_sums), c(0, 0))
  expect_identical(names(age_sums), c("0-4", "5-9"))
})

test_that("rates_to_derivative conserves total population for infection-only SEIR rates", {
  ages <- derivative_test_ages()
  state <- data.frame(
    compartment = rep(c("S", "E", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 4),
    value = c(90, 180, 5, 15, 10, 20, 0, 0),
    stringsAsFactors = FALSE
  )
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  rates <- transition_rates(
    state = state,
    model = SEIRModel(sigma = 0.3, gamma = 0.2),
    age_structure = ages,
    contact_matrix = contact_matrix
  )
  derivative <- rates_to_derivative(rates, c("S", "E", "I", "R"), ages)

  age_sums <- tapply(derivative$derivative, derivative$age_group, sum)
  expect_equal(as.numeric(age_sums), c(0, 0))
  expect_equal(sum(derivative$derivative), 0)
})


test_that("rates_to_derivative handles zero rates", {
  rates <- data.frame(
    from = c("S", "I"),
    to = c("I", "R"),
    age_group = c("0-4", "5-9"),
    rate = c(0, 0),
    stringsAsFactors = FALSE
  )

  derivative <- rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages())

  expect_equal(derivative$derivative, rep(0, 6))
  expect_type(derivative$derivative, "double")
})

test_that("rates_to_derivative rejects missing required columns", {
  rates <- derivative_test_rates()
  rates$to <- NULL

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "missing required column"
  )
})

test_that("rates_to_derivative rejects non-data-frame transition rate input", {
  expect_error(
    rates_to_derivative(list(from = "S"), c("S", "I", "R"), derivative_test_ages()),
    "must be a data frame"
  )
})

test_that("rates_to_derivative rejects unknown source compartments", {
  rates <- derivative_test_rates()
  rates$from[1] <- "E"

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "unknown source compartment"
  )
})

test_that("rates_to_derivative rejects unknown destination compartments", {
  rates <- derivative_test_rates()
  rates$to[1] <- "E"

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "unknown destination compartment"
  )
})

test_that("rates_to_derivative rejects unknown age groups", {
  rates <- derivative_test_rates()
  rates$age_group[1] <- "10-14"

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "unknown age_group"
  )
})

test_that("rates_to_derivative rejects non-numeric rates", {
  rates <- derivative_test_rates()
  rates$rate <- as.character(rates$rate)

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "rate must be numeric"
  )
})

test_that("rates_to_derivative rejects missing rates", {
  rates <- derivative_test_rates()
  rates$rate[1] <- NA_real_

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "rate cannot contain missing"
  )
})

test_that("rates_to_derivative rejects negative rates", {
  rates <- derivative_test_rates()
  rates$rate[1] <- -1

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "rate cannot contain negative"
  )
})

test_that("rates_to_derivative rejects invalid compartments vectors", {
  expect_error(
    rates_to_derivative(derivative_test_rates(), c(1, 2, 3), derivative_test_ages()),
    "compartments must be a character vector"
  )
  expect_error(
    rates_to_derivative(derivative_test_rates(), character(), derivative_test_ages()),
    "compartments must contain at least one"
  )
})

test_that("rates_to_derivative rejects duplicate compartment names", {
  expect_error(
    rates_to_derivative(derivative_test_rates(), c("S", "I", "I"), derivative_test_ages()),
    "compartments must be unique"
  )
})

test_that("rates_to_derivative rejects invalid age structures", {
  ages <- derivative_test_ages()
  ages$n_age_groups <- 3

  expect_error(
    rates_to_derivative(derivative_test_rates(), c("S", "I", "R"), ages),
    "n_age_groups"
  )
})

test_that("rates_to_derivative rejects exact duplicate transition rows", {
  rates <- rbind(derivative_test_rates(), derivative_test_rates()[1, ])

  expect_error(
    rates_to_derivative(rates, c("S", "I", "R"), derivative_test_ages()),
    "duplicate transition row"
  )
})
