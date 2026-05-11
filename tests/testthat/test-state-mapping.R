test_that("long-format state converts to a compartment-major vector", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")
  state_long <- data.frame(
    age_group = c("0-4", "0-4", "0-4", "5-9", "5-9", "5-9"),
    compartment = c("S", "I", "R", "S", "I", "R"),
    value = c(100, 2, 0, 120, 3, 1)
  )

  state_vector <- state_long_to_vector(state_long, ages, compartments)

  expect_equal(state_vector, c(100, 120, 2, 3, 0, 1), ignore_attr = TRUE)
  expect_identical(names(state_vector), c("S_0-4", "S_5-9", "I_0-4", "I_5-9", "R_0-4", "R_5-9"))
})

test_that("state vectors convert back to long format in the same ordering", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")

  state_long <- state_vector_to_long(c(100, 120, 2, 3, 0, 1), ages, compartments)

  expect_identical(state_long$compartment, c("S", "S", "I", "I", "R", "R"))
  expect_identical(state_long$age_group, c("0-4", "5-9", "0-4", "5-9", "0-4", "5-9"))
  expect_equal(state_long$value, c(100, 120, 2, 3, 0, 1))
})

test_that("round trip preserves totals by age and compartment", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")
  state_long <- data.frame(
    age_group = c("5-9", "0-4", "5-9", "0-4", "5-9", "0-4"),
    compartment = c("R", "S", "I", "R", "S", "I"),
    value = c(1, 100, 3, 0, 120, 2)
  )

  round_trip <- state_vector_to_long(
    state_long_to_vector(state_long, ages, compartments),
    ages,
    compartments
  )

  expect_equal(
    tapply(round_trip$value, round_trip$age_group, sum),
    tapply(state_long$value, state_long$age_group, sum)
  )
  expect_equal(
    tapply(round_trip$value, round_trip$compartment, sum),
    tapply(state_long$value, state_long$compartment, sum)
  )
})

test_that("round trip works with more than two compartments and age groups", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-14"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, 14)
  )
  compartments <- c("S", "E", "I", "R")
  state_long <- expand.grid(
    age_group = c("10-14", "0-4", "5-9"),
    compartment = c("I", "S", "R", "E"),
    stringsAsFactors = FALSE
  )
  state_long$value <- seq_len(nrow(state_long))

  state_vector <- state_long_to_vector(state_long, ages, compartments)
  round_trip <- state_vector_to_long(state_vector, ages, compartments)

  expect_identical(
    names(state_vector),
    c(
      "S_0-4", "S_5-9", "S_10-14",
      "E_0-4", "E_5-9", "E_10-14",
      "I_0-4", "I_5-9", "I_10-14",
      "R_0-4", "R_5-9", "R_10-14"
    )
  )
  expect_equal(sum(round_trip$value), sum(state_long$value))
  expect_identical(round_trip$compartment, rep(compartments, each = 3))
  expect_identical(round_trip$age_group, rep(ages$age_groups, times = 4))
})

test_that("state mapping requires exactly one row per compartment-age pair", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I")

  missing_row <- data.frame(
    age_group = c("0-4", "5-9", "0-4"),
    compartment = c("S", "S", "I"),
    value = c(1, 2, 3)
  )
  duplicate_row <- data.frame(
    age_group = c("0-4", "5-9", "0-4", "0-4"),
    compartment = c("S", "S", "I", "I"),
    value = c(1, 2, 3, 4)
  )

  expect_error(state_long_to_vector(missing_row, ages, compartments), "missing")
  expect_error(state_long_to_vector(duplicate_row, ages, compartments), "duplicate")
})

test_that("state mapping rejects non-numeric and missing state values", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I")

  non_numeric <- data.frame(
    age_group = c("0-4", "5-9", "0-4", "5-9"),
    compartment = c("S", "S", "I", "I"),
    value = c("1", "2", "3", "4")
  )
  missing_value <- data.frame(
    age_group = c("0-4", "5-9", "0-4", "5-9"),
    compartment = c("S", "S", "I", "I"),
    value = c(1, 2, NA, 4)
  )

  expect_error(state_long_to_vector(non_numeric, ages, compartments), "numeric")
  expect_error(state_long_to_vector(missing_value, ages, compartments), "missing")
})

test_that("state vectors must match the expected length", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    state_vector_to_long(c(1, 2, 3), ages, c("S", "I")),
    "length"
  )
})

test_that("state vector names are ignored when converting to long format", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  state_vector <- c(foo = 100, bar = 120, baz = 2, qux = 3)

  state_long <- state_vector_to_long(state_vector, ages, c("S", "I"))

  expect_identical(state_long$age_group, c("0-4", "5-9", "0-4", "5-9"))
  expect_identical(state_long$compartment, c("S", "S", "I", "I"))
  expect_equal(state_long$value, c(100, 120, 2, 3))
})
