summary_test_output <- function() {
  data.frame(
    time = rep(c(0, 1), each = 6),
    compartment = rep(c("M", "M", "N", "N", "X", "X"), times = 2),
    age_group = rep(c("young", "old", "young", "old", "young", "old"), times = 2),
    value = c(1, 2, 3, 4, 5, 6, 10, 20, 30, 40, 50, 60),
    stringsAsFactors = FALSE
  )
}

summary_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

summary_test_state <- function() {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(90, 180, 10, 20, 0, 0),
    stringsAsFactors = FALSE
  )
}

summary_test_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

test_that("compartment_totals sums across age groups", {
  totals <- compartment_totals(summary_test_output())

  expect_identical(names(totals), c("time", "compartment", "value"))
  expect_identical(totals$time, c(0, 0, 0, 1, 1, 1))
  expect_identical(totals$compartment, c("M", "N", "X", "M", "N", "X"))
  expect_equal(totals$value, c(3, 7, 11, 30, 70, 110))
})

test_that("age_group_totals sums across compartments", {
  totals <- age_group_totals(summary_test_output())

  expect_identical(names(totals), c("time", "age_group", "value"))
  expect_identical(totals$time, c(0, 0, 1, 1))
  expect_identical(totals$age_group, c("young", "old", "young", "old"))
  expect_equal(totals$value, c(9, 12, 90, 120))
})

test_that("total_population sums across compartments and age groups", {
  totals <- total_population(summary_test_output())

  expect_identical(names(totals), c("time", "value"))
  expect_identical(totals$time, c(0, 1))
  expect_equal(totals$value, c(21, 210))
})

test_that("summary helpers preserve first-seen time ordering", {
  output <- summary_test_output()
  output <- output[output$time == 1 | output$time == 0, ]
  output <- rbind(output[output$time == 1, ], output[output$time == 0, ])
  row.names(output) <- NULL

  expect_identical(total_population(output)$time, c(1, 0))
  expect_identical(age_group_totals(output)$time, c(1, 1, 0, 0))
})

test_that("summary helpers validate required columns", {
  output <- summary_test_output()

  expect_error(
    compartment_totals(output[, setdiff(names(output), "age_group")]),
    "missing required column"
  )
})

test_that("summary helpers validate value column", {
  output <- summary_test_output()

  non_numeric <- output
  non_numeric$value <- as.character(non_numeric$value)
  expect_error(age_group_totals(non_numeric), "value must be numeric")

  missing_value <- output
  missing_value$value[1] <- NA_real_
  expect_error(total_population(missing_value), "missing values")

  infinite_value <- output
  infinite_value$value[1] <- Inf
  expect_error(compartment_totals(infinite_value), "non-finite values")

  negative_value <- output
  negative_value$value[1] <- -1
  expect_error(age_group_totals(negative_value), "cannot be negative")
})

test_that("summary helpers work with simulate_deterministic output", {
  output <- simulate_deterministic(
    initial_state = summary_test_state(),
    times = seq(0, 0.3, by = 0.1),
    model = SIRModel(gamma = 0.2),
    age_structure = summary_test_ages(),
    contact_matrix = summary_test_contacts()
  )

  by_age <- age_group_totals(output)
  total <- total_population(output)

  expect_equal(by_age$value[by_age$age_group == "0-4"], rep(100, 4))
  expect_equal(by_age$value[by_age$age_group == "5-9"], rep(200, 4))
  expect_equal(total$value, rep(300, 4))
})
