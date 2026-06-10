calibration_test_trajectory <- function() {
  data.frame(
    time = rep(c(0, 1), each = 4),
    compartment = rep(c("S", "S", "I", "I"), times = 2),
    age_group = rep(c("0-4", "5+"), times = 4),
    value = c(90, 180, 10, 20, 85, 170, 15, 30),
    stringsAsFactors = FALSE
  )
}

calibration_test_cumulative <- function() {
  data.frame(
    time = rep(c(0, 1), each = 2),
    cumulative_name = "infections",
    transition_id = "infection:S->I",
    from = "S",
    to = "I",
    age_group = rep(c("0-4", "5+"), times = 2),
    value = c(0, 0, 5, 10),
    stringsAsFactors = FALSE
  )
}

test_that("trajectory calibration target evaluates Gaussian likelihood", {
  target <- CalibrationTarget(
    name = "infectious_by_age",
    data = data.frame(
      time = c(1, 1),
      age_group = c("0-4", "5+"),
      value = c(14, 31),
      sd = c(2, 2),
      stringsAsFactors = FALSE
    ),
    type = "trajectory",
    observation_model = "gaussian",
    compartment = "I"
  )

  result <- evaluate_calibration_target(target, calibration_test_trajectory())

  expected <- sum(stats::dnorm(c(14, 31), mean = c(15, 30), sd = c(2, 2), log = TRUE))
  expect_equal(result$target, "infectious_by_age")
  expect_equal(result$log_likelihood, expected)
  expect_equal(result$objective, -expected)
  expect_equal(result$n_used, 2)
  expect_equal(result$n_missing, 0)
  expect_equal(result$details$predicted, c(15, 30))
})

test_that("trajectory calibration target evaluates lognormal likelihood", {
  target <- CalibrationTarget(
    name = "infectious_total",
    data = data.frame(time = 1, value = 44, sdlog = 0.1),
    type = "trajectory",
    observation_model = "lognormal",
    compartment = "I"
  )

  result <- evaluate_calibration_target(target, calibration_test_trajectory())

  expected <- stats::dlnorm(44, meanlog = log(45), sdlog = 0.1, log = TRUE)
  expect_equal(result$log_likelihood, expected)
  expect_equal(result$details$predicted, 45)
})

test_that("cumulative-flow calibration target evaluates Poisson likelihood", {
  target <- CalibrationTarget(
    name = "infections",
    data = data.frame(
      time = c(1, 1),
      age_group = c("0-4", "5+"),
      value = c(4, 12),
      stringsAsFactors = FALSE
    ),
    type = "cumulative",
    observation_model = "poisson",
    cumulative_name = "infections"
  )

  result <- evaluate_calibration_target(target, calibration_test_cumulative())

  expected <- sum(stats::dpois(c(4, 12), lambda = c(5, 10), log = TRUE))
  expect_equal(result$log_likelihood, expected)
  expect_equal(result$details$predicted, c(5, 10))
})

test_that("multiple calibration targets sum into a total objective", {
  trajectory_target <- CalibrationTarget(
    name = "infectious",
    data = data.frame(time = 1, value = 44, sd = 5),
    type = "trajectory",
    observation_model = "gaussian",
    compartment = "I"
  )
  cumulative_target <- CalibrationTarget(
    name = "infections",
    data = data.frame(time = 1, value = 16),
    type = "cumulative",
    observation_model = "poisson",
    cumulative_name = "infections"
  )
  output <- list(
    trajectory = calibration_test_trajectory(),
    cumulative = calibration_test_cumulative()
  )

  result <- evaluate_calibration_objective(list(trajectory_target, cumulative_target), output)

  expected_log_likelihood <- stats::dnorm(44, mean = 45, sd = 5, log = TRUE) +
    stats::dpois(16, lambda = 15, log = TRUE)
  expect_equal(result$log_likelihood, expected_log_likelihood)
  expect_equal(result$objective, -expected_log_likelihood)
  expect_equal(result$n_used, 2)
  expect_length(result$targets, 2)
})

test_that("calibration target validation catches missing required columns", {
  expect_error(
    CalibrationTarget(
      name = "bad",
      data = data.frame(value = 1),
      type = "trajectory",
      observation_model = "gaussian",
      compartment = "I"
    ),
    "missing required column"
  )
})

test_that("calibration target validation catches unknown observation models", {
  expect_error(
    CalibrationTarget(
      name = "bad",
      data = data.frame(time = 1, value = 1),
      type = "trajectory",
      observation_model = "mystery",
      compartment = "I"
    ),
    "unsupported observation model"
  )
})

test_that("Poisson calibration rejects negative observed counts", {
  target <- CalibrationTarget(
    name = "counts",
    data = data.frame(time = 1, value = -1),
    type = "cumulative",
    observation_model = "poisson",
    cumulative_name = "infections"
  )

  expect_error(
    evaluate_calibration_target(target, calibration_test_cumulative()),
    "negative observed counts"
  )
})

test_that("calibration evaluation reports missing model outputs for requested rows", {
  target <- CalibrationTarget(
    name = "missing_time",
    data = data.frame(time = c(1, 2), value = c(44, 50), sd = c(5, 5)),
    type = "trajectory",
    observation_model = "gaussian",
    compartment = "I"
  )

  result <- evaluate_calibration_target(target, calibration_test_trajectory())

  expect_equal(result$n_observations, 2)
  expect_equal(result$n_used, 1)
  expect_equal(result$n_missing, 1)
  expect_equal(result$missing$time, 2)
})

test_that("calibration evaluation errors when selected model outputs are absent", {
  target <- CalibrationTarget(
    name = "absent_compartment",
    data = data.frame(time = 1, value = 1, sd = 1),
    type = "trajectory",
    observation_model = "gaussian",
    compartment = "R"
  )

  expect_error(
    evaluate_calibration_target(target, calibration_test_trajectory()),
    "no model outputs matched"
  )
})

test_that("calibration evaluation supports deterministic output lists with trajectory", {
  target <- CalibrationTarget(
    name = "infectious",
    data = data.frame(time = 1, value = 45, sd = 1),
    type = "trajectory",
    observation_model = "gaussian",
    compartment = "I"
  )

  result <- evaluate_calibration_target(
    target,
    list(trajectory = calibration_test_trajectory())
  )

  expect_equal(result$n_used, 1)
  expect_equal(result$details$predicted, 45)
})

test_that("calibration evaluation supports deterministic output lists with cumulative", {
  target <- CalibrationTarget(
    name = "infections",
    data = data.frame(time = 1, value = 15),
    type = "cumulative",
    observation_model = "poisson",
    cumulative_name = "infections"
  )

  result <- evaluate_calibration_target(
    target,
    list(cumulative = calibration_test_cumulative())
  )

  expect_equal(result$n_used, 1)
  expect_equal(result$details$predicted, 15)
})

test_that("calibration evaluation errors clearly for malformed outputs", {
  target <- CalibrationTarget(
    name = "infectious",
    data = data.frame(time = 1, value = 45, sd = 1),
    type = "trajectory",
    observation_model = "gaussian",
    compartment = "I"
  )

  expect_error(
    evaluate_calibration_target(target, list(events = data.frame())),
    "list with \\$trajectory"
  )
  expect_error(
    evaluate_calibration_target(target, calibration_test_trajectory()[, c("time", "compartment", "value")]),
    "missing required column"
  )

  bad_output <- calibration_test_trajectory()
  bad_output$value[1] <- Inf
  expect_error(
    evaluate_calibration_target(target, bad_output),
    "finite and non-missing"
  )
})

test_that("calibration evaluation errors clearly for unsupported target output types", {
  target <- CalibrationTarget(
    name = "infectious",
    data = data.frame(time = 1, value = 45, sd = 1),
    type = "trajectory",
    observation_model = "gaussian",
    compartment = "I"
  )
  target$type <- "events"

  expect_error(
    evaluate_calibration_target(target, calibration_test_trajectory()),
    "unsupported calibration output type"
  )
})
