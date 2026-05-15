demography_validation_age_structure <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
}

demography_validation_observed <- function(population = c(100, 200, 300, 110, 210, 330)) {
  Demography(
    data.frame(
      time = rep(c(0, 1), each = 3),
      age_group = rep(c("0-4", "5-9", "10+"), times = 2),
      population = population,
      stringsAsFactors = FALSE
    ),
    demography_validation_age_structure()
  )
}

demography_validation_simulated <- function(population = c(100, 200, 300, 110, 210, 330)) {
  data.frame(
    time = rep(c(0, 1), each = 3),
    age_group = rep(c("0-4", "5-9", "10+"), times = 2),
    population = population,
    stringsAsFactors = FALSE
  )
}

test_that("compare_demography_to_observed returns zero errors for a perfect match", {
  observed <- demography_validation_observed()
  simulated <- demography_validation_simulated()

  comparison <- compare_demography_to_observed(simulated, observed)

  expect_equal(comparison$absolute_error, rep(0, 6))
  expect_equal(comparison$relative_error, rep(0, 6))
  expect_equal(comparison$total_absolute_error, rep(0, 6))
  expect_equal(comparison$total_relative_error, rep(0, 6))
  expect_equal(comparison$age_share_error, rep(0, 6))
})

test_that("compare_demography_to_observed computes hand-checkable mismatches", {
  observed <- Demography(
    data.frame(
      time = c(0, 0),
      age_group = c("0-4", "5-9"),
      population = c(100, 300),
      stringsAsFactors = FALSE
    ),
    AgeStructure(c("0-4", "5-9"), c(0, 5), c(4, 9))
  )
  simulated <- data.frame(
    time = c(0, 0),
    age_group = c("0-4", "5-9"),
    population = c(110, 270),
    stringsAsFactors = FALSE
  )

  comparison <- compare_demography_to_observed(simulated, observed)

  expect_equal(comparison$absolute_error, c(10, -30))
  expect_equal(comparison$relative_error, c(0.1, -0.1))
  expect_equal(comparison$simulated_total_population, c(380, 380))
  expect_equal(comparison$observed_total_population, c(400, 400))
  expect_equal(comparison$total_absolute_error, c(-20, -20))
  expect_equal(comparison$total_relative_error, c(-0.05, -0.05))
  expect_equal(comparison$simulated_age_share, c(110 / 380, 270 / 380))
  expect_equal(comparison$observed_age_share, c(100 / 400, 300 / 400))
  expect_equal(
    comparison$age_share_error,
    c(110 / 380 - 100 / 400, 270 / 380 - 300 / 400)
  )
})

test_that("compare_demography_to_observed compares only shared exact times", {
  observed <- demography_validation_observed()
  simulated <- rbind(
    demography_validation_simulated(),
    data.frame(
      time = 2,
      age_group = c("0-4", "5-9", "10+"),
      population = c(1, 2, 3),
      stringsAsFactors = FALSE
    )
  )

  comparison <- compare_demography_to_observed(simulated, observed)

  expect_identical(unique(comparison$time), c(0, 1))
  expect_equal(nrow(comparison), 6)

  simulated_no_shared_time <- transform(simulated, time = time + 10)
  expect_error(
    compare_demography_to_observed(simulated_no_shared_time, observed),
    "no shared exact time"
  )
})

test_that("compare_demography_to_observed errors clearly for age-group mismatches", {
  observed <- demography_validation_observed()
  simulated <- demography_validation_simulated()
  simulated$age_group[simulated$age_group == "10+"] <- "10-14"

  expect_error(
    compare_demography_to_observed(simulated, observed),
    "age groups must match exactly"
  )

  simulated <- demography_validation_simulated()
  simulated <- simulated[-3, ]
  expect_error(
    compare_demography_to_observed(simulated, observed),
    "age groups must match exactly"
  )
})

test_that("compare_demography_to_observed uses NA relative error for zero observed population", {
  observed <- demography_validation_observed(
    population = c(0, 200, 300, 110, 210, 330)
  )
  simulated <- demography_validation_simulated(
    population = c(10, 200, 300, 110, 210, 330)
  )

  comparison <- compare_demography_to_observed(simulated, observed)

  expect_equal(comparison$absolute_error[1], 10)
  expect_true(is.na(comparison$relative_error[1]))
})

test_that("summarise_demography_comparison returns one row per time with summaries", {
  observed <- demography_validation_observed()
  simulated <- demography_validation_simulated(
    population = c(110, 190, 300, 100, 220, 360)
  )

  comparison <- compare_demography_to_observed(simulated, observed)
  summary <- summarise_demography_comparison(comparison)

  expect_identical(summary$time, c(0, 1))
  expect_equal(summary$simulated_total_population, c(600, 680))
  expect_equal(summary$observed_total_population, c(600, 650))
  expect_equal(summary$total_absolute_error, c(0, 30))
  expect_equal(summary$total_relative_error, c(0, 30 / 650))
  expect_equal(summary$mean_absolute_age_error, c(mean(c(10, 10, 0)), mean(c(10, 10, 30))))
  expect_equal(summary$max_absolute_age_error, c(10, 30))
  expect_equal(
    summary$mean_absolute_age_share_error,
    as.numeric(tapply(abs(comparison$age_share_error), comparison$time, mean))
  )
  expect_equal(
    summary$max_absolute_age_share_error,
    as.numeric(tapply(abs(comparison$age_share_error), comparison$time, max))
  )
})

test_that("compare_demography_to_observed output is ordered by time and observed age order", {
  observed <- demography_validation_observed()
  simulated <- demography_validation_simulated()
  simulated <- simulated[c(6, 4, 5, 3, 1, 2), ]

  comparison <- compare_demography_to_observed(simulated, observed)

  expect_identical(comparison$time, rep(c(0, 1), each = 3))
  expect_identical(comparison$age_group, rep(c("0-4", "5-9", "10+"), times = 2))
})
