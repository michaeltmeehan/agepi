test_that("rate_from_epiparameter converts mean delay to a Markov rate", {
  skip_if_not_installed("epiparameter")

  incubation <- epiparameter::epiparameter_db(
    disease = "COVID-19",
    epi_name = "incubation period",
    single_epiparameter = TRUE
  )

  expect_equal(rate_from_epiparameter(incubation), 1 / mean(incubation))
})

test_that("rate_from_epiparameter rejects non-epiparameter inputs", {
  skip_if_not_installed("epiparameter")

  expect_error(
    rate_from_epiparameter(5),
    "x must be an <epiparameter> object"
  )
})

test_that("rate_from_epiparameter rejects missing or invalid means", {
  skip_if_not_installed("epiparameter")

  fake_epiparameter <- structure(list(), class = c("fake_epiparameter", "epiparameter"))

  mean.fake_epiparameter <- function(x, ...) NA_real_
  expect_error(
    rate_from_epiparameter(fake_epiparameter),
    "mean delay must be a finite numeric scalar"
  )

  mean.fake_epiparameter <- function(x, ...) Inf
  expect_error(
    rate_from_epiparameter(fake_epiparameter),
    "mean delay must be a finite numeric scalar"
  )

  mean.fake_epiparameter <- function(x, ...) 0
  expect_error(
    rate_from_epiparameter(fake_epiparameter),
    "mean delay must be positive"
  )

  mean.fake_epiparameter <- function(x, ...) -1
  expect_error(
    rate_from_epiparameter(fake_epiparameter),
    "mean delay must be positive"
  )
})
