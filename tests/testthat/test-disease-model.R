test_that("SIRModel constructs a valid inspectable model", {
  model <- SIRModel(gamma = 0.2)

  expect_s3_class(model, "DiseaseModel")
  expect_identical(model$model_type, "SIR")
  expect_silent(validate_disease_model(model))
})

test_that("SIRModel stores S, I, R compartments", {
  model <- SIRModel(gamma = 0.2)

  expect_identical(model$compartments, c("S", "I", "R"))
})

test_that("SIRModel stores S to I and I to R transitions", {
  model <- SIRModel(gamma = 0.2)

  expect_identical(model$transitions$from, c("S", "I"))
  expect_identical(model$transitions$to, c("I", "R"))
})

test_that("SIRModel stores gamma", {
  model <- SIRModel(gamma = 0.2)

  expect_identical(model$gamma, 0.2)
})

test_that("SIRModel rejects invalid gamma", {
  expect_error(SIRModel(), "gamma is required")
  expect_error(SIRModel("0.2"), "gamma must be a finite numeric scalar")
  expect_error(SIRModel(c(0.1, 0.2)), "gamma must be a finite numeric scalar")
  expect_error(SIRModel(NA_real_), "gamma must be a finite numeric scalar")
  expect_error(SIRModel(Inf), "gamma must be a finite numeric scalar")
  expect_error(SIRModel(-0.1), "gamma cannot be negative")
})

test_that("validate_disease_model rejects malformed model objects", {
  expect_error(validate_disease_model("SIR"), "model must be a list")

  missing_gamma <- SIRModel(0.2)
  missing_gamma$gamma <- NULL
  expect_error(validate_disease_model(missing_gamma), "missing required field")

  bad_type <- SIRModel(0.2)
  bad_type$model_type <- "SEIR"
  expect_error(validate_disease_model(bad_type), "unsupported disease model type")

  bad_compartments <- SIRModel(0.2)
  bad_compartments$compartments <- c("S", "E", "I", "R")
  expect_error(validate_disease_model(bad_compartments), "compartments must be S, I, R")

  bad_transitions <- SIRModel(0.2)
  bad_transitions$transitions$to[1] <- "R"
  expect_error(validate_disease_model(bad_transitions), "transitions must be S -> I and I -> R")

  bad_gamma <- SIRModel(0.2)
  bad_gamma$gamma <- NA_real_
  expect_error(validate_disease_model(bad_gamma), "gamma must be a finite numeric scalar")
})
