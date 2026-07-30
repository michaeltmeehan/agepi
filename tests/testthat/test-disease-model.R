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

test_that("SEIRModel constructs a valid inspectable model", {
  model <- SEIRModel(sigma = 0.3, gamma = 0.2)

  expect_s3_class(model, "DiseaseModel")
  expect_identical(model$model_type, "SEIR")
  expect_identical(model$compartments, c("S", "E", "I", "R"))
  expect_identical(model$transitions$from, c("S", "E", "I"))
  expect_identical(model$transitions$to, c("E", "I", "R"))
  expect_identical(model$sigma, 0.3)
  expect_identical(model$gamma, 0.2)
  expect_silent(validate_disease_model(model))
})

test_that("SEIRModel rejects invalid parameters", {
  expect_error(SEIRModel(gamma = 0.2), "sigma is required")
  expect_error(SEIRModel(sigma = 0.3), "gamma is required")
  expect_error(SEIRModel("0.3", 0.2), "sigma must be a finite numeric scalar")
  expect_error(SEIRModel(c(0.1, 0.2), 0.2), "sigma must be a finite numeric scalar")
  expect_error(SEIRModel(NA_real_, 0.2), "sigma must be a finite numeric scalar")
  expect_error(SEIRModel(Inf, 0.2), "sigma must be a finite numeric scalar")
  expect_error(SEIRModel(-0.1, 0.2), "sigma cannot be negative")
  expect_error(SEIRModel(0.3, "0.2"), "gamma must be a finite numeric scalar")
  expect_error(SEIRModel(0.3, NA_real_), "gamma must be a finite numeric scalar")
  expect_error(SEIRModel(0.3, -0.1), "gamma cannot be negative")
})

test_that("validate_disease_model validates stored beta values", {
  model <- SIRModel(gamma = 0.2, beta = 0.4)
  expect_silent(validate_disease_model(model))

  zero_beta <- model
  zero_beta$beta <- 0
  expect_silent(validate_disease_model(zero_beta))

  null_beta <- model
  null_beta$beta <- NULL
  expect_silent(validate_disease_model(null_beta))

  bad_betas <- list(-1, NA_real_, NaN, Inf, c(0.1, 0.2), "0.1")
  for (bad_beta in bad_betas) {
    malformed <- model
    malformed$beta <- bad_beta
    if (identical(bad_beta, -1)) {
      expect_error(validate_disease_model(malformed), "beta cannot be negative")
    } else {
      expect_error(validate_disease_model(malformed), "beta must be a finite numeric scalar")
    }
  }
})

test_that("validate_disease_model rejects malformed model objects", {
  expect_error(validate_disease_model("SIR"), "model must be a list")

  missing_gamma <- SIRModel(0.2)
  missing_gamma$gamma <- NULL
  expect_error(validate_disease_model(missing_gamma), "missing required field")

  bad_type <- SIRModel(0.2)
  bad_type$model_type <- "SIRS"
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

test_that("validate_disease_model rejects malformed SEIR model objects", {
  missing_sigma <- SEIRModel(0.3, 0.2)
  missing_sigma$sigma <- NULL
  expect_error(validate_disease_model(missing_sigma), "missing required field")

  bad_compartments <- SEIRModel(0.3, 0.2)
  bad_compartments$compartments <- c("S", "I", "R")
  expect_error(validate_disease_model(bad_compartments), "SEIR model compartments must be S, E, I, R")

  bad_transitions <- SEIRModel(0.3, 0.2)
  bad_transitions$transitions$to[1] <- "I"
  expect_error(validate_disease_model(bad_transitions), "SEIR model transitions")

  bad_sigma <- SEIRModel(0.3, 0.2)
  bad_sigma$sigma <- NA_real_
  expect_error(validate_disease_model(bad_sigma), "sigma must be a finite numeric scalar")
})
