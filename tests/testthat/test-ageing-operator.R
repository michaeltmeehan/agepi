test_that("AgeingOperator works on default one-year WPP age structure", {
  ages <- wpp_age_structure_1year()
  ageing <- AgeingOperator(ages)

  expect_s3_class(ageing, "agepi_ageing_operator")
  expect_identical(ageing$age_structure, ages)
  expect_identical(ageing$age_groups, ages$age_groups)
  expect_identical(ageing$n_age_groups, ages$n_age_groups)
  expect_identical(ageing$lower_bounds, ages$lower_bounds)
  expect_identical(ageing$upper_bounds, ages$upper_bounds)
  expect_equal(ageing$widths[1:100], rep(1, 100))
  expect_true(is.infinite(ageing$widths[101]))
  expect_equal(ageing$departure_rate[1:100], rep(1, 100))
  expect_identical(ageing$departure_rate[101], 0)
  expect_identical(ageing$destination_index[1:100], 2:101)
  expect_true(is.na(ageing$destination_index[101]))
  expect_silent(validate_ageing_operator(ageing))
})

test_that("AgeingOperator uses inverse widths on five-year age structure", {
  ageing <- AgeingOperator(wpp_age_structure_5year())

  expect_equal(ageing$widths[1:20], rep(5, 20))
  expect_true(is.infinite(ageing$widths[21]))
  expect_equal(ageing$departure_rate[1:20], rep(1 / 5, 20))
  expect_identical(ageing$departure_rate[21], 0)
  expect_identical(ageing$destination_index[1:20], 2:21)
  expect_true(is.na(ageing$destination_index[21]))
})

test_that("AgeingOperator rejects non-contiguous age structures", {
  ages <- AgeStructure(
    age_groups = c("0-4", "10-14", "15+"),
    lower_bounds = c(0, 10, 15),
    upper_bounds = c(4, 14, Inf)
  )

  expect_error(AgeingOperator(ages), "contiguous")
})

test_that("AgeingOperator requires an open-ended final bin", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(AgeingOperator(ages), "open-ended final")
})

test_that("validate_ageing_operator rejects corrupted departure rates", {
  ageing <- AgeingOperator(wpp_age_structure_5year())

  ageing$departure_rate[1] <- -1
  expect_error(validate_ageing_operator(ageing), "negative")

  ageing <- AgeingOperator(wpp_age_structure_5year())
  ageing$departure_rate[1] <- Inf
  expect_error(validate_ageing_operator(ageing), "finite numeric")

  ageing <- AgeingOperator(wpp_age_structure_5year())
  ageing$departure_rate[21] <- 1
  expect_error(validate_ageing_operator(ageing), "final bin")
})

test_that("validate_ageing_operator rejects corrupted destination indices", {
  ageing <- AgeingOperator(wpp_age_structure_5year())

  ageing$destination_index[1] <- 3L
  expect_error(validate_ageing_operator(ageing), "destination_index")

  ageing <- AgeingOperator(wpp_age_structure_5year())
  ageing$destination_index[21] <- 21L
  expect_error(validate_ageing_operator(ageing), "destination_index")
})

test_that("validate_ageing_operator rejects inconsistent structure", {
  ageing <- AgeingOperator(wpp_age_structure_5year())

  ageing$age_groups <- ageing$age_groups[-1]
  expect_error(validate_ageing_operator(ageing), "length")

  ageing <- AgeingOperator(wpp_age_structure_5year())
  class(ageing) <- "list"
  expect_error(validate_ageing_operator(ageing), "agepi_ageing_operator")
})
