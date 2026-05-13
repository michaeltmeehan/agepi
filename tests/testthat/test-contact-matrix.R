test_that("valid contact matrices pass validation", {
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  expect_silent(validate_contact_matrix(contact_matrix))
  expect_identical(validate_contact_matrix(contact_matrix), contact_matrix)
})

test_that("contact matrices validate against age structures", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_silent(validate_contact_matrix(diag(2), ages))
  expect_error(validate_contact_matrix(diag(3), ages), "dimensions")
})

test_that("contact matrices require numeric matrix input", {
  expect_error(
    validate_contact_matrix(data.frame(a = c(1, 0), b = c(0, 1))),
    "contact_matrix must be a numeric matrix"
  )
  expect_error(
    validate_contact_matrix(matrix(c("1", "0", "0", "1"), nrow = 2)),
    "contact_matrix must be a numeric matrix"
  )
})

test_that("contact matrices reject missing and non-finite values", {
  expect_error(
    validate_contact_matrix(matrix(c(1, NA, 0, 1), nrow = 2)),
    "missing or non-finite"
  )
  expect_error(
    validate_contact_matrix(matrix(c(1, Inf, 0, 1), nrow = 2)),
    "missing or non-finite"
  )
})

test_that("contact matrices reject negative entries", {
  contact_matrix <- matrix(c(
    1, -1,
    0, 1
  ), nrow = 2, byrow = TRUE)

  expect_error(
    validate_contact_matrix(contact_matrix),
    "contact_matrix cannot contain negative"
  )
})

test_that("contact matrices must be square", {
  expect_error(
    validate_contact_matrix(matrix(1, nrow = 2, ncol = 3)),
    "contact_matrix must be square"
  )
})

test_that("contact matrix validation validates supplied age structures", {
  ages <- list(
    age_groups = c("0-4", "5-9"),
    n_age_groups = 3,
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    validate_contact_matrix(diag(2), ages),
    "n_age_groups"
  )
})
