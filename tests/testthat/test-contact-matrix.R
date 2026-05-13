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

test_that("plain numeric matrix is returned unchanged under recipient_source orientation", {
  contact_matrix <- matrix(c(
    1, 2,
    3, 4
  ), nrow = 2, byrow = TRUE)

  expect_equal(as_agepi_contact_matrix(contact_matrix), contact_matrix)
})

test_that("source_recipient orientation transposes the matrix", {
  contact_matrix <- matrix(c(
    1, 2,
    3, 4
  ), nrow = 2, byrow = TRUE)

  expect_equal(
    as_agepi_contact_matrix(contact_matrix, orientation = "source_recipient"),
    t(contact_matrix)
  )
})

test_that("transpose is applied after orientation handling", {
  contact_matrix <- matrix(c(
    1, 2,
    3, 4
  ), nrow = 2, byrow = TRUE)

  expect_equal(
    as_agepi_contact_matrix(
      contact_matrix,
      orientation = "source_recipient",
      transpose = TRUE
    ),
    contact_matrix
  )
})

test_that("numeric data frame matrix input works", {
  contact_data <- data.frame(
    young = c(1, 3),
    adult = c(2, 4)
  )

  expect_equal(
    as_agepi_contact_matrix(contact_data),
    matrix(c(1, 3, 2, 4), nrow = 2)
  )
})

test_that("socialmixr-like list with matrix works", {
  contact_matrix <- matrix(c(
    1, 2,
    3, 4
  ), nrow = 2, byrow = TRUE)
  socialmixr_like <- list(matrix = contact_matrix, participants = data.frame())

  expect_equal(as_agepi_contact_matrix(socialmixr_like), contact_matrix)
})

test_that("contact_matrix_from_socialmixr extracts a plain matrix", {
  contact_matrix <- matrix(c(
    1, 2,
    3, 4
  ), nrow = 2, byrow = TRUE)

  expect_equal(contact_matrix_from_socialmixr(contact_matrix), contact_matrix)
})

test_that("contact_matrix_from_socialmixr extracts list matrix element", {
  contact_matrix <- matrix(c(
    1, 2,
    3, 4
  ), nrow = 2, byrow = TRUE)
  socialmixr_like <- structure(
    list(matrix = contact_matrix, participants = data.frame()),
    class = "contact_matrix"
  )

  expect_equal(contact_matrix_from_socialmixr(socialmixr_like), contact_matrix)
})

test_that("contact_matrix_from_socialmixr handles explicit orientation", {
  contact_matrix <- matrix(c(
    1, 2,
    3, 4
  ), nrow = 2, byrow = TRUE)
  socialmixr_like <- list(matrix = contact_matrix)

  expect_equal(
    contact_matrix_from_socialmixr(socialmixr_like, orientation = "source_recipient"),
    t(contact_matrix)
  )
})

test_that("contact_matrix_from_socialmixr applies and checks age_structure labels", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  contact_matrix <- contact_matrix_from_socialmixr(list(matrix = diag(2)), age_structure = ages)

  expect_identical(rownames(contact_matrix), ages$age_groups)
  expect_identical(colnames(contact_matrix), ages$age_groups)

  labelled_matrix <- diag(2)
  dimnames(labelled_matrix) <- list(c("5-9", "0-4"), c("5-9", "0-4"))

  expect_error(
    contact_matrix_from_socialmixr(list(matrix = labelled_matrix), age_structure = ages),
    "rownames"
  )
})

test_that("contact_matrix_from_socialmixr rejects invalid matrices", {
  expect_error(
    contact_matrix_from_socialmixr(list(matrix = matrix(1, nrow = 2, ncol = 3))),
    "square"
  )

  expect_error(
    contact_matrix_from_socialmixr(list(matrix = matrix(c(1, Inf, 0, 1), nrow = 2))),
    "non-finite"
  )

  expect_error(
    contact_matrix_from_socialmixr(list(matrix = matrix(c(1, NA, 0, 1), nrow = 2))),
    "missing"
  )

  expect_error(
    contact_matrix_from_socialmixr(list(participants = data.frame())),
    "matrix element"
  )
})

test_that("contact_matrix_from_socialmixr example skips cleanly without socialmixr", {
  testthat::skip_if_not_installed("socialmixr")

  expect_true(requireNamespace("socialmixr", quietly = TRUE))
})

test_that("contact_matrix_from_socialmixr accepts a small socialmixr result", {
  testthat::skip_if_not_installed("socialmixr")

  socialmixr_datasets <- utils::data(package = "socialmixr")$results[, "Item"]
  testthat::skip_if_not("polymod" %in% socialmixr_datasets)

  data_environment <- new.env(parent = emptyenv())
  utils::data("polymod", package = "socialmixr", envir = data_environment)

  socialmixr_matrix <- suppressWarnings(
    socialmixr::contact_matrix(
      get("polymod", envir = data_environment),
      countries = "United Kingdom",
      age_limits = c(0, 5, 10),
      symmetric = FALSE
    )
  )

  contact_matrix <- contact_matrix_from_socialmixr(socialmixr_matrix)

  expect_type(contact_matrix, "double")
  expect_equal(dim(contact_matrix), c(3L, 3L))
  expect_false(anyNA(contact_matrix))
})

test_that("conmat-style long table is converted correctly", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9"),
    contacts = c(1, 2, 3, 4)
  )

  expected <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("0-4", "5-9"), c("0-4", "5-9"))
  )

  expect_equal(as_agepi_contact_matrix(conmat_long), expected)
})

test_that("contact_matrix_from_conmat converts conmat-style long table", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9"),
    contacts = c(1, 2, 3, 4)
  )

  expected <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("0-4", "5-9"), c("0-4", "5-9"))
  )

  expect_equal(contact_matrix_from_conmat(conmat_long), expected)
})

test_that("contact_matrix_from_conmat delegates orientation handling", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9"),
    contacts = c(1, 2, 3, 4)
  )

  expected_recipient_source <- contact_matrix_from_conmat(conmat_long)

  expect_equal(
    contact_matrix_from_conmat(conmat_long, orientation = "source_recipient"),
    t(expected_recipient_source)
  )
})

test_that("contact_matrix_from_conmat rejects missing required long-table columns", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9"),
    age_group_to = c("0-4", "0-4")
  )

  expect_error(
    contact_matrix_from_conmat(conmat_long),
    "age_group_from, age_group_to, contacts"
  )
})

test_that("contact_matrix_from_conmat rejects duplicate age-pair records", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "0-4"),
    age_group_to = c("5-9", "5-9"),
    contacts = c(1, 2)
  )

  expect_error(
    contact_matrix_from_conmat(conmat_long),
    "duplicate"
  )
})

test_that("contact_matrix_from_conmat checks age labels against age_structure", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  conmat_long <- data.frame(
    age_group_from = c("0-4", "10-14", "0-4", "10-14"),
    age_group_to = c("0-4", "0-4", "10-14", "10-14"),
    contacts = c(1, 2, 3, 4)
  )

  expect_error(
    contact_matrix_from_conmat(conmat_long, age_structure = ages),
    "age groups must match"
  )
})

test_that("contact_matrix_from_conmat rejects invalid contacts values", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9"),
    contacts = c(1, 2, 3, 4)
  )

  nonnumeric_contacts <- conmat_long
  nonnumeric_contacts$contacts <- as.character(nonnumeric_contacts$contacts)
  expect_error(
    contact_matrix_from_conmat(nonnumeric_contacts),
    "contacts column must be numeric"
  )

  missing_contacts <- conmat_long
  missing_contacts$contacts[1] <- NA_real_
  expect_error(
    contact_matrix_from_conmat(missing_contacts),
    "missing"
  )

  infinite_contacts <- conmat_long
  infinite_contacts$contacts[1] <- Inf
  expect_error(
    contact_matrix_from_conmat(infinite_contacts),
    "non-finite"
  )

  negative_contacts <- conmat_long
  negative_contacts$contacts[1] <- -1
  expect_error(
    contact_matrix_from_conmat(negative_contacts),
    "negative"
  )
})

test_that("contact_matrix_from_conmat example skips cleanly without conmat", {
  testthat::skip_if_not_installed("conmat")

  expect_true(requireNamespace("conmat", quietly = TRUE))
})

test_that("contact_matrix_from_conmat accepts a small conmat-generated long table", {
  testthat::skip_if_not_installed("conmat")

  source_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("0-4", "5-9"), c("0-4", "5-9"))
  )

  conmat_long <- suppressWarnings(conmat::matrix_to_predictions(source_matrix))
  contact_matrix <- contact_matrix_from_conmat(conmat_long)

  expect_equal(contact_matrix, source_matrix)
})

test_that("conmat-style long table uses age_structure order", {
  ages <- AgeStructure(
    age_groups = c("5-9", "0-4"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  ages$age_groups <- c("5-9", "0-4")

  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9"),
    contacts = c(1, 2, 3, 4)
  )

  expected <- matrix(
    c(4, 3, 2, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("5-9", "0-4"), c("5-9", "0-4"))
  )

  expect_equal(as_agepi_contact_matrix(conmat_long, age_structure = ages), expected)
})

test_that("conmat-style duplicate pairs are rejected", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "0-4"),
    age_group_to = c("5-9", "5-9"),
    contacts = c(1, 2)
  )

  expect_error(
    as_agepi_contact_matrix(conmat_long),
    "duplicate"
  )
})

test_that("conmat-style missing matrix combinations are rejected", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9", "0-4"),
    age_group_to = c("0-4", "0-4", "5-9"),
    contacts = c(1, 2, 3)
  )

  expect_error(
    as_agepi_contact_matrix(conmat_long),
    "every age_group_from/age_group_to combination"
  )
})

test_that("nonnumeric input is rejected", {
  expect_error(
    as_agepi_contact_matrix(matrix(c("1", "0", "0", "1"), nrow = 2)),
    "numeric"
  )

  expect_error(
    as_agepi_contact_matrix(data.frame(a = c("1", "0"), b = c("0", "1"))),
    "numeric"
  )
})

test_that("missing, infinite, and negative entries are rejected", {
  expect_error(
    as_agepi_contact_matrix(matrix(c(1, NA, 0, 1), nrow = 2)),
    "missing"
  )

  expect_error(
    as_agepi_contact_matrix(matrix(c(1, Inf, 0, 1), nrow = 2)),
    "non-finite"
  )

  expect_error(
    as_agepi_contact_matrix(matrix(c(1, -1, 0, 1), nrow = 2)),
    "negative"
  )
})

test_that("nonsquare matrix is rejected", {
  expect_error(
    as_agepi_contact_matrix(matrix(1, nrow = 2, ncol = 3)),
    "square"
  )
})

test_that("age_structure dimension mismatch is rejected", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    as_agepi_contact_matrix(diag(3), age_structure = ages),
    "dimensions"
  )
})

test_that("age_structure names are applied when absent", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  contact_matrix <- as_agepi_contact_matrix(diag(2), age_structure = ages)

  expect_identical(rownames(contact_matrix), ages$age_groups)
  expect_identical(colnames(contact_matrix), ages$age_groups)
})

test_that("age_structure names are checked when present", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  contact_matrix <- diag(2)
  dimnames(contact_matrix) <- list(c("5-9", "0-4"), c("5-9", "0-4"))

  expect_error(
    as_agepi_contact_matrix(contact_matrix, age_structure = ages),
    "rownames"
  )
})

test_that("returned matrix works with force_of_infection", {
  conmat_long <- data.frame(
    age_group_from = c("0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9"),
    contacts = c(2, 1, 3, 4)
  )

  contact_matrix <- as_agepi_contact_matrix(conmat_long)
  lambda <- force_of_infection(
    infectious = c(10, 20),
    population = c(100, 200),
    contact_matrix = contact_matrix
  )

  expect_equal(lambda, c(0.3, 0.7), ignore_attr = TRUE)
})
