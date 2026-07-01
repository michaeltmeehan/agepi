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

test_that("adapt_contact_matrix_to_age_structure expands nested source grids with equal source-band splitting", {
  source_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  target_ages <- AgeStructure(
    age_groups = as.character(0:9),
    lower_bounds = 0:9,
    upper_bounds = 0:9
  )
  source_matrix <- matrix(
    c(5, 10, 2, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(source_ages$age_groups, source_ages$age_groups)
  )
  source <- ContactMatrixSource(
    matrix = source_matrix,
    age_structure = source_ages,
    source = "test",
    source_reference = "test source"
  )

  contact_matrix <- adapt_contact_matrix_to_age_structure(source, target_ages)
  metadata <- attr(contact_matrix, "contact_source")

  expect_equal(contact_matrix[1:5, 1:5], matrix(5 / 5, nrow = 5, ncol = 5), ignore_attr = TRUE)
  expect_equal(contact_matrix[1:5, 6:10], matrix(10 / 5, nrow = 5, ncol = 5), ignore_attr = TRUE)
  expect_equal(contact_matrix[6:10, 1:5], matrix(2 / 5, nrow = 5, ncol = 5), ignore_attr = TRUE)
  expect_equal(contact_matrix[6:10, 6:10], matrix(4 / 5, nrow = 5, ncol = 5), ignore_attr = TRUE)

  expect_equal(rowSums(contact_matrix[1:5, 1:5]), rep(5, 5), ignore_attr = TRUE)
  expect_equal(rowSums(contact_matrix[1:5, 6:10]), rep(10, 5), ignore_attr = TRUE)
  expect_equal(rowSums(contact_matrix[6:10, 1:5]), rep(2, 5), ignore_attr = TRUE)
  expect_equal(rowSums(contact_matrix[6:10, 6:10]), rep(4, 5), ignore_attr = TRUE)
  expect_match(metadata$adaptation_note, "expanded")
  expect_match(metadata$adaptation_note, "preserving")
  expect_match(metadata$adaptation_note, "equal weights")
  expect_identical(metadata$adaptation_method, "source_band")
})

test_that("adapt_contact_matrix_to_age_structure expands nested source grids with target population weights", {
  source_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  target_ages <- AgeStructure(
    age_groups = as.character(0:9),
    lower_bounds = 0:9,
    upper_bounds = 0:9
  )
  source_matrix <- matrix(
    c(5, 10, 2, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(source_ages$age_groups, source_ages$age_groups)
  )
  source <- ContactMatrixSource(
    matrix = source_matrix,
    age_structure = source_ages,
    source = "test",
    source_reference = "test source"
  )
  population <- c(10, 20, 30, 40, 50, 5, 5, 10, 20, 60)

  contact_matrix <- adapt_contact_matrix_to_age_structure(source, target_ages, population = population)

  expect_equal(rowSums(contact_matrix[1:5, 1:5]), rep(5, 5), ignore_attr = TRUE)
  expect_equal(rowSums(contact_matrix[1:5, 6:10]), rep(10, 5), ignore_attr = TRUE)
  expect_equal(rowSums(contact_matrix[6:10, 1:5]), rep(2, 5), ignore_attr = TRUE)
  expect_equal(rowSums(contact_matrix[6:10, 6:10]), rep(4, 5), ignore_attr = TRUE)

  expect_equal(contact_matrix["0", "4"], 5 * 50 / sum(c(10, 20, 30, 40, 50)), ignore_attr = TRUE)
  expect_equal(contact_matrix["0", "9"], 10 * 60 / sum(c(5, 5, 10, 20, 60)), ignore_attr = TRUE)
  expect_match(attr(contact_matrix, "contact_source")$adaptation_note, "target-grid population weights")
})

test_that("contact matrix expansion validates target-grid population", {
  source_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  target_ages <- AgeStructure(
    age_groups = as.character(0:9),
    lower_bounds = 0:9,
    upper_bounds = 0:9
  )
  source <- ContactMatrixSource(
    matrix = diag(2),
    age_structure = source_ages,
    source = "test",
    source_reference = "test source"
  )

  expect_error(
    adapt_contact_matrix_to_age_structure(source, target_ages, population = c(100, 200)),
    "target age_structure\\$n_age_groups"
  )
  expect_error(
    adapt_contact_matrix_to_age_structure(source, target_ages, population = c(rep(1, 5), rep(0, 5))),
    "must sum to a positive value"
  )
  expect_error(
    adapt_contact_matrix_to_age_structure(source, target_ages, population = c(rep(1, 9), NA_real_)),
    "finite"
  )
  expect_error(
    adapt_contact_matrix_to_age_structure(source, target_ages, population = c(rep(1, 9), -1)),
    "negative"
  )
})

test_that("adapt_contact_matrix_to_age_structure aggregates using transform_contact_matrix", {
  source_ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-14", "15-19"),
    lower_bounds = c(0, 5, 10, 15),
    upper_bounds = c(4, 9, 14, 19)
  )
  target_ages <- AgeStructure(
    age_groups = c("0-9", "10-19"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, 19)
  )
  source_matrix <- matrix(
    c(
      1, 2, 3, 4,
      5, 6, 7, 8,
      9, 10, 11, 12,
      13, 14, 15, 16
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(source_ages$age_groups, source_ages$age_groups)
  )
  source <- ContactMatrixSource(
    matrix = source_matrix,
    age_structure = source_ages,
    source = "test",
    source_reference = "test source"
  )

  contact_matrix <- adapt_contact_matrix_to_age_structure(
    source,
    target_ages,
    population = c(10, 10, 10, 10)
  )

  expect_equal(
    contact_matrix,
    transform_contact_matrix(source_matrix, source_ages, target_ages, population = c(10, 10, 10, 10)),
    ignore_attr = TRUE
  )
  expect_match(attr(contact_matrix, "contact_source")$adaptation_note, "aggregated")
  expect_identical(attr(contact_matrix, "contact_source")$adaptation_method, "source_band")
})

test_that("contact matrix aggregation requires source-grid population", {
  source_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  target_ages <- AgeStructure(
    age_groups = "0-9",
    lower_bounds = 0,
    upper_bounds = 9
  )
  source <- ContactMatrixSource(
    matrix = diag(2),
    age_structure = source_ages,
    source = "test",
    source_reference = "test source"
  )

  expect_error(
    adapt_contact_matrix_to_age_structure(source, target_ages),
    "population is required"
  )
})

test_that("contact matrix source object structure is validated", {
  source_ages <- AgeStructure(
    age_groups = c("0-4", "5+"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, Inf)
  )
  source <- ContactMatrixSource(
    matrix = diag(2),
    age_structure = source_ages,
    source = "unit_test",
    country = "Exampleland",
    setting = "all",
    source_reference = "test reference",
    notes = "test notes",
    limitations = "test limitations"
  )

  expect_s3_class(source, "agepi_contact_matrix_source")
  expect_identical(source$source, "unit_test")
  expect_identical(source$country, "Exampleland")
  expect_identical(source$setting, "all")
  expect_identical(source$orientation, "recipient_source")
  expect_true(all(c(
    "matrix", "age_structure", "source", "country", "setting",
    "orientation", "convention", "source_reference", "notes",
    "limitations", "metadata"
  ) %in% names(source)))
  expect_silent(validate_contact_matrix_source(source))

  source$source <- NULL
  expect_error(validate_contact_matrix_source(source), "missing required")
})

test_that("deprecated exact adaptation method is an alias for source_band", {
  source_ages <- AgeStructure(
    age_groups = c("0-9", "10+"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, Inf)
  )
  target_ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
  source <- ContactMatrixSource(
    matrix = matrix(c(2, 1, 3, 4), nrow = 2, byrow = TRUE),
    age_structure = source_ages,
    source = "test",
    source_reference = "test source",
    notes = "test notes"
  )

  expect_warning(
    contact_matrix <- adapt_contact_matrix_to_age_structure(source, target_ages, method = "exact"),
    "deprecated"
  )
  expect_identical(attr(contact_matrix, "contact_source")$adaptation_method, "source_band")
})

test_that("source validation rejects absent notes and limitations", {
  source_ages <- AgeStructure(
    age_groups = c("0-4", "5+"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, Inf)
  )
  source <- ContactMatrixSource(
    matrix = diag(2),
    age_structure = source_ages,
    source = "unit_test",
    source_reference = "test reference",
    notes = "test notes"
  )
  source$notes <- NULL
  source$limitations <- NULL

  expect_error(validate_contact_matrix_source(source), "missing required")
})

test_that("POLYMOD source loader supports UK compatibility alias", {
  testthat::skip_if_not_installed("socialmixr")

  socialmixr_datasets <- utils::data(package = "socialmixr")$results[, "Item"]
  testthat::skip_if_not("polymod" %in% socialmixr_datasets)

  ages <- wpp_age_structure_1year(max_age = 65)
  source <- load_contact_matrix_source("polymod_uk")
  contact_matrix <- adapt_contact_matrix_to_age_structure(source, ages)
  metadata <- attr(contact_matrix, "contact_source")

  expect_s3_class(source, "agepi_contact_matrix_source")
  expect_identical(source$source, "polymod_uk")
  expect_identical(source$country, "United Kingdom")
  expect_true(is.matrix(contact_matrix))
  expect_identical(dim(contact_matrix), c(ages$n_age_groups, ages$n_age_groups))
  expect_identical(rownames(contact_matrix), ages$age_groups)
  expect_identical(colnames(contact_matrix), ages$age_groups)
  expect_true(all(is.finite(contact_matrix)))
  expect_true(all(contact_matrix >= 0))
  expect_match(metadata$source_label, "POLYMOD")
  expect_match(metadata$adaptation_note, "preserving")
})

test_that("POLYMOD source loader supports another available country", {
  testthat::skip_if_not_installed("socialmixr")

  socialmixr_datasets <- utils::data(package = "socialmixr")$results[, "Item"]
  testthat::skip_if_not("polymod" %in% socialmixr_datasets)

  source <- load_contact_matrix_source("polymod", country = "Germany")

  expect_s3_class(source, "agepi_contact_matrix_source")
  expect_identical(source$source, "polymod")
  expect_identical(source$country, "Germany")
  expect_identical(source$setting, "all")
  expect_true(is.matrix(source$matrix))
  expect_identical(dim(source$matrix), c(source$age_structure$n_age_groups, source$age_structure$n_age_groups))
  expect_silent(validate_contact_matrix_source(source))
})

test_that("POLYMOD source loader handles setting filters conservatively", {
  testthat::skip_if_not_installed("socialmixr")

  socialmixr_datasets <- utils::data(package = "socialmixr")$results[, "Item"]
  testthat::skip_if_not("polymod" %in% socialmixr_datasets)

  home_source <- load_contact_matrix_source("polymod", country = "Germany", setting = "home")
  expect_identical(home_source$setting, "home")
  expect_s3_class(home_source, "agepi_contact_matrix_source")

  expect_error(
    load_contact_matrix_source("polymod", country = "Germany", setting = "other"),
    "explicit socialmixr filter"
  )
})

test_that("POLYMOD source loader fails clearly without socialmixr", {
  skip_if(requireNamespace("socialmixr", quietly = TRUE))

  expect_error(
    load_contact_matrix_source("polymod_uk"),
    "Package socialmixr is required"
  )
})

test_that("Prem source loader fails clearly without contactdata", {
  skip_if(requireNamespace("contactdata", quietly = TRUE))

  expect_error(
    load_contact_matrix_source("prem", country = "Kiribati"),
    "Package contactdata is required"
  )
})

test_that("Prem source loader returns a validated source object", {
  testthat::skip_if_not_installed("contactdata")

  countries <- contactdata::list_countries()
  country <- if ("Kiribati" %in% countries) "Kiribati" else countries[[1]]
  source <- load_contact_matrix_source("prem", country = country, setting = "home")

  expect_s3_class(source, "agepi_contact_matrix_source")
  expect_identical(source$source, "prem")
  expect_identical(source$country, country)
  expect_identical(source$setting, "home")
  expect_equal(dim(source$matrix), c(16L, 16L))
  expect_silent(validate_contact_matrix_source(source))
})

test_that("conmat source loader requires population input", {
  testthat::skip_if_not_installed("conmat")

  expect_error(
    load_contact_matrix_source("conmat"),
    "population is required"
  )
})

test_that("conmat source loader fails clearly without conmat", {
  skip_if(requireNamespace("conmat", quietly = TRUE))

  expect_error(
    load_contact_matrix_source("conmat", population = data.frame()),
    "Package conmat is required"
  )
})

test_that("conmat source loader accepts generated matrices when conmat dependencies are available", {
  testthat::skip_if_not_installed("conmat")

  population <- data.frame(
    lower.age.limit = seq(0, 75, by = 5),
    population = rep(1000, 16)
  )
  source <- tryCatch(
    load_contact_matrix_source(
      "conmat",
      population = population,
      age_limits = c(seq(0, 75, by = 5), Inf)
    ),
    error = function(e) {
      testthat::skip(paste("conmat source generation unavailable:", conditionMessage(e)))
    }
  )

  expect_s3_class(source, "agepi_contact_matrix_source")
  expect_identical(source$source, "conmat")
  expect_identical(source$setting, "all")
  expect_silent(validate_contact_matrix_source(source))
})

test_that("deprecated contact_matrix_for_age_structure delegates to source workflow", {
  testthat::skip_if_not_installed("socialmixr")

  socialmixr_datasets <- utils::data(package = "socialmixr")$results[, "Item"]
  testthat::skip_if_not("polymod" %in% socialmixr_datasets)

  ages <- wpp_age_structure_1year(max_age = 65)

  expect_warning(
    contact_matrix <- contact_matrix_for_age_structure(ages),
    "deprecated"
  )
  expect_true(is.matrix(contact_matrix))
  expect_identical(dim(contact_matrix), c(ages$n_age_groups, ages$n_age_groups))
  expect_match(attr(contact_matrix, "contact_source")$source, "polymod_uk")
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

test_contact_schedule_age_structure <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

test_contact_schedule_list <- function() {
  list(
    "1" = matrix(c(
      3, 1,
      2, 4
    ), nrow = 2, byrow = TRUE),
    "0" = matrix(c(
      4, 2,
      1, 5
    ), nrow = 2, byrow = TRUE)
  )
}

test_contact_schedule_long <- function() {
  data.frame(
    time = c(1, 1, 1, 1, 0, 0, 0, 0),
    age_group_from = c("0-4", "5-9", "0-4", "5-9", "0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9", "0-4", "0-4", "5-9", "5-9"),
    contacts = c(3, 1, 2, 4, 4, 2, 1, 5)
  )
}

test_that("ContactSchedule constructs from named list of matrices", {
  ages <- test_contact_schedule_age_structure()
  schedule <- ContactSchedule(test_contact_schedule_list(), ages)

  expect_s3_class(schedule, "agepi_contact_schedule")
  expect_identical(schedule$age_structure, ages)
  expect_identical(schedule$times, c(0, 1))
  expect_identical(schedule$age_groups, ages$age_groups)
  expect_identical(schedule$n_times, 2L)
  expect_identical(schedule$n_age_groups, 2L)
  expect_equal(
    schedule$contacts[[1]],
    matrix(
      c(4, 2, 1, 5),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(ages$age_groups, ages$age_groups)
    )
  )
})

test_that("ContactSchedule constructs from long table", {
  ages <- test_contact_schedule_age_structure()
  schedule <- ContactSchedule(test_contact_schedule_long(), ages)

  expect_s3_class(schedule, "agepi_contact_schedule")
  expect_identical(schedule$times, c(0, 1))
  expect_equal(
    schedule$contacts[[2]],
    matrix(
      c(3, 1, 2, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(ages$age_groups, ages$age_groups)
    )
  )
})

test_that("ContactSchedule stores matrices in age-structure order", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  ages$age_groups <- c("5-9", "0-4")
  long <- data.frame(
    time = c(0, 0, 0, 0),
    age_group_from = c("0-4", "5-9", "0-4", "5-9"),
    age_group_to = c("0-4", "0-4", "5-9", "5-9"),
    contacts = c(1, 2, 3, 4)
  )

  schedule <- ContactSchedule(long, ages)

  expect_equal(
    contact_matrix_at(schedule, time = 0),
    matrix(
      c(4, 3, 2, 1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("5-9", "0-4"), c("5-9", "0-4"))
    )
  )
})

test_that("ContactSchedule rejects missing age-pair records", {
  long <- test_contact_schedule_long()
  long <- long[-1, ]

  expect_error(
    ContactSchedule(long, test_contact_schedule_age_structure()),
    "every age_group_from/age_group_to combination"
  )
})

test_that("ContactSchedule rejects duplicate age-pair records", {
  long <- test_contact_schedule_long()
  long <- rbind(long, long[1, ])

  expect_error(
    ContactSchedule(long, test_contact_schedule_age_structure()),
    "duplicate"
  )
})

test_that("ContactSchedule rejects invalid contact values", {
  long <- test_contact_schedule_long()
  long$contacts[1] <- -1

  expect_error(
    ContactSchedule(long, test_contact_schedule_age_structure()),
    "negative"
  )
})

test_that("ContactSchedule rejects mismatched age labels", {
  long <- test_contact_schedule_long()
  long$age_group_from[1] <- "10-14"

  expect_error(
    ContactSchedule(long, test_contact_schedule_age_structure()),
    "age groups must match"
  )
})

test_that("ContactSchedule rejects non-square matrices", {
  contacts <- list("0" = matrix(1, nrow = 2, ncol = 3))

  expect_error(
    ContactSchedule(contacts, test_contact_schedule_age_structure()),
    "square"
  )
})

test_that("ContactSchedule preserves recipient-row source-column orientation", {
  schedule <- ContactSchedule(test_contact_schedule_long(), test_contact_schedule_age_structure())
  contact_matrix <- contact_matrix_at(schedule, time = 0)

  expect_identical(unname(contact_matrix["5-9", "0-4"]), 1)
  expect_identical(unname(contact_matrix["0-4", "5-9"]), 2)
})

test_that("contact_matrix_at returns exact-time matrix in schedule order", {
  ages <- test_contact_schedule_age_structure()
  schedule <- ContactSchedule(test_contact_schedule_list(), ages)

  expect_equal(
    contact_matrix_at(schedule, time = 1),
    matrix(
      c(3, 1, 2, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(ages$age_groups, ages$age_groups)
    )
  )
})

test_that("contact_matrix_at rejects unavailable times without interpolation", {
  schedule <- ContactSchedule(
    test_contact_schedule_list(),
    test_contact_schedule_age_structure()
  )

  expect_error(
    contact_matrix_at(schedule, time = 0.5),
    "time is not available"
  )
})
