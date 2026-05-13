#' Coerce an object to an agepi contact matrix
#'
#' Converts common contact-matrix inputs to the agepi convention used by
#' [force_of_infection()]: rows are recipient age groups, columns are source age
#' groups, and `contact_matrix[a, b]` gives contacts made by recipient group
#' `a` with source group `b`.
#'
#' Supported inputs are numeric matrices, data frames that can be converted to
#' numeric matrices, socialmixr-like list objects with a numeric matrix at
#' `x$matrix`, and conmat-style long data frames with `age_group_from`,
#' `age_group_to`, and `contacts` columns. For conmat-style long data,
#' `age_group_to` is used for recipient rows and `age_group_from` is used for
#' source columns.
#'
#' Orientation is handled before the explicit `transpose` argument. That is,
#' `orientation = "source_recipient"` first transposes the coerced matrix into
#' recipient-source orientation, then `transpose = TRUE` transposes that result.
#'
#' @param x A supported contact-matrix input.
#' @param age_structure Optional valid age structure used to validate dimensions
#'   and age-group names.
#' @param orientation Orientation of matrix-like inputs. Use
#'   `"recipient_source"` when rows are recipients and columns are sources. Use
#'   `"source_recipient"` when rows are sources and columns are recipients.
#' @param transpose Logical scalar. If `TRUE`, transpose the matrix after
#'   `orientation` handling.
#'
#' @return A numeric square matrix in agepi recipient-source orientation.
#' @export
as_agepi_contact_matrix <- function(
  x,
  age_structure = NULL,
  orientation = c("recipient_source", "source_recipient"),
  transpose = FALSE
) {
  orientation <- match.arg(orientation)

  if (!is.logical(transpose) || length(transpose) != 1 || anyNA(transpose)) {
    stop("transpose must be a non-missing logical scalar.", call. = FALSE)
  }

  if (!is.null(age_structure)) {
    validate_age_structure(age_structure)
  }

  contact_matrix <- coerce_contact_matrix_input(x, age_structure)

  if (orientation == "source_recipient") {
    contact_matrix <- t(contact_matrix)
  }

  if (transpose) {
    contact_matrix <- t(contact_matrix)
  }

  contact_matrix <- apply_contact_matrix_age_structure(contact_matrix, age_structure)
  validate_contact_matrix(contact_matrix, age_structure)

  contact_matrix
}

coerce_contact_matrix_input <- function(x, age_structure) {
  if (is_conmat_long_data_frame(x)) {
    return(conmat_long_to_contact_matrix(x, age_structure))
  }

  if (is.list(x) && !is.data.frame(x) && !is.null(x$matrix)) {
    x <- x$matrix
  }

  if (is.matrix(x)) {
    return(coerce_matrix_like_contact_matrix(x))
  }

  if (is.data.frame(x)) {
    contact_matrix <- coerce_matrix_like_contact_matrix(as.matrix(x))
    dimnames(contact_matrix) <- NULL
    return(contact_matrix)
  }

  stop(
    "x must be a numeric matrix, numeric data frame, socialmixr-like list with x$matrix, ",
    "or conmat-style long data frame.",
    call. = FALSE
  )
}

coerce_matrix_like_contact_matrix <- function(x) {
  if (!is.numeric(x)) {
    stop("contact matrix input must be numeric.", call. = FALSE)
  }

  matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
}

is_conmat_long_data_frame <- function(x) {
  is.data.frame(x) &&
    all(c("age_group_from", "age_group_to", "contacts") %in% names(x))
}

conmat_long_to_contact_matrix <- function(x, age_structure) {
  age_group_from <- as.character(x$age_group_from)
  age_group_to <- as.character(x$age_group_to)
  contacts <- x$contacts

  if (!is.numeric(contacts)) {
    stop("conmat-style contacts column must be numeric.", call. = FALSE)
  }

  if (anyNA(age_group_from) || anyNA(age_group_to)) {
    stop("conmat-style age_group_from and age_group_to cannot contain missing values.", call. = FALSE)
  }

  pair_key <- paste(age_group_from, age_group_to, sep = "\r")
  if (any(duplicated(pair_key))) {
    stop("conmat-style input cannot contain duplicate age_group_from/age_group_to pairs.", call. = FALSE)
  }

  age_groups <- conmat_contact_matrix_age_groups(age_group_from, age_group_to, age_structure)
  contact_matrix <- matrix(
    NA_real_,
    nrow = length(age_groups),
    ncol = length(age_groups),
    dimnames = list(age_groups, age_groups)
  )

  row_index <- match(age_group_to, age_groups)
  col_index <- match(age_group_from, age_groups)

  if (anyNA(row_index) || anyNA(col_index)) {
    stop("conmat-style age groups must match age_structure$age_groups.", call. = FALSE)
  }

  contact_matrix[cbind(row_index, col_index)] <- contacts

  if (anyNA(contact_matrix)) {
    stop("conmat-style input must contain every age_group_from/age_group_to combination.", call. = FALSE)
  }

  contact_matrix
}

conmat_contact_matrix_age_groups <- function(age_group_from, age_group_to, age_structure) {
  if (!is.null(age_structure)) {
    return(age_structure$age_groups)
  }

  unique(c(age_group_to, age_group_from))
}

apply_contact_matrix_age_structure <- function(contact_matrix, age_structure) {
  if (is.null(age_structure)) {
    return(contact_matrix)
  }

  validate_contact_matrix_length(contact_matrix, age_structure$n_age_groups)
  age_groups <- age_structure$age_groups

  if (is.null(rownames(contact_matrix))) {
    rownames(contact_matrix) <- age_groups
  } else if (!identical(rownames(contact_matrix), age_groups)) {
    stop("contact_matrix rownames must match age_structure$age_groups exactly.", call. = FALSE)
  }

  if (is.null(colnames(contact_matrix))) {
    colnames(contact_matrix) <- age_groups
  } else if (!identical(colnames(contact_matrix), age_groups)) {
    stop("contact_matrix colnames must match age_structure$age_groups exactly.", call. = FALSE)
  }

  contact_matrix
}
