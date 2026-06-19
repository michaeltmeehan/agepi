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

#' Convert a socialmixr contact matrix output to an agepi contact matrix
#'
#' Optional, dependency-free adapter for objects shaped like the output of
#' `socialmixr::contact_matrix()`. The `socialmixr` package is not required for
#' core agepi functionality. This adapter accepts either a numeric matrix or a
#' list-like object with a numeric `matrix` element, then delegates validation,
#' orientation handling, optional transposition, and age-structure labelling to
#' [as_agepi_contact_matrix()].
#'
#' Returned matrices use the agepi convention: rows are recipient age groups and
#' columns are source age groups. No interpolation, reciprocity correction,
#' population balancing, or age-bin splitting is applied.
#'
#' @param x A numeric matrix or socialmixr contact-matrix-like list object with
#'   a numeric `matrix` element.
#' @param age_structure Optional valid age structure used to validate dimensions
#'   and age-group names.
#' @param orientation Orientation of the extracted matrix. Use
#'   `"recipient_source"` when rows are recipients and columns are sources. Use
#'   `"source_recipient"` when rows are sources and columns are recipients.
#' @param transpose Logical scalar. If `TRUE`, transpose the matrix after
#'   `orientation` handling.
#'
#' @return A numeric square matrix in agepi recipient-source orientation.
#'
#' @examples
#' age_structure <- AgeStructure(
#'   age_groups = c("0-4", "5-9"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, 9)
#' )
#'
#' socialmixr_like <- list(
#'   matrix = matrix(
#'     c(4, 2,
#'       1, 5),
#'     nrow = 2,
#'     byrow = TRUE
#'   )
#' )
#'
#' contact_matrix_from_socialmixr(
#'   socialmixr_like,
#'   age_structure = age_structure
#' )
#' @export
contact_matrix_from_socialmixr <- function(
  x,
  age_structure = NULL,
  orientation = c("recipient_source", "source_recipient"),
  transpose = FALSE
) {
  matrix_input <- x

  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(x$matrix)) {
      stop("socialmixr contact matrix output must contain a matrix element.", call. = FALSE)
    }
    matrix_input <- x$matrix
  }

  as_agepi_contact_matrix(
    matrix_input,
    age_structure = age_structure,
    orientation = orientation,
    transpose = transpose
  )
}

#' Build a published proxy contact matrix for a target age structure
#'
#' Currently supports a POLYMOD United Kingdom proxy through the optional
#' `socialmixr` package. The source matrix is built on a six-band age grid and
#' expanded to the target age structure by assigning each target age group to
#' the source band containing its lower bound. No reciprocity correction,
#' population balancing, or Kiribati-specific calibration is applied.
#'
#' @param age_structure Target valid age structure.
#' @param source Contact source. Currently only `"polymod_uk"` is supported.
#'
#' @return A numeric contact matrix in agepi recipient-source orientation. A
#'   metadata list is attached as attribute `"contact_source"`.
#' @export
contact_matrix_for_age_structure <- function(
  age_structure,
  source = c("polymod_uk")
) {
  source <- match.arg(source)
  validate_age_structure(age_structure)

  switch(
    source,
    polymod_uk = polymod_uk_contact_matrix_for_age_structure(age_structure)
  )
}

#' Convert conmat-style output to an agepi contact matrix
#'
#' Optional, dependency-free adapter for conmat-style long tables. The `conmat`
#' package is not required for core agepi functionality. This adapter expects
#' `age_group_from`, `age_group_to`, and `contacts` columns, then delegates
#' conversion, validation, orientation handling, optional transposition, and
#' age-structure labelling to [as_agepi_contact_matrix()].
#'
#' In conmat-style long input, `age_group_to` becomes recipient rows and
#' `age_group_from` becomes source columns. No interpolation, reciprocity
#' correction, population balancing, demographic prediction, or age-bin
#' splitting is applied.
#'
#' @param x A conmat-style long table or object coercible with `as.data.frame()`.
#' @param age_structure Optional valid age structure used to validate dimensions
#'   and age-group names.
#' @param orientation Orientation of the converted matrix. Use
#'   `"recipient_source"` when rows are recipients and columns are sources. Use
#'   `"source_recipient"` when rows are sources and columns are recipients.
#' @param transpose Logical scalar. If `TRUE`, transpose the matrix after
#'   `orientation` handling.
#'
#' @return A numeric square matrix in agepi recipient-source orientation.
#'
#' @examples
#' age_structure <- AgeStructure(
#'   age_groups = c("0-4", "5-9"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, 9)
#' )
#'
#' conmat_like <- data.frame(
#'   age_group_from = c("0-4", "5-9", "0-4", "5-9"),
#'   age_group_to = c("0-4", "0-4", "5-9", "5-9"),
#'   contacts = c(4, 2, 1, 5)
#' )
#'
#' contact_matrix_from_conmat(
#'   conmat_like,
#'   age_structure = age_structure
#' )
#' @export
contact_matrix_from_conmat <- function(
  x,
  age_structure = NULL,
  orientation = c("recipient_source", "source_recipient"),
  transpose = FALSE
) {
  orientation <- match.arg(orientation)

  conmat_data <- as.data.frame(x)
  required_columns <- c("age_group_from", "age_group_to", "contacts")
  missing_columns <- setdiff(required_columns, names(conmat_data))

  if (length(missing_columns) > 0) {
    stop(
      "conmat-style input must contain columns: ",
      paste(required_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  as_agepi_contact_matrix(
    conmat_data,
    age_structure = age_structure,
    orientation = orientation,
    transpose = transpose
  )
}

#' Construct a contact schedule
#'
#' Stores externally supplied contact matrices indexed by time. Matrices use the
#' agepi force-of-infection convention: rows are recipient age groups and
#' columns are source age groups. No interpolation, reciprocity correction, or
#' population balancing is applied.
#'
#' @param data Either a named list of contact matrices with names coercible to
#'   numeric times, or a long data frame with columns `time`, `age_group_from`,
#'   `age_group_to`, and `contacts`.
#' @param age_structure Age-structure object validated by
#'   [validate_age_structure()].
#'
#' @return An `agepi_contact_schedule` list.
#'
#' @examples
#' age_structure <- AgeStructure(
#'   age_groups = c("0-4", "5-9"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, 9)
#' )
#'
#' contacts <- list(
#'   "0" = matrix(c(4, 2, 1, 5), nrow = 2, byrow = TRUE),
#'   "1" = matrix(c(3, 2, 2, 4), nrow = 2, byrow = TRUE)
#' )
#'
#' schedule <- ContactSchedule(contacts, age_structure)
#' contact_matrix_at(schedule, time = 1)
#' @export
ContactSchedule <- function(data, age_structure) {
  validate_age_structure(age_structure)

  if (is_contact_schedule_long_data_frame(data)) {
    contacts <- contact_schedule_from_long_table(data, age_structure)
  } else if (is.list(data) && !is.data.frame(data)) {
    contacts <- contact_schedule_from_named_list(data, age_structure)
  } else {
    stop(
      "data must be a named list of contact matrices or a long data frame with ",
      "time, age_group_from, age_group_to, and contacts columns.",
      call. = FALSE
    )
  }

  times <- as.numeric(names(contacts))
  time_order <- order(times)
  times <- times[time_order]
  contacts <- contacts[time_order]
  names(contacts) <- as.character(times)

  contact_schedule <- list(
    contacts = contacts,
    age_structure = age_structure,
    times = times,
    age_groups = age_structure$age_groups,
    n_times = length(times),
    n_age_groups = age_structure$n_age_groups
  )

  class(contact_schedule) <- c("agepi_contact_schedule", "list")
  contact_schedule
}

#' Get a contact matrix at an exact time
#'
#' Retrieves an externally supplied contact matrix at one exact available time
#' point. No interpolation or nearest-time lookup is applied.
#'
#' @param contact_schedule An `agepi_contact_schedule` object.
#' @param time Exact time point to retrieve.
#'
#' @return Numeric contact matrix ordered by `contact_schedule$age_groups`.
#' @export
contact_matrix_at <- function(contact_schedule, time) {
  validate_agepi_contact_schedule(contact_schedule)
  validate_contact_schedule_time(time, contact_schedule$times)

  time_index <- match(time, contact_schedule$times)
  contact_schedule$contacts[[time_index]]
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

  if (anyNA(contacts) || any(!is.finite(contacts))) {
    stop("conmat-style contacts column cannot contain missing or non-finite values.", call. = FALSE)
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

  validate_contact_matrix(contact_matrix)
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

is_contact_schedule_long_data_frame <- function(x) {
  is.data.frame(x) &&
    all(c("time", "age_group_from", "age_group_to", "contacts") %in% names(x))
}

contact_schedule_from_named_list <- function(data, age_structure) {
  if (length(data) == 0) {
    stop("contact schedule list cannot be empty.", call. = FALSE)
  }

  if (is.null(names(data)) || any(names(data) == "")) {
    stop("contact schedule list must have one name per time point.", call. = FALSE)
  }

  times <- suppressWarnings(as.numeric(names(data)))
  if (anyNA(times) || any(!is.finite(times))) {
    stop("contact schedule list names must be finite numeric times.", call. = FALSE)
  }

  if (any(duplicated(times))) {
    stop("contact schedule cannot contain duplicate time points.", call. = FALSE)
  }

  contacts <- lapply(data, as_agepi_contact_matrix, age_structure = age_structure)
  names(contacts) <- as.character(times)
  contacts
}

contact_schedule_from_long_table <- function(data, age_structure) {
  required_columns <- c("time", "age_group_from", "age_group_to", "contacts")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "contact schedule data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(data$time)) {
    stop("contact schedule time must be numeric.", call. = FALSE)
  }

  if (anyNA(data$time) || any(!is.finite(data$time))) {
    stop("contact schedule time must be finite and non-missing.", call. = FALSE)
  }

  if (nrow(data) == 0) {
    stop("contact schedule long table cannot be empty.", call. = FALSE)
  }

  times <- sort(unique(data$time))
  contacts <- vector("list", length(times))
  names(contacts) <- as.character(times)

  for (i in seq_along(times)) {
    this_time <- times[i]
    time_data <- data[data$time == this_time, required_columns[-1], drop = FALSE]
    contacts[[i]] <- as_agepi_contact_matrix(time_data, age_structure = age_structure)
  }

  contacts
}

validate_agepi_contact_schedule <- function(contact_schedule) {
  if (!inherits(contact_schedule, "agepi_contact_schedule")) {
    stop("contact_schedule must be an agepi_contact_schedule object.", call. = FALSE)
  }

  required_fields <- c(
    "contacts",
    "age_structure",
    "times",
    "age_groups",
    "n_times",
    "n_age_groups"
  )
  missing_fields <- setdiff(required_fields, names(contact_schedule))
  if (length(missing_fields) > 0) {
    stop(
      "contact_schedule is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(contact_schedule)
}

validate_contact_schedule_time <- function(time, available_times) {
  if (!is.numeric(time) || length(time) != 1 || anyNA(time) || !is.finite(time)) {
    stop("time must be a finite numeric scalar.", call. = FALSE)
  }

  if (!time %in% available_times) {
    stop(
      "time is not available in contact_schedule: ",
      time,
      ". Available time point(s): ",
      paste(available_times, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(time)
}

polymod_uk_contact_matrix_for_age_structure <- function(age_structure) {
  if (!requireNamespace("socialmixr", quietly = TRUE)) {
    stop(
      "Package socialmixr is required for source = 'polymod_uk'.",
      call. = FALSE
    )
  }

  source_age_structure <- AgeStructure(
    age_groups = c("0-4", "5-14", "15-24", "25-44", "45-64", "65+"),
    lower_bounds = c(0, 5, 15, 25, 45, 65),
    upper_bounds = c(4, 14, 24, 44, 64, Inf)
  )

  data_environment <- new.env(parent = emptyenv())
  utils::data("polymod", package = "socialmixr", envir = data_environment)
  polymod <- get("polymod", envir = data_environment, inherits = FALSE)
  socialmixr_matrix <- suppressWarnings(socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age_limits = source_age_structure$lower_bounds,
    symmetric = FALSE
  ))

  source_matrix <- socialmixr_matrix$matrix
  dimnames(source_matrix) <- list(
    source_age_structure$age_groups,
    source_age_structure$age_groups
  )
  source_matrix <- contact_matrix_from_socialmixr(
    source_matrix,
    age_structure = source_age_structure
  )

  age_band <- target_age_groups_to_source_bands(age_structure, source_age_structure)
  contact_matrix <- source_matrix[age_band, age_band]
  dimnames(contact_matrix) <- list(age_structure$age_groups, age_structure$age_groups)
  validate_contact_matrix(contact_matrix, age_structure)

  attr(contact_matrix, "contact_source") <- list(
    source_label = paste(
      "POLYMOD United Kingdom empirical social-contact survey matrix from",
      "socialmixr::contact_matrix(); used as a published proxy, not a",
      "Kiribati-specific matrix"
    ),
    source_reference = paste(
      "Mossong et al. 2008 / POLYMOD, distributed through socialmixr's",
      "bundled polymod dataset"
    ),
    source_age_grid = paste(source_age_structure$age_groups, collapse = ", "),
    model_age_grid = paste(age_structure$age_groups, collapse = ", "),
    expansion_note = paste(
      "The validated six-age-band source matrix is expanded to the target",
      "model grid by assigning each target age group to its source age band",
      "and using constant contacts within each source band."
    ),
    limitation = paste(
      "This proxy is European and pre-pandemic; it is not calibrated to",
      "Kiribati household structure, crowding, school attendance, or",
      "TB-relevant prolonged indoor exposure."
    )
  )

  contact_matrix
}

target_age_groups_to_source_bands <- function(target_age_structure, source_age_structure) {
  vapply(
    target_age_structure$lower_bounds,
    function(lower_bound) {
      source_index <- which(
        source_age_structure$lower_bounds <= lower_bound &
          source_age_structure$upper_bounds >= lower_bound
      )
      if (length(source_index) != 1) {
        stop("target age group lower bounds must map to exactly one source contact age band.", call. = FALSE)
      }
      source_age_structure$age_groups[source_index]
    },
    character(1)
  )
}
