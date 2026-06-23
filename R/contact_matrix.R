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

#' Load a contact matrix source
#'
#' Loads a contact matrix with source metadata on its native/source age grid.
#' The returned matrix is not adapted to a model age grid; use
#' [adapt_contact_matrix_to_age_structure()] for age-grid adaptation.
#'
#' Supported sources are:
#' - `"polymod"`: empirical POLYMOD matrices through optional `socialmixr`;
#' - `"polymod_uk"`: compatibility alias for `source = "polymod"` and
#'   `country = "United Kingdom"`;
#' - `"prem"`: Prem et al. synthetic country matrices through optional
#'   `contactdata`;
#' - `"conmat"`: conmat-generated matrices from caller-supplied population data
#'   through optional `conmat`.
#'
#' @param source Contact source.
#' @param country Country/location name, when applicable.
#' @param setting Contact setting. For `"prem"` and `"conmat"`, one of
#'   `"all"`, `"home"`, `"school"`, `"work"`, or `"other"`. For `"polymod"`,
#'   `"all"` applies no contact filter; `"home"`, `"school"`, and `"work"`
#'   filter the POLYMOD contact records when the setting variables are
#'   available. POLYMOD `"other"` spans multiple columns and requires an
#'   explicit `filter`.
#' @param age_limits Lower age limits for source bins. Defaults to six broad
#'   reporting bands for POLYMOD and five-year bins for Prem/conmat.
#' @param geographic_setting Prem/contactdata geographic setting: `"all"`,
#'   `"rural"`, or `"urban"`.
#' @param data_source Prem/contactdata source vintage: `"2020"` or `"2017"`.
#' @param filter Optional socialmixr filter list for POLYMOD. If supplied, this
#'   overrides `setting`.
#' @param symmetric Whether to request a symmetric socialmixr POLYMOD matrix.
#' @param population Population input required for `source = "conmat"`.
#' @param per_capita_household_size Optional conmat household-size adjustment.
#'
#' @return An `agepi_contact_matrix_source` list with `matrix`,
#'   `age_structure`, `source`, optional `country`, optional `setting`,
#'   `orientation`, `convention`, `source_reference`, `notes`, `limitations`,
#'   and source-specific `metadata`.
#' @export
load_contact_matrix_source <- function(source = c("polymod_uk", "polymod", "prem", "conmat"),
                                       country = NULL,
                                       setting = "all",
                                       age_limits = NULL,
                                       geographic_setting = "all",
                                       data_source = "2020",
                                       filter = NULL,
                                       symmetric = FALSE,
                                       population = NULL,
                                       per_capita_household_size = NULL) {
  source <- match.arg(source)

  switch(
    source,
    polymod_uk = load_polymod_contact_matrix_source(
      country = "United Kingdom",
      setting = setting,
      age_limits = age_limits,
      filter = filter,
      symmetric = symmetric,
      alias = "polymod_uk"
    ),
    polymod = load_polymod_contact_matrix_source(
      country = country,
      setting = setting,
      age_limits = age_limits,
      filter = filter,
      symmetric = symmetric,
      alias = NULL
    ),
    prem = load_prem_contact_matrix_source(
      country = country,
      setting = setting,
      age_limits = age_limits,
      geographic_setting = geographic_setting,
      data_source = data_source
    ),
    conmat = load_conmat_contact_matrix_source(
      population = population,
      setting = setting,
      age_limits = age_limits,
      per_capita_household_size = per_capita_household_size
    )
  )
}

#' Adapt a contact matrix source to a target age structure
#'
#' Adapts a loaded source contact matrix to a requested age grid using
#' source-band assumptions. Fine-to-coarse aggregation uses
#' [transform_contact_matrix()] with recipient-population weighting and
#' therefore requires a source-grid population vector. Coarse-to-fine expansion
#' assigns each target age group to the source age band that fully contains it
#' and uses constant contacts within the source band.
#'
#' @param source_contact_matrix An `agepi_contact_matrix_source` object returned
#'   by [load_contact_matrix_source()].
#' @param age_structure Target valid age structure.
#' @param population Optional source-grid population vector required for
#'   fine-to-coarse aggregation.
#' @param method Adaptation method. Currently only `"source_band"` is
#'   supported. It uses recipient-population weighting for fine-to-coarse
#'   aggregation and constant contacts within each source age band for
#'   coarse-to-fine expansion. `"exact"` is accepted as a deprecated alias for
#'   `"source_band"`; it is exact only under the same source-band/piecewise
#'   constant assumption.
#'
#' @return A numeric contact matrix in agepi recipient-source orientation. A
#'   metadata list is attached as attribute `"contact_source"`.
#' @export
adapt_contact_matrix_to_age_structure <- function(source_contact_matrix,
                                                  age_structure,
                                                  population = NULL,
                                                  method = c("source_band", "exact")) {
  method <- match.arg(method)
  if (identical(method, "exact")) {
    warning(
      "method = 'exact' is deprecated; use method = 'source_band'. ",
      "'exact' meant exact only under the source-band/piecewise-constant assumption.",
      call. = FALSE
    )
    method <- "source_band"
  }
  validate_contact_matrix_source(source_contact_matrix)
  validate_age_structure(age_structure)

  source_age_structure <- source_contact_matrix$age_structure
  source_matrix <- validate_contact_matrix(source_contact_matrix$matrix, source_age_structure)
  metadata <- source_contact_matrix$metadata

  if (identical(source_age_structure$age_groups, age_structure$age_groups) &&
      identical(source_age_structure$lower_bounds, age_structure$lower_bounds) &&
      identical(source_age_structure$upper_bounds, age_structure$upper_bounds)) {
    contact_matrix <- as_agepi_contact_matrix(source_matrix, age_structure = age_structure)
    metadata$adaptation_method <- method
    metadata$adaptation_note <- "Source contact matrix already matches the target age grid."
    attr(contact_matrix, "contact_source") <- contact_source_metadata(source_contact_matrix, metadata)
    return(contact_matrix)
  }

  mapping <- AgeGridMapping(source_age_structure, age_structure, open_ended = "include")
  if (isTRUE(mapping$can_aggregate)) {
    if (is.null(population)) {
      stop("population is required when aggregating a source contact matrix to a coarser age grid.", call. = FALSE)
    }
    contact_matrix <- transform_contact_matrix(
      source_matrix,
      from_age_structure = source_age_structure,
      to_age_structure = age_structure,
      population = population
    )
    metadata$adaptation_note <- paste(
      "The source matrix was aggregated to the target age grid using",
      "recipient-population weighting."
    )
  } else if (isTRUE(mapping$can_expand)) {
    source_index <- mapping$target_to_source_index
    contact_matrix <- source_matrix[source_index, source_index, drop = FALSE]
    dimnames(contact_matrix) <- list(age_structure$age_groups, age_structure$age_groups)
    contact_matrix <- validate_contact_matrix(contact_matrix, age_structure)
    metadata$adaptation_note <- paste(
      "The source matrix was expanded to the target age grid by assigning",
      "each target age group to its containing source age band and using",
      "constant contacts within each source band."
    )
  } else {
    stop("source and target age grids are not compatible for source-band contact-matrix adaptation.", call. = FALSE)
  }

  metadata$adaptation_method <- method
  metadata$model_age_grid <- paste(age_structure$age_groups, collapse = ", ")
  attr(contact_matrix, "contact_source") <- contact_source_metadata(source_contact_matrix, metadata)
  contact_matrix
}

#' Build a published proxy contact matrix for a target age structure
#'
#' Deprecated compatibility wrapper. Prefer [load_contact_matrix_source()]
#' followed by [adapt_contact_matrix_to_age_structure()] so source loading and
#' age-grid adaptation are explicit.
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
  .Deprecated("adapt_contact_matrix_to_age_structure")
  source_contact_matrix <- load_contact_matrix_source(source = source)
  adapt_contact_matrix_to_age_structure(source_contact_matrix, age_structure)
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

load_polymod_contact_matrix_source <- function(country,
                                               setting,
                                               age_limits,
                                               filter,
                                               symmetric,
                                               alias) {
  if (!requireNamespace("socialmixr", quietly = TRUE)) {
    stop(
      "Package socialmixr is required for source = 'polymod'. Install socialmixr or choose another source.",
      call. = FALSE
    )
  }

  if (is.null(country)) {
    country <- "United Kingdom"
  }
  country <- validate_contact_source_scalar(country, "country")
  setting <- validate_contact_source_scalar(setting, "setting")
  if (!setting %in% c("all", "home", "school", "work", "other")) {
    stop("setting must be one of: all, home, school, work, other.", call. = FALSE)
  }
  if (!is.logical(symmetric) || length(symmetric) != 1 || anyNA(symmetric)) {
    stop("symmetric must be TRUE or FALSE.", call. = FALSE)
  }

  data_environment <- new.env(parent = emptyenv())
  utils::data("polymod", package = "socialmixr", envir = data_environment)
  polymod <- get("polymod", envir = data_environment, inherits = FALSE)

  available_countries <- sort(unique(as.character(polymod$participants$country)))
  if (!country %in% available_countries) {
    stop(
      "POLYMOD country is not available: ",
      country,
      ". Available countries: ",
      paste(available_countries, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(age_limits)) {
    age_limits <- c(0, 5, 15, 25, 45, 65)
  }
  source_age_structure <- contact_source_age_structure_from_limits(age_limits)

  socialmixr_filter <- filter
  setting_note <- NULL
  if (is.null(socialmixr_filter) && !identical(setting, "all")) {
    if (all(c("cnt_home", "cnt_school", "cnt_work", "cnt_transport", "cnt_leisure", "cnt_otherplace") %in% names(polymod$contacts))) {
      if (identical(setting, "other")) {
        stop(
          "POLYMOD setting = 'other' spans multiple contact-location columns; ",
          "supply an explicit socialmixr filter instead.",
          call. = FALSE
        )
      }
      setting_columns <- switch(
        setting,
        home = "cnt_home",
        school = "cnt_school",
        work = "cnt_work"
      )
      socialmixr_filter <- stats::setNames(as.list(rep(1, length(setting_columns))), setting_columns)
      setting_note <- paste("POLYMOD contacts were filtered using column(s):", paste(setting_columns, collapse = ", "))
    } else {
      stop("POLYMOD setting filters are not available in this socialmixr dataset.", call. = FALSE)
    }
  }

  socialmixr_matrix <- suppressWarnings(socialmixr::contact_matrix(
    polymod,
    countries = country,
    age_limits = source_age_structure$lower_bounds,
    filter = socialmixr_filter,
    symmetric = symmetric
  ))

  source_matrix <- socialmixr_matrix$matrix
  source_matrix <- contact_source_matrix_on_age_structure(source_matrix, source_age_structure)

  notes <- c(
    paste(
      "POLYMOD empirical social-contact survey matrix from",
      "socialmixr::contact_matrix()."
    ),
    if (!is.null(setting_note)) setting_note else NULL,
    "No reciprocity correction or age-grid adaptation is applied by the source loader."
  )

  ContactMatrixSource(
    matrix = source_matrix,
    age_structure = source_age_structure,
    source = if (is.null(alias)) "polymod" else alias,
    country = country,
    setting = setting,
    orientation = "recipient_source",
    source_reference = paste(
      "Mossong et al. 2008 / POLYMOD, distributed through socialmixr's",
      "bundled polymod dataset"
    ),
    notes = notes,
    limitations = paste(
      "POLYMOD countries are European and pre-pandemic; matrices are not",
      "calibrated to non-survey populations unless the caller makes an",
      "explicit proxy assumption."
    ),
    metadata = list(
      source_label = paste(
        "POLYMOD",
        country,
        "empirical social-contact survey matrix from socialmixr::contact_matrix()"
      ),
      symmetric = symmetric,
      filter = socialmixr_filter
    )
  )
}

load_prem_contact_matrix_source <- function(country,
                                            setting,
                                            age_limits,
                                            geographic_setting,
                                            data_source) {
  if (!requireNamespace("contactdata", quietly = TRUE)) {
    stop(
      "Package contactdata is required for source = 'prem'. Install contactdata or choose another source.",
      call. = FALSE
    )
  }

  country <- validate_contact_source_scalar(country, "country")
  setting <- match.arg(setting, c("all", "home", "school", "work", "other"))
  geographic_setting <- match.arg(geographic_setting, c("all", "rural", "urban"))
  data_source <- match.arg(as.character(data_source), c("2020", "2017"))

  available_countries <- contactdata::list_countries(
    geographic_setting = geographic_setting,
    data_source = data_source
  )
  if (!country %in% available_countries) {
    stop(
      "Prem/contactdata country is not available: ",
      country,
      ".",
      call. = FALSE
    )
  }

  source_matrix <- contactdata::contact_matrix(
    country = country,
    location = setting,
    geographic_setting = geographic_setting,
    data_source = data_source
  )

  if (is.null(age_limits)) {
    age_limits <- seq(0, 75, by = 5)
  }
  source_age_structure <- contact_source_age_structure_from_limits(age_limits)
  source_matrix <- as_agepi_contact_matrix(
    contact_source_matrix_on_age_structure(source_matrix, source_age_structure),
    age_structure = source_age_structure
  )

  ContactMatrixSource(
    matrix = source_matrix,
    age_structure = source_age_structure,
    source = "prem",
    country = country,
    setting = setting,
    orientation = "recipient_source",
    source_reference = paste(
      "Prem et al. 2017 and Prem et al. 2021 synthetic contact matrices,",
      "accessed through the optional contactdata package."
    ),
    notes = paste(
      "Prem/contactdata matrix loaded on its native five-year age grid;",
      "no age-grid adaptation, reciprocity correction, or calibration is",
      "applied by agepi."
    ),
    limitations = paste(
      "Prem matrices are synthetic country-level projections and may not",
      "represent subnational, temporal, or disease-specific contact structure."
    ),
    metadata = list(
      source_label = paste("Prem/contactdata", country, setting, "contact matrix"),
      geographic_setting = geographic_setting,
      data_source = data_source
    )
  )
}

load_conmat_contact_matrix_source <- function(population,
                                              setting,
                                              age_limits,
                                              per_capita_household_size) {
  if (!requireNamespace("conmat", quietly = TRUE)) {
    stop(
      "Package conmat is required for source = 'conmat'. Install conmat or choose another source.",
      call. = FALSE
    )
  }
  if (is.null(population)) {
    stop("population is required for source = 'conmat'.", call. = FALSE)
  }
  setting <- match.arg(setting, c("all", "home", "school", "work", "other"))
  if (is.null(age_limits)) {
    age_limits <- c(seq(0, 75, by = 5), Inf)
  }
  source_age_structure <- contact_source_age_structure_from_limits(age_limits)

  conmat_result <- tryCatch(
    conmat::extrapolate_polymod(
      population = population,
      age_breaks = age_limits,
      per_capita_household_size = per_capita_household_size
    ),
    error = function(e) {
      stop("conmat::extrapolate_polymod() failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  source_matrix <- extract_conmat_source_matrix(conmat_result, setting, source_age_structure)

  ContactMatrixSource(
    matrix = source_matrix,
    age_structure = source_age_structure,
    source = "conmat",
    country = NULL,
    setting = setting,
    orientation = "recipient_source",
    source_reference = paste(
      "conmat::extrapolate_polymod(), based on the conmat POLYMOD-fitted",
      "setting models and caller-supplied population data."
    ),
    notes = paste(
      "conmat source generation requires caller-supplied population data;",
      "agepi only validates and stores the generated matrix."
    ),
    limitations = paste(
      "conmat outputs depend on the supplied population, model defaults,",
      "and optional household-size adjustment. agepi does not refit models",
      "or apply additional reciprocity correction."
    ),
    metadata = list(
      source_label = paste("conmat POLYMOD-extrapolated", setting, "contact matrix"),
      per_capita_household_size = per_capita_household_size
    )
  )
}

ContactMatrixSource <- function(matrix,
                                age_structure,
                                source,
                                country = NULL,
                                setting = NULL,
                                orientation = "recipient_source",
                                source_reference,
                                notes = NULL,
                                limitations = NULL,
                                metadata = list()) {
  validate_age_structure(age_structure)
  matrix <- validate_contact_matrix(matrix, age_structure)
  source <- validate_contact_source_scalar(source, "source")
  orientation <- validate_contact_source_scalar(orientation, "orientation")
  source_reference <- validate_contact_source_scalar(source_reference, "source_reference")

  if (!is.null(country)) {
    country <- validate_contact_source_scalar(country, "country")
  }
  if (!is.null(setting)) {
    setting <- validate_contact_source_scalar(setting, "setting")
  }
  if (is.null(notes) && is.null(limitations)) {
    notes <- "No source notes supplied."
  }
  if (!is.null(notes)) {
    notes <- as.character(notes)
  }
  if (!is.null(limitations)) {
    limitations <- as.character(limitations)
  }
  if (!is.list(metadata) || is.data.frame(metadata)) {
    stop("metadata must be a list.", call. = FALSE)
  }

  result <- list(
    matrix = matrix,
    age_structure = age_structure,
    source = source,
    country = country,
    setting = setting,
    orientation = orientation,
    convention = "rows are recipients/contact-makers; columns are sources/contacts",
    source_reference = source_reference,
    notes = notes,
    limitations = limitations,
    metadata = metadata
  )
  class(result) <- c("agepi_contact_matrix_source", "list")
  validate_contact_matrix_source(result)
}

#' Validate a contact matrix source object
#'
#' Checks that a contact source object has the common agepi source-layer
#' structure: `matrix`, `age_structure`, `source`, optional `country`, optional
#' `setting`, orientation/convention metadata, `source_reference`, notes or
#' limitations, source-specific `metadata`, and class
#' `agepi_contact_matrix_source`.
#'
#' @param source_contact_matrix Object to validate.
#'
#' @return Invisibly returns `source_contact_matrix` when valid.
#' @export
validate_contact_matrix_source <- function(source_contact_matrix) {
  if (!inherits(source_contact_matrix, "agepi_contact_matrix_source")) {
    stop("source_contact_matrix must be an agepi_contact_matrix_source object.", call. = FALSE)
  }

  required_fields <- c(
    "matrix",
    "age_structure",
    "source",
    "country",
    "setting",
    "orientation",
    "convention",
    "source_reference",
    "notes",
    "limitations",
    "metadata"
  )
  missing_fields <- setdiff(required_fields, names(source_contact_matrix))
  if (length(missing_fields) > 0) {
    stop(
      "source_contact_matrix is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  validate_age_structure(source_contact_matrix$age_structure)
  validate_contact_matrix(source_contact_matrix$matrix, source_contact_matrix$age_structure)
  validate_contact_source_scalar(source_contact_matrix$source, "source")
  validate_contact_source_scalar(source_contact_matrix$orientation, "orientation")
  if (!source_contact_matrix$orientation %in% c("recipient_source", "source_recipient")) {
    stop("source_contact_matrix$orientation must be 'recipient_source' or 'source_recipient'.", call. = FALSE)
  }
  validate_contact_source_scalar(source_contact_matrix$convention, "convention")
  validate_contact_source_scalar(source_contact_matrix$source_reference, "source_reference")
  if (!is.null(source_contact_matrix$country)) {
    validate_contact_source_scalar(source_contact_matrix$country, "country")
  }
  if (!is.null(source_contact_matrix$setting)) {
    validate_contact_source_scalar(source_contact_matrix$setting, "setting")
  }
  if (!is.null(source_contact_matrix$notes) && !is.character(source_contact_matrix$notes)) {
    stop("source_contact_matrix$notes must be NULL or character.", call. = FALSE)
  }
  if (!is.null(source_contact_matrix$limitations) && !is.character(source_contact_matrix$limitations)) {
    stop("source_contact_matrix$limitations must be NULL or character.", call. = FALSE)
  }
  if (is.null(source_contact_matrix$notes) && is.null(source_contact_matrix$limitations)) {
    stop("source_contact_matrix must include notes or limitations.", call. = FALSE)
  }
  if (!is.list(source_contact_matrix$metadata) || is.data.frame(source_contact_matrix$metadata)) {
    stop("source_contact_matrix$metadata must be a list.", call. = FALSE)
  }
  invisible(source_contact_matrix)
}

contact_source_metadata <- function(source_contact_matrix, adaptation_metadata = list()) {
  metadata <- c(
    list(
      source = source_contact_matrix$source,
      country = source_contact_matrix$country,
      setting = source_contact_matrix$setting,
      orientation = source_contact_matrix$orientation,
      convention = source_contact_matrix$convention,
      source_reference = source_contact_matrix$source_reference,
      source_age_grid = paste(source_contact_matrix$age_structure$age_groups, collapse = ", "),
      notes = paste(source_contact_matrix$notes, collapse = " "),
      limitations = paste(source_contact_matrix$limitations, collapse = " ")
    ),
    source_contact_matrix$metadata,
    adaptation_metadata
  )

  metadata[!vapply(metadata, is.null, logical(1))]
}

contact_source_age_structure_from_limits <- function(age_limits) {
  if (!is.numeric(age_limits) || length(age_limits) < 2 || anyNA(age_limits)) {
    stop("age_limits must be a numeric vector of at least two lower age limits.", call. = FALSE)
  }
  if (any(!is.finite(age_limits[-length(age_limits)]))) {
    stop("Only the final age limit can be Inf.", call. = FALSE)
  }
  if (is.unsorted(age_limits, strictly = TRUE) || any(duplicated(age_limits))) {
    stop("age_limits must be strictly increasing.", call. = FALSE)
  }

  if (is.infinite(age_limits[length(age_limits)])) {
    lower_bounds <- age_limits[-length(age_limits)]
    upper_bounds <- c(age_limits[-c(1, length(age_limits))] - 1, Inf)
  } else {
    lower_bounds <- age_limits
    upper_bounds <- c(age_limits[-1] - 1, Inf)
  }

  age_groups <- ifelse(
    is.infinite(upper_bounds),
    paste0(lower_bounds, "+"),
    ifelse(lower_bounds == upper_bounds, as.character(lower_bounds), paste0(lower_bounds, "-", upper_bounds))
  )

  AgeStructure(age_groups, lower_bounds, upper_bounds)
}

validate_contact_source_scalar <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop(name, " must be a single non-empty string.", call. = FALSE)
  }
  x
}

extract_conmat_source_matrix <- function(conmat_result, setting, age_structure) {
  matrix_candidate <- conmat_result
  if (is.list(conmat_result) && !is.data.frame(conmat_result) && !is.matrix(conmat_result)) {
    if (setting %in% names(conmat_result)) {
      matrix_candidate <- conmat_result[[setting]]
    } else if (setting == "all" && "all" %in% names(conmat_result)) {
      matrix_candidate <- conmat_result$all
    } else if (setting == "all" && "contact_matrix" %in% names(conmat_result)) {
      matrix_candidate <- conmat_result$contact_matrix
    } else {
      stop(
        "conmat result does not contain requested setting '",
        setting,
        "'. Available component(s): ",
        paste(names(conmat_result), collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (is.data.frame(matrix_candidate)) {
    return(contact_matrix_from_conmat(matrix_candidate, age_structure = age_structure))
  }

  as_agepi_contact_matrix(
    contact_source_matrix_on_age_structure(matrix_candidate, age_structure),
    age_structure = age_structure
  )
}

contact_source_matrix_on_age_structure <- function(contact_matrix, age_structure) {
  if (!is.matrix(contact_matrix) || !is.numeric(contact_matrix)) {
    stop("source contact matrix must be a numeric matrix.", call. = FALSE)
  }
  if (!identical(dim(contact_matrix), c(age_structure$n_age_groups, age_structure$n_age_groups))) {
    stop("source contact matrix dimensions must match source age_structure.", call. = FALSE)
  }
  contact_matrix <- matrix(
    as.numeric(contact_matrix),
    nrow = age_structure$n_age_groups,
    ncol = age_structure$n_age_groups,
    dimnames = list(age_structure$age_groups, age_structure$age_groups)
  )
  validate_contact_matrix(contact_matrix, age_structure)
}
