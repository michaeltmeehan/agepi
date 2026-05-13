#' Validate an age-specific demography table
#'
#' Checks a tabular age-specific population trajectory against an
#' age structure. The table must contain exactly one population value
#' for every age group at every time point.
#'
#' @param demography Data frame with `time`, `age_group`, and `population`
#'   columns.
#' @param age_structure Age-structure object validated by
#'   [validate_age_structure()].
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_demography_table <- function(demography, age_structure) {
  validate_age_structure(age_structure)

  if (!is.data.frame(demography)) {
    stop("demography must be a data frame.", call. = FALSE)
  }

  required_columns <- c("time", "age_group", "population")
  missing_columns <- setdiff(required_columns, names(demography))
  if (length(missing_columns) > 0) {
    stop(
      "demography is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  time <- demography$time
  age_group <- demography$age_group
  population <- demography$population

  if (!is.numeric(time)) {
    stop("demography time must be numeric.", call. = FALSE)
  }

  if (anyNA(time) || any(!is.finite(time))) {
    stop("demography time must be finite and non-missing.", call. = FALSE)
  }

  if (!is.numeric(population)) {
    stop("demography population must be numeric.", call. = FALSE)
  }

  if (anyNA(population) || any(!is.finite(population))) {
    stop("demography population must be finite and non-missing.", call. = FALSE)
  }

  if (any(population < 0)) {
    stop("demography population cannot be negative.", call. = FALSE)
  }

  if (anyNA(age_group)) {
    stop("demography age_group cannot contain missing values.", call. = FALSE)
  }

  age_group <- as.character(age_group)
  expected_age_groups <- age_structure$age_groups
  extra_age_groups <- setdiff(unique(age_group), expected_age_groups)
  if (length(extra_age_groups) > 0) {
    stop(
      "demography contains age_group value(s) not in age_structure: ",
      paste(extra_age_groups, collapse = ", "),
      call. = FALSE
    )
  }

  duplicate_rows <- duplicated(data.frame(time = time, age_group = age_group))
  if (any(duplicate_rows)) {
    stop("demography contains duplicate time-age_group rows.", call. = FALSE)
  }

  times <- unique(time)
  for (this_time in times) {
    observed_age_groups <- age_group[time == this_time]
    missing_age_groups <- setdiff(expected_age_groups, observed_age_groups)
    if (length(missing_age_groups) > 0) {
      stop(
        "demography is missing age_group value(s) at time ",
        this_time,
        ": ",
        paste(missing_age_groups, collapse = ", "),
        call. = FALSE
      )
    }

    if (length(observed_age_groups) != length(expected_age_groups)) {
      stop(
        "demography must contain exactly one row per age group at each time point.",
        call. = FALSE
      )
    }
  }

  invisible(demography)
}

#' Construct a demography object
#'
#' Creates a simple internal representation of age-specific population
#' trajectories from a validated tabular input. Rows are ordered by
#' increasing `time`, then by the order of `age_structure$age_groups`.
#'
#' @param demography Data frame with `time`, `age_group`, and `population`
#'   columns.
#' @param age_structure Age-structure object validated by
#'   [validate_age_structure()].
#'
#' @return An `agepi_demography` list.
#' @export
Demography <- function(demography, age_structure) {
  validate_demography_table(demography, age_structure)

  demography <- demography[, c("time", "age_group", "population")]
  demography$age_group <- as.character(demography$age_group)
  row_order <- order(demography$time, match(demography$age_group, age_structure$age_groups))
  demography <- demography[row_order, ]
  row.names(demography) <- NULL

  times <- sort(unique(demography$time))

  demography_object <- list(
    demography = demography,
    age_structure = age_structure,
    times = times,
    n_times = length(times),
    age_groups = age_structure$age_groups,
    n_age_groups = age_structure$n_age_groups
  )

  class(demography_object) <- "agepi_demography"
  demography_object
}

#' Get demography time points
#'
#' @param demography An `agepi_demography` object.
#'
#' @return Numeric vector of available time points.
#' @export
demography_times <- function(demography) {
  validate_agepi_demography(demography)
  demography$times
}

#' Get an age-specific population vector
#'
#' Returns population values for one exact time point. Values are ordered
#' and named by the demography object's age groups.
#'
#' @param demography An `agepi_demography` object.
#' @param time Exact time point to retrieve.
#'
#' @return Numeric vector of population values.
#' @export
demography_population_vector <- function(demography, time) {
  validate_agepi_demography(demography)
  validate_demography_time(time, demography$times)

  population_table <- demography_population_table(demography, time = time)
  population <- population_table$population
  names(population) <- demography$age_groups
  population
}

#' Get the demography table
#'
#' Returns the internal tidy demography table, optionally filtered to one
#' exact time point.
#'
#' @param demography An `agepi_demography` object.
#' @param time Optional exact time point to retrieve.
#'
#' @return Data frame with `time`, `age_group`, and `population` columns.
#' @export
demography_population_table <- function(demography, time = NULL) {
  validate_agepi_demography(demography)

  if (is.null(time)) {
    return(demography$demography)
  }

  validate_demography_time(time, demography$times)
  population_table <- demography$demography[demography$demography$time == time, ]
  row.names(population_table) <- NULL
  population_table
}

validate_agepi_demography <- function(demography) {
  if (!inherits(demography, "agepi_demography")) {
    stop("demography must be an agepi_demography object.", call. = FALSE)
  }

  required_fields <- c(
    "demography",
    "age_structure",
    "times",
    "n_times",
    "age_groups",
    "n_age_groups"
  )
  missing_fields <- setdiff(required_fields, names(demography))
  if (length(missing_fields) > 0) {
    stop(
      "demography is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(demography)
}

validate_demography_time <- function(time, available_times) {
  if (!is.numeric(time) || length(time) != 1 || anyNA(time) || !is.finite(time)) {
    stop("time must be a finite numeric scalar.", call. = FALSE)
  }

  if (!time %in% available_times) {
    stop(
      "time is not available in demography: ",
      time,
      ". Available time point(s): ",
      paste(available_times, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(time)
}
