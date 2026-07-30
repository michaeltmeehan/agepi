#' Convert long-format state data to a numeric state vector
#'
#' Converts one row per compartment-age cell into deterministic
#' compartment-major, age-group-minor order. For compartments `c("S", "I",
#' "R")` and age groups `c("0-4", "5-9")`, the vector order is `S_0-4`,
#' `S_5-9`, `I_0-4`, `I_5-9`, `R_0-4`, `R_5-9`.
#'
#' @param state_long Data frame with `age_group`, `compartment`, and `value`.
#' @param age_structure Valid age structure.
#' @param compartments Character vector defining compartment order.
#'
#' @return Numeric vector in compartment-major order.
#' @export
state_long_to_vector <- function(state_long, age_structure, compartments) {
  validate_age_structure(age_structure)
  validate_compartments(compartments)
  validate_state_long(state_long, age_structure, compartments)

  ordered_state <- merge(
    state_order(age_structure, compartments),
    state_long,
    by = c("compartment", "age_group"),
    sort = FALSE
  )
  ordered_state <- ordered_state[order(ordered_state$.order), ]

  values <- ordered_state$value
  names(values) <- state_vector_names(age_structure, compartments)
  values
}

#' Convert a numeric state vector to long-format state data
#'
#' Interprets `state_vector` using compartment-major, age-group-minor order.
#' Any names on `state_vector` are ignored; ordering is determined only by
#' `age_structure` and `compartments`.
#'
#' @param state_vector Numeric state vector.
#' @param age_structure Valid age structure.
#' @param compartments Character vector defining compartment order.
#'
#' @return Data frame with `age_group`, `compartment`, and `value`.
#' @export
state_vector_to_long <- function(state_vector, age_structure, compartments) {
  validate_age_structure(age_structure)
  validate_compartments(compartments)
  validate_state_vector(state_vector, age_structure, compartments)

  data.frame(
    age_group = rep(age_structure$age_groups, times = length(compartments)),
    compartment = rep(compartments, each = age_structure$n_age_groups),
    value = as.numeric(state_vector),
    stringsAsFactors = FALSE
  )
}

#' Initialise compartment counts from scalar or age-specific proportions
#'
#' Allocates an age-specific population vector into compartments by multiplying
#' supplied proportions by the population. Each supplied compartment proportion
#' may be a scalar, an unnamed full-length vector in model age order, or a
#' named full-length vector that is realigned to `age_structure$age_groups`.
#' Any valid non-residual compartment omitted from `proportions` defaults to
#' zero in every age group. The residual population is assigned to
#' `residual_compartment`.
#'
#' @param population Numeric age-specific population vector.
#' @param proportions Named list of numeric scalars or numeric vectors, one per
#'   non-residual compartment.
#' @param residual_compartment Compartment receiving population not allocated
#'   by `proportions`.
#' @param compartments Optional full compartment order. Defaults to
#'   `residual_compartment` followed by `names(proportions)`.
#' @param age_structure Optional age structure used to validate and name ages.
#'
#' @return Long-format state data with columns `compartment`, `age_group`, and
#'   `value`.
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("child", "adult"),
#'   lower_bounds = c(0, 18),
#'   upper_bounds = c(17, Inf)
#' )
#' population <- c(child = 1000, adult = 2000)
#'
#' initialise_compartments_from_proportions(
#'   population = population,
#'   proportions = list(I = 0.01),
#'   residual_compartment = "S",
#'   compartments = c("S", "I", "R"),
#'   age_structure = ages
#' )
#'
#' initialise_compartments_from_proportions(
#'   population = population,
#'   proportions = list(
#'     I = c(child = 0.001, adult = 0.01)
#'   ),
#'   residual_compartment = "S",
#'   compartments = c("S", "I", "R"),
#'   age_structure = ages
#' )
#' @export
initialise_compartments_from_proportions <- function(population,
                                                     proportions,
                                                     residual_compartment,
                                                     compartments = NULL,
                                                     age_structure = NULL) {
  if (!is.null(age_structure)) {
    validate_age_structure(age_structure)
  }

  age_groups <- validate_population_for_initialisation(population, age_structure)
  if (is.null(age_structure)) {
    age_structure <- list(
      age_groups = age_groups,
      n_age_groups = length(age_groups)
    )
  }

  validate_initialisation_proportions(proportions)

  if (!is.character(residual_compartment) || length(residual_compartment) != 1 ||
      is.na(residual_compartment) || !nzchar(residual_compartment)) {
    stop("residual_compartment must be a single non-empty string.", call. = FALSE)
  }

  if (is.null(compartments)) {
    compartments <- c(residual_compartment, names(proportions))
  }
  validate_compartments(compartments)

  if (!residual_compartment %in% compartments) {
    stop("residual_compartment must be included in compartments.", call. = FALSE)
  }

  validate_initialisation_proportion_names(
    proportions = proportions,
    compartments = compartments,
    residual_compartment = residual_compartment
  )

  n_age_groups <- age_structure$n_age_groups
  normalised_proportions <- vector("list", length(compartments))
  names(normalised_proportions) <- compartments
  non_residual_compartments <- compartments[compartments != residual_compartment]

  for (compartment in non_residual_compartments) {
    x <- if (compartment %in% names(proportions)) proportions[[compartment]] else 0
    normalised_proportions[[compartment]] <- normalise_initial_proportion(
      x = x,
      compartment = compartment,
      age_structure = age_structure
    )
  }

  if (length(non_residual_compartments) > 0) {
    specified_proportions <- do.call(cbind, normalised_proportions[non_residual_compartments])
    residual_proportion <- 1 - rowSums(specified_proportions)
  } else {
    residual_proportion <- rep(1, n_age_groups)
  }

  tolerance <- sqrt(.Machine$double.eps)
  if (any(residual_proportion < -tolerance)) {
    bad_age_group <- age_groups[which(residual_proportion < -tolerance)[1]]
    stop(
      "Specified proportions exceed 1 for age group `",
      bad_age_group,
      "`.",
      call. = FALSE
    )
  }
  residual_proportion[residual_proportion < 0] <- 0
  normalised_proportions[[residual_compartment]] <- residual_proportion

  population <- as.numeric(population)
  values_by_compartment <- lapply(
    normalised_proportions[compartments],
    function(proportion) population * proportion
  )

  state <- data.frame(
    compartment = rep(compartments, each = length(population)),
    age_group = rep(age_groups, times = length(compartments)),
    value = unlist(values_by_compartment[compartments], use.names = FALSE),
    stringsAsFactors = FALSE
  )

  if (!is.null(age_structure)) {
    validate_state_long(state, age_structure, compartments)
  }
  state
}

validate_compartments <- function(compartments) {
  if (!is.character(compartments)) {
    stop("compartments must be a character vector.", call. = FALSE)
  }

  if (length(compartments) == 0) {
    stop("compartments must contain at least one compartment.", call. = FALSE)
  }

  if (anyNA(compartments) || any(compartments == "")) {
    stop("compartments cannot contain missing or empty values.", call. = FALSE)
  }

  duplicated_compartments <- unique(compartments[duplicated(compartments)])
  if (length(duplicated_compartments) > 0) {
    stop(
      "compartments must be unique; duplicate compartment(s): ",
      paste(duplicated_compartments, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(compartments)
}

validate_state_long <- function(state_long, age_structure, compartments) {
  if (!is.data.frame(state_long)) {
    stop("state_long must be a data frame.", call. = FALSE)
  }

  required_columns <- c("age_group", "compartment", "value")
  missing_columns <- setdiff(required_columns, names(state_long))
  if (length(missing_columns) > 0) {
    stop(
      "state_long is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(state_long$age_group) || anyNA(state_long$compartment)) {
    stop("state_long age_group and compartment cannot contain missing values.", call. = FALSE)
  }

  if (!is.numeric(state_long$value) || anyNA(state_long$value)) {
    stop("state_long value must be numeric and cannot contain missing values.", call. = FALSE)
  }

  unknown_ages <- setdiff(unique(state_long$age_group), age_structure$age_groups)
  if (length(unknown_ages) > 0) {
    stop("state_long contains unknown age_group value(s): ", paste(unknown_ages, collapse = ", "), call. = FALSE)
  }

  unknown_compartments <- setdiff(unique(state_long$compartment), compartments)
  if (length(unknown_compartments) > 0) {
    stop(
      "state_long contains unknown compartment value(s): ",
      paste(unknown_compartments, collapse = ", "),
      call. = FALSE
    )
  }

  expected <- state_order(age_structure, compartments)
  observed <- state_long[, c("compartment", "age_group")]

  duplicate_rows <- duplicated(observed)
  if (any(duplicate_rows)) {
    duplicated_keys <- observed[duplicate_rows, , drop = FALSE]
    stop(
      "state_long contains duplicate compartment-age rows, including: ",
      duplicated_keys$compartment[1],
      "/",
      duplicated_keys$age_group[1],
      call. = FALSE
    )
  }

  missing_rows <- merge(
    expected[, c("compartment", "age_group")],
    observed,
    by = c("compartment", "age_group"),
    all.x = TRUE
  )
  missing_rows <- missing_rows[is.na(match(
    paste(missing_rows$compartment, missing_rows$age_group),
    paste(observed$compartment, observed$age_group)
  )), , drop = FALSE]

  if (nrow(missing_rows) > 0) {
    stop(
      "state_long is missing compartment-age row: ",
      missing_rows$compartment[1],
      "/",
      missing_rows$age_group[1],
      call. = FALSE
    )
  }

  invisible(state_long)
}

validate_state_vector <- function(state_vector, age_structure, compartments) {
  if (!is.numeric(state_vector)) {
    stop("state_vector must be numeric.", call. = FALSE)
  }

  expected_length <- age_structure$n_age_groups * length(compartments)
  if (length(state_vector) != expected_length) {
    stop(
      "state_vector length must equal number of compartments times number of age groups: ",
      expected_length,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(state_vector)) {
    stop("state_vector cannot contain missing values.", call. = FALSE)
  }

  invisible(state_vector)
}

# Deterministic compartment-major ordering:
# S_1, S_2, ..., S_A, I_1, I_2, ..., I_A, etc.,
# where compartments are supplied order and ages follow age_structure$age_groups.
state_order <- function(age_structure, compartments) {
  data.frame(
    compartment = rep(compartments, each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = length(compartments)),
    .order = seq_len(age_structure$n_age_groups * length(compartments)),
    stringsAsFactors = FALSE
  )
}

state_vector_names <- function(age_structure, compartments) {
  paste(
    rep(compartments, each = age_structure$n_age_groups),
    rep(age_structure$age_groups, times = length(compartments)),
    sep = "_"
  )
}

validate_population_for_initialisation <- function(population, age_structure) {
  if (!is.numeric(population)) {
    stop("population must be numeric.", call. = FALSE)
  }

  if (length(population) == 0 || anyNA(population) || any(!is.finite(population))) {
    stop("population must contain finite non-missing values.", call. = FALSE)
  }

  if (any(population < 0)) {
    stop("population cannot be negative.", call. = FALSE)
  }

  if (!is.null(age_structure)) {
    if (length(population) != age_structure$n_age_groups) {
      stop("population length must match age_structure$n_age_groups.", call. = FALSE)
    }
    return(age_structure$age_groups)
  }

  if (is.null(names(population)) || any(names(population) == "")) {
    stop("population must be named when age_structure is not supplied.", call. = FALSE)
  }

  names(population)
}

validate_initialisation_proportions <- function(proportions) {
  if (!is.list(proportions) || is.data.frame(proportions)) {
    stop("proportions must be a named list.", call. = FALSE)
  }

  proportion_names <- names(proportions)
  if (length(proportions) == 0 && is.null(proportion_names)) {
    stop("proportions must be a named list.", call. = FALSE)
  }

  if (length(proportions) > 0 && is.null(proportion_names)) {
    stop("proportions must be a named list.", call. = FALSE)
  }

  if (length(proportions) > 0 && (anyNA(proportion_names) || any(proportion_names == ""))) {
    stop("proportions must not contain missing or empty compartment names.", call. = FALSE)
  }

  if (length(proportions) > 0) {
    duplicated_compartments <- unique(proportion_names[duplicated(proportion_names)])
    if (length(duplicated_compartments) > 0) {
      stop(
        "Duplicate compartment name(s) in `proportions`: ",
        paste(duplicated_compartments, collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(proportions)
}

validate_initialisation_proportion_names <- function(proportions,
                                                     compartments,
                                                     residual_compartment) {
  if (length(proportions) == 0) {
    return(invisible(proportions))
  }

  proportion_names <- names(proportions)

  if (residual_compartment %in% proportion_names) {
    stop(
      "The residual compartment `",
      residual_compartment,
      "` must not be supplied in `proportions`.",
      call. = FALSE
    )
  }

  unknown_compartments <- setdiff(proportion_names, compartments)
  if (length(unknown_compartments) > 0) {
    stop(
      "Unknown compartment in `proportions`: ",
      paste(unknown_compartments, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(proportions)
}

normalise_initial_proportion <- function(x, compartment, age_structure) {
  age_groups <- age_structure$age_groups
  n_age_groups <- age_structure$n_age_groups

  if (is.matrix(x) || !is.numeric(x)) {
    stop(
      "Proportion for compartment `",
      compartment,
      "` must be numeric.",
      call. = FALSE
    )
  }

  if (length(x) != 1 && length(x) != n_age_groups) {
    stop(
      "Proportion for compartment `",
      compartment,
      "` must have length 1 or ",
      n_age_groups,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop(
      "Proportion for compartment `",
      compartment,
      "` must be finite and non-missing.",
      call. = FALSE
    )
  }

  values <- as.numeric(x)

  if (length(values) == 1) {
    values <- rep(values, n_age_groups)
  } else {
    value_names <- names(x)
    if (!is.null(value_names) && length(value_names) > 0) {
      if (length(value_names) != n_age_groups ||
          anyNA(value_names) ||
          any(value_names == "") ||
          anyDuplicated(value_names)) {
        stop(
          "Named proportion vector for compartment `",
          compartment,
          "` must contain each model age group exactly once.",
          call. = FALSE
        )
      }

      unknown_age_groups <- setdiff(value_names, age_groups)
      if (length(unknown_age_groups) > 0) {
        stop(
          "Named proportion vector for compartment `",
          compartment,
          "` must contain each model age group exactly once.",
          call. = FALSE
        )
      }

      values <- values[match(age_groups, value_names)]
    }
  }

  if (any(values < 0) || any(values > 1)) {
    stop(
      "Proportion for compartment `",
      compartment,
      "` must be between 0 and 1.",
      call. = FALSE
    )
  }

  stats::setNames(values, age_groups)
}
