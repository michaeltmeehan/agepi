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

#' Initialise compartment counts from age-specific proportions
#'
#' Allocates an age-specific population vector into compartments by multiplying
#' supplied age-specific proportions by the population. The residual population
#' is assigned to `residual_compartment`.
#'
#' @param population Numeric age-specific population vector.
#' @param proportions Named list of numeric vectors, one per non-residual
#'   compartment.
#' @param residual_compartment Compartment receiving population not allocated
#'   by `proportions`.
#' @param compartments Optional full compartment order. Defaults to
#'   `residual_compartment` followed by `names(proportions)`.
#' @param age_structure Optional age structure used to validate and name ages.
#'
#' @return Long-format state data with columns `compartment`, `age_group`, and
#'   `value`.
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
  validate_initialisation_proportions(proportions, length(population))

  if (!is.character(residual_compartment) || length(residual_compartment) != 1 ||
      is.na(residual_compartment) || !nzchar(residual_compartment)) {
    stop("residual_compartment must be a single non-empty string.", call. = FALSE)
  }

  if (is.null(compartments)) {
    compartments <- c(residual_compartment, names(proportions))
  }
  validate_compartments(compartments)

  unknown_proportion_names <- setdiff(names(proportions), compartments)
  if (length(unknown_proportion_names) > 0) {
    stop(
      "proportions contains compartment(s) not in compartments: ",
      paste(unknown_proportion_names, collapse = ", "),
      call. = FALSE
    )
  }

  if (!residual_compartment %in% compartments) {
    stop("residual_compartment must be included in compartments.", call. = FALSE)
  }

  allocated <- rep(0, length(population))
  names(allocated) <- age_groups
  values_by_compartment <- vector("list", length(compartments))
  names(values_by_compartment) <- compartments

  for (compartment in names(proportions)) {
    values <- as.numeric(population) * as.numeric(proportions[[compartment]])
    values_by_compartment[[compartment]] <- values
    allocated <- allocated + values
  }

  residual <- as.numeric(population) - allocated
  if (any(residual < -1e-8)) {
    stop("proportions allocate more than the population for at least one age group.", call. = FALSE)
  }
  residual[residual < 0] <- 0
  values_by_compartment[[residual_compartment]] <- residual

  for (compartment in compartments) {
    if (is.null(values_by_compartment[[compartment]])) {
      values_by_compartment[[compartment]] <- rep(0, length(population))
    }
  }

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

validate_initialisation_proportions <- function(proportions, expected_length) {
  if (!is.list(proportions) || is.data.frame(proportions)) {
    stop("proportions must be a named list.", call. = FALSE)
  }

  if (length(proportions) == 0) {
    stop("proportions must contain at least one compartment.", call. = FALSE)
  }

  if (is.null(names(proportions)) || any(names(proportions) == "")) {
    stop("proportions must be a named list.", call. = FALSE)
  }

  for (compartment in names(proportions)) {
    values <- proportions[[compartment]]
    if (!is.numeric(values) || length(values) != expected_length ||
        anyNA(values) || any(!is.finite(values))) {
      stop(
        "proportion vector for compartment ",
        compartment,
        " must be finite numeric and match population length.",
        call. = FALSE
      )
    }

    if (any(values < 0 | values > 1)) {
      stop("proportions must be between 0 and 1.", call. = FALSE)
    }
  }

  invisible(proportions)
}
