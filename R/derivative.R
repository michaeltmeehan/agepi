#' Convert transition rates to deterministic derivatives
#'
#' Converts long-form transition rates into one derivative per
#' compartment-age cell. Each transition rate is subtracted from its source
#' compartment and added to its destination compartment within the same age
#' group.
#'
#' Output rows use deterministic compartment-major ordering: compartments are
#' outermost, and age groups are innermost following `age_structure$age_groups`.
#' Exact duplicate transition rows with the same `from`, `to`, and `age_group`
#' are rejected as ambiguous.
#'
#' @param transition_rate_table Data frame with columns `from`, `to`,
#'   `age_group`, and `rate`.
#' @param compartments Character vector defining compartment order.
#' @param age_structure Valid age structure.
#'
#' @return Data frame with columns `compartment`, `age_group`, and
#'   `derivative`.
#' @export
rates_to_derivative <- function(
  transition_rate_table,
  compartments,
  age_structure
) {
  validate_age_structure(age_structure)
  validate_compartments(compartments)
  validate_transition_rate_table(
    transition_rate_table = transition_rate_table,
    compartments = compartments,
    age_structure = age_structure
  )

  derivative <- state_order(age_structure, compartments)
  derivative$derivative <- 0

  for (i in seq_len(nrow(transition_rate_table))) {
    from_index <- match_derivative_cell(
      derivative,
      transition_rate_table$from[i],
      transition_rate_table$age_group[i]
    )
    to_index <- match_derivative_cell(
      derivative,
      transition_rate_table$to[i],
      transition_rate_table$age_group[i]
    )

    derivative$derivative[from_index] <- derivative$derivative[from_index] -
      transition_rate_table$rate[i]
    derivative$derivative[to_index] <- derivative$derivative[to_index] +
      transition_rate_table$rate[i]
  }

  derivative[, c("compartment", "age_group", "derivative")]
}

validate_transition_rate_table <- function(
  transition_rate_table,
  compartments,
  age_structure
) {
  if (!is.data.frame(transition_rate_table)) {
    stop("transition_rate_table must be a data frame.", call. = FALSE)
  }

  required_columns <- c("from", "to", "age_group", "rate")
  missing_columns <- setdiff(required_columns, names(transition_rate_table))
  if (length(missing_columns) > 0) {
    stop(
      "transition_rate_table is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(transition_rate_table$from) ||
      anyNA(transition_rate_table$to) ||
      anyNA(transition_rate_table$age_group)) {
    stop("transition_rate_table from, to, and age_group cannot contain missing values.", call. = FALSE)
  }

  if (!is.numeric(transition_rate_table$rate)) {
    stop("transition_rate_table rate must be numeric.", call. = FALSE)
  }

  if (anyNA(transition_rate_table$rate)) {
    stop("transition_rate_table rate cannot contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(transition_rate_table$rate))) {
    stop("transition_rate_table rate cannot contain non-finite values.", call. = FALSE)
  }

  if (any(transition_rate_table$rate < 0)) {
    stop("transition_rate_table rate cannot contain negative values.", call. = FALSE)
  }

  unknown_from <- setdiff(unique(transition_rate_table$from), compartments)
  if (length(unknown_from) > 0) {
    stop(
      "transition_rate_table contains unknown source compartment value(s): ",
      paste(unknown_from, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_to <- setdiff(unique(transition_rate_table$to), compartments)
  if (length(unknown_to) > 0) {
    stop(
      "transition_rate_table contains unknown destination compartment value(s): ",
      paste(unknown_to, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_ages <- setdiff(unique(transition_rate_table$age_group), age_structure$age_groups)
  if (length(unknown_ages) > 0) {
    stop(
      "transition_rate_table contains unknown age_group value(s): ",
      paste(unknown_ages, collapse = ", "),
      call. = FALSE
    )
  }

  transition_keys <- transition_rate_table[, c("from", "to", "age_group")]
  duplicate_rows <- duplicated(transition_keys)
  if (any(duplicate_rows)) {
    duplicated_key <- transition_keys[which(duplicate_rows)[1], , drop = FALSE]
    stop(
      "transition_rate_table contains duplicate transition row: ",
      duplicated_key$from,
      "->",
      duplicated_key$to,
      "/",
      duplicated_key$age_group,
      call. = FALSE
    )
  }

  invisible(transition_rate_table)
}

match_derivative_cell <- function(derivative, compartment, age_group) {
  match(
    paste(compartment, age_group),
    paste(derivative$compartment, derivative$age_group)
  )
}
