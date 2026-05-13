#' Summarise deterministic simulation output by compartment
#'
#' Totals tidy deterministic simulation output across age groups for each time
#' and compartment.
#'
#' @param simulation_output Data frame produced by [simulate_deterministic()],
#'   with columns `time`, `compartment`, `age_group`, and `value`.
#'
#' @return Data frame with columns `time`, `compartment`, and `value`.
#' @export
compartment_totals <- function(simulation_output) {
  validate_simulation_output_summary_input(simulation_output)
  simulation_output_sum_by(simulation_output, c("time", "compartment"))
}

#' Summarise deterministic simulation output by age group
#'
#' Totals tidy deterministic simulation output across compartments for each time
#' and age group.
#'
#' @param simulation_output Data frame produced by [simulate_deterministic()],
#'   with columns `time`, `compartment`, `age_group`, and `value`.
#'
#' @return Data frame with columns `time`, `age_group`, and `value`.
#' @export
age_group_totals <- function(simulation_output) {
  validate_simulation_output_summary_input(simulation_output)
  simulation_output_sum_by(simulation_output, c("time", "age_group"))
}

#' Summarise total population from deterministic simulation output
#'
#' Totals tidy deterministic simulation output across all compartments and age
#' groups for each time.
#'
#' @param simulation_output Data frame produced by [simulate_deterministic()],
#'   with columns `time`, `compartment`, `age_group`, and `value`.
#'
#' @return Data frame with columns `time` and `value`.
#' @export
total_population <- function(simulation_output) {
  validate_simulation_output_summary_input(simulation_output)
  simulation_output_sum_by(simulation_output, "time")
}

validate_simulation_output_summary_input <- function(simulation_output) {
  if (!is.data.frame(simulation_output)) {
    stop("simulation_output must be a data frame.", call. = FALSE)
  }

  required_columns <- c("time", "compartment", "age_group", "value")
  missing_columns <- setdiff(required_columns, names(simulation_output))
  if (length(missing_columns) > 0) {
    stop(
      "simulation_output is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(simulation_output$value)) {
    stop("simulation_output value must be numeric.", call. = FALSE)
  }

  if (anyNA(simulation_output$value)) {
    stop("simulation_output value cannot contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(simulation_output$value))) {
    stop("simulation_output value cannot contain non-finite values.", call. = FALSE)
  }

  if (any(simulation_output$value < 0)) {
    stop("simulation_output value cannot be negative.", call. = FALSE)
  }

  invisible(simulation_output)
}

simulation_output_sum_by <- function(simulation_output, by) {
  grouping_key <- do.call(
    paste,
    c(simulation_output[by], list(sep = "\r"))
  )
  group_levels <- unique(grouping_key)
  group_index <- match(grouping_key, group_levels)

  first_rows <- match(group_levels, grouping_key)
  totals <- simulation_output[first_rows, by, drop = FALSE]
  totals$value <- as.numeric(tapply(
    simulation_output$value,
    factor(group_index, levels = seq_along(group_levels)),
    sum
  ))

  row.names(totals) <- NULL
  totals
}
