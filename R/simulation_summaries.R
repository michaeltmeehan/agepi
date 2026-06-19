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

#' Summarise cumulative-flow output across age groups and transitions
#'
#' @param cumulative_flows Data frame returned in the `cumulative` component of
#'   [simulate_deterministic()] or [simulate_stochastic()].
#'
#' @return Data frame with columns `time`, `cumulative_name`, and `value`.
#' @export
cumulative_flow_totals <- function(cumulative_flows) {
  validate_cumulative_flow_summary_input(cumulative_flows)
  simulation_output_sum_by(cumulative_flows, c("time", "cumulative_name"))
}

#' Reshape cumulative-flow totals to one row per time
#'
#' @inheritParams cumulative_flow_totals
#' @param prefix Prefix for cumulative-flow value columns.
#'
#' @return A wide data frame with one row per time.
#' @export
cumulative_flow_totals_wide <- function(cumulative_flows, prefix = "cumulative_") {
  if (!is.character(prefix) || length(prefix) != 1 || is.na(prefix)) {
    stop("prefix must be a single string.", call. = FALSE)
  }

  totals <- cumulative_flow_totals(cumulative_flows)
  wide <- stats::reshape(
    totals,
    idvar = "time",
    timevar = "cumulative_name",
    direction = "wide"
  )
  names(wide) <- sub("^value[.]", prefix, names(wide))
  row.names(wide) <- NULL
  wide
}

#' Calculate interval increments from cumulative-flow totals
#'
#' @inheritParams cumulative_flow_totals
#' @param times Optional numeric vector of cumulative times to retain before
#'   differencing.
#' @param value_col Name for the increment column.
#'
#' @return Data frame with cumulative totals and interval increment values.
#' @export
cumulative_flow_increments <- function(cumulative_flows,
                                       times = NULL,
                                       value_col = "increment_value") {
  if (!is.character(value_col) || length(value_col) != 1 || is.na(value_col) || !nzchar(value_col)) {
    stop("value_col must be a single non-empty string.", call. = FALSE)
  }

  totals <- cumulative_flow_totals(cumulative_flows)
  if (!is.null(times)) {
    if (!is.numeric(times) || anyNA(times) || any(!is.finite(times))) {
      stop("times must be a numeric vector of finite non-missing values.", call. = FALSE)
    }
    totals <- totals[totals$time %in% times, , drop = FALSE]
  }

  totals <- totals[order(totals$cumulative_name, totals$time), ]
  totals[[value_col]] <- stats::ave(
    totals$value,
    totals$cumulative_name,
    FUN = function(x) c(x[1], diff(x))
  )
  row.names(totals) <- NULL
  totals
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

validate_cumulative_flow_summary_input <- function(cumulative_flows) {
  if (!is.data.frame(cumulative_flows)) {
    stop("cumulative_flows must be a data frame.", call. = FALSE)
  }

  required_columns <- c("time", "cumulative_name", "age_group", "value")
  missing_columns <- setdiff(required_columns, names(cumulative_flows))
  if (length(missing_columns) > 0) {
    stop(
      "cumulative_flows is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(cumulative_flows$value)) {
    stop("cumulative_flows value must be numeric.", call. = FALSE)
  }

  if (anyNA(cumulative_flows$value)) {
    stop("cumulative_flows value cannot contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(cumulative_flows$value))) {
    stop("cumulative_flows value cannot contain non-finite values.", call. = FALSE)
  }

  if (any(cumulative_flows$value < 0)) {
    stop("cumulative_flows value cannot be negative.", call. = FALSE)
  }

  invisible(cumulative_flows)
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
