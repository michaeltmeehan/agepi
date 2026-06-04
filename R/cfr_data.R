#' Convert cumulative flow output for CFR or severity analysis
#'
#' Converts selected cumulative flow outputs from [simulate_deterministic()] or
#' [simulate_stochastic()] into a simple time series with cases and deaths. This
#' helper does not import or require `cfr`; it only prepares an ordinary
#' `data.frame` that can be passed to downstream tools.
#'
#' @param simulation_output Output returned by [simulate_deterministic()] or
#'   [simulate_stochastic()] when `cumulative_flows` was supplied.
#' @param cases Flow selector for cases. A character scalar can match either
#'   `cumulative_name` or `transition_id`. A list can specify `cumulative_name`,
#'   `transition_id`, or both.
#' @param deaths Flow selector for deaths, using the same format as `cases`.
#' @param cumulative Logical. If `TRUE`, selected cumulative flow values are
#'   converted to interval increments. If `FALSE`, values are returned as-is.
#' @param by_age_group Logical. If `TRUE`, return one row per time and age group.
#'   If `FALSE`, aggregate selected flows across age groups.
#'
#' @return A data frame with columns `time`, `cases`, and `deaths`. When
#'   `by_age_group = TRUE`, an `age_group` column is also included.
#' @export
#'
#' @examples
#' ages <- AgeStructure(c("0-4", "5+"), c(0, 5), c(4, Inf))
#' initial_state <- data.frame(
#'   compartment = rep(c("S", "I", "R"), each = 2),
#'   age_group = rep(ages$age_groups, times = 3),
#'   value = c(99, 198, 1, 2, 0, 0)
#' )
#' output <- simulate_deterministic(
#'   initial_state = initial_state,
#'   times = 0:2,
#'   model = SIRModel(gamma = 0.2),
#'   age_structure = ages,
#'   contact_matrix = diag(2),
#'   beta = 0.05,
#'   method = "euler",
#'   cumulative_flows = list(
#'     cases = list(from = "S", to = "I"),
#'     deaths = list(from = "I", to = "R")
#'   )
#' )
#' as_cfr_data(output, cases = "cases", deaths = "deaths")
as_cfr_data <- function(
  simulation_output,
  cases,
  deaths,
  cumulative = TRUE,
  by_age_group = FALSE
) {
  cumulative_output <- extract_cumulative_output(simulation_output)
  validate_cfr_cumulative_output(cumulative_output)

  case_rows <- select_cfr_flow(cumulative_output, cases, "cases")
  death_rows <- select_cfr_flow(cumulative_output, deaths, "deaths")

  group_columns <- if (by_age_group) c("time", "age_group") else "time"
  case_series <- cfr_flow_series(case_rows, group_columns, "cases")
  death_series <- cfr_flow_series(death_rows, group_columns, "deaths")

  result <- merge(case_series, death_series, by = group_columns, all = TRUE, sort = FALSE)
  result <- result[do.call(order, result[group_columns]), , drop = FALSE]
  result$cases[is.na(result$cases)] <- 0
  result$deaths[is.na(result$deaths)] <- 0

  if (cumulative) {
    result$cases <- cfr_interval_increments(result, "cases", group_columns)
    result$deaths <- cfr_interval_increments(result, "deaths", group_columns)
  }

  row.names(result) <- NULL
  result
}

extract_cumulative_output <- function(simulation_output) {
  if (is.list(simulation_output) && !is.data.frame(simulation_output) &&
      !is.null(simulation_output$cumulative)) {
    return(simulation_output$cumulative)
  }

  if (is.data.frame(simulation_output)) {
    return(simulation_output)
  }

  stop(
    "simulation_output must be a cumulative-flow simulation result or a cumulative data frame.",
    call. = FALSE
  )
}

validate_cfr_cumulative_output <- function(cumulative_output) {
  required_columns <- c("time", "cumulative_name", "transition_id", "age_group", "value")
  missing_columns <- setdiff(required_columns, names(cumulative_output))
  if (length(missing_columns) > 0) {
    stop(
      "cumulative output is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(cumulative_output$value)) {
    stop("cumulative output value must be numeric.", call. = FALSE)
  }

  if (anyNA(cumulative_output$value) || any(!is.finite(cumulative_output$value))) {
    stop("cumulative output value must be finite and non-missing.", call. = FALSE)
  }

  invisible(cumulative_output)
}

select_cfr_flow <- function(cumulative_output, selector, role) {
  selector <- normalize_cfr_flow_selector(selector, role)

  if (!is.null(selector$value)) {
    matches <- cumulative_output$cumulative_name == selector$value |
      cumulative_output$transition_id == selector$value
  } else {
    matches <- rep(TRUE, nrow(cumulative_output))
    if (!is.null(selector$cumulative_name)) {
      matches <- matches & cumulative_output$cumulative_name == selector$cumulative_name
    }
    if (!is.null(selector$transition_id)) {
      matches <- matches & cumulative_output$transition_id == selector$transition_id
    }
  }

  selected <- cumulative_output[matches, , drop = FALSE]
  if (nrow(selected) == 0) {
    stop(role, " flow selector did not match any cumulative output rows.", call. = FALSE)
  }

  flow_keys <- unique(selected[, c("cumulative_name", "transition_id"), drop = FALSE])
  if (nrow(flow_keys) != 1) {
    stop(
      role,
      " flow selector is ambiguous; specify cumulative_name and/or transition_id.",
      call. = FALSE
    )
  }

  selected
}

normalize_cfr_flow_selector <- function(selector, role) {
  if (is.character(selector) && length(selector) == 1 && !is.na(selector) && selector != "") {
    return(list(value = selector, cumulative_name = NULL, transition_id = NULL))
  }

  if (!is.list(selector) || is.data.frame(selector)) {
    stop(role, " flow selector must be a character scalar or a list.", call. = FALSE)
  }

  allowed_names <- c("cumulative_name", "transition_id")
  unknown_names <- setdiff(names(selector), allowed_names)
  if (length(unknown_names) > 0) {
    stop(
      role,
      " flow selector has unknown field(s): ",
      paste(unknown_names, collapse = ", "),
      call. = FALSE
    )
  }

  normalized <- list(
    cumulative_name = selector$cumulative_name,
    transition_id = selector$transition_id
  )

  if (is.null(normalized$cumulative_name) && is.null(normalized$transition_id)) {
    stop(role, " flow selector must specify cumulative_name or transition_id.", call. = FALSE)
  }

  for (field in allowed_names) {
    value <- normalized[[field]]
    if (!is.null(value) && (!is.character(value) || length(value) != 1 || is.na(value) || value == "")) {
      stop(role, " flow selector ", field, " must be a non-empty character scalar.", call. = FALSE)
    }
  }

  normalized
}

cfr_flow_series <- function(flow_rows, group_columns, value_name) {
  grouping_key <- do.call(paste, c(flow_rows[group_columns], list(sep = "\r")))
  group_levels <- unique(grouping_key)
  group_index <- match(grouping_key, group_levels)
  first_rows <- match(group_levels, grouping_key)

  series <- flow_rows[first_rows, group_columns, drop = FALSE]
  series[[value_name]] <- as.numeric(tapply(
    flow_rows$value,
    factor(group_index, levels = seq_along(group_levels)),
    sum
  ))
  row.names(series) <- NULL
  series
}

cfr_interval_increments <- function(result, value_column, group_columns) {
  if (!"age_group" %in% group_columns) {
    return(c(result[[value_column]][1], diff(result[[value_column]])))
  }

  increments <- numeric(nrow(result))
  for (age_group in unique(result$age_group)) {
    rows <- which(result$age_group == age_group)
    values <- result[[value_column]][rows]
    increments[rows] <- c(values[1], diff(values))
  }

  increments
}
