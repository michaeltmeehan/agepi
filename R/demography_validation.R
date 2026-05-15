#' Compare simulated demography to observed population trajectories
#'
#' Diagnostic helper for comparing output from [simulate_demography()] against
#' externally supplied population trajectories stored as an [Demography()] object,
#' including WPP-style objects created by [demography_from_wpp()].
#'
#' The comparison is diagnostic only. It does not change simulation dynamics, add
#' residual forcing, interpolate between time points, transform age groups, or
#' imply that [simulate_demography()] reproduces WPP projections. Only exact
#' shared time points and exact matching age groups are compared.
#'
#' If `simulated` contains extra time points, they are ignored. If `observed`
#' contains extra time points, they are ignored. If there are no shared exact
#' time points, the function errors. Age structures must match exactly for this
#' first-pass diagnostic.
#'
#' @param simulated Data frame returned by [simulate_demography()] with columns
#'   `time`, `age_group`, and `population`.
#' @param observed An `agepi_demography` object from [Demography()] or
#'   [demography_from_wpp()].
#'
#' @return A data frame with one row per shared exact time and age group, ordered
#'   by time and then by the age-group order in `observed`. Columns are `time`,
#'   `age_group`, `simulated_population`, `observed_population`,
#'   `absolute_error`, `relative_error`, `simulated_total_population`,
#'   `observed_total_population`, `total_absolute_error`,
#'   `total_relative_error`, `simulated_age_share`, `observed_age_share`, and
#'   `age_share_error`.
#'
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("0-4", "5+"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, Inf)
#' )
#'
#' observed <- Demography(
#'   data.frame(
#'     time = c(2020, 2020, 2021, 2021),
#'     age_group = c("0-4", "5+", "0-4", "5+"),
#'     population = c(100, 200, 95, 210)
#'   ),
#'   ages
#' )
#'
#' simulated <- data.frame(
#'   time = c(2020, 2020, 2021, 2021),
#'   age_group = c("0-4", "5+", "0-4", "5+"),
#'   population = c(100, 200, 97, 208)
#' )
#'
#' compare_demography_to_observed(simulated, observed)
#' @export
compare_demography_to_observed <- function(simulated, observed) {
  validate_simulated_demography_output(simulated)
  validate_agepi_demography(observed)

  simulated$age_group <- as.character(simulated$age_group)
  simulated_age_groups <- unique(simulated$age_group)
  observed_age_groups <- observed$age_groups

  if (
    length(simulated_age_groups) != length(observed_age_groups) ||
      !setequal(simulated_age_groups, observed_age_groups)
  ) {
    stop(
      "simulated and observed age groups must match exactly; ",
      "no age transformation is performed.",
      call. = FALSE
    )
  }

  shared_times <- sort(intersect(unique(simulated$time), observed$times))
  if (length(shared_times) == 0) {
    stop(
      "simulated and observed have no shared exact time points; no interpolation is performed.",
      call. = FALSE
    )
  }

  simulated <- simulated[simulated$time %in% shared_times, c("time", "age_group", "population")]
  for (this_time in shared_times) {
    age_groups_at_time <- simulated$age_group[simulated$time == this_time]
    if (
      length(age_groups_at_time) != length(observed_age_groups) ||
        !setequal(age_groups_at_time, observed_age_groups)
    ) {
      stop(
        "simulated and observed age groups must match exactly at each shared time; ",
        "no age transformation is performed.",
        call. = FALSE
      )
    }
  }

  observed_table <- observed$demography[
    observed$demography$time %in% shared_times,
    c("time", "age_group", "population")
  ]

  row_order <- order(
    simulated$time,
    match(simulated$age_group, observed_age_groups)
  )
  simulated <- simulated[row_order, ]
  row.names(simulated) <- NULL

  row_order <- order(
    observed_table$time,
    match(observed_table$age_group, observed_age_groups)
  )
  observed_table <- observed_table[row_order, ]
  row.names(observed_table) <- NULL

  comparison <- data.frame(
    time = observed_table$time,
    age_group = observed_table$age_group,
    simulated_population = simulated$population,
    observed_population = observed_table$population,
    stringsAsFactors = FALSE
  )

  comparison$absolute_error <- comparison$simulated_population - comparison$observed_population
  comparison$relative_error <- ifelse(
    comparison$observed_population > 0,
    comparison$absolute_error / comparison$observed_population,
    NA_real_
  )

  totals <- stats::aggregate(
    cbind(simulated_population, observed_population) ~ time,
    data = comparison,
    FUN = sum
  )
  names(totals) <- c(
    "time",
    "simulated_total_population",
    "observed_total_population"
  )
  totals$total_absolute_error <-
    totals$simulated_total_population - totals$observed_total_population
  totals$total_relative_error <- ifelse(
    totals$observed_total_population > 0,
    totals$total_absolute_error / totals$observed_total_population,
    NA_real_
  )

  comparison <- merge(comparison, totals, by = "time", sort = FALSE)
  comparison$simulated_age_share <- ifelse(
    comparison$simulated_total_population > 0,
    comparison$simulated_population / comparison$simulated_total_population,
    NA_real_
  )
  comparison$observed_age_share <- ifelse(
    comparison$observed_total_population > 0,
    comparison$observed_population / comparison$observed_total_population,
    NA_real_
  )
  comparison$age_share_error <- comparison$simulated_age_share - comparison$observed_age_share

  comparison <- comparison[
    order(comparison$time, match(comparison$age_group, observed_age_groups)),
    c(
      "time",
      "age_group",
      "simulated_population",
      "observed_population",
      "absolute_error",
      "relative_error",
      "simulated_total_population",
      "observed_total_population",
      "total_absolute_error",
      "total_relative_error",
      "simulated_age_share",
      "observed_age_share",
      "age_share_error"
    )
  ]
  row.names(comparison) <- NULL

  comparison
}

#' Summarise a demographic comparison by time
#'
#' Summarises output from [compare_demography_to_observed()] to one row per time
#' point. This is a diagnostic summary only; it does not interpolate, transform
#' age groups, or modify simulation dynamics.
#'
#' @param comparison Data frame returned by [compare_demography_to_observed()].
#'
#' @return A data frame with one row per time and columns `time`,
#'   `simulated_total_population`, `observed_total_population`,
#'   `total_absolute_error`, `total_relative_error`,
#'   `mean_absolute_age_error`, `max_absolute_age_error`,
#'   `mean_absolute_age_share_error`, and `max_absolute_age_share_error`.
#'
#' @examples
#' comparison <- data.frame(
#'   time = c(2020, 2020),
#'   age_group = c("0-4", "5+"),
#'   simulated_population = c(105, 195),
#'   observed_population = c(100, 200),
#'   absolute_error = c(5, -5),
#'   relative_error = c(0.05, -0.025),
#'   simulated_total_population = c(300, 300),
#'   observed_total_population = c(300, 300),
#'   total_absolute_error = c(0, 0),
#'   total_relative_error = c(0, 0),
#'   simulated_age_share = c(0.35, 0.65),
#'   observed_age_share = c(1 / 3, 2 / 3),
#'   age_share_error = c(0.35 - 1 / 3, 0.65 - 2 / 3)
#' )
#' summarise_demography_comparison(comparison)
#' @export
summarise_demography_comparison <- function(comparison) {
  validate_demography_comparison(comparison)

  times <- sort(unique(comparison$time))
  summaries <- vector("list", length(times))

  for (i in seq_along(times)) {
    rows <- comparison$time == times[i]
    time_comparison <- comparison[rows, ]

    summaries[[i]] <- data.frame(
      time = times[i],
      simulated_total_population = time_comparison$simulated_total_population[1],
      observed_total_population = time_comparison$observed_total_population[1],
      total_absolute_error = time_comparison$total_absolute_error[1],
      total_relative_error = time_comparison$total_relative_error[1],
      mean_absolute_age_error = mean(abs(time_comparison$absolute_error)),
      max_absolute_age_error = max(abs(time_comparison$absolute_error)),
      mean_absolute_age_share_error = mean(abs(time_comparison$age_share_error), na.rm = TRUE),
      max_absolute_age_share_error = max(abs(time_comparison$age_share_error), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, summaries)
  row.names(summary) <- NULL
  summary
}

validate_simulated_demography_output <- function(simulated) {
  if (!is.data.frame(simulated)) {
    stop("simulated must be a data frame.", call. = FALSE)
  }

  required_columns <- c("time", "age_group", "population")
  missing_columns <- setdiff(required_columns, names(simulated))
  if (length(missing_columns) > 0) {
    stop(
      "simulated is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(simulated$time)) {
    stop("simulated time must be numeric.", call. = FALSE)
  }
  if (anyNA(simulated$time) || any(!is.finite(simulated$time))) {
    stop("simulated time must be finite and non-missing.", call. = FALSE)
  }
  if (!is.numeric(simulated$population)) {
    stop("simulated population must be numeric.", call. = FALSE)
  }
  if (anyNA(simulated$population) || any(!is.finite(simulated$population))) {
    stop("simulated population must be finite and non-missing.", call. = FALSE)
  }
  if (any(simulated$population < 0)) {
    stop("simulated population cannot be negative.", call. = FALSE)
  }
  if (anyNA(simulated$age_group)) {
    stop("simulated age_group cannot contain missing values.", call. = FALSE)
  }

  duplicate_rows <- duplicated(data.frame(
    time = simulated$time,
    age_group = as.character(simulated$age_group)
  ))
  if (any(duplicate_rows)) {
    stop("simulated contains duplicate time-age_group rows.", call. = FALSE)
  }

  invisible(simulated)
}

validate_demography_comparison <- function(comparison) {
  if (!is.data.frame(comparison)) {
    stop("comparison must be a data frame.", call. = FALSE)
  }

  required_columns <- c(
    "time",
    "age_group",
    "simulated_population",
    "observed_population",
    "absolute_error",
    "relative_error",
    "simulated_total_population",
    "observed_total_population",
    "total_absolute_error",
    "total_relative_error",
    "simulated_age_share",
    "observed_age_share",
    "age_share_error"
  )
  missing_columns <- setdiff(required_columns, names(comparison))
  if (length(missing_columns) > 0) {
    stop(
      "comparison is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(comparison)
}
