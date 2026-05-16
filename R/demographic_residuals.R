#' Compute implied demographic residual flows
#'
#' Diagnostic helper for one-step Euler accounting against an observed
#' demographic trajectory. For each observed interval, the function starts from
#' the observed population at `time_start`, evaluates [demographic_derivative()]
#' at that exact time using `process`, predicts the next population with one
#' Euler step, and reports the residual needed to match the observed population
#' at `time_end`.
#'
#' The residual is diagnostic only. It is not automatically applied to a
#' process, converted to migration, used as residual forcing, or fed back into
#' simulation dynamics. No interpolation or age transformation is performed.
#' Residuals may be interpreted as implied unmodelled net flows, such as omitted
#' migration, data mismatch, or rate approximation error. If `process` already
#' includes migration, the residual is the remaining gap relative to that full
#' process.
#'
#' For each age group and interval:
#'
#' `predicted_end = observed_start + dt * model_derivative`
#'
#' `residual_count = observed_end - predicted_end`
#'
#' `residual_rate = residual_count / (dt * observed_start)`
#'
#' `residual_rate` is `NA_real_` where `observed_start` is zero.
#'
#' @param observed An `agepi_demography` object from [Demography()] or
#'   [demography_from_wpp()].
#' @param process Valid `agepi_demographic_process` object.
#' @param type Which residual columns to include: `"count"`, `"rate"`, or
#'   `"both"`.
#'
#' @return A data frame with one row per observed interval and age group,
#'   ordered by `time_start` and the process age-group order. Common columns are
#'   `time_start`, `time_end`, `dt`, `age_group`, `observed_start`,
#'   `observed_end`, `predicted_end`, and `model_derivative`. Residual columns
#'   are controlled by `type`.
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("0-4", "5+"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, Inf)
#' )
#' process <- build_demographic_process(ages)
#' observed <- Demography(
#'   data.frame(
#'     time = c(0, 0, 1, 1),
#'     age_group = c("0-4", "5+", "0-4", "5+"),
#'     population = c(100, 200, 85, 215)
#'   ),
#'   ages
#' )
#' implied_demographic_residual(observed, process)
#' @export
implied_demographic_residual <- function(observed,
                                         process,
                                         type = c("both", "count", "rate")) {
  type <- match.arg(type)
  validate_implied_demographic_residual_inputs(observed, process)

  times <- observed$times
  age_groups <- process$age_structure$age_groups
  n_intervals <- length(times) - 1
  residuals <- vector("list", n_intervals)

  for (i in seq_len(n_intervals)) {
    time_start <- times[i]
    time_end <- times[i + 1]
    dt <- time_end - time_start

    observed_start <- demography_population_vector(observed, time_start)
    observed_end <- demography_population_vector(observed, time_end)
    model_derivative <- demographic_derivative(
      state = observed_start,
      time = time_start,
      process = process
    )
    predicted_end <- as.numeric(observed_start) + dt * as.numeric(model_derivative)
    residual_count <- as.numeric(observed_end) - predicted_end
    residual_rate <- ifelse(
      as.numeric(observed_start) > 0,
      residual_count / (dt * as.numeric(observed_start)),
      NA_real_
    )

    residuals[[i]] <- data.frame(
      time_start = rep(time_start, length(age_groups)),
      time_end = rep(time_end, length(age_groups)),
      dt = rep(dt, length(age_groups)),
      age_group = age_groups,
      observed_start = as.numeric(observed_start),
      observed_end = as.numeric(observed_end),
      predicted_end = predicted_end,
      residual_count = residual_count,
      residual_rate = residual_rate,
      model_derivative = as.numeric(model_derivative),
      stringsAsFactors = FALSE
    )
  }

  residual <- do.call(rbind, residuals)
  row.names(residual) <- NULL

  keep_columns <- c(
    "time_start",
    "time_end",
    "dt",
    "age_group",
    "observed_start",
    "observed_end",
    "predicted_end"
  )
  if (type %in% c("count", "both")) {
    keep_columns <- c(keep_columns, "residual_count")
  }
  if (type %in% c("rate", "both")) {
    keep_columns <- c(keep_columns, "residual_rate")
  }
  keep_columns <- c(keep_columns, "model_derivative")

  residual[, keep_columns, drop = FALSE]
}

validate_implied_demographic_residual_inputs <- function(observed, process) {
  validate_agepi_demography(observed)
  validate_demographic_process(process)

  if (length(observed$times) < 2) {
    stop("observed must contain at least two exact time points.", call. = FALSE)
  }

  validate_same_age_structure(
    process$age_structure,
    observed$age_structure,
    "observed"
  )

  invisible(NULL)
}
