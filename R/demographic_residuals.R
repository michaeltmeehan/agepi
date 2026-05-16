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

#' Convert implied residual diagnostics to a migration schedule
#'
#' Converts output from [implied_demographic_residual()] into a
#' [MigrationSchedule()] using an explicit Euler interval convention. This is a
#' conversion helper only: it does not automatically force simulations, modify a
#' [DemographicProcess()], interpolate, or apply residuals to any model.
#'
#' `residual_count` is interpreted as the total interval gap over
#' `[time_start, time_end]`. For `use = "count"`, this interval total is
#' converted to the per-time additive flow used by [demographic_derivative()]:
#' `migration_count = residual_count / dt`. For `use = "rate"`,
#' `residual_rate` is used directly as `migration_rate`. Schedule times are the
#' interval start times, `time_start`, matching agepi's exact-time Euler
#' schedule lookup semantics.
#'
#' Missing residual rates are not replaced. If `use = "rate"` and any
#' `residual_rate` value is `NA_real_`, the function errors.
#'
#' @param residual Data frame returned by [implied_demographic_residual()], or a
#'   data frame with required columns `time_start`, `time_end`, `dt`,
#'   `age_group`, `residual_count`, and `residual_rate`.
#' @param age_structure Age structure validated by [validate_age_structure()].
#'   Residual age groups must match this age structure exactly.
#' @param use Which residual convention to convert: `"count"` converts
#'   interval residual counts to per-time `migration_count` flows via
#'   `residual_count / dt`; `"rate"` uses `residual_rate` directly as
#'   `migration_rate`.
#'
#' @return An `agepi_migration_schedule` object.
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("0-4", "5+"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, Inf)
#' )
#' residual <- data.frame(
#'   time_start = c(0, 0),
#'   time_end = c(2, 2),
#'   dt = c(2, 2),
#'   age_group = c("0-4", "5+"),
#'   residual_count = c(10, -4),
#'   residual_rate = c(0.05, -0.01)
#' )
#' residual_to_migration_schedule(residual, ages, use = "count")
#' @export
residual_to_migration_schedule <- function(residual,
                                           age_structure,
                                           use = c("count", "rate")) {
  use <- match.arg(use)
  validate_residual_to_migration_inputs(residual, age_structure, use)

  if (use == "count") {
    migration <- data.frame(
      time = residual$time_start,
      age_group = as.character(residual$age_group),
      migration_count = residual$residual_count / residual$dt,
      stringsAsFactors = FALSE
    )
  } else {
    migration <- data.frame(
      time = residual$time_start,
      age_group = as.character(residual$age_group),
      migration_rate = residual$residual_rate,
      stringsAsFactors = FALSE
    )
  }

  MigrationSchedule(migration, age_structure)
}

validate_residual_to_migration_inputs <- function(residual, age_structure, use) {
  validate_age_structure(age_structure)

  if (!is.data.frame(residual)) {
    stop("residual must be a data frame.", call. = FALSE)
  }

  required_columns <- c(
    "time_start",
    "time_end",
    "dt",
    "age_group",
    "residual_count",
    "residual_rate"
  )
  missing_columns <- setdiff(required_columns, names(residual))
  if (length(missing_columns) > 0) {
    stop(
      "residual is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(residual$time_start) ||
      anyNA(residual$time_start) ||
      any(!is.finite(residual$time_start))) {
    stop("residual time_start must contain finite non-missing numeric values.", call. = FALSE)
  }
  if (!is.numeric(residual$time_end) ||
      anyNA(residual$time_end) ||
      any(!is.finite(residual$time_end))) {
    stop("residual time_end must contain finite non-missing numeric values.", call. = FALSE)
  }
  if (!is.numeric(residual$dt) || anyNA(residual$dt) || any(!is.finite(residual$dt))) {
    stop("residual dt must contain finite non-missing numeric values.", call. = FALSE)
  }
  if (any(residual$dt <= 0)) {
    stop("residual dt must be positive.", call. = FALSE)
  }

  residual_age_groups <- as.character(residual$age_group)
  if (anyNA(residual_age_groups)) {
    stop("residual age_group cannot contain missing values.", call. = FALSE)
  }
  unknown_ages <- setdiff(unique(residual_age_groups), age_structure$age_groups)
  if (length(unknown_ages) > 0) {
    stop(
      "residual contains age_group value(s) not in age_structure: ",
      paste(unknown_ages, collapse = ", "),
      call. = FALSE
    )
  }

  for (this_time in unique(residual$time_start)) {
    observed_age_groups <- residual_age_groups[residual$time_start == this_time]
    missing_age_groups <- setdiff(age_structure$age_groups, observed_age_groups)
    if (length(missing_age_groups) > 0) {
      stop(
        "residual is missing age_group value(s) at time_start ",
        this_time,
        ": ",
        paste(missing_age_groups, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.numeric(residual$residual_count) ||
      anyNA(residual$residual_count) ||
      any(!is.finite(residual$residual_count))) {
    stop("residual_count must contain finite non-missing numeric values.", call. = FALSE)
  }

  if (!is.numeric(residual$residual_rate) || any(!is.finite(residual$residual_rate) & !is.na(residual$residual_rate))) {
    stop("residual_rate must contain numeric finite values or NA.", call. = FALSE)
  }
  if (use == "rate" && anyNA(residual$residual_rate)) {
    stop("residual_rate cannot contain NA values when use = \"rate\".", call. = FALSE)
  }

  invisible(NULL)
}
