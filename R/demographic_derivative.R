#' Compute a demographic-only population derivative
#'
#' Computes `dN/dt` for a one-compartment age-structured population using
#' ageing, optional fertility, optional mortality, and optional migration from
#' an [DemographicProcess()] object. This is a demographic-only derivative: it
#' does not simulate trajectories, couple to infection compartments, or
#' implement any infection dynamics.
#'
#' Schedule values use exact-time lookup by default. If `time_policy = "step"`,
#' schedule values use left-continuous interval-start lookup: for evaluation
#' time `t`, the greatest schedule time less than or equal to `t` is used. With
#' `time_policy = "linear"`, rate-like schedules are linearly interpolated
#' between bracketing schedule times without extrapolation. Fertility currently
#' uses the total state in each age group as the exposure population, as a
#' first-pass convention, scaled by `process$fertility_exposure_fraction`, and
#' births enter the youngest age group only.
#'
#' @param state Numeric vector of non-negative population counts ordered by the
#'   age groups in `process`.
#' @param time Single finite numeric time at which to evaluate supplied
#'   schedules.
#' @param process Valid `agepi_demographic_process` object.
#' @param time_policy Schedule lookup policy. `"exact"` requires `time` to be
#'   present in each supplied schedule. `"step"` uses the greatest schedule time
#'   less than or equal to `time`. `"linear"` interpolates rate-like schedules
#'   between bracketing schedule times. Neither `"step"` nor `"linear"`
#'   extrapolates outside the schedule range.
#'
#' @return A named numeric vector of derivatives ordered by the age groups in
#'   `process`.
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("0-4", "5+"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, Inf)
#' )
#' process <- DemographicProcess(age_structure = ages)
#' demographic_derivative(c("0-4" = 100, "5+" = 200), time = 2020, process)
#' @export
demographic_derivative <- function(state,
                                   time,
                                   process,
                                   time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)
  validate_demographic_derivative_inputs(state, time, process)

  state <- as.numeric(state)
  age_groups <- process$age_structure$age_groups
  n_age_groups <- process$age_structure$n_age_groups
  ageing <- process$ageing_operator

  dNdt <- numeric(n_age_groups)

  ageing_out <- ageing$departure_rate * state
  dNdt <- dNdt - ageing_out
  has_destination <- !is.na(ageing$destination_index)
  dNdt[ageing$destination_index[has_destination]] <-
    dNdt[ageing$destination_index[has_destination]] + ageing_out[has_destination]

  fertility_rates <- fertility_rates_at(process$fertility_schedule, time, age_groups, time_policy)
  births <- sum(fertility_rates * process$fertility_exposure_fraction * state)
  dNdt[1] <- dNdt[1] + births

  mortality_rates <- mortality_rates_at(process$mortality_schedule, time, age_groups, time_policy)
  dNdt <- dNdt - mortality_rates * state

  migration <- migration_values_at(process$migration_schedule, time, state, age_groups, time_policy)
  dNdt <- dNdt + migration

  if (any(!is.finite(dNdt))) {
    stop("demographic_derivative result must contain only finite values.", call. = FALSE)
  }

  stats::setNames(dNdt, age_groups)
}

validate_demographic_time_policy <- function(time_policy) {
  if (!is.character(time_policy) || anyNA(time_policy) || any(time_policy == "")) {
    stop("time_policy must be non-missing character value(s).", call. = FALSE)
  }

  time_policy <- tryCatch(
    match.arg(time_policy, c("exact", "step", "linear")),
    error = function(e) {
      stop("unsupported time_policy: ", paste(time_policy, collapse = ", "), call. = FALSE)
    }
  )

  time_policy
}

validate_demographic_derivative_inputs <- function(state, time, process) {
  validate_demographic_process(process)

  n_age_groups <- process$age_structure$n_age_groups
  if (!is.numeric(state) || length(state) != n_age_groups) {
    stop("state must be a numeric vector with length equal to the number of age groups in process.", call. = FALSE)
  }

  if (anyNA(state) || any(!is.finite(state))) {
    stop("state values must be finite and non-missing.", call. = FALSE)
  }

  if (any(state < 0)) {
    stop("state values must be non-negative.", call. = FALSE)
  }

  if (!is.numeric(time) || length(time) != 1 || anyNA(time) || !is.finite(time)) {
    stop("time must be a single finite numeric value.", call. = FALSE)
  }

  invisible(NULL)
}

schedule_values_at <- function(schedule,
                               time,
                               value_column,
                               age_groups,
                               fill_value = NA_real_,
                               time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)
  values <- rep(fill_value, length(age_groups))

  if (is.null(schedule)) {
    return(values)
  }

  if (time_policy == "linear") {
    return(interpolate_demographic_schedule_values(
      schedule = schedule,
      time = time,
      value_column = value_column,
      age_groups = age_groups,
      fill_value = fill_value
    ))
  }

  schedule_time <- lookup_demographic_schedule_time(schedule, time, time_policy)
  values <- demographic_schedule_values_for_time(
    schedule = schedule,
    schedule_time = schedule_time,
    value_column = value_column,
    age_groups = age_groups,
    fill_value = fill_value
  )
  values
}

lookup_demographic_schedule_time <- function(schedule, time, time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)

  if (time_policy == "exact" || (time_policy == "linear" && time %in% schedule$times)) {
    if (!time %in% schedule$times) {
      stop(
        "Exact time ",
        time,
        " is not available in the supplied schedule; exact-time lookup only, no interpolation is performed.",
        call. = FALSE
      )
    }

    return(time)
  }

  if (time_policy == "linear" && length(schedule$times) < 2) {
    stop(
      "linear time_policy requires at least two schedule times for off-grid interpolation.",
      call. = FALSE
    )
  }

  if (time < min(schedule$times)) {
    stop(
      "Time ",
      time,
      " is before the first available schedule time ",
      min(schedule$times),
      "; ",
      time_policy,
      " lookup cannot extrapolate before the schedule starts.",
      call. = FALSE
    )
  }

  if (time > max(schedule$times)) {
    stop(
      "Time ",
      time,
      " is after the final available schedule time ",
      max(schedule$times),
      "; ",
      time_policy,
      " lookup does not silently extend schedules beyond their final time.",
      call. = FALSE
    )
  }

  if (time_policy == "linear") {
    return(c(
      lower = max(schedule$times[schedule$times < time]),
      upper = min(schedule$times[schedule$times > time])
    ))
  }

  max(schedule$times[schedule$times <= time])
}

demographic_schedule_values_for_time <- function(schedule,
                                                 schedule_time,
                                                 value_column,
                                                 age_groups,
                                                 fill_value) {
  values <- rep(fill_value, length(age_groups))
  rows <- schedule$data$time == schedule_time
  age_index <- match(schedule$data$age_group[rows], age_groups)
  values[age_index] <- schedule$data[[value_column]][rows]
  values
}

interpolate_demographic_schedule_values <- function(schedule,
                                                    time,
                                                    value_column,
                                                    age_groups,
                                                    fill_value) {
  schedule_times <- lookup_demographic_schedule_time(schedule, time, "linear")

  if (length(schedule_times) == 1) {
    return(demographic_schedule_values_for_time(
      schedule = schedule,
      schedule_time = schedule_times,
      value_column = value_column,
      age_groups = age_groups,
      fill_value = fill_value
    ))
  }

  lower_time <- unname(schedule_times["lower"])
  upper_time <- unname(schedule_times["upper"])
  lower_rows <- schedule$data$time == lower_time
  upper_rows <- schedule$data$time == upper_time
  lower_age_groups <- schedule$data$age_group[lower_rows]
  upper_age_groups <- schedule$data$age_group[upper_rows]

  if (!setequal(lower_age_groups, upper_age_groups)) {
    stop(
      "linear time_policy requires consistent age-group coverage across schedule times.",
      call. = FALSE
    )
  }

  lower_values <- demographic_schedule_values_for_time(
    schedule = schedule,
    schedule_time = lower_time,
    value_column = value_column,
    age_groups = age_groups,
    fill_value = fill_value
  )
  upper_values <- demographic_schedule_values_for_time(
    schedule = schedule,
    schedule_time = upper_time,
    value_column = value_column,
    age_groups = age_groups,
    fill_value = fill_value
  )
  weight <- (time - lower_time) / (upper_time - lower_time)
  lower_values + weight * (upper_values - lower_values)
}

fertility_rates_at <- function(schedule, time, age_groups, time_policy = c("exact", "step", "linear")) {
  schedule_values_at(
    schedule = schedule,
    time = time,
    value_column = "fertility_rate",
    age_groups = age_groups,
    fill_value = 0,
    time_policy = time_policy
  )
}

mortality_rates_at <- function(schedule, time, age_groups, time_policy = c("exact", "step", "linear")) {
  schedule_values_at(
    schedule = schedule,
    time = time,
    value_column = "mortality_rate",
    age_groups = age_groups,
    fill_value = 0,
    time_policy = time_policy
  )
}

migration_values_at <- function(schedule, time, state, age_groups, time_policy = c("exact", "step", "linear")) {
  if (is.null(schedule)) {
    return(rep(0, length(age_groups)))
  }

  value_column <- paste0("migration_", schedule$migration_type)
  values <- schedule_values_at(
    schedule = schedule,
    time = time,
    value_column = value_column,
    age_groups = age_groups,
    fill_value = 0,
    time_policy = time_policy
  )

  if (identical(schedule$migration_type, "rate")) {
    return(values * state)
  }

  values
}
