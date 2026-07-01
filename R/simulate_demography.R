#' Simulate demographic-only one-compartment population dynamics
#'
#' Simulates one-compartment age-structured population dynamics from an
#' [DemographicProcess()] using [demographic_derivative()]. This is a
#' demographic-only wrapper: it does not implement infection dynamics,
#' infection-demography coupling, residual forcing, WPP validation diagnostics,
#' or WPP parsing.
#'
#' Schedule lookup is exact-time only by default. With `time_policy = "exact"`,
#' derivative evaluation times must exactly match required schedule times. With
#' `time_policy = "step"`, schedules use left-continuous interval-start lookup:
#' a row at schedule time `t_i` applies over `[t_i, t_{i + 1})`. Euler
#' intervals are evaluated at their left endpoints, `times[i]`. With
#' `time_policy = "linear"`, fertility, mortality, and migration schedules are
#' linearly interpolated without extrapolation. The optional `method =
#' "deSolve"` and its alias `method = "ode"` use `deSolve::ode()` when the
#' suggested `deSolve` package is installed.
#'
#' Euler updates are not truncated: if an update would produce a negative
#' population value, simulation stops with an error.
#'
#' @param process Valid `agepi_demographic_process` object.
#' @param initial_state Numeric vector of non-negative population counts ordered
#'   by `process$age_structure$age_groups`. Names are ignored on input and
#'   restored to the process age groups internally.
#' @param times Numeric vector of finite, non-missing, strictly increasing time
#'   points. Must have length at least two.
#' @param method Simulation method. `NULL` preserves the demographic wrapper's
#'   explicit Euler default. `"deSolve"` and `"ode"` request the optional
#'   `deSolve::ode()` backend. `"euler"` requests explicit Euler time steps.
#' @param time_policy Schedule lookup policy. `"exact"` requires exact schedule
#'   times and remains the default. `"step"` uses interval-start stepwise lookup.
#'   `"linear"` interpolates rate-like demographic schedules only.
#' @param ageing_policy Ageing implementation. `"exponential"` preserves the
#'   existing derivative-based ageing implementation. `"annual_cohort"` applies
#'   `annual_cohort_demographic_step()` once per annual interval and requires a
#'   complete 1-year age grid ending in an open-ended age group.
#' @param ... Reserved for future method-specific arguments. Currently unused.
#'
#' @return Data frame with columns `time`, `age_group`, and `population`,
#'   ordered by time outermost and process age groups innermost.
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("0-4", "5+"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, Inf)
#' )
#' process <- DemographicProcess(age_structure = ages)
#' simulate_demography(
#'   process = process,
#'   initial_state = c(100, 50),
#'   times = c(2020, 2021)
#' )
#' @export
simulate_demography <- function(
  process,
  initial_state,
  times,
  method = NULL,
  time_policy = c("exact", "step", "linear"),
  ageing_policy = c("exponential", "annual_cohort"),
  ...
) {
  method <- if (is.null(method)) "euler" else validate_simulation_method(method)
  validate_demographic_process(process)
  validate_simulation_times(times)
  time_policy <- validate_demographic_time_policy(time_policy)
  ageing_policy <- validate_demography_ageing_policy(ageing_policy)
  check_dots_empty(...)

  state <- validate_demography_initial_state(initial_state, process)

  if (identical(ageing_policy, "annual_cohort")) {
    return(simulate_demography_annual_cohort(
      process = process,
      initial_state = state,
      times = times,
      method = method,
      time_policy = time_policy
    ))
  }

  simulate_demography_integrated(
    process = process,
    initial_state = state,
    times = times,
    method = method,
    time_policy = time_policy
  )
}

simulate_demography_annual_cohort <- function(process,
                                              initial_state,
                                              times,
                                              method,
                                              time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)
  if (!identical(method, "euler")) {
    stop("ageing_policy = \"annual_cohort\" requires method = NULL or method = \"euler\".", call. = FALSE)
  }

  validate_annual_cohort_simulation_times(times)
  validate_simulate_demography_annual_cohort_age_structure(process$age_structure)
  validate_demography_schedule_coverage(
    process,
    times,
    time_policy,
    include_output_times = FALSE
  )

  age_groups <- process$age_structure$age_groups
  output <- vector("list", length(times))
  state <- stats::setNames(as.numeric(initial_state), age_groups)
  output[[1]] <- demographic_state_output(state, time = times[1], process = process)

  for (time_index in seq_len(length(times) - 1)) {
    interval_start <- times[time_index]
    fertility <- fertility_rates_at(
      process$fertility_schedule,
      time = interval_start,
      age_groups = age_groups,
      time_policy = time_policy
    )
    mortality <- mortality_rates_at(
      process$mortality_schedule,
      time = interval_start,
      age_groups = age_groups,
      time_policy = time_policy
    )
    migration <- annual_cohort_migration_values_at(
      process$migration_schedule,
      time = interval_start,
      age_groups = age_groups,
      time_policy = time_policy
    )

    step <- annual_cohort_demographic_step(
      population = state,
      age_structure = process$age_structure,
      fertility = stats::setNames(fertility, age_groups),
      mortality = stats::setNames(mortality, age_groups),
      migration = migration$values,
      fertility_exposure_fraction = process$fertility_exposure_fraction,
      migration_type = migration$type
    )

    state <- stats::setNames(step$population, age_groups)
    output[[time_index + 1]] <- demographic_state_output(
      state,
      time = times[time_index + 1],
      process = process
    )
  }

  do.call(rbind, output)
}

annual_cohort_migration_values_at <- function(schedule,
                                              time,
                                              age_groups,
                                              time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)
  if (is.null(schedule)) {
    return(list(
      type = NULL,
      values = NULL
    ))
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

  list(
    type = schedule$migration_type,
    values = stats::setNames(values, age_groups)
  )
}

validate_demography_ageing_policy <- function(ageing_policy = c("exponential", "annual_cohort")) {
  match.arg(ageing_policy)
}

validate_annual_cohort_simulation_times <- function(times,
                                                    tolerance = sqrt(.Machine$double.eps)) {
  time_steps <- diff(times)
  if (any(abs(time_steps - 1) > tolerance)) {
    stop("ageing_policy = \"annual_cohort\" requires annual time steps of exactly 1 year.", call. = FALSE)
  }

  invisible(times)
}

validate_simulate_demography_annual_cohort_age_structure <- function(age_structure) {
  validate_annual_cohort_age_structure(age_structure)

  if (!identical(age_structure$lower_bounds[1], 0)) {
    stop("ageing_policy = \"annual_cohort\" requires the first age group to start at age 0.", call. = FALSE)
  }

  if (!is.infinite(age_structure$upper_bounds[age_structure$n_age_groups])) {
    stop("ageing_policy = \"annual_cohort\" requires a terminal open-ended age group.", call. = FALSE)
  }

  invisible(age_structure)
}

simulate_demography_euler <- function(process, initial_state, times, time_policy = c("exact", "step", "linear")) {
  simulate_demography_integrated(
    process = process,
    initial_state = initial_state,
    times = times,
    method = "euler",
    time_policy = time_policy
  )
}

simulate_demography_integrated <- function(process,
                                           initial_state,
                                           times,
                                           method,
                                           time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)
  validate_demography_schedule_coverage(
    process,
    times,
    time_policy,
    include_output_times = TRUE
  )

  integrate_state_trajectory(
    initial_state = initial_state,
    times = times,
    method = method,
    derivative = function(time, state) {
      demographic_derivative(
        state,
        time = time,
        process = process,
        time_policy = time_policy
      )
    },
    output = function(state, time) {
      demographic_state_output(
        stats::setNames(as.numeric(state), process$age_structure$age_groups),
        time = time,
        process = process
      )
    },
    non_negative = validate_non_negative_demography_euler_state,
    tcrit = desolve_schedule_tcrit(process, times),
    desolve_error = paste(
      "method = \"deSolve\" requires the deSolve package.",
      "Install deSolve or use method = \"euler\"."
    )
  )
}

validate_demography_schedule_coverage <- function(process,
                                                  times,
                                                  time_policy = c("exact", "step", "linear"),
                                                  include_output_times = FALSE) {
  time_policy <- validate_demographic_time_policy(time_policy)
  lookup_times <- times[-length(times)]

  if (include_output_times && time_policy == "linear") {
    lookup_times <- times
  }

  validate_one_demography_schedule_coverage(process$fertility_schedule, lookup_times, time_policy)
  validate_one_demography_schedule_coverage(process$mortality_schedule, lookup_times, time_policy)
  validate_one_demography_schedule_coverage(process$migration_schedule, lookup_times, time_policy)

  invisible(NULL)
}

validate_one_demography_schedule_coverage <- function(schedule, times, time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)

  if (is.null(schedule)) {
    return(invisible(NULL))
  }

  for (time in times) {
    lookup_demographic_schedule_time(schedule, time, time_policy)
  }

  invisible(NULL)
}

demographic_state_output <- function(state, time, process) {
  data.frame(
    time = rep(time, process$age_structure$n_age_groups),
    age_group = process$age_structure$age_groups,
    population = as.numeric(state),
    stringsAsFactors = FALSE
  )
}

validate_demography_initial_state <- function(initial_state, process) {
  n_age_groups <- process$age_structure$n_age_groups

  if (!is.numeric(initial_state) || is.matrix(initial_state) || is.data.frame(initial_state)) {
    stop("initial_state must be a numeric vector.", call. = FALSE)
  }

  if (length(initial_state) != n_age_groups) {
    stop("initial_state length must equal the number of age groups in process.", call. = FALSE)
  }

  if (anyNA(initial_state) || any(!is.finite(initial_state))) {
    stop("initial_state values must be finite and non-missing.", call. = FALSE)
  }

  if (any(initial_state < 0)) {
    stop("initial_state values must be non-negative.", call. = FALSE)
  }

  stats::setNames(as.numeric(initial_state), process$age_structure$age_groups)
}

validate_demography_simulation_method <- function(method) {
  validate_simulation_method(method)
}

validate_non_negative_demography_euler_state <- function(state, time) {
  if (any(!is.finite(state))) {
    stop("Euler step produced non-finite population value at time ", time, ".", call. = FALSE)
  }

  if (any(state < 0)) {
    stop("Euler step produced negative population value at time ", time, ".", call. = FALSE)
  }

  invisible(state)
}

check_dots_empty <- function(...) {
  dots <- list(...)
  if (length(dots) > 0) {
    stop("unused arguments in ... are not currently supported.", call. = FALSE)
  }

  invisible(NULL)
}
