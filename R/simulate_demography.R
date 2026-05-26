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
#' a row at schedule time `t_i` applies over `[t_i, t_{i + 1})`, and Euler
#' intervals are evaluated at their left endpoints, `times[i]`. With
#' `time_policy = "linear"`, fertility, mortality, and migration schedules are
#' linearly interpolated at Euler left endpoints without extrapolation. The
#' optional `method = "deSolve"` and its alias `method = "ode"` are reserved for
#' future support and currently error clearly.
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
#' @param method Simulation method. `"euler"` is the default. `"deSolve"` and
#'   `"ode"` are accepted for future compatibility but currently error for this
#'   demographic wrapper.
#' @param time_policy Schedule lookup policy. `"exact"` requires exact schedule
#'   times and remains the default. `"step"` uses interval-start stepwise lookup.
#'   `"linear"` interpolates rate-like demographic schedules only.
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
  method = c("euler", "deSolve", "ode"),
  time_policy = c("exact", "step", "linear"),
  ...
) {
  if (missing(method)) {
    method <- "euler"
  } else {
    method <- validate_demography_simulation_method(method)
  }
  validate_demographic_process(process)
  validate_simulation_times(times)
  time_policy <- validate_demographic_time_policy(time_policy)
  check_dots_empty(...)

  state <- validate_demography_initial_state(initial_state, process)

  if (method == "deSolve") {
    stop(
      "method = \"deSolve\" is not yet supported by simulate_demography(). ",
      "Demographic schedules use fixed Euler evaluation points only in this wrapper, ",
      "so adaptive deSolve solver time points cannot be evaluated safely. Use method = \"euler\".",
      call. = FALSE
    )
  }

  simulate_demography_euler(
    process = process,
    initial_state = state,
    times = times,
    time_policy = time_policy
  )
}

simulate_demography_euler <- function(process, initial_state, times, time_policy = c("exact", "step", "linear")) {
  time_policy <- validate_demographic_time_policy(time_policy)
  validate_demography_schedule_coverage(
    process,
    times,
    time_policy,
    include_output_times = TRUE
  )

  output <- vector("list", length(times))
  output[[1]] <- demographic_state_output(initial_state, time = times[1], process = process)

  current_state <- initial_state
  for (i in seq_len(length(times) - 1)) {
    dt <- times[i + 1] - times[i]
    derivative <- demographic_derivative(
      current_state,
      time = times[i],
      process = process,
      time_policy = time_policy
    )
    next_state <- as.numeric(current_state) + dt * as.numeric(derivative)
    validate_non_negative_demography_euler_state(next_state, time = times[i + 1])

    current_state <- stats::setNames(next_state, process$age_structure$age_groups)
    output[[i + 1]] <- demographic_state_output(current_state, time = times[i + 1], process = process)
  }

  result <- do.call(rbind, output)
  row.names(result) <- NULL
  result
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
  if (!is.character(method) || length(method) != 1 || anyNA(method) || method == "") {
    stop("method must be a non-missing character scalar.", call. = FALSE)
  }

  if (!method %in% c("euler", "deSolve", "ode")) {
    stop("unsupported simulation method: ", method, call. = FALSE)
  }

  if (method == "ode") {
    return("deSolve")
  }

  method
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
