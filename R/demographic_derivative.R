#' Compute a demographic-only population derivative
#'
#' Computes `dN/dt` for a one-compartment age-structured population using
#' ageing, optional fertility, optional mortality, and optional migration from
#' an [DemographicProcess()] object. This is a demographic-only derivative: it
#' does not simulate trajectories, interpolate schedules, couple to infection
#' compartments, or implement any infection dynamics.
#'
#' Schedule values use exact-time lookup only. If a supplied schedule does not
#' contain `time`, an error is raised; no interpolation is performed. Fertility
#' currently uses the total state in each age group as the exposure population,
#' as a first-pass convention, and births enter the youngest age group only.
#'
#' @param state Numeric vector of non-negative population counts ordered by the
#'   age groups in `process`.
#' @param time Single finite numeric time at which to evaluate supplied
#'   schedules.
#' @param process Valid `agepi_demographic_process` object.
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
demographic_derivative <- function(state, time, process) {
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

  fertility_rates <- fertility_rates_at(process$fertility_schedule, time, age_groups)
  births <- sum(fertility_rates * state)
  dNdt[1] <- dNdt[1] + births

  mortality_rates <- mortality_rates_at(process$mortality_schedule, time, age_groups)
  dNdt <- dNdt - mortality_rates * state

  migration <- migration_values_at(process$migration_schedule, time, state, age_groups)
  dNdt <- dNdt + migration

  if (any(!is.finite(dNdt))) {
    stop("demographic_derivative result must contain only finite values.", call. = FALSE)
  }

  stats::setNames(dNdt, age_groups)
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

schedule_values_at <- function(schedule, time, value_column, age_groups, fill_value = NA_real_) {
  values <- rep(fill_value, length(age_groups))

  if (is.null(schedule)) {
    return(values)
  }

  if (!time %in% schedule$times) {
    stop(
      "Exact time ",
      time,
      " is not available in the supplied schedule; exact-time lookup only, no interpolation is performed.",
      call. = FALSE
    )
  }

  rows <- schedule$data$time == time
  age_index <- match(schedule$data$age_group[rows], age_groups)
  values[age_index] <- schedule$data[[value_column]][rows]
  values
}

fertility_rates_at <- function(schedule, time, age_groups) {
  schedule_values_at(
    schedule = schedule,
    time = time,
    value_column = "fertility_rate",
    age_groups = age_groups,
    fill_value = 0
  )
}

mortality_rates_at <- function(schedule, time, age_groups) {
  schedule_values_at(
    schedule = schedule,
    time = time,
    value_column = "mortality_rate",
    age_groups = age_groups,
    fill_value = 0
  )
}

migration_values_at <- function(schedule, time, state, age_groups) {
  if (is.null(schedule)) {
    return(rep(0, length(age_groups)))
  }

  value_column <- paste0("migration_", schedule$migration_type)
  values <- schedule_values_at(
    schedule = schedule,
    time = time,
    value_column = value_column,
    age_groups = age_groups,
    fill_value = 0
  )

  if (identical(schedule$migration_type, "rate")) {
    return(values * state)
  }

  values
}
