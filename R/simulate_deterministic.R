#' Simulate a deterministic SIR or SEIR model
#'
#' Runs a simple deterministic prototype simulation. The default
#' `method = "euler"` uses explicit Euler time steps. Optional
#' `method = "deSolve"` uses `deSolve::ode()` when the suggested `deSolve`
#' package is installed.
#'
#' The initial state may be supplied either as long-form state data or as a
#' numeric vector in the existing compartment-major, age-group-minor ordering.
#' Names on numeric state vectors are ignored, matching the state-mapping
#' helpers. The initial state is included in the returned output.
#'
#' Supported combinations are:
#'
#' - SIR + Euler, infection-only;
#' - SEIR + Euler, infection-only;
#' - SIR + deSolve, infection-only;
#' - SEIR + deSolve, infection-only;
#' - SIR + demography + Euler;
#' - SEIR + demography + Euler;
#' - SIR + demography + deSolve;
#' - SEIR + demography + deSolve.
#'
#' Both solver backends use the same transition-rate and derivative pathway:
#' static contact matrix, static `beta`, static susceptibility, and static
#' infectiousness. The deSolve backend also supports demographic coupling via
#' the existing demographic derivative path.
#'
#' A first-pass SIR/SEIR-demography coupling is available via
#' `demographic_process`. In this coupled mode, births enter the youngest
#' susceptible age group, background mortality applies independently to every
#' disease compartment, ageing moves each compartment independently using the
#' process [AgeingOperator()] convention, and net migration is applied to
#' susceptible compartments only. Fertility and migration-rate exposure use the
#' current total infection-state population by age group: `S + I + R` for SIR
#' and `S + E + I + R` for SEIR. This is not WPP projection matching, and
#' `S`-only migration is an allocation convention for age-total net migration
#' inputs rather than a mechanistic movement model. For adaptive deSolve runs,
#' `time_policy = "linear"` is generally recommended for demographic schedules;
#' `time_policy = "step"` gives piecewise-constant rates, and
#' `time_policy = "exact"` will usually fail unless every solver evaluation
#' time is present in the schedules.
#'
#' Euler updates are intentionally not truncated: if a step would produce
#' negative compartment values, simulation stops with an error.
#'
#' @param initial_state Long-form state data frame with columns `compartment`,
#'   `age_group`, and `value`, or numeric state vector.
#' @param times Numeric vector of finite, non-missing, strictly increasing time
#'   points. Must have length at least two.
#' @param model Disease model. `SIRModel()` and `SEIRModel()` output are
#'   supported for infection-only deterministic simulation.
#' @param age_structure Valid age structure.
#' @param contact_matrix Numeric contact matrix with rows as recipient age
#'   groups and columns as source age groups.
#' @param beta Non-negative finite transmission scaling parameter.
#' @param susceptibility Optional non-negative numeric vector by recipient age
#'   group.
#' @param infectiousness Optional non-negative numeric vector by source age
#'   group.
#' @param method Simulation method. `"euler"` is the default. `"deSolve"` and
#'   `"ode"` request the optional `deSolve::ode()` backend.
#' @param demographic_process Optional [DemographicProcess()] object for
#'   first-pass deterministic SIR/SEIR-demography coupling. Defaults to `NULL`,
#'   which preserves infection-only simulation.
#' @param time_policy Demographic schedule lookup policy used only when
#'   `demographic_process` is supplied. `"exact"` requires exact schedule times;
#'   `"step"` uses left-continuous interval-start lookup; `"linear"`
#'   interpolates rate-like demographic schedules through the same demographic
#'   derivative path.
#'
#' @return Data frame with columns `time`, `compartment`, `age_group`, and
#'   `value`, ordered by time outermost, compartment next, and age group
#'   innermost.
#' @export
simulate_deterministic <- function(
  initial_state,
  times,
  model,
  age_structure,
  contact_matrix,
  beta = 1,
  susceptibility = NULL,
  infectiousness = NULL,
  method = "euler",
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear")
) {
  method <- validate_simulation_method(method)
  validate_simulation_times(times)
  validate_disease_model(model)
  validate_age_structure(age_structure)
  time_policy <- validate_simulation_demography_inputs(
    demographic_process = demographic_process,
    time_policy = time_policy,
    method = method,
    model = model,
    age_structure = age_structure,
    times = times
  )

  state_vector <- simulation_state_to_vector(
    initial_state,
    age_structure,
    model$compartments
  )
  validate_non_negative_simulation_state(state_vector, "initial_state")

  if (method == "euler") {
    return(simulate_deterministic_euler(
      state_vector = state_vector,
      times = times,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness,
      demographic_process = demographic_process,
      time_policy = time_policy
    ))
  }

  simulate_deterministic_desolve(
    state_vector = state_vector,
    times = times,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    demographic_process = demographic_process,
    time_policy = time_policy
  )
}

simulate_deterministic_euler <- function(
  state_vector,
  times,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear")
) {
  output <- vector("list", length(times))
  output[[1]] <- simulation_state_output(
    state_vector,
    time = times[1],
    age_structure = age_structure,
    compartments = model$compartments
  )

  current_state <- state_vector
  for (i in seq_len(length(times) - 1)) {
    dt <- times[i + 1] - times[i]
    derivative <- deterministic_euler_derivative(
      state_vector = current_state,
      time = times[i],
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness,
      demographic_process = demographic_process,
      time_policy = time_policy
    )

    next_state <- as.numeric(current_state) + dt * derivative
    validate_non_negative_euler_state(next_state, time = times[i + 1])

    current_state <- next_state
    output[[i + 1]] <- simulation_state_output(
      current_state,
      time = times[i + 1],
      age_structure = age_structure,
      compartments = model$compartments
    )
  }

  result <- do.call(rbind, output)
  row.names(result) <- NULL
  result
}

deterministic_euler_derivative <- function(
  state_vector,
  time,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear")
) {
  rates <- transition_rates(
    state = state_vector,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness
  )
  derivative <- rates_to_derivative(
    transition_rate_table = rates,
    compartments = model$compartments,
    age_structure = age_structure
  )
  infection_derivative <- derivative$derivative

  if (is.null(demographic_process)) {
    return(infection_derivative)
  }

  infection_derivative + compartment_demographic_derivative(
    state_vector = state_vector,
    time = time,
    model = model,
    age_structure = age_structure,
    demographic_process = demographic_process,
    time_policy = time_policy
  )
}

compartment_demographic_derivative <- function(
  state_vector,
  time,
  model,
  age_structure,
  demographic_process,
  time_policy = c("exact", "step", "linear")
) {
  time_policy <- validate_demographic_time_policy(time_policy)
  validate_disease_model(model)
  validate_age_structure(age_structure)
  validate_demographic_process(demographic_process)

  if (!"S" %in% model$compartments) {
    stop("demographic_process coupling requires an S compartment.", call. = FALSE)
  }

  state_long <- state_vector_to_long(state_vector, age_structure, model$compartments)
  age_groups <- age_structure$age_groups
  n_age_groups <- age_structure$n_age_groups

  compartment_indices <- lapply(seq_along(model$compartments), function(compartment_position) {
    ((compartment_position - 1) * n_age_groups) + seq_len(n_age_groups)
  })
  names(compartment_indices) <- model$compartments

  population <- numeric(n_age_groups)
  for (compartment in model$compartments) {
    population <- population +
      transition_compartment_values(state_long, age_structure, compartment)
  }

  derivative <- numeric(length(state_vector))
  S_index <- compartment_indices[["S"]]

  ageing <- demographic_process$ageing_operator
  for (index in compartment_indices) {
    compartment_state <- as.numeric(state_vector[index])
    ageing_out <- ageing$departure_rate * compartment_state
    derivative[index] <- derivative[index] - ageing_out

    has_destination <- !is.na(ageing$destination_index)
    derivative[index[ageing$destination_index[has_destination]]] <-
      derivative[index[ageing$destination_index[has_destination]]] +
      ageing_out[has_destination]
  }

  fertility_rates <- fertility_rates_at(
    demographic_process$fertility_schedule,
    time,
    age_groups,
    time_policy
  )
  derivative[S_index[1]] <- derivative[S_index[1]] + sum(fertility_rates * population)

  mortality_rates <- mortality_rates_at(
    demographic_process$mortality_schedule,
    time,
    age_groups,
    time_policy
  )
  for (index in compartment_indices) {
    derivative[index] <- derivative[index] - mortality_rates * as.numeric(state_vector[index])
  }

  migration <- migration_values_at(
    demographic_process$migration_schedule,
    time,
    population,
    age_groups,
    time_policy
  )
  derivative[S_index] <- derivative[S_index] + migration

  if (any(!is.finite(derivative))) {
    stop("compartment_demographic_derivative result must contain only finite values.", call. = FALSE)
  }

  derivative
}

simulate_deterministic_desolve <- function(
  state_vector,
  times,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear")
) {
  if (!desolve_is_available()) {
    stop(
      "method = \"deSolve\" requires the optional deSolve package. ",
      "Install deSolve or use method = \"euler\".",
      call. = FALSE
    )
  }

  solved <- deSolve::ode(
    y = as.numeric(state_vector),
    times = times,
    func = deterministic_derivative_vector,
    parms = list(
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness,
      demographic_process = demographic_process,
      time_policy = time_policy
    )
  )

  output <- vector("list", length(times))
  state_columns <- seq_len(length(state_vector)) + 1
  for (i in seq_along(times)) {
    output[[i]] <- simulation_state_output(
      as.numeric(solved[i, state_columns]),
      time = solved[i, "time"],
      age_structure = age_structure,
      compartments = model$compartments
    )
  }

  result <- do.call(rbind, output)
  row.names(result) <- NULL
  result
}

desolve_is_available <- function() {
  requireNamespace("deSolve", quietly = TRUE)
}

deterministic_derivative_vector <- function(time, state, parms) {
  rates <- transition_rates(
    state = as.numeric(state),
    model = parms$model,
    age_structure = parms$age_structure,
    contact_matrix = parms$contact_matrix,
    beta = parms$beta,
    susceptibility = parms$susceptibility,
    infectiousness = parms$infectiousness
  )
  derivative <- rates_to_derivative(
    transition_rate_table = rates,
    compartments = parms$model$compartments,
    age_structure = parms$age_structure
  )

  infection_derivative <- as.numeric(derivative$derivative)

  if (is.null(parms$demographic_process)) {
    return(list(infection_derivative))
  }

  list(infection_derivative + compartment_demographic_derivative(
    state_vector = as.numeric(state),
    time = time,
    model = parms$model,
    age_structure = parms$age_structure,
    demographic_process = parms$demographic_process,
    time_policy = parms$time_policy
  ))
}

validate_simulation_method <- function(method) {
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

validate_simulation_demography_inputs <- function(
  demographic_process,
  time_policy,
  method,
  model,
  age_structure,
  times
) {
  if (is.null(demographic_process)) {
    return(validate_demographic_time_policy(time_policy))
  }

  time_policy <- validate_demographic_time_policy(time_policy)
  validate_demographic_process(demographic_process)

  if (!model$model_type %in% c("SIR", "SEIR")) {
    stop("demographic_process coupling currently supports only SIR and SEIR models.", call. = FALSE)
  }

  validate_same_age_structure(
    age_structure,
    demographic_process$age_structure,
    "demographic_process"
  )
  validate_demography_schedule_coverage(demographic_process, times, time_policy)

  time_policy
}

validate_simulation_times <- function(times) {
  if (!is.numeric(times) || is.matrix(times) || is.data.frame(times)) {
    stop("times must be a numeric vector.", call. = FALSE)
  }

  if (length(times) < 2) {
    stop("times must contain at least two time points.", call. = FALSE)
  }

  if (anyNA(times)) {
    stop("times cannot contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(times))) {
    stop("times cannot contain non-finite values.", call. = FALSE)
  }

  if (any(diff(times) <= 0)) {
    stop("times must be strictly increasing.", call. = FALSE)
  }

  invisible(times)
}

simulation_state_to_vector <- function(initial_state, age_structure, compartments) {
  if (is.data.frame(initial_state)) {
    return(state_long_to_vector(initial_state, age_structure, compartments))
  }

  if (is.numeric(initial_state) && !is.matrix(initial_state)) {
    validate_state_vector(initial_state, age_structure, compartments)
    return(as.numeric(initial_state))
  }

  stop(
    "initial_state must be a long-form data frame or a numeric vector.",
    call. = FALSE
  )
}

validate_non_negative_simulation_state <- function(state_vector, name) {
  if (any(!is.finite(state_vector))) {
    stop(name, " values cannot contain non-finite values.", call. = FALSE)
  }

  if (any(state_vector < 0)) {
    stop(name, " values cannot be negative.", call. = FALSE)
  }

  invisible(state_vector)
}

validate_non_negative_euler_state <- function(state_vector, time) {
  if (any(state_vector < 0)) {
    first_negative <- which(state_vector < 0)[1]
    stop(
      "Euler step produced negative compartment value at time ",
      time,
      " for state index ",
      first_negative,
      ".",
      call. = FALSE
    )
  }

  invisible(state_vector)
}

simulation_state_output <- function(state_vector, time, age_structure, compartments) {
  state_long <- state_vector_to_long(state_vector, age_structure, compartments)
  data.frame(
    time = rep(time, nrow(state_long)),
    compartment = state_long$compartment,
    age_group = state_long$age_group,
    value = state_long$value,
    stringsAsFactors = FALSE
  )
}
