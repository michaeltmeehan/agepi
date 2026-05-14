#' Simulate a deterministic SIR model
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
#' The deSolve backend currently supports the same static deterministic SIR
#' model as the Euler backend: static contact matrix, static `beta`, static
#' susceptibility, and static infectiousness. It does not add time-varying
#' demography, time-varying contacts, interpolation, or demographic dynamics.
#'
#' Euler updates are intentionally not truncated: if a step would produce
#' negative compartment values, simulation stops with an error.
#'
#' @param initial_state Long-form state data frame with columns `compartment`,
#'   `age_group`, and `value`, or numeric state vector.
#' @param times Numeric vector of finite, non-missing, strictly increasing time
#'   points. Must have length at least two.
#' @param model Disease model. Currently only `SIRModel()` output is supported.
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
  method = "euler"
) {
  method <- validate_simulation_method(method)
  validate_simulation_times(times)
  validate_disease_model(model)
  validate_age_structure(age_structure)

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
      infectiousness = infectiousness
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
    infectiousness = infectiousness
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
  infectiousness
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
    rates <- transition_rates(
      state = current_state,
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

    next_state <- as.numeric(current_state) + dt * derivative$derivative
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

simulate_deterministic_desolve <- function(
  state_vector,
  times,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness
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
      infectiousness = infectiousness
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

  list(as.numeric(derivative$derivative))
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
