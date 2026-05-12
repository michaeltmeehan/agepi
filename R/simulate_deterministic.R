#' Simulate a deterministic SIR model
#'
#' Runs a simple deterministic prototype simulation using explicit Euler time
#' steps. The initial state may be supplied either as long-form state data or as
#' a numeric vector in the existing compartment-major, age-group-minor ordering.
#' Names on numeric state vectors are ignored, matching the state-mapping
#' helpers.
#'
#' Currently only `method = "euler"` is supported. The initial state is included
#' in the returned output. Euler updates are intentionally not truncated: if a
#' step would produce negative compartment values, simulation stops with an
#' error. This function is a transparent prototype layer, not a final production
#' ODE backend.
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
#' @param method Simulation method. Currently only `"euler"` is supported.
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
  validate_simulation_method(method)
  validate_simulation_times(times)
  validate_disease_model(model)
  validate_age_structure(age_structure)

  state_vector <- simulation_state_to_vector(
    initial_state,
    age_structure,
    model$compartments
  )
  validate_non_negative_simulation_state(state_vector, "initial_state")

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

validate_simulation_method <- function(method) {
  if (!is.character(method) || length(method) != 1 || anyNA(method) || method == "") {
    stop("method must be a non-missing character scalar.", call. = FALSE)
  }

  if (method != "euler") {
    stop("unsupported simulation method: ", method, call. = FALSE)
  }

  invisible(method)
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
