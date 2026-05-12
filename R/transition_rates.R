#' Compute disease transition rates
#'
#' Computes transition hazards/flows for a disease model at the supplied state.
#' The returned rates are not updated compartment counts. Currently only SIR
#' disease models are supported.
#'
#' `state` may be either a long-form data frame with columns `compartment`,
#' `age_group`, and `value`, or a numeric vector using the existing
#' compartment-major, age-group-minor ordering. Names on numeric state vectors
#' are ignored.
#'
#' @param state Long-form state data frame or numeric state vector.
#' @param model Disease model. Currently only `SIRModel()` output is supported.
#' @param age_structure Valid age structure.
#' @param contact_matrix Numeric contact matrix with rows as recipient age
#'   groups and columns as source age groups.
#' @param beta Non-negative finite transmission scaling parameter.
#' @param susceptibility Optional non-negative numeric vector by recipient age
#'   group.
#' @param infectiousness Optional non-negative numeric vector by source age
#'   group.
#'
#' @return Data frame with columns `from`, `to`, `age_group`, and `rate`,
#'   ordered by age group outermost and transition innermost.
#' @export
transition_rates <- function(
  state,
  model,
  age_structure,
  contact_matrix,
  beta = 1,
  susceptibility = NULL,
  infectiousness = NULL
) {
  validate_disease_model(model)
  validate_age_structure(age_structure)

  if (model$model_type != "SIR") {
    stop("transition_rates() currently supports only SIR models.", call. = FALSE)
  }

  state_long <- transition_state_to_long(state, age_structure, model$compartments)
  validate_non_negative_state_values(state_long)

  S <- transition_compartment_values(state_long, age_structure, "S")
  I <- transition_compartment_values(state_long, age_structure, "I")
  R <- transition_compartment_values(state_long, age_structure, "R")
  population <- S + I + R
  validate_positive_age_populations(population, age_structure)

  lambda <- force_of_infection(
    infectious = I,
    population = population,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    age_structure = age_structure
  )

  infection_rates <- as.numeric(lambda) * S
  recovery_rates <- model$gamma * I

  data.frame(
    from = rep(model$transitions$from, times = age_structure$n_age_groups),
    to = rep(model$transitions$to, times = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, each = nrow(model$transitions)),
    rate = as.numeric(rbind(infection_rates, recovery_rates)),
    stringsAsFactors = FALSE
  )
}

transition_state_to_long <- function(state, age_structure, compartments) {
  if (is.data.frame(state)) {
    validate_state_long(state, age_structure, compartments)
    return(state)
  }

  if (is.numeric(state) && !is.matrix(state)) {
    return(state_vector_to_long(state, age_structure, compartments))
  }

  stop(
    "state must be a long-form data frame or a numeric vector.",
    call. = FALSE
  )
}

validate_non_negative_state_values <- function(state_long) {
  if (any(!is.finite(state_long$value))) {
    stop("state values cannot contain non-finite values.", call. = FALSE)
  }

  if (any(state_long$value < 0)) {
    stop("state values cannot be negative.", call. = FALSE)
  }

  invisible(state_long)
}

validate_positive_age_populations <- function(population, age_structure) {
  zero_population <- population == 0
  if (any(zero_population)) {
    stop(
      "state population must be positive in every age group; zero population in age_group: ",
      age_structure$age_groups[which(zero_population)[1]],
      call. = FALSE
    )
  }

  invisible(population)
}

transition_compartment_values <- function(state_long, age_structure, compartment) {
  values <- state_long$value[
    match(
      paste(compartment, age_structure$age_groups),
      paste(state_long$compartment, state_long$age_group)
    )
  ]

  as.numeric(values)
}
