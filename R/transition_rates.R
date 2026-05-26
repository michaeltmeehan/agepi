#' Compute disease transition rates
#'
#' Computes transition hazards/flows for a disease model at the supplied state.
#' The returned rates are not updated compartment counts. SIR and SEIR disease
#' models and generic `CompartmentModel()` models are supported.
#'
#' `state` may be either a long-form data frame with columns `compartment`,
#' `age_group`, and `value`, or a numeric vector using the existing
#' compartment-major, age-group-minor ordering. Names on numeric state vectors
#' are ignored.
#'
#' @param state Long-form state data frame or numeric state vector.
#' @param model Disease model. `SIRModel()`, `SEIRModel()`, and
#'   `CompartmentModel()` output are supported.
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

  state_long <- transition_state_to_long(state, age_structure, model$compartments)
  validate_non_negative_state_values(state_long)

  if (model$model_type == "CompartmentModel") {
    return(generic_transition_rates(
      state_long = state_long,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness
    ))
  }

  S <- transition_compartment_values(state_long, age_structure, "S")
  I <- transition_compartment_values(state_long, age_structure, "I")
  population <- transition_population_by_age(state_long, age_structure, model$compartments)
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

  if (model$model_type == "SEIR") {
    E <- transition_compartment_values(state_long, age_structure, "E")
    progression_rates <- model$sigma * E
    recovery_rates <- model$gamma * I

    return(data.frame(
      from = rep(model$transitions$from, times = age_structure$n_age_groups),
      to = rep(model$transitions$to, times = age_structure$n_age_groups),
      age_group = rep(age_structure$age_groups, each = nrow(model$transitions)),
      rate = as.numeric(rbind(infection_rates, progression_rates, recovery_rates)),
      stringsAsFactors = FALSE
    ))
  }

  recovery_rates <- model$gamma * I

  data.frame(
    from = rep(model$transitions$from, times = age_structure$n_age_groups),
    to = rep(model$transitions$to, times = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, each = nrow(model$transitions)),
    rate = as.numeric(rbind(infection_rates, recovery_rates)),
    stringsAsFactors = FALSE
  )
}

generic_transition_rates <- function(
  state_long,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness
) {
  population <- transition_population_by_age(state_long, age_structure, model$compartments)
  validate_positive_age_populations(population, age_structure)

  infection_rates <- generic_infection_rates(
    state_long = state_long,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    population = population
  )
  per_capita_rates <- generic_per_capita_transition_rates(
    state_long = state_long,
    model = model,
    age_structure = age_structure
  )

  rates <- rbind(infection_rates, per_capita_rates)
  rates$.transition_order <- seq_len(nrow(rates))
  rates <- rates[order(rates$age_index, rates$.transition_order), ]
  row.names(rates) <- NULL

  rates[, c("from", "to", "age_group", "rate")]
}

generic_infection_rates <- function(
  state_long,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  population
) {
  if (nrow(model$infection_transitions) == 0) {
    return(generic_empty_rate_table())
  }

  infectious <- numeric(age_structure$n_age_groups)
  for (i in seq_along(model$infectious_compartments)) {
    infectious <- infectious +
      model$infectiousness_weights[i] *
      transition_compartment_values(
        state_long,
        age_structure,
        model$infectious_compartments[i]
      )
  }

  lambda <- force_of_infection(
    infectious = infectious,
    population = population,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    age_structure = age_structure
  )

  rows <- vector("list", nrow(model$infection_transitions))
  for (i in seq_len(nrow(model$infection_transitions))) {
    from_values <- transition_compartment_values(
      state_long,
      age_structure,
      model$infection_transitions$from[i]
    )
    rows[[i]] <- data.frame(
      from = model$infection_transitions$from[i],
      to = model$infection_transitions$to[i],
      age_group = age_structure$age_groups,
      rate = as.numeric(lambda) * from_values,
      age_index = seq_len(age_structure$n_age_groups),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

generic_per_capita_transition_rates <- function(state_long, model, age_structure) {
  if (nrow(model$transitions) == 0) {
    return(generic_empty_rate_table())
  }

  rows <- vector("list", nrow(model$transitions))
  for (i in seq_len(nrow(model$transitions))) {
    from_values <- transition_compartment_values(
      state_long,
      age_structure,
      model$transitions$from[i]
    )
    per_capita_rate <- validate_transition_rate_for_age(
      model$transitions$rate[[i]],
      age_structure,
      paste0("transition rate for ", model$transitions$from[i], "->", model$transitions$to[i])
    )
    rows[[i]] <- data.frame(
      from = model$transitions$from[i],
      to = model$transitions$to[i],
      age_group = age_structure$age_groups,
      rate = per_capita_rate * from_values,
      age_index = seq_len(age_structure$n_age_groups),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

generic_empty_rate_table <- function() {
  data.frame(
    from = character(),
    to = character(),
    age_group = character(),
    rate = numeric(),
    age_index = integer(),
    stringsAsFactors = FALSE
  )
}

validate_transition_rate_for_age <- function(rate, age_structure, name) {
  if (!is.numeric(rate) || is.matrix(rate) || is.data.frame(rate) ||
      length(rate) == 0 || anyNA(rate) || any(!is.finite(rate))) {
    stop(name, " must contain finite non-missing numeric value(s).", call. = FALSE)
  }

  if (!length(rate) %in% c(1, age_structure$n_age_groups)) {
    stop(
      name,
      " length must be 1 or match the number of age groups: ",
      age_structure$n_age_groups,
      ".",
      call. = FALSE
    )
  }

  if (any(rate < 0)) {
    stop(name, " cannot contain negative values.", call. = FALSE)
  }

  if (length(rate) == 1) {
    return(rep(as.numeric(rate), age_structure$n_age_groups))
  }

  as.numeric(rate)
}

transition_population_by_age <- function(state_long, age_structure, compartments) {
  population <- numeric(age_structure$n_age_groups)
  for (compartment in compartments) {
    population <- population +
      transition_compartment_values(state_long, age_structure, compartment)
  }

  population
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
