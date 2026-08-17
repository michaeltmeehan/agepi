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
#' @param beta Optional non-negative finite transmission scaling parameter.
#'   When omitted, infectious models use `model$beta` if present and error if
#'   no default is available. Models without infection transitions do not
#'   require beta; an explicit beta is validated and otherwise ignored.
#'
#' @return Data frame with columns `from`, `to`, `age_group`, `rate`, and
#'   `transition_id`, ordered by age group outermost and transition innermost.
#'   `transition_id` identifies the logical transition and is shared by all age
#'   rows for that transition.
#' @export
transition_rates <- function(
  state,
  model,
  age_structure,
  contact_matrix,
  beta = NULL
) {
  validate_disease_model(model)
  validate_age_structure(age_structure)
  if (model_has_infection_process(model)) {
    beta <- resolve_transmission_beta(model, beta)
  } else if (!is.null(beta)) {
    beta <- validate_force_beta(beta)
  }

  state_long <- transition_state_to_long(state, age_structure, model$compartments)
  validate_non_negative_state_values(state_long)

  if (model$model_type == "CompartmentModel") {
    return(generic_transition_rates(
      state_long = state_long,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta
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
      transition_id = rep(
        transition_identifiers(
          from = model$transitions$from,
          to = model$transitions$to,
          transition_type = c("infection", "transition", "transition")
        ),
        times = age_structure$n_age_groups
      ),
      transition_label = rep(NA_character_, nrow(model$transitions) * age_structure$n_age_groups),
      transition_type = rep(c("infection", "transition", "transition"), times = age_structure$n_age_groups),
      stringsAsFactors = FALSE
    ))
  }

  recovery_rates <- model$gamma * I

  data.frame(
    from = rep(model$transitions$from, times = age_structure$n_age_groups),
    to = rep(model$transitions$to, times = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, each = nrow(model$transitions)),
    rate = as.numeric(rbind(infection_rates, recovery_rates)),
    transition_id = rep(
      transition_identifiers(
        from = model$transitions$from,
        to = model$transitions$to,
        transition_type = c("infection", "transition")
      ),
      times = age_structure$n_age_groups
    ),
    transition_label = rep(NA_character_, nrow(model$transitions) * age_structure$n_age_groups),
    transition_type = rep(c("infection", "transition"), times = age_structure$n_age_groups),
    stringsAsFactors = FALSE
  )
}

generic_transition_rates <- function(
  state_long,
  model,
  age_structure,
  contact_matrix,
  beta
) {
  population <- transition_population_by_age(state_long, age_structure, model$compartments)

  infection_rates <- generic_infection_rates(
    state_long = state_long,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
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

  rates[, c("from", "to", "age_group", "rate", "transition_id", "transition_label", "transition_type")]
}

generic_infection_rates <- function(
  state_long,
  model,
  age_structure,
  contact_matrix,
  beta,
  population
) {
  if (nrow(model$infection_transitions) == 0) {
    return(generic_empty_rate_table())
  }

  susceptibility_matrix <- infection_transition_susceptibility_matrix(
    model = model,
    age_structure = age_structure
  )

  infectiousness_matrix <- generic_infectiousness_weight_matrix(
    model = model,
    age_structure = age_structure
  )
  infectious_state_matrix <- generic_infectious_state_matrix(
    state_long = state_long,
    model = model,
    age_structure = age_structure
  )
  infectious <- if (nrow(infectious_state_matrix) == 0) {
    numeric(age_structure$n_age_groups)
  } else {
    colSums(infectious_state_matrix * infectiousness_matrix)
  }

  lambda <- force_of_infection(
    infectious = infectious,
    population = population,
    contact_matrix = contact_matrix,
    beta = beta,
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
      from = rep(model$infection_transitions$from[i], age_structure$n_age_groups),
      to = rep(model$infection_transitions$to[i], age_structure$n_age_groups),
      age_group = age_structure$age_groups,
      rate = as.numeric(lambda) * susceptibility_matrix[i, ] * from_values,
      transition_id = rep(
        transition_identifiers(
        from = model$infection_transitions$from[i],
        to = model$infection_transitions$to[i],
        transition_type = "infection"
        ),
        age_structure$n_age_groups
      ),
      transition_label = rep(model$infection_transitions$transition_label[i], age_structure$n_age_groups),
      transition_type = rep("infection", age_structure$n_age_groups),
      age_index = seq_len(age_structure$n_age_groups),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

generic_infectious_state_matrix <- function(state_long, model, age_structure) {
  if (length(model$infectious_compartments) == 0) {
    return(matrix(
      numeric(),
      nrow = 0,
      ncol = age_structure$n_age_groups,
      dimnames = list(character(), age_structure$age_groups)
    ))
  }

  rows <- lapply(model$infectious_compartments, function(compartment) {
    transition_compartment_values(state_long, age_structure, compartment)
  })
  state_matrix <- do.call(rbind, rows)
  row.names(state_matrix) <- model$infectious_compartments
  colnames(state_matrix) <- age_structure$age_groups
  state_matrix
}

generic_infectiousness_weight_matrix <- function(model, age_structure) {
  validate_age_structure(age_structure)

  if (length(model$infectious_compartments) == 0) {
    return(matrix(
      numeric(),
      nrow = 0,
      ncol = age_structure$n_age_groups,
      dimnames = list(character(), age_structure$age_groups)
    ))
  }

  weights <- validate_generic_infectiousness_weights(
    model$infectiousness_weights,
    model$infectious_compartments
  )

  rows <- lapply(seq_along(weights), function(i) {
    normalize_generic_infectiousness_weight_for_age(
      susceptibility = weights[[i]],
      age_structure = age_structure,
      name = paste0("infectiousness_weights for infectious compartment ", names(weights)[i])
    )
  })

  weight_matrix <- do.call(rbind, rows)
  row.names(weight_matrix) <- names(weights)
  colnames(weight_matrix) <- age_structure$age_groups
  weight_matrix
}

infection_transition_susceptibility_matrix <- function(model, age_structure) {
  validate_age_structure(age_structure)

  if (nrow(model$infection_transitions) == 0) {
    return(matrix(
      numeric(),
      nrow = 0,
      ncol = age_structure$n_age_groups,
      dimnames = list(character(), age_structure$age_groups)
    ))
  }

  if ("susceptibility" %in% names(model$infection_transitions)) {
    susceptibility_specs <- model$infection_transitions$susceptibility
  } else {
    susceptibility_specs <- rep(list(1), nrow(model$infection_transitions))
  }

  rows <- vector("list", nrow(model$infection_transitions))
  for (i in seq_len(nrow(model$infection_transitions))) {
    rows[[i]] <- normalize_generic_infectiousness_weight_for_age(
      susceptibility = susceptibility_specs[[i]],
      age_structure = age_structure,
      name = paste0(
        "susceptibility for infection transition ",
        model$infection_transitions$from[i],
        "->",
        model$infection_transitions$to[i]
      )
    )
  }

  susceptibility_matrix <- do.call(rbind, rows)
  row.names(susceptibility_matrix) <- transition_identifiers(
    from = model$infection_transitions$from,
    to = model$infection_transitions$to,
    transition_type = "infection"
  )
  susceptibility_matrix
}

normalize_generic_infectiousness_weight_for_age <- function(
  susceptibility,
  age_structure,
  name
) {
  if (is.null(susceptibility)) {
    return(rep(1, age_structure$n_age_groups))
  }

  if (is.data.frame(susceptibility) || is.matrix(susceptibility) || is.list(susceptibility)) {
    stop(name, " must be numeric.", call. = FALSE)
  }

  if (!is.numeric(susceptibility) || length(susceptibility) == 0 ||
      anyNA(susceptibility) || any(!is.finite(susceptibility))) {
    stop(name, " must be finite numeric value(s).", call. = FALSE)
  }

  if (any(susceptibility < 0)) {
    stop(name, " cannot contain negative values.", call. = FALSE)
  }

  if (length(susceptibility) == 1) {
    return(rep(as.numeric(susceptibility), age_structure$n_age_groups))
  }

  if (length(susceptibility) != age_structure$n_age_groups) {
    stop(
      name,
      " length must be 1 or match the number of age groups: ",
      age_structure$n_age_groups,
      ".",
      call. = FALSE
    )
  }

  susceptibility_names <- names(susceptibility)
  if (is.null(susceptibility_names)) {
    return(as.numeric(susceptibility))
  }

  if (anyNA(susceptibility_names) || any(susceptibility_names == "")) {
    stop(name, " age-specific vector names must be non-empty.", call. = FALSE)
  }

  duplicated_names <- unique(susceptibility_names[duplicated(susceptibility_names)])
  if (length(duplicated_names) > 0) {
    stop(
      name,
      " age-specific vector names must be unique; duplicate age group(s): ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_names <- setdiff(susceptibility_names, age_structure$age_groups)
  if (length(unknown_names) > 0) {
    stop(
      name,
      " age-specific vector contains unknown age_group value(s): ",
      paste(unknown_names, collapse = ", "),
      call. = FALSE
    )
  }

  missing_names <- setdiff(age_structure$age_groups, susceptibility_names)
  if (length(missing_names) > 0) {
    stop(
      name,
      " age-specific vector is missing age_group value(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  as.numeric(susceptibility[age_structure$age_groups])
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
      paste0("transition rate for ", transition_row_label(model$transitions, i))
    )
    transition_id <- if ("transition_id" %in% names(model$transitions)) {
      model$transitions$transition_id[i]
    } else {
      transition_row_identifier(model$transitions, i)
    }
    transition_type <- if ("transition_type" %in% names(model$transitions)) {
      model$transitions$transition_type[i]
    } else ifelse(is.na(model$transitions$to[i]), "outflow", "transition")
    transition_type[transition_type == "internal"] <- "transition"
    rows[[i]] <- data.frame(
      from = model$transitions$from[i],
      to = model$transitions$to[i],
      age_group = age_structure$age_groups,
      rate = per_capita_rate * from_values,
      transition_id = transition_id,
      transition_label = if ("transition_label" %in% names(model$transitions)) model$transitions$transition_label[i] else NA_character_,
      transition_type = transition_type,
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
    transition_id = character(),
    transition_label = character(),
    transition_type = character(),
    age_index = integer(),
    stringsAsFactors = FALSE
  )
}

transition_identifiers <- function(from, to, transition_type, id = NULL, label = NULL) {
  transition_type <- as.character(transition_type)
  from <- as.character(from)
  if (!is.null(id) && !is.character(id)) {
    id <- as.character(id)
  }
  if (!is.null(label) && !is.character(label)) {
    label <- as.character(label)
  }
  out <- character(length(from))

  outflow <- transition_type == "outflow"
  if (any(outflow)) {
    fallback <- from[outflow]
    if (!is.null(label)) {
      fallback <- ifelse(!is.na(label[outflow]), label[outflow], fallback)
    }
    if (!is.null(id)) {
      fallback <- ifelse(!is.na(id[outflow]), id[outflow], fallback)
    }
    out[outflow] <- paste0(
      "outflow:",
      fallback
    )
  }

  if (any(!outflow)) {
    out[!outflow] <- paste0(
      transition_type[!outflow],
      ":",
      from[!outflow],
      "->",
      to[!outflow]
    )
  }

  out
}

transition_row_label <- function(transitions, i) {
  if ("transition_label" %in% names(transitions) && !is.na(transitions$transition_label[i]) && transitions$transition_label[i] != "") {
    return(transitions$transition_label[i])
  }

  if (!"transition_type" %in% names(transitions)) {
    if (is.na(transitions$to[i])) {
      return(paste0(transitions$from[i], "->outside"))
    }
    return(paste0(transitions$from[i], "->", transitions$to[i]))
  }

  if (transitions$transition_type[i] == "outflow") {
    return(paste0(transitions$from[i], "->outside"))
  }

  paste0(transitions$from[i], "->", transitions$to[i])
}

transition_row_identifier <- function(transitions, i) {
  if ("transition_id" %in% names(transitions) && !is.na(transitions$transition_id[i])) {
    return(transitions$transition_id[i])
  }

  if (!"transition_type" %in% names(transitions) && is.na(transitions$to[i])) {
    return(transition_identifiers(
      from = transitions$from[i],
      to = transitions$to[i],
      transition_type = "outflow",
      label = if ("transition_label" %in% names(transitions)) transitions$transition_label[i] else NULL
    ))
  }

  if ("transition_type" %in% names(transitions) && transitions$transition_type[i] == "outflow") {
    return(transition_identifiers(
      from = transitions$from[i],
      to = transitions$to[i],
      transition_type = "outflow",
      label = if ("transition_label" %in% names(transitions)) transitions$transition_label[i] else NULL
    ))
  }

  transition_identifiers(
    from = transitions$from[i],
    to = transitions$to[i],
    transition_type = "transition"
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

  rate_names <- names(rate)
  if (is.null(rate_names) || anyNA(rate_names) || any(rate_names == "")) {
    stop(name, " age-specific vector must be named by age group.", call. = FALSE)
  }

  duplicated_names <- unique(rate_names[duplicated(rate_names)])
  if (length(duplicated_names) > 0) {
    stop(
      name,
      " age-specific vector names must be unique; duplicate age group(s): ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_names <- setdiff(rate_names, age_structure$age_groups)
  if (length(unknown_names) > 0) {
    stop(
      name,
      " age-specific vector contains unknown age_group value(s): ",
      paste(unknown_names, collapse = ", "),
      call. = FALSE
    )
  }

  missing_names <- setdiff(age_structure$age_groups, rate_names)
  if (length(missing_names) > 0) {
    stop(
      name,
      " age-specific vector is missing age_group value(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  as.numeric(rate[age_structure$age_groups])
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
