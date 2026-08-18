prepare_transition_rate_context <- function(model, age_structure, contact_matrix, beta = NULL) {
  validate_disease_model(model)
  validate_age_structure(age_structure)
  prepare_transition_rate_context_validated(
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta
  )
}

prepare_transition_rate_context_validated <- function(model, age_structure, contact_matrix, beta = NULL) {
  validate_contact_matrix(contact_matrix, age_structure)
  beta <- resolve_transmission_beta_validated(model, beta)

  context <- list(
    model = model,
    age_structure = age_structure,
    age_groups = age_structure$age_groups,
    n_age_groups = age_structure$n_age_groups,
    compartments = model$compartments,
    contact_matrix = contact_matrix,
    beta = beta,
    has_infection_process = model_has_infection_process(model),
    state_order = state_order(age_structure, model$compartments),
    stochastic_transition_order = stochastic_transition_order(model),
    specialized_transition_ids = if (nrow(model$transitions) > 0) {
      if (identical(model$model_type, "CompartmentModel")) {
        transition_order_from_compartment_transitions(model$transitions)$transition_id
      } else {
        transition_order_from_specialized_transitions(model$transitions)$transition_id
      }
    } else {
      character()
    }
  )

  if (identical(model$model_type, "CompartmentModel")) {
    context$infection_susceptibility_matrix <- infection_transition_susceptibility_matrix_validated(
      model = model,
      age_groups = context$age_groups,
      n_age_groups = context$n_age_groups
    )
    context$infectiousness_matrix <- generic_infectiousness_weight_matrix_validated(
      model = model,
      age_groups = context$age_groups,
      n_age_groups = context$n_age_groups
    )
    context$transition_rate_vectors <- generic_transition_rate_vectors_validated(
      model = model,
      age_groups = context$age_groups,
      n_age_groups = context$n_age_groups
    )
  }

  context
}

resolve_transmission_beta_validated <- function(model, beta = NULL) {
  if (!model_has_infection_process(model)) {
    if (is.null(beta)) {
      return(NULL)
    }
    return(validate_force_beta(beta))
  }

  if (!is.null(beta)) {
    return(validate_force_beta(beta))
  }

  if ("beta" %in% names(model) && !is.null(model$beta)) {
    return(validate_force_beta(model$beta))
  }

  stop(
    "beta is required when the model contains infection transitions. ",
    "Store a default in model$beta or supply beta at simulation/rate-evaluation time.",
    call. = FALSE
  )
}

transition_rates_from_state_long <- function(state_long, context) {
  validate_non_negative_state_values(state_long)
  transition_rates_from_state_long_unchecked(state_long, context)
}

transition_rates_from_state_vector <- function(state_vector, context) {
  state_long <- state_vector_to_long_unchecked(
    state_vector = state_vector,
    age_structure = context$age_structure,
    compartments = context$compartments
  )
  transition_rates_from_state_long_unchecked(state_long, context)
}

transition_rates_from_state_long_unchecked <- function(state_long, context) {
  model <- context$model
  age_structure <- context$age_structure

  if (identical(model$model_type, "CompartmentModel")) {
    return(generic_transition_rates_from_state_long(state_long, context))
  }

  S <- transition_compartment_values(state_long, age_structure, "S")
  I <- transition_compartment_values(state_long, age_structure, "I")
  population <- transition_population_by_age(state_long, age_structure, model$compartments)
  validate_positive_age_populations(population, age_structure)

  lambda <- as.numeric(context$beta * (context$contact_matrix %*% (I / population)))
  infection_rates <- lambda * S

  if (identical(model$model_type, "SEIR")) {
    E <- transition_compartment_values(state_long, age_structure, "E")
    progression_rates <- model$sigma * E
    recovery_rates <- model$gamma * I

    return(data.frame(
      from = rep(model$transitions$from, times = age_structure$n_age_groups),
      to = rep(model$transitions$to, times = age_structure$n_age_groups),
      age_group = rep(age_structure$age_groups, each = nrow(model$transitions)),
      rate = as.numeric(rbind(infection_rates, progression_rates, recovery_rates)),
      transition_id = rep(
        context$specialized_transition_ids,
        times = age_structure$n_age_groups
      ),
      transition_label = rep(
        if ("transition_label" %in% names(model$transitions)) model$transitions$transition_label else NA_character_,
        times = age_structure$n_age_groups
      ),
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
      context$specialized_transition_ids,
      times = age_structure$n_age_groups
    ),
    transition_label = rep(
      if ("transition_label" %in% names(model$transitions)) model$transitions$transition_label else NA_character_,
      times = age_structure$n_age_groups
    ),
    transition_type = rep(c("infection", "transition"), times = age_structure$n_age_groups),
    stringsAsFactors = FALSE
  )
}

generic_transition_rates_from_state_long <- function(state_long, context) {
  model <- context$model
  age_structure <- context$age_structure
  population <- transition_population_by_age(state_long, age_structure, model$compartments)

  infection_rates <- generic_infection_rates_from_state_long(state_long, context, population)
  per_capita_rates <- generic_per_capita_transition_rates_from_state_long(state_long, context)

  rates <- rbind(infection_rates, per_capita_rates)
  rates$.transition_order <- seq_len(nrow(rates))
  rates <- rates[order(rates$age_index, rates$.transition_order), ]
  row.names(rates) <- NULL

  rates[, c("from", "to", "age_group", "rate", "transition_id", "transition_label", "transition_type")]
}

generic_infection_rates_from_state_long <- function(state_long, context, population) {
  model <- context$model
  age_structure <- context$age_structure

  if (nrow(model$infection_transitions) == 0) {
    return(generic_empty_rate_table())
  }

  infectious_state_matrix <- generic_infectious_state_matrix(
    state_long = state_long,
    model = model,
    age_structure = age_structure
  )
  infectious <- if (nrow(infectious_state_matrix) == 0) {
    numeric(context$n_age_groups)
  } else {
    colSums(infectious_state_matrix * context$infectiousness_matrix)
  }

  validate_positive_age_populations(population, age_structure)
  lambda <- as.numeric(context$beta * (context$contact_matrix %*% (infectious / population)))

  rows <- vector("list", nrow(model$infection_transitions))
  for (i in seq_len(nrow(model$infection_transitions))) {
    from_values <- transition_compartment_values(
      state_long,
      age_structure,
      model$infection_transitions$from[i]
    )
    rows[[i]] <- data.frame(
      from = rep(model$infection_transitions$from[i], context$n_age_groups),
      to = rep(model$infection_transitions$to[i], context$n_age_groups),
      age_group = context$age_groups,
      rate = lambda * context$infection_susceptibility_matrix[i, ] * from_values,
      transition_id = rep(
        row.names(context$infection_susceptibility_matrix)[i],
        context$n_age_groups
      ),
      transition_label = rep(model$infection_transitions$transition_label[i], context$n_age_groups),
      transition_type = rep("infection", context$n_age_groups),
      age_index = seq_len(context$n_age_groups),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

generic_per_capita_transition_rates_from_state_long <- function(state_long, context) {
  model <- context$model
  age_structure <- context$age_structure

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
    rows[[i]] <- data.frame(
      from = model$transitions$from[i],
      to = model$transitions$to[i],
      age_group = context$age_groups,
      rate = context$transition_rate_vectors[[i]] * from_values,
      transition_id = model$transitions$transition_id[i],
      transition_label = if ("transition_label" %in% names(model$transitions)) model$transitions$transition_label[i] else NA_character_,
      transition_type = if ("transition_type" %in% names(model$transitions)) {
        transition_type <- model$transitions$transition_type[i]
        ifelse(transition_type == "internal", "transition", transition_type)
      } else {
        ifelse(is.na(model$transitions$to[i]), "outflow", "transition")
      },
      age_index = seq_len(context$n_age_groups),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

generic_transition_rate_vectors_validated <- function(model, age_groups, n_age_groups) {
  if (nrow(model$transitions) == 0) {
    return(list())
  }

  fake_age_structure <- list(
    age_groups = age_groups,
    n_age_groups = n_age_groups
  )

  lapply(seq_len(nrow(model$transitions)), function(i) {
    validate_transition_rate_for_age(
      model$transitions$rate[[i]],
      fake_age_structure,
      paste0("transition rate for ", transition_row_label(model$transitions, i))
    )
  })
}

generic_infectiousness_weight_matrix_validated <- function(model, age_groups, n_age_groups) {
  if (length(model$infectious_compartments) == 0) {
    return(matrix(
      numeric(),
      nrow = 0,
      ncol = n_age_groups,
      dimnames = list(character(), age_groups)
    ))
  }

  fake_age_structure <- list(
    age_groups = age_groups,
    n_age_groups = n_age_groups
  )

  weights <- validate_generic_infectiousness_weights(
    model$infectiousness_weights,
    model$infectious_compartments
  )

  rows <- lapply(seq_along(weights), function(i) {
    normalize_generic_infectiousness_weight_for_age(
      susceptibility = weights[[i]],
      age_structure = fake_age_structure,
      name = paste0("infectiousness_weights for infectious compartment ", names(weights)[i])
    )
  })

  weight_matrix <- do.call(rbind, rows)
  row.names(weight_matrix) <- names(weights)
  colnames(weight_matrix) <- age_groups
  weight_matrix
}

infection_transition_susceptibility_matrix_validated <- function(model, age_groups, n_age_groups) {
  if (nrow(model$infection_transitions) == 0) {
    return(matrix(
      numeric(),
      nrow = 0,
      ncol = n_age_groups,
      dimnames = list(character(), age_groups)
    ))
  }

  fake_age_structure <- list(
    age_groups = age_groups,
    n_age_groups = n_age_groups
  )

  if ("susceptibility" %in% names(model$infection_transitions)) {
    susceptibility_specs <- model$infection_transitions$susceptibility
  } else {
    susceptibility_specs <- rep(list(1), nrow(model$infection_transitions))
  }

  rows <- vector("list", nrow(model$infection_transitions))
  for (i in seq_len(nrow(model$infection_transitions))) {
    rows[[i]] <- normalize_generic_infectiousness_weight_for_age(
      susceptibility = susceptibility_specs[[i]],
      age_structure = fake_age_structure,
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

rates_to_derivative_from_rates <- function(transition_rate_table, context) {
  derivative <- context$state_order
  derivative$derivative <- 0

  for (i in seq_len(nrow(transition_rate_table))) {
    from_index <- match_derivative_cell(
      derivative,
      transition_rate_table$from[i],
      transition_rate_table$age_group[i]
    )

    derivative$derivative[from_index] <- derivative$derivative[from_index] - transition_rate_table$rate[i]
    if (!is.na(transition_rate_table$to[i])) {
      to_index <- match_derivative_cell(
        derivative,
        transition_rate_table$to[i],
        transition_rate_table$age_group[i]
      )
      derivative$derivative[to_index] <- derivative$derivative[to_index] + transition_rate_table$rate[i]
    }
  }

  derivative[, c("compartment", "age_group", "derivative")]
}

stochastic_event_table_from_state_vector <- function(state_vector, context) {
  rates <- transition_rates_from_state_vector(state_vector, context)
  transition_order <- context$stochastic_transition_order
  rates$.transition_index <- match(rates$transition_id, transition_order$transition_id)
  rates$age_index <- match(rates$age_group, context$age_groups)
  if (!"transition_type" %in% names(rates)) {
    rates$transition_type <- transition_type_from_ids(rates$transition_id)
  }
  rates$event <- stochastic_event_labels(rates, context$model)
  rates <- rates[order(rates$.transition_index, rates$age_index), ]
  row.names(rates) <- NULL

  rates[, c("event", "transition_type", "transition_id", "age_group", "age_index", "from", "to", "rate")]
}
