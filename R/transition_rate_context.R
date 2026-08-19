prepare_transition_rate_context <- function(
  model,
  age_structure,
  contact_matrix,
  beta = NULL,
  include_public_template = TRUE,
  include_transition_id = FALSE
) {
  validate_disease_model(model)
  validate_age_structure(age_structure)
  prepare_transition_rate_context_validated(
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    include_public_template = include_public_template,
    include_transition_id = include_transition_id
  )
}

prepare_transition_rate_context_validated <- function(
  model,
  age_structure,
  contact_matrix,
  beta = NULL,
  include_public_template = TRUE,
  include_transition_id = FALSE
) {
  validate_contact_matrix(contact_matrix, age_structure)
  beta <- resolve_transmission_beta_validated(model, beta)

  state_order <- state_order(age_structure, model$compartments)

  context <- list(
    model = model,
    age_structure = age_structure,
    age_groups = age_structure$age_groups,
    n_age_groups = age_structure$n_age_groups,
    compartments = model$compartments,
    compartment_index = stats::setNames(seq_along(model$compartments), model$compartments),
    contact_matrix = contact_matrix,
    beta = beta,
    has_infection_process = model_has_infection_process(model),
    state_order = state_order,
    state_output_template = state_order[, c("compartment", "age_group"), drop = FALSE],
    state_order_key = paste(
      rep(model$compartments, each = age_structure$n_age_groups),
      rep(age_structure$age_groups, times = length(model$compartments))
    ),
    transition_rate_structure = prepare_transition_rate_structure(
      model,
      age_structure,
      include_transition_id = include_transition_id
    )
  )
  context$derivative_index_lookup <- prepare_derivative_index_lookup(
    transition_rate_table = context$transition_rate_structure,
    state_order_key = context$state_order_key
  )
  if (include_public_template) {
    context$stochastic_transition_order <- stochastic_transition_order(model)
    context$specialized_transition_ids <- if (nrow(model$transitions) > 0) {
      if (identical(model$model_type, "CompartmentModel")) {
        transition_order_from_compartment_transitions(model$transitions)$transition_id
      } else {
        transition_order_from_specialized_transitions(model$transitions)$transition_id
      }
    } else {
      character()
    }
    context$transition_rate_template <- prepare_transition_rate_template(model, age_structure)
    context$stochastic_event_order <- order(
      match(
        context$transition_rate_template$transition_id,
        context$stochastic_transition_order$transition_id
      ),
      context$transition_rate_template$age_index
    )
    context$stochastic_event_template <- context$transition_rate_template[context$stochastic_event_order, c(
      "event",
      "transition_type",
      "transition_id",
      "age_group",
      "age_index",
      "from",
      "to"
    ), drop = FALSE]
  }

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
    context$generic_infectious_compartment_indices <- if (length(model$infectious_compartments) > 0) {
      unname(context$compartment_index[model$infectious_compartments])
    } else {
      integer()
    }
    context$generic_infectiousness_matrix_by_age <- t(context$infectiousness_matrix)
    context$generic_infection_from_indices <- if (nrow(model$infection_transitions) > 0) {
      unname(context$compartment_index[model$infection_transitions$from])
    } else {
      integer()
    }
    context$generic_transition_from_indices <- if (nrow(model$transitions) > 0) {
      unname(context$compartment_index[model$transitions$from])
    } else {
      integer()
    }
    context$generic_transition_rate_matrix <- if (nrow(model$transitions) > 0) {
      do.call(rbind, context$transition_rate_vectors)
    } else {
      matrix(numeric(), nrow = 0, ncol = context$n_age_groups)
    }
  }

  context
}

prepare_transition_rate_structure <- function(model, age_structure, include_transition_id = FALSE) {
  if (identical(model$model_type, "CompartmentModel")) {
    return(prepare_generic_transition_rate_structure(
      model,
      age_structure,
      include_transition_id = include_transition_id
    ))
  }

  prepare_specialized_transition_rate_structure(
    model,
    age_structure,
    include_transition_id = include_transition_id
  )
}

prepare_specialized_transition_rate_structure <- function(
  model,
  age_structure,
  include_transition_id = FALSE
) {
  n_transitions <- nrow(model$transitions)
  if (n_transitions == 0) {
    return(empty_transition_rate_structure(include_transition_id = include_transition_id))
  }

  structure <- data.frame(
    from = rep(model$transitions$from, times = age_structure$n_age_groups),
    to = rep(model$transitions$to, times = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, each = n_transitions),
    stringsAsFactors = FALSE
  )
  if (include_transition_id) {
    structure$transition_id <- rep(
      transition_order_from_specialized_transitions(model$transitions)$transition_id,
      times = age_structure$n_age_groups
    )
    structure <- structure[, c("from", "to", "age_group", "transition_id"), drop = FALSE]
  }
  structure
}

prepare_generic_transition_rate_structure <- function(
  model,
  age_structure,
  include_transition_id = FALSE
) {
  if (nrow(model$infection_transitions) == 0 && nrow(model$transitions) == 0) {
    return(empty_transition_rate_structure(include_transition_id = include_transition_id))
  }

  infection_transition_ids <- if (nrow(model$infection_transitions) > 0) {
    transition_order_from_infection_transitions(model$infection_transitions)$transition_id
  } else {
    character()
  }

  from_values <- c(
    if (length(infection_transition_ids) > 0) model$infection_transitions$from else character(),
    if (nrow(model$transitions) > 0) model$transitions$from else character()
  )
  to_values <- c(
    if (length(infection_transition_ids) > 0) model$infection_transitions$to else character(),
    if (nrow(model$transitions) > 0) model$transitions$to else character()
  )
  structure <- data.frame(
    from = rep(from_values, times = age_structure$n_age_groups),
    to = rep(to_values, times = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, each = length(from_values)),
    stringsAsFactors = FALSE
  )
  if (include_transition_id) {
    transition_ids <- c(
      if (nrow(model$infection_transitions) > 0) {
        transition_order_from_infection_transitions(model$infection_transitions)$transition_id
      } else {
        character()
      },
      if (nrow(model$transitions) > 0) {
        if ("transition_id" %in% names(model$transitions)) {
          model$transitions$transition_id
        } else {
          transition_order_from_compartment_transitions(model$transitions)$transition_id
        }
      } else {
        character()
      }
    )
    structure$transition_id <- rep(transition_ids, times = age_structure$n_age_groups)
    structure <- structure[, c("from", "to", "age_group", "transition_id"), drop = FALSE]
  }
  structure
}

empty_transition_rate_structure <- function(include_transition_id = FALSE) {
  if (include_transition_id) {
    return(data.frame(
      from = character(),
      to = character(),
      age_group = character(),
      transition_id = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    from = character(),
    to = character(),
    age_group = character(),
    stringsAsFactors = FALSE
  )
}

prepare_transition_rate_template <- function(model, age_structure) {
  if (identical(model$model_type, "CompartmentModel")) {
    return(prepare_generic_transition_rate_template(model, age_structure))
  }

  prepare_specialized_transition_rate_template(model, age_structure)
}

prepare_specialized_transition_rate_template <- function(model, age_structure) {
  n_transitions <- nrow(model$transitions)
  if (n_transitions == 0) {
    return(empty_transition_rate_template())
  }

  transition_ids <- transition_order_from_specialized_transitions(model$transitions)$transition_id
  transition_labels <- if ("transition_label" %in% names(model$transitions)) {
    model$transitions$transition_label
  } else {
    rep(NA_character_, n_transitions)
  }
  transition_types <- specialized_transition_types(model)
  event_labels <- specialized_transition_event_labels(model)

  rows <- vector("list", age_structure$n_age_groups)
  for (age_index in seq_len(age_structure$n_age_groups)) {
    rows[[age_index]] <- data.frame(
      from = model$transitions$from,
      to = model$transitions$to,
      age_group = age_structure$age_groups[age_index],
      transition_id = transition_ids,
      transition_label = transition_labels,
      transition_type = transition_types,
      age_index = age_index,
      event = event_labels,
      stringsAsFactors = FALSE
    )
  }

  template <- do.call(rbind, rows)
  row.names(template) <- NULL
  template
}

prepare_generic_transition_rate_template <- function(model, age_structure) {
  if (nrow(model$infection_transitions) == 0 && nrow(model$transitions) == 0) {
    return(empty_transition_rate_template())
  }

  infection_transition_ids <- if (nrow(model$infection_transitions) > 0) {
    transition_order_from_infection_transitions(model$infection_transitions)$transition_id
  } else {
    character()
  }
  infection_transition_labels <- if (nrow(model$infection_transitions) > 0) {
    if ("transition_label" %in% names(model$infection_transitions)) {
      model$infection_transitions$transition_label
    } else {
      rep(NA_character_, nrow(model$infection_transitions))
    }
  } else {
    character()
  }
  infection_transition_types <- rep("infection", length(infection_transition_ids))
  infection_event_labels <- rep("infection", length(infection_transition_ids))

  per_capita_transition_ids <- if (nrow(model$transitions) > 0) {
    if ("transition_id" %in% names(model$transitions)) {
      model$transitions$transition_id
    } else {
      transition_order_from_compartment_transitions(model$transitions)$transition_id
    }
  } else {
    character()
  }
  per_capita_transition_labels <- if (nrow(model$transitions) > 0) {
    if ("transition_label" %in% names(model$transitions)) {
      model$transitions$transition_label
    } else {
      rep(NA_character_, nrow(model$transitions))
    }
  } else {
    character()
  }
  per_capita_transition_types <- if (nrow(model$transitions) > 0) {
    if ("transition_type" %in% names(model$transitions)) {
      transition_type <- model$transitions$transition_type
      ifelse(transition_type == "internal", "transition", transition_type)
    } else {
      ifelse(is.na(model$transitions$to), "outflow", "transition")
    }
  } else {
    character()
  }
  per_capita_event_labels <- if (nrow(model$transitions) > 0) {
    stochastic_event_labels_from_metadata(
      model = model,
      transition_id = per_capita_transition_ids,
      from = model$transitions$from,
      to = model$transitions$to,
      transition_type = per_capita_transition_types
    )
  } else {
    character()
  }

  rows <- vector("list", age_structure$n_age_groups)
  for (age_index in seq_len(age_structure$n_age_groups)) {
    rows[[age_index]] <- data.frame(
      from = c(
        if (length(infection_transition_ids) > 0) model$infection_transitions$from else character(),
        if (nrow(model$transitions) > 0) model$transitions$from else character()
      ),
      to = c(
        if (length(infection_transition_ids) > 0) model$infection_transitions$to else character(),
        if (nrow(model$transitions) > 0) model$transitions$to else character()
      ),
      age_group = age_structure$age_groups[age_index],
      transition_id = c(infection_transition_ids, per_capita_transition_ids),
      transition_label = c(infection_transition_labels, per_capita_transition_labels),
      transition_type = c(infection_transition_types, per_capita_transition_types),
      age_index = age_index,
      event = c(infection_event_labels, per_capita_event_labels),
      stringsAsFactors = FALSE
    )
  }

  template <- do.call(rbind, rows)
  row.names(template) <- NULL
  template
}

empty_transition_rate_template <- function() {
  data.frame(
    from = character(),
    to = character(),
    age_group = character(),
    transition_id = character(),
    transition_label = character(),
    transition_type = character(),
    age_index = integer(),
    event = character(),
    stringsAsFactors = FALSE
  )
}

specialized_transition_types <- function(model) {
  if (identical(model$model_type, "SIR")) {
    return(c("infection", "transition"))
  }

  c("infection", "transition", "transition")
}

specialized_transition_event_labels <- function(model) {
  if (identical(model$model_type, "SIR")) {
    return(c("infection", "recovery"))
  }

  if (identical(model$model_type, "SEIR")) {
    return(c("infection", "progression", "recovery"))
  }

  rep("transition", nrow(model$transitions))
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
  state_matrix <- transition_state_matrix_from_vector(state_vector, context)
  if (!is.null(context$transition_rate_template)) {
    return(transition_rates_from_state_matrix(state_matrix, context))
  }

  transition_rates_from_state_structure(state_matrix, context)
}

transition_rates_from_state_long_unchecked <- function(state_long, context) {
  state_matrix <- transition_state_matrix_from_long(state_long, context)
  transition_rates_from_state_matrix(state_matrix, context)
}

transition_rates_from_state_matrix <- function(state_matrix, context) {
  model <- context$model

  rate_values <- if (identical(model$model_type, "CompartmentModel")) {
    generic_transition_rate_values_from_state_matrix(state_matrix, context)
  } else {
    specialized_transition_rate_values_from_state_matrix(state_matrix, context)
  }

  rates <- context$transition_rate_template
  rates$rate <- rate_values
  rates[, c("from", "to", "age_group", "rate", "transition_id", "transition_label", "transition_type")]
}

generic_transition_rates_from_state_long <- function(state_long, context) {
  state_matrix <- transition_state_matrix_from_long(state_long, context)
  generic_transition_rates_from_state_matrix(state_matrix, context)
}

generic_transition_rates_from_state_matrix <- function(state_matrix, context) {
  rate_values <- generic_transition_rate_values_from_state_matrix(state_matrix, context)
  rates <- context$transition_rate_template
  rates$rate <- rate_values
  rates[, c("from", "to", "age_group", "rate", "transition_id", "transition_label", "transition_type")]
}

transition_rates_from_state_structure <- function(state_matrix, context) {
  model <- context$model

  rate_values <- if (identical(model$model_type, "CompartmentModel")) {
    generic_transition_rate_values_from_state_matrix(state_matrix, context)
  } else {
    specialized_transition_rate_values_from_state_matrix(state_matrix, context)
  }

  rates <- context$transition_rate_structure
  rates$rate <- rate_values
  rates
}

transition_state_matrix_from_long <- function(state_long, context) {
  ordered_values <- state_long$value[match(
    context$state_order_key,
    paste(state_long$compartment, state_long$age_group)
  )]

  matrix(
    ordered_values,
    nrow = context$n_age_groups,
    ncol = length(context$compartments),
    dimnames = list(context$age_groups, context$compartments)
  )
}

transition_state_matrix_from_vector <- function(state_vector, context) {
  matrix(
    as.numeric(state_vector),
    nrow = context$n_age_groups,
    ncol = length(context$compartments),
    dimnames = list(context$age_groups, context$compartments)
  )
}

transition_compartment_values_from_state_matrix <- function(state_matrix, compartment, context) {
  state_matrix[, context$compartment_index[[compartment]]]
}

transition_population_by_age_from_state_matrix <- function(state_matrix) {
  rowSums(state_matrix)
}

specialized_transition_rate_values_from_state_matrix <- function(state_matrix, context) {
  model <- context$model
  age_structure <- context$age_structure

  S <- transition_compartment_values_from_state_matrix(state_matrix, "S", context)
  I <- transition_compartment_values_from_state_matrix(state_matrix, "I", context)
  population <- transition_population_by_age_from_state_matrix(state_matrix)
  validate_positive_age_populations(population, age_structure)

  lambda <- as.numeric(context$beta * (context$contact_matrix %*% (I / population)))
  infection_rates <- lambda * S

  if (identical(model$model_type, "SEIR")) {
    E <- transition_compartment_values_from_state_matrix(state_matrix, "E", context)
    progression_rates <- model$sigma * E
    recovery_rates <- model$gamma * I
    n_age_groups <- context$n_age_groups
    rates <- numeric(3L * n_age_groups)
    infection_index <- seq.int(from = 1L, by = 3L, length.out = n_age_groups)
    rates[infection_index] <- infection_rates
    rates[infection_index + 1L] <- progression_rates
    rates[infection_index + 2L] <- recovery_rates
    return(rates)
  }

  recovery_rates <- model$gamma * I
  n_age_groups <- context$n_age_groups
  rates <- numeric(2L * n_age_groups)
  infection_index <- seq.int(from = 1L, by = 2L, length.out = n_age_groups)
  rates[infection_index] <- infection_rates
  rates[infection_index + 1L] <- recovery_rates
  rates
}

stochastic_event_labels_from_metadata <- function(
  model,
  transition_id,
  from,
  to,
  transition_type
) {
  labels <- ifelse(is.na(transition_type), transition_type_from_ids(transition_id), transition_type)

  if (identical(model$model_type, "SIR")) {
    labels[transition_id == "infection:S->I"] <- "infection"
    labels[transition_id == "transition:I->R"] <- "recovery"
    return(labels)
  }

  if (identical(model$model_type, "SEIR")) {
    labels[transition_id == "infection:S->E"] <- "infection"
    labels[transition_id == "transition:E->I"] <- "progression"
    labels[transition_id == "transition:I->R"] <- "recovery"
    return(labels)
  }

  if (identical(model$model_type, "CompartmentModel")) {
    labels[transition_type == "outflow"] <- "outflow"
    labels[transition_type == "infection"] <- "infection"
    labels[transition_type == "transition" & from == "E" & to == "I"] <- "progression"
    labels[transition_type == "transition" & to == "R"] <- "recovery"
    return(labels)
  }

  labels
}

generic_transition_rate_values_from_state_matrix <- function(state_matrix, context) {
  infection_rates <- generic_infection_rate_matrix_from_state_matrix(state_matrix, context)
  per_capita_rates <- generic_per_capita_rate_matrix_from_state_matrix(state_matrix, context)

  if (nrow(infection_rates) == 0 && nrow(per_capita_rates) == 0) {
    return(numeric())
  }

  n_infection_rates <- nrow(infection_rates)
  n_per_capita_rates <- nrow(per_capita_rates)
  n_age_groups <- context$n_age_groups
  rates <- matrix(
    numeric((n_infection_rates + n_per_capita_rates) * n_age_groups),
    nrow = n_infection_rates + n_per_capita_rates,
    ncol = n_age_groups
  )

  if (n_infection_rates > 0) {
    rates[seq_len(n_infection_rates), ] <- infection_rates
  }
  if (n_per_capita_rates > 0) {
    rates[n_infection_rates + seq_len(n_per_capita_rates), ] <- per_capita_rates
  }

  as.numeric(rates)
}

generic_infection_rate_matrix_from_state_matrix <- function(state_matrix, context) {
  model <- context$model
  if (nrow(model$infection_transitions) == 0) {
    return(matrix(numeric(), nrow = 0, ncol = context$n_age_groups))
  }

  population <- transition_population_by_age_from_state_matrix(state_matrix)
  validate_positive_age_populations(population, context$age_structure)

  infectious <- numeric(context$n_age_groups)
  if (length(context$generic_infectious_compartment_indices) > 0) {
    for (i in seq_along(context$generic_infectious_compartment_indices)) {
      infectious <- infectious +
        state_matrix[, context$generic_infectious_compartment_indices[i]] *
        context$generic_infectiousness_matrix_by_age[, i]
    }
  }

  lambda <- as.numeric(context$beta * (context$contact_matrix %*% (infectious / population)))

  n_infection_transitions <- nrow(model$infection_transitions)
  rates <- matrix(numeric(n_infection_transitions * context$n_age_groups), nrow = n_infection_transitions, ncol = context$n_age_groups)

  for (i in seq_len(n_infection_transitions)) {
    from_values <- state_matrix[, context$generic_infection_from_indices[i]]
    rates[i, ] <- lambda * context$infection_susceptibility_matrix[i, ] * from_values
  }

  rates
}

generic_per_capita_rate_matrix_from_state_matrix <- function(state_matrix, context) {
  model <- context$model
  if (nrow(model$transitions) == 0) {
    return(matrix(numeric(), nrow = 0, ncol = context$n_age_groups))
  }

  n_transitions <- nrow(model$transitions)
  rates <- matrix(numeric(n_transitions * context$n_age_groups), nrow = n_transitions, ncol = context$n_age_groups)

  for (i in seq_len(n_transitions)) {
    from_values <- state_matrix[, context$generic_transition_from_indices[i]]
    rates[i, ] <- context$generic_transition_rate_matrix[i, ] * from_values
  }

  rates
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
  derivative_index_lookup <- context$derivative_index_lookup
  if (is.null(derivative_index_lookup) ||
      length(derivative_index_lookup$from_state_index) != nrow(transition_rate_table)) {
    derivative_index_lookup <- prepare_derivative_index_lookup(
      transition_rate_table = transition_rate_table,
      state_order_key = context$state_order_key
    )
  }

  derivative$derivative <- accumulate_transition_derivative(
    rate = transition_rate_table$rate,
    from_state_index = derivative_index_lookup$from_state_index,
    to_state_index = derivative_index_lookup$to_state_index,
    state_length = nrow(derivative)
  )

  derivative[, c("compartment", "age_group", "derivative")]
}

stochastic_event_table_from_state_vector <- function(state_vector, context) {
  state_matrix <- transition_state_matrix_from_vector(state_vector, context)
  rate_values <- if (identical(context$model$model_type, "CompartmentModel")) {
    generic_transition_rate_values_from_state_matrix(state_matrix, context)
  } else {
    specialized_transition_rate_values_from_state_matrix(state_matrix, context)
  }

  event_table <- context$stochastic_event_template
  event_table$rate <- rate_values[context$stochastic_event_order]
  event_table
}
