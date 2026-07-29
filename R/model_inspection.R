#' Inspect logical disease-model transitions
#'
#' Returns a tidy, structural description of the logical transitions stored in a
#' supported disease model. This is a model-level inspection only: no state,
#' contact matrix, or transmission parameter is required.
#'
#' The returned table preserves the model's deterministic transition ordering and
#' covers built-in SIR/SEIR transitions, generic infection transitions, generic
#' per-capita transitions, and explicit outflows from the represented
#' compartment system.
#'
#' @param model Disease model created by [SIRModel()], [SEIRModel()], or
#'   [CompartmentModel()].
#'
#' @return A data frame with one row per logical transition and columns:
#'   `transition_id`, `transition_type`, `source`, `destination`,
#'   `rate_definition`, `rate_definition_type`, `rate_source`, and
#'   `is_infection_driven`.
#' @export
#' @examples
#' model <- CompartmentModel(
#'   compartments = c("S", "I", "R"),
#'   infection_transitions = data.frame(from = "S", to = "I"),
#'   transitions = data.frame(from = "I", to = "R", rate = 0.2)
#' )
#' inspect_transitions(model)
inspect_transitions <- function(model) {
  validate_disease_model(model)
  normalise_model_transitions(model, strict = TRUE)
}

#' Inspect inflows and outflows for model compartments
#'
#' Summarises how each compartment participates in the model transition graph,
#' or, when `compartment` is supplied, returns the transition-level inflow and
#' outflow table for that compartment.
#'
#' @param model Disease model created by [SIRModel()], [SEIRModel()], or
#'   [CompartmentModel()].
#' @param compartment Optional compartment name to inspect in detail.
#'
#' @return When `compartment = NULL`, a data frame with one row per compartment
#'   and columns `compartment`, `n_inflows`, `n_outflows`, `inflow_transition_ids`,
#'   `outflow_transition_ids`, `inflow_sources`, and `outflow_destinations`.
#'   When `compartment` is supplied, a transition-level data frame with columns
#'   `direction`, `transition_id`, `source`, `destination`, `transition_type`,
#'   and `rate_definition`.
#' @export
#' @examples
#' model <- CompartmentModel(
#'   compartments = c("S", "I", "R"),
#'   infection_transitions = data.frame(from = "S", to = "I"),
#'   transitions = data.frame(from = "I", to = "R", rate = 0.2)
#' )
#' inspect_compartment_flows(model)
#' inspect_compartment_flows(model, "I")
inspect_compartment_flows <- function(model, compartment = NULL) {
  transitions <- inspect_transitions(model)
  compartments <- validate_compartment_flow_model(model)

  if (is.null(compartment)) {
    rows <- lapply(compartments, function(compartment_name) {
      inflow_rows <- transitions[
        !is.na(transitions$destination) & transitions$destination == compartment_name,
        ,
        drop = FALSE
      ]
      outflow_rows <- transitions[
        transitions$source == compartment_name &
          (is.na(transitions$destination) | transitions$destination != compartment_name),
        ,
        drop = FALSE
      ]

      data.frame(
        compartment = compartment_name,
        n_inflows = nrow(inflow_rows),
        n_outflows = nrow(outflow_rows),
        stringsAsFactors = FALSE
      )
    })

    summary <- do.call(rbind, rows)
    summary$inflow_transition_ids <- I(lapply(compartments, function(compartment_name) {
      transitions$transition_id[
        !is.na(transitions$destination) & transitions$destination == compartment_name
      ]
    }))
    summary$outflow_transition_ids <- I(lapply(compartments, function(compartment_name) {
      transitions$transition_id[
        transitions$source == compartment_name &
          (is.na(transitions$destination) | transitions$destination != compartment_name)
      ]
    }))
    summary$inflow_sources <- I(lapply(compartments, function(compartment_name) {
      transitions$source[
        !is.na(transitions$destination) & transitions$destination == compartment_name
      ]
    }))
    summary$outflow_destinations <- I(lapply(compartments, function(compartment_name) {
      transitions$destination[
        transitions$source == compartment_name &
          (is.na(transitions$destination) | transitions$destination != compartment_name)
      ]
    }))
    row.names(summary) <- NULL
    return(summary)
  }

  if (!is.character(compartment) || length(compartment) != 1 || is.na(compartment) || !nzchar(compartment)) {
    stop("compartment must be a single non-empty string.", call. = FALSE)
  }
  if (!compartment %in% compartments) {
    stop("compartment is not present in model compartments: ", compartment, call. = FALSE)
  }

  inflow_rows <- transitions[
    !is.na(transitions$destination) & transitions$destination == compartment,
    ,
    drop = FALSE
  ]
  outflow_rows <- transitions[
    transitions$source == compartment &
    (is.na(transitions$destination) | transitions$destination != compartment),
    ,
    drop = FALSE
  ]

  if (nrow(inflow_rows) > 0) {
    inflow_rows$direction <- "inflow"
  }
  if (nrow(outflow_rows) > 0) {
    outflow_rows$direction <- "outflow"
  }

  selected <- rbind(inflow_rows, outflow_rows)
  if (nrow(selected) == 0) {
    return(data.frame(
      direction = character(),
      transition_id = character(),
      source = character(),
      destination = character(),
      transition_type = character(),
      rate_definition = I(list()),
      stringsAsFactors = FALSE
    ))
  }

  selected <- selected[, c(
    "direction",
    "transition_id",
    "source",
    "destination",
    "transition_type",
    "rate_definition"
  ), drop = FALSE]
  row.names(selected) <- NULL
  selected
}

#' Diagnose the structural integrity of a disease model
#'
#' Produces a structured report about logical transition duplicates,
#' compartment connectivity, and optional reachability from an initial state.
#' This is a structural diagnostic only; it does not evaluate transition rates.
#'
#' @param model Disease model created by [SIRModel()], [SEIRModel()], or
#'   [CompartmentModel()].
#' @param initial_state Optional initial state used only to identify occupied
#'   compartments for reachability checks. The state may be a numeric vector or a
#'   data frame with at least `compartment` and `value` columns.
#'
#' @return A data frame with columns `severity`, `check`, `compartment`,
#'   `transition_id`, `source`, `destination`, and `details`.
#' @export
#' @examples
#' model <- CompartmentModel(
#'   compartments = c("S", "I", "R"),
#'   infection_transitions = data.frame(from = "S", to = "I"),
#'   transitions = data.frame(from = "I", to = "R", rate = 0.2)
#' )
#' initial_state <- data.frame(
#'   compartment = rep(c("S", "I", "R"), each = 2),
#'   age_group = rep(c("0-4", "5-9"), times = 3),
#'   value = c(95, 90, 5, 10, 0, 0)
#' )
#' diagnose_model_structure(model, initial_state)
diagnose_model_structure <- function(model, initial_state = NULL) {
  validate_supported_model_object(model)
  transitions <- diagnostic_model_transitions(model)
  compartments <- validate_compartment_flow_model(model)
  report <- list()

  add_row <- function(severity, check, compartment = NA_character_, transition_id = NA_character_,
                      source = NA_character_, destination = NA_character_, details) {
    report[[length(report) + 1L]] <<- data.frame(
      severity = severity,
      check = check,
      compartment = compartment,
      transition_id = transition_id,
      source = source,
      destination = destination,
      details = details,
      stringsAsFactors = FALSE
    )
  }

  duplicate_ids <- transitions$transition_id[duplicated(transitions$transition_id)]
  if (length(duplicate_ids) > 0) {
    for (transition_id in unique(duplicate_ids)) {
      rows <- transitions[transitions$transition_id == transition_id, , drop = FALSE]
      add_row(
        "error",
        "duplicate logical transition id",
        transition_id = transition_id,
        source = rows$source[1],
        destination = rows$destination[1],
        details = paste("transition_id appears", nrow(rows), "times.")
      )
    }
  } else {
    add_row("pass", "duplicate logical transition id", details = "No duplicate transition IDs detected.")
  }

  duplicate_pairs <- duplicated(data.frame(
    source = transitions$source,
    destination = transitions$destination,
    transition_type = transitions$transition_type,
    stringsAsFactors = FALSE
  ))
  if (any(duplicate_pairs)) {
    repeated <- unique(transitions[duplicate_pairs, c("source", "destination", "transition_type"), drop = FALSE])
    for (i in seq_len(nrow(repeated))) {
      rows <- transitions[
        transitions$source == repeated$source[i] &
          transitions$destination == repeated$destination[i] &
          transitions$transition_type == repeated$transition_type[i],
        ,
        drop = FALSE
      ]
      add_row(
        "error",
        "duplicate logical transition pair",
        source = repeated$source[i],
        destination = repeated$destination[i],
        details = paste(
          "transition pair appears",
          nrow(rows),
          "times for transition_type",
          repeated$transition_type[i],
          "."
        )
      )
    }
  } else {
    add_row("pass", "duplicate logical transition pair", details = "No duplicate source-destination pairs detected.")
  }

  flow_summary <- diagnostic_compartment_flow_summary(transitions, compartments)
  for (i in seq_len(nrow(flow_summary))) {
    compartment_name <- flow_summary$compartment[i]
    n_inflows <- flow_summary$n_inflows[i]
    n_outflows <- flow_summary$n_outflows[i]
    if (n_inflows == 0 && n_outflows == 0) {
      add_row(
        "warning",
        "isolated compartment",
        compartment = compartment_name,
        details = "No inflows and no outflows."
      )
    } else if (n_inflows == 0) {
      add_row(
        "warning",
        "source-only compartment",
        compartment = compartment_name,
        details = "Outflows exist but no inflows were found."
      )
    } else if (n_outflows == 0) {
      add_row(
        "note",
        "sink-only compartment",
        compartment = compartment_name,
        details = "Inflows exist but no outflows were found."
      )
    } else {
      add_row(
        "pass",
        "compartment connectivity",
        compartment = compartment_name,
        details = "At least one inflow and one outflow were found."
      )
    }
  }

  if (is.null(initial_state)) {
    add_row(
      "note",
      "reachability",
      details = "Reachability was not evaluated because initial_state was not supplied."
    )
  } else {
    occupied <- initially_occupied_compartments(initial_state, compartments)
    reachable <- reachable_compartments_from_initial_state(transitions, occupied)
    unreachable_compartments <- setdiff(compartments, reachable)

    add_row(
      "pass",
      "reachability",
      details = paste(
        "Occupied compartments:",
        if (length(occupied) == 0) "<none>" else paste(occupied, collapse = ", ")
      )
    )

    if (length(unreachable_compartments) == 0) {
      add_row("pass", "reachable compartments", details = "All compartments are reachable from the initial state.")
    } else {
      for (compartment_name in unreachable_compartments) {
        add_row(
          "warning",
          "unreachable compartment",
          compartment = compartment_name,
          details = "No directed path exists from the occupied initial compartments."
        )
      }
    }

    unreachable_transition_rows <- transitions[!is.na(transitions$source) & !(transitions$source %in% reachable), , drop = FALSE]
    if (nrow(unreachable_transition_rows) == 0) {
      add_row("pass", "unreachable transitions", details = "Every transition source is reachable.")
    } else {
      for (i in seq_len(nrow(unreachable_transition_rows))) {
        add_row(
          "note",
          "unreachable transition",
          transition_id = unreachable_transition_rows$transition_id[i],
          source = unreachable_transition_rows$source[i],
          destination = unreachable_transition_rows$destination[i],
          details = "The source compartment is unreachable from the initial state."
        )
      }
    }
  }

  do.call(rbind, report)
}

#' Inspect evaluated transition rates for a supplied state
#'
#' A user-facing wrapper around [transition_rates()] that adds source
#' populations, per-capita source rates, and simple diagnostic flags.
#'
#' @param state Long-form state data frame or numeric state vector.
#' @param model Disease model created by [SIRModel()], [SEIRModel()], or
#'   [CompartmentModel()].
#' @param age_structure Valid age structure.
#' @param contact_matrix Numeric contact matrix with rows as recipient age
#'   groups and columns as source age groups.
#' @param ... Additional arguments forwarded to [transition_rates()].
#'
#' @return A data frame containing the evaluated transition table plus source
#'   population, per-capita source rate, and diagnostic columns.
#' @export
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("0-4", "5-9"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, 9)
#' )
#' model <- SIRModel(gamma = 0.2)
#' state <- data.frame(
#'   compartment = rep(c("S", "I", "R"), each = 2),
#'   age_group = rep(ages$age_groups, times = 3),
#'   value = c(95, 90, 5, 10, 0, 0)
#' )
#' inspect_transition_rates(state, model, ages, diag(2), beta = 0.1)
inspect_transition_rates <- function(state,
                                     model,
                                     age_structure,
                                     contact_matrix,
                                     ...) {
  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    ...
  )
  transitions <- inspect_transitions(model)
  state_long <- transition_state_to_long(state, age_structure, model$compartments)
  state_key <- paste(state_long$compartment, state_long$age_group, sep = "\r")
  source_key <- paste(rates$from, rates$age_group, sep = "\r")
  source_population <- state_long$value[match(source_key, state_key)]

  structural <- transitions[match(rates$transition_id, transitions$transition_id), , drop = FALSE]

  out <- cbind(
    rates,
    transition_type = structural$transition_type,
    rate_definition = I(structural$rate_definition),
    rate_definition_type = structural$rate_definition_type,
    rate_source = structural$rate_source,
    is_infection_driven = structural$is_infection_driven,
    source_population = as.numeric(source_population)
  )
  out$per_capita_source_rate <- ifelse(
    out$source_population > 0,
    out$rate / out$source_population,
    NA_real_
  )
  out$is_finite <- is.finite(out$rate)
  out$is_negative <- out$rate < 0
  out$source_empty <- out$source_population == 0
  out$nonzero_from_empty_source <- out$source_empty & abs(out$rate) > 0
  out[, c(
    "from",
    "to",
    "age_group",
    "transition_id",
    "transition_type",
    "rate_definition",
    "rate_definition_type",
    "rate_source",
    "is_infection_driven",
    "rate",
    "source_population",
    "per_capita_source_rate",
    "is_finite",
    "is_negative",
    "source_empty",
    "nonzero_from_empty_source"
  )]
}

#' Check population balance from evaluated disease-model transitions
#'
#' Evaluates transition rates, reconstructs compartment derivatives, and
#' compares the total derivative against the net external contribution implied
#' by the evaluated transition table.
#'
#' This check covers disease-model transitions only. It does not include a
#' separate demographic process, ageing operator, or other simulation-side
#' adjustments. External inflows are only counted when they are represented in
#' the evaluated transition table; the current disease-model transition
#' machinery primarily represents external outflows via `to = NA`.
#'
#' @param state Long-form state data frame or numeric state vector.
#' @param model Disease model created by [SIRModel()], [SEIRModel()], or
#'   [CompartmentModel()].
#' @param age_structure Valid age structure.
#' @param contact_matrix Numeric contact matrix with rows as recipient age
#'   groups and columns as source age groups.
#' @param tolerance Numeric scalar tolerance for the balance check.
#' @param ... Additional arguments forwarded to [transition_rates()].
#'
#' @return A list with `summary` and `by_age` data frames.
#' @export
#' @examples
#' ages <- AgeStructure(
#'   age_groups = c("0-4", "5-9"),
#'   lower_bounds = c(0, 5),
#'   upper_bounds = c(4, 9)
#' )
#' model <- SIRModel(gamma = 0.2)
#' state <- data.frame(
#'   compartment = rep(c("S", "I", "R"), each = 2),
#'   age_group = rep(ages$age_groups, times = 3),
#'   value = c(95, 90, 5, 10, 0, 0)
#' )
#' check_population_balance(state, model, ages, diag(2), beta = 0.1)
check_population_balance <- function(state,
                                     model,
                                     age_structure,
                                     contact_matrix,
                                     ...,
                                     tolerance = sqrt(.Machine$double.eps)) {
  if (!is.numeric(tolerance) || length(tolerance) != 1 || is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a single finite non-negative numeric value.", call. = FALSE)
  }

  rates <- transition_rates(
    state = state,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    ...
  )
  derivative <- rates_to_derivative(
    transition_rate_table = rates,
    compartments = model$compartments,
    age_structure = age_structure
  )

  internal_rows <- !is.na(rates$to)
  external_outflow_rows <- is.na(rates$to)
  external_inflow_rows <- is.na(rates$from)

  by_age_rows <- lapply(age_structure$age_groups, function(age_group) {
    age_rates <- rates[rates$age_group == age_group, , drop = FALSE]
    age_derivative <- derivative[derivative$age_group == age_group, , drop = FALSE]
    internal_inflow <- sum(age_rates$rate[!is.na(age_rates$to)])
    internal_outflow <- sum(age_rates$rate[!is.na(age_rates$to)])
    external_inflow <- sum(age_rates$rate[is.na(age_rates$from)])
    external_outflow <- sum(age_rates$rate[is.na(age_rates$to)])
    total_derivative <- sum(age_derivative$derivative)
    expected_net_population_change <- external_inflow - external_outflow
    residual_balance_error <- total_derivative - expected_net_population_change

    data.frame(
      age_group = age_group,
      total_internal_inflow = internal_inflow,
      total_internal_outflow = internal_outflow,
      net_internal_transfer = internal_inflow - internal_outflow,
      total_external_inflow = external_inflow,
      total_external_outflow = external_outflow,
      total_derivative = total_derivative,
      expected_net_population_change = expected_net_population_change,
      residual_balance_error = residual_balance_error,
      passes_balance_check = abs(residual_balance_error) <= tolerance,
      stringsAsFactors = FALSE
    )
  })
  by_age <- do.call(rbind, by_age_rows)
  row.names(by_age) <- NULL

  summary <- data.frame(
    total_internal_inflow = sum(rates$rate[internal_rows]),
    total_internal_outflow = sum(rates$rate[internal_rows]),
    net_internal_transfer = sum(rates$rate[internal_rows]) - sum(rates$rate[internal_rows]),
    total_external_inflow = sum(rates$rate[external_inflow_rows]),
    total_external_outflow = sum(rates$rate[external_outflow_rows]),
    total_derivative = sum(derivative$derivative),
    expected_net_population_change = sum(rates$rate[external_inflow_rows]) - sum(rates$rate[external_outflow_rows]),
    residual_balance_error = sum(derivative$derivative) -
      (sum(rates$rate[external_inflow_rows]) - sum(rates$rate[external_outflow_rows])),
    passes_balance_check = abs(
      sum(derivative$derivative) -
        (sum(rates$rate[external_inflow_rows]) - sum(rates$rate[external_outflow_rows]))
    ) <= tolerance,
    tolerance = tolerance,
    stringsAsFactors = FALSE
  )

  list(summary = summary, by_age = by_age)
}

validate_supported_model_object <- function(model) {
  if (!is.list(model)) {
    stop("model must be a supported disease-model object.", call. = FALSE)
  }
  if (is.null(model$model_type) || !is.character(model$model_type) || length(model$model_type) != 1 || is.na(model$model_type)) {
    stop("model must define a supported disease-model type.", call. = FALSE)
  }
  if (!model$model_type %in% c("SIR", "SEIR", "CompartmentModel")) {
    stop("model must be a supported disease-model object.", call. = FALSE)
  }
  invisible(model)
}

validate_compartment_flow_model <- function(model) {
  validate_supported_model_object(model)
  validate_compartments(model$compartments)
  model$compartments
}

normalise_model_transitions <- function(model, strict = TRUE) {
  validate_supported_model_object(model)
  if (strict) {
    validate_disease_model(model)
  }

  if (model$model_type == "SIR") {
    return(normalise_sir_transitions(model))
  }

  if (model$model_type == "SEIR") {
    return(normalise_seir_transitions(model))
  }

  normalise_compartment_model_transitions(model)
}

normalise_sir_transitions <- function(model) {
  transitions <- data.frame(
    transition_id = c("infection:S->I", "transition:I->R"),
    transition_type = c("infection", "transition"),
    source = c("S", "I"),
    destination = c("I", "R"),
    rate_definition_type = c("derived", "scalar"),
    rate_source = c("force_of_infection()", "model$gamma"),
    is_infection_driven = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  transitions$rate_definition <- I(list(NA_real_, model$gamma))
  transitions
}

normalise_seir_transitions <- function(model) {
  transitions <- data.frame(
    transition_id = c("infection:S->E", "transition:E->I", "transition:I->R"),
    transition_type = c("infection", "transition", "transition"),
    source = c("S", "E", "I"),
    destination = c("E", "I", "R"),
    rate_definition_type = c("derived", "scalar", "scalar"),
    rate_source = c("force_of_infection()", "model$sigma", "model$gamma"),
    is_infection_driven = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  transitions$rate_definition <- I(list(NA_real_, model$sigma, model$gamma))
  transitions
}

normalise_compartment_model_transitions <- function(model) {
  rows <- list()
  row_index <- 0L

  if (!is.null(model$infection_transitions) && nrow(model$infection_transitions) > 0) {
    for (i in seq_len(nrow(model$infection_transitions))) {
      row_index <- row_index + 1L
      susceptibility <- model$infection_transitions$susceptibility[[i]]
      rows[[row_index]] <- data.frame(
        transition_id = transition_identifiers(
          from = model$infection_transitions$from[i],
          to = model$infection_transitions$to[i],
          transition_type = "infection"
        ),
        transition_type = "infection",
        source = model$infection_transitions$from[i],
        destination = model$infection_transitions$to[i],
        rate_definition_type = transition_rate_definition_type(susceptibility),
        rate_definition = I(list(susceptibility)),
        rate_source = "infection_transitions$susceptibility",
        is_infection_driven = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!is.null(model$transitions) && nrow(model$transitions) > 0) {
    for (i in seq_len(nrow(model$transitions))) {
      row_index <- row_index + 1L
      transition_type <- if ("transition_type" %in% names(model$transitions)) {
        model$transitions$transition_type[i]
      } else if (is.na(model$transitions$to[i])) {
        "outflow"
      } else {
        "transition"
      }
      if (identical(transition_type, "internal")) {
        transition_type <- "transition"
      }
      rate_definition <- model$transitions$rate[[i]]
      rows[[row_index]] <- data.frame(
        transition_id = if ("transition_id" %in% names(model$transitions)) {
          model$transitions$transition_id[i]
        } else {
          transition_identifiers(
            from = model$transitions$from[i],
            to = model$transitions$to[i],
            transition_type = if (transition_type == "outflow") "outflow" else "transition"
          )
        },
        transition_type = transition_type,
        source = model$transitions$from[i],
        destination = model$transitions$to[i],
        rate_definition_type = transition_rate_definition_type(rate_definition),
        rate_definition = I(list(rate_definition)),
        rate_source = "transitions$rate",
        is_infection_driven = FALSE,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(data.frame(
      transition_id = character(),
      transition_type = character(),
      source = character(),
      destination = character(),
      rate_definition = I(list()),
      rate_definition_type = character(),
      rate_source = character(),
      is_infection_driven = logical(),
      stringsAsFactors = FALSE
    ))
  }

  transitions <- do.call(rbind, rows)
  row.names(transitions) <- NULL
  transitions
}

transition_rate_definition_type <- function(x) {
  if (is.null(x)) {
    return("missing")
  }
  if (is.atomic(x) && length(x) == 1 && is.na(x)) {
    return("derived")
  }
  if (is.atomic(x) && length(x) == 1) {
    return("scalar")
  }
  if (is.atomic(x) && length(x) > 1) {
    return("age_specific")
  }
  if (is.list(x)) {
    return("list")
  }
  "other"
}

diagnostic_model_transitions <- function(model) {
  if (model$model_type == "SIR") {
    return(normalise_sir_transitions(model))
  }

  if (model$model_type == "SEIR") {
    return(normalise_seir_transitions(model))
  }

  infection_transitions <- model$infection_transitions
  if (is.null(infection_transitions) || !is.data.frame(infection_transitions) || nrow(infection_transitions) == 0) {
    infection_rows <- empty_transition_inspection_table()
  } else {
    susceptibility <- if ("susceptibility" %in% names(infection_transitions)) {
      infection_transitions$susceptibility
    } else {
      rep(list(1), nrow(infection_transitions))
    }
    if (!is.list(susceptibility)) {
      susceptibility <- as.list(susceptibility)
    }
    infection_rows <- data.frame(
      transition_id = if ("transition_id" %in% names(infection_transitions)) {
        infection_transitions$transition_id
      } else {
        transition_identifiers(
          from = infection_transitions$from,
          to = infection_transitions$to,
          transition_type = "infection"
        )
      },
      transition_type = "infection",
      source = infection_transitions$from,
      destination = infection_transitions$to,
      rate_definition = I(susceptibility),
      rate_definition_type = vapply(
        susceptibility,
        transition_rate_definition_type,
        character(1)
      ),
      rate_source = "infection_transitions$susceptibility",
      is_infection_driven = TRUE,
      stringsAsFactors = FALSE
    )
  }

  transitions <- model$transitions
  if (is.null(transitions) || !is.data.frame(transitions) || nrow(transitions) == 0) {
    transition_rows <- empty_transition_inspection_table()
  } else {
    rate_definition <- if ("rate" %in% names(transitions)) {
      transitions$rate
    } else {
      rep(list(NA_real_), nrow(transitions))
    }
    if (!is.list(rate_definition)) {
      rate_definition <- as.list(rate_definition)
    }
    transition_type <- if ("transition_type" %in% names(transitions)) {
      as.character(transitions$transition_type)
    } else {
      ifelse(is.na(transitions$to), "outflow", "transition")
    }
    transition_type[transition_type == "internal"] <- "transition"
    transition_type[is.na(transition_type) & is.na(transitions$to)] <- "outflow"
    transition_type[is.na(transition_type)] <- "transition"
    transition_rows <- data.frame(
      transition_id = if ("transition_id" %in% names(transitions)) {
        transitions$transition_id
      } else {
        transition_identifiers(
          from = transitions$from,
          to = transitions$to,
          transition_type = ifelse(is.na(transitions$to), "outflow", "transition")
        )
      },
      transition_type = transition_type,
      source = transitions$from,
      destination = transitions$to,
      rate_definition = I(rate_definition),
      rate_definition_type = vapply(rate_definition, transition_rate_definition_type, character(1)),
      rate_source = "transitions$rate",
      is_infection_driven = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (nrow(infection_rows) == 0) {
    transitions <- transition_rows
  } else if (nrow(transition_rows) == 0) {
    transitions <- infection_rows
  } else {
    transitions <- rbind(infection_rows, transition_rows)
  }

  row.names(transitions) <- NULL
  transitions
}

empty_transition_inspection_table <- function() {
  data.frame(
    transition_id = character(),
    transition_type = character(),
    source = character(),
    destination = character(),
    rate_definition = I(list()),
    rate_definition_type = character(),
    rate_source = character(),
    is_infection_driven = logical(),
    stringsAsFactors = FALSE
  )
}

diagnostic_compartment_flow_summary <- function(transitions, compartments) {
  rows <- lapply(compartments, function(compartment_name) {
    inflow_rows <- transitions[
      !is.na(transitions$destination) & transitions$destination == compartment_name,
      ,
      drop = FALSE
    ]
    outflow_rows <- transitions[
      transitions$source == compartment_name &
        (is.na(transitions$destination) | transitions$destination != compartment_name),
      ,
      drop = FALSE
    ]

    data.frame(
      compartment = compartment_name,
      n_inflows = nrow(inflow_rows),
      n_outflows = nrow(outflow_rows),
      stringsAsFactors = FALSE
    )
  })

  summary <- do.call(rbind, rows)
  summary$inflow_transition_ids <- I(lapply(compartments, function(compartment_name) {
    transitions$transition_id[
      !is.na(transitions$destination) & transitions$destination == compartment_name
    ]
  }))
  summary$outflow_transition_ids <- I(lapply(compartments, function(compartment_name) {
    transitions$transition_id[
      transitions$source == compartment_name &
        (is.na(transitions$destination) | transitions$destination != compartment_name)
    ]
  }))
  summary$inflow_sources <- I(lapply(compartments, function(compartment_name) {
    transitions$source[
      !is.na(transitions$destination) & transitions$destination == compartment_name
    ]
  }))
  summary$outflow_destinations <- I(lapply(compartments, function(compartment_name) {
    transitions$destination[
      transitions$source == compartment_name &
        (is.na(transitions$destination) | transitions$destination != compartment_name)
    ]
  }))
  row.names(summary) <- NULL
  summary
}

initially_occupied_compartments <- function(initial_state, compartments) {
  if (is.data.frame(initial_state)) {
    if (!"compartment" %in% names(initial_state) || !"value" %in% names(initial_state)) {
      stop("initial_state data frame must contain compartment and value columns.", call. = FALSE)
    }
    if (!is.numeric(initial_state$value) || anyNA(initial_state$value) || any(!is.finite(initial_state$value))) {
      stop("initial_state value must be finite numeric values.", call. = FALSE)
    }
    if (any(initial_state$value < 0)) {
      stop("initial_state value cannot contain negative values.", call. = FALSE)
    }
    if (anyNA(initial_state$compartment)) {
      stop("initial_state compartment cannot contain missing values.", call. = FALSE)
    }
    occupied <- unique(as.character(initial_state$compartment[initial_state$value > 0]))
    unknown <- setdiff(occupied, compartments)
    if (length(unknown) > 0) {
      stop("initial_state contains unknown compartment value(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    return(intersect(compartments, occupied))
  }

  if (is.numeric(initial_state) && !is.matrix(initial_state)) {
    if (length(initial_state) != length(compartments)) {
      stop("initial_state length must match the number of model compartments when supplied as a numeric vector.", call. = FALSE)
    }
    if (anyNA(initial_state) || any(!is.finite(initial_state))) {
      stop("initial_state must contain finite non-missing values.", call. = FALSE)
    }
    if (any(initial_state < 0)) {
      stop("initial_state cannot contain negative values.", call. = FALSE)
    }
    if (!is.null(names(initial_state)) && any(names(initial_state) != "")) {
      unknown <- setdiff(names(initial_state), compartments)
      if (length(unknown) > 0) {
        stop("initial_state contains unknown compartment name(s): ", paste(unknown, collapse = ", "), call. = FALSE)
      }
      return(compartments[compartments %in% names(initial_state)[initial_state > 0]])
    }
    return(compartments[initial_state > 0])
  }

  stop("initial_state must be a numeric vector or a data frame.", call. = FALSE)
}

reachable_compartments_from_initial_state <- function(transitions, occupied) {
  reachable <- occupied
  queue <- occupied
  adjacency <- split(
    transitions$destination[!is.na(transitions$destination)],
    transitions$source[!is.na(transitions$destination)]
  )

  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    next_compartments <- unique(adjacency[[current]])
    next_compartments <- next_compartments[!is.na(next_compartments)]
    new_compartments <- setdiff(next_compartments, reachable)
    if (length(new_compartments) > 0) {
      reachable <- c(reachable, new_compartments)
      queue <- c(queue, new_compartments)
    }
  }

  unique(reachable)
}
