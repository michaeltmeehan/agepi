#' Simulate a stochastic compartment model with Gillespie's direct method
#'
#' Runs a narrow stochastic simulation for supported age-structured
#' compartmental models: fixed population, no demography, no ageing, no
#' fertility, no mortality, no migration, and no tau-leaping.
#'
#' The simulator uses Gillespie's continuous-time direct method. Infection
#' propensities reuse the package transition-rate interface. Each transition
#' rate row becomes one individual event type within an age group. Infection
#' propensities therefore use the existing force-of-infection convention: rows
#' of `contact_matrix` are recipient age groups and columns are infectious
#' source age groups.
#'
#' The initial state may be supplied either as long-form state data or as a
#' numeric vector in the existing compartment-major, age-group-minor ordering.
#' Names on numeric state vectors are ignored, matching the state-mapping
#' helpers. Output rows are aligned to the requested `times`; if no event occurs
#' before an output time, the current state is carried forward.
#'
#' If `seed` is supplied, the R random-number state is restored on exit when
#' possible, so the caller's RNG stream is not permanently advanced by this
#' function.
#'
#' @param initial_state Long-form state data frame with columns `compartment`,
#'   `age_group`, and `value`, or numeric state vector.
#' @param times Numeric vector of finite, non-missing, strictly increasing time
#'   points. Must have length at least two.
#' @param model Disease model. `SIRModel()`, `SEIRModel()`, and
#'   `CompartmentModel()` outputs with supported fixed-population transition
#'   structures expressible through `transition_rates()` are currently
#'   supported.
#' @param age_structure Valid age structure.
#' @param contact_matrix Numeric contact matrix with rows as recipient age
#'   groups and columns as source age groups.
#' @param beta Non-negative finite transmission scaling parameter.
#' @param susceptibility Optional non-negative numeric vector by recipient age
#'   group.
#' @param infectiousness Optional non-negative numeric vector by source age
#'   group.
#' @param population Optional fixed positive population vector by age group.
#'   When supplied, it must match the age-specific totals in `initial_state`.
#'   When omitted, age-specific totals from `initial_state` are used.
#' @param method Simulation method. Only `"gillespie"` is currently supported.
#' @param seed Optional integer seed for reproducible simulation.
#' @param return_events Logical scalar. If `FALSE`, return only the trajectory.
#'   If `TRUE`, return a list with `trajectory` and `events`.
#' @param demographic_process Unsupported. Stochastic demography is not included
#'   and must be `NULL`.
#' @param cumulative_flows Optional named list or data frame specifying disease
#'   transition flows to summarise from realised stochastic events. Named list
#'   entries must contain `from` and `to` fields; data frames must contain
#'   `name`, `from`, and `to` columns. Cumulative flows are derived from the
#'   event log and do not add stochastic state variables or propensities.
#'
#' @return If `return_events = FALSE`, a data frame with columns `time`,
#'   `compartment`, `age_group`, and `value`, ordered by time outermost,
#'   compartment next, and age group innermost. If `return_events = TRUE`, a
#'   list with `trajectory` and an event log containing `time`, `event`,
#'   `transition_id`, `age_group`, `from`, `to`, and `rate`. When
#'   `cumulative_flows` is supplied, a list is returned with `trajectory` and
#'   `cumulative`, and also `events` when `return_events = TRUE`. `cumulative`
#'   contains columns `time`, `cumulative_name`, `transition_id`, `from`, `to`,
#'   `age_group`, and `value`.
#' @export
simulate_stochastic <- function(
  initial_state,
  times,
  model,
  age_structure,
  contact_matrix,
  beta = 1,
  susceptibility = NULL,
  infectiousness = NULL,
  population = NULL,
  method = "gillespie",
  seed = NULL,
  return_events = FALSE,
  demographic_process = NULL,
  cumulative_flows = NULL
) {
  validate_stochastic_method(method)
  validate_stochastic_return_events(return_events)
  validate_stochastic_seed(seed)
  validate_simulation_times(times)
  validate_disease_model(model)
  validate_age_structure(age_structure)
  validate_stochastic_model(model)
  validate_stochastic_no_demography(demographic_process)
  validate_contact_matrix(contact_matrix, age_structure)
  beta <- validate_force_beta(beta)
  susceptibility <- validate_optional_force_vector(
    susceptibility,
    age_structure$n_age_groups,
    "susceptibility"
  )
  infectiousness <- validate_optional_force_vector(
    infectiousness,
    age_structure$n_age_groups,
    "infectiousness"
  )

  state_vector <- simulation_state_to_vector(
    initial_state,
    age_structure,
    model$compartments
  )
  validate_non_negative_simulation_state(state_vector, "initial_state")
  validate_stochastic_integer_state(state_vector)

  population <- validate_stochastic_population(
    population = population,
    state_vector = state_vector,
    age_structure = age_structure,
    compartments = model$compartments
  )

  cumulative_spec <- NULL
  if (!is.null(cumulative_flows)) {
    cumulative_spec <- prepare_stochastic_cumulative_flows(
      cumulative_flows = cumulative_flows,
      state_vector = state_vector,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness
    )
  }

  with_stochastic_seed(seed, {
    stochastic_gillespie_fixed_population(
      state_vector = state_vector,
      times = times,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness,
      population = population,
      return_events = return_events,
      cumulative_spec = cumulative_spec
    )
  })
}

stochastic_gillespie_fixed_population <- function(
  state_vector,
  times,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  population,
  return_events,
  cumulative_spec = NULL
) {
  current_state <- as.numeric(state_vector)
  current_time <- times[1]
  output_states <- vector("list", length(times))
  output_states[[1]] <- simulation_state_output(
    current_state,
    time = times[1],
    age_structure = age_structure,
    compartments = model$compartments
  )
  event_rows <- list()
  event_count <- 0L
  output_index <- 2L

  while (output_index <= length(times)) {
    propensities <- stochastic_propensities(
      state_vector = current_state,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness,
      population = population
    )
    total_rate <- sum(propensities$rate)

    if (total_rate <= 0) {
      while (output_index <= length(times)) {
        output_states[[output_index]] <- simulation_state_output(
          current_state,
          time = times[output_index],
          age_structure = age_structure,
          compartments = model$compartments
        )
        output_index <- output_index + 1L
      }
      break
    }

    next_event_time <- current_time + stats::rexp(1, rate = total_rate)

    while (output_index <= length(times) && times[output_index] < next_event_time) {
      output_states[[output_index]] <- simulation_state_output(
        current_state,
        time = times[output_index],
        age_structure = age_structure,
        compartments = model$compartments
      )
      output_index <- output_index + 1L
    }

    if (output_index > length(times)) {
      break
    }

    event_index <- stochastic_sample_event(propensities$rate, total_rate)
    event <- propensities[event_index, , drop = FALSE]
    current_state <- stochastic_apply_event(
      state_vector = current_state,
      event = event,
      age_structure = age_structure,
      compartments = model$compartments
    )
    current_time <- next_event_time

    if (return_events || !is.null(cumulative_spec)) {
      event_count <- event_count + 1L
      event_rows[[event_count]] <- data.frame(
        time = current_time,
        event = event$event,
        transition_id = event$transition_id,
        age_group = event$age_group,
        from = event$from,
        to = event$to,
        rate = event$rate,
        stringsAsFactors = FALSE
      )
    }

    while (output_index <= length(times) && times[output_index] == current_time) {
      output_states[[output_index]] <- simulation_state_output(
        current_state,
        time = times[output_index],
        age_structure = age_structure,
        compartments = model$compartments
      )
      output_index <- output_index + 1L
    }
  }

  trajectory <- do.call(rbind, output_states)
  row.names(trajectory) <- NULL

  if (!return_events && is.null(cumulative_spec)) {
    return(trajectory)
  }

  events <- stochastic_event_log(event_rows)
  if (is.null(cumulative_spec)) {
    return(list(trajectory = trajectory, events = events))
  }

  result <- list(
    trajectory = trajectory,
    cumulative = stochastic_cumulative_output(
      events = events,
      times = times,
      cumulative_spec = cumulative_spec
    )
  )
  if (return_events) {
    result$events <- events
    result <- result[c("trajectory", "events", "cumulative")]
  }
  result
}

stochastic_propensities <- function(
  state_vector,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  population
) {
  stochastic_event_table(
    state_vector = state_vector,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    population = population
  )
}

stochastic_event_table <- function(
  state_vector,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  population
) {
  rates <- transition_rates(
    state = as.numeric(state_vector),
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness
  )
  validate_stochastic_transition_rates(rates, model, age_structure)

  transition_order <- stochastic_transition_order(model)
  rates$.transition_index <- match(
    paste(rates$from, rates$to, sep = "->"),
    paste(transition_order$from, transition_order$to, sep = "->")
  )
  rates$age_index <- match(rates$age_group, age_structure$age_groups)
  rates$event <- stochastic_event_labels(rates, model)
  rates <- rates[order(rates$.transition_index, rates$age_index), ]
  row.names(rates) <- NULL

  rates[, c("event", "transition_id", "age_group", "age_index", "from", "to", "rate")]
}

stochastic_transition_order <- function(model) {
  if (identical(model$model_type, "CompartmentModel")) {
    transitions <- rbind(
      model$infection_transitions[, c("from", "to"), drop = FALSE],
      model$transitions[, c("from", "to"), drop = FALSE]
    )
    row.names(transitions) <- NULL
    return(transitions)
  }

  model$transitions[, c("from", "to"), drop = FALSE]
}

stochastic_event_labels <- function(rates, model) {
  labels <- paste(rates$from, rates$to, sep = "->")

  if (identical(model$model_type, "SIR")) {
    labels[labels == "S->I"] <- "infection"
    labels[labels == "I->R"] <- "recovery"
    return(labels)
  }

  if (identical(model$model_type, "SEIR")) {
    labels[labels == "S->E"] <- "infection"
    labels[labels == "E->I"] <- "progression"
    labels[labels == "I->R"] <- "recovery"
    return(labels)
  }

  infection_keys <- paste(
    model$infection_transitions$from,
    model$infection_transitions$to,
    sep = "->"
  )
  labels[labels %in% infection_keys] <- "infection"
  labels[labels == "E->I"] <- "progression"
  labels[labels == "I->R"] <- "recovery"
  labels
}

validate_stochastic_transition_rates <- function(rates, model, age_structure) {
  validate_transition_rate_table(
    transition_rate_table = rates,
    compartments = model$compartments,
    age_structure = age_structure
  )

  transition_order <- stochastic_transition_order(model)
  expected_keys <- paste(transition_order$from, transition_order$to, sep = "->")
  observed_keys <- paste(rates$from, rates$to, sep = "->")
  unknown_keys <- setdiff(unique(observed_keys), expected_keys)
  if (length(unknown_keys) > 0) {
    stop(
      "stochastic transition rates contain transition(s) not declared by the model: ",
      paste(unknown_keys, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(rates)
}

stochastic_sample_event <- function(rates, total_rate) {
  threshold <- stats::runif(1, min = 0, max = total_rate)
  which(cumsum(rates) > threshold)[1]
}

stochastic_apply_event <- function(state_vector, event, age_structure, compartments) {
  from_index <- stochastic_state_index(event$from, event$age_index, age_structure, compartments)
  to_index <- stochastic_state_index(event$to, event$age_index, age_structure, compartments)

  if (state_vector[from_index] < 1) {
    stop("selected stochastic event would make the state negative.", call. = FALSE)
  }

  state_vector[from_index] <- state_vector[from_index] - 1
  state_vector[to_index] <- state_vector[to_index] + 1
  state_vector
}

stochastic_state_index <- function(compartment, age_index, age_structure, compartments) {
  compartment_index <- match(compartment, compartments)
  ((compartment_index - 1L) * age_structure$n_age_groups) + age_index
}

stochastic_event_log <- function(event_rows) {
  if (length(event_rows) == 0) {
    return(data.frame(
      time = numeric(),
      event = character(),
      transition_id = character(),
      age_group = character(),
      from = character(),
      to = character(),
      rate = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  events <- do.call(rbind, event_rows)
  row.names(events) <- NULL
  events
}

prepare_stochastic_cumulative_flows <- function(
  cumulative_flows,
  state_vector,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness
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
  flows <- validate_cumulative_flows(cumulative_flows, rates)

  output_order <- merge(
    flows,
    data.frame(
      age_group = age_structure$age_groups,
      .age_order = seq_along(age_structure$age_groups),
      stringsAsFactors = FALSE
    ),
    all = TRUE,
    sort = FALSE
  )
  output_order$.flow_order <- match(output_order$cumulative_name, flows$cumulative_name)
  output_order <- output_order[order(output_order$.flow_order, output_order$.age_order), ]
  row.names(output_order) <- NULL

  list(
    flows = flows,
    output_order = output_order[, c(
      "cumulative_name",
      "transition_id",
      "from",
      "to",
      "age_group"
    )]
  )
}

stochastic_cumulative_output <- function(events, times, cumulative_spec) {
  output_order <- cumulative_spec$output_order
  rows <- vector("list", length(times))

  for (i in seq_along(times)) {
    cumulative <- output_order
    cumulative$time <- times[i]
    cumulative$value <- numeric(nrow(cumulative))

    for (j in seq_len(nrow(cumulative))) {
      cumulative$value[j] <- sum(
        events$time <= times[i] &
          events$transition_id == cumulative$transition_id[j] &
          events$age_group == cumulative$age_group[j]
      )
    }

    rows[[i]] <- cumulative[, c(
      "time",
      "cumulative_name",
      "transition_id",
      "from",
      "to",
      "age_group",
      "value"
    )]
  }

  cumulative <- do.call(rbind, rows)
  row.names(cumulative) <- NULL
  cumulative
}

validate_stochastic_method <- function(method) {
  if (!is.character(method) || length(method) != 1 || anyNA(method) || method == "") {
    stop("method must be a non-missing character scalar.", call. = FALSE)
  }

  if (method != "gillespie") {
    stop("unsupported stochastic simulation method: ", method, call. = FALSE)
  }

  invisible(method)
}

validate_stochastic_return_events <- function(return_events) {
  if (!is.logical(return_events) || length(return_events) != 1 || anyNA(return_events)) {
    stop("return_events must be a non-missing logical scalar.", call. = FALSE)
  }

  invisible(return_events)
}

validate_stochastic_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(seed))
  }

  if (!is.numeric(seed) || length(seed) != 1 || anyNA(seed) || !is.finite(seed)) {
    stop("seed must be a finite numeric scalar.", call. = FALSE)
  }

  invisible(seed)
}

validate_stochastic_model <- function(model) {
  if (!model$model_type %in% c("SIR", "SEIR", "CompartmentModel")) {
    stop(
      "simulate_stochastic() currently supports SIRModel(), SEIRModel(), and supported CompartmentModel() models.",
      call. = FALSE
    )
  }

  invisible(model)
}

validate_stochastic_no_demography <- function(demographic_process) {
  if (!is.null(demographic_process)) {
    stop(
      "simulate_stochastic() does not include stochastic demography in this milestone; ",
      "demographic_process must be NULL.",
      call. = FALSE
    )
  }

  invisible(demographic_process)
}

validate_stochastic_integer_state <- function(state_vector) {
  if (any(state_vector != floor(state_vector))) {
    stop("initial_state values must be whole-number counts for stochastic simulation.", call. = FALSE)
  }

  invisible(state_vector)
}

validate_stochastic_population <- function(population, state_vector, age_structure, compartments) {
  initial_population <- transition_population_by_age(
    state_vector_to_long(state_vector, age_structure, compartments),
    age_structure,
    compartments
  )

  if (is.null(population)) {
    population <- initial_population
  } else {
    population <- validate_force_vector(
      population,
      expected_length = age_structure$n_age_groups,
      name = "population",
      allow_zero = FALSE
    )
    if (!isTRUE(all.equal(as.numeric(population), as.numeric(initial_population)))) {
      stop("population must match initial age-specific state totals.", call. = FALSE)
    }
  }

  validate_positive_age_populations(population, age_structure)
  as.numeric(population)
}

with_stochastic_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}
