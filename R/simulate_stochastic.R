#' Simulate a stochastic SIR or SEIR model with Gillespie's direct method
#'
#' Runs a narrow stochastic simulation for age-structured SIR and SEIR models:
#' fixed population, no demography, no ageing, no fertility, no mortality, no
#' migration, and no tau-leaping.
#'
#' The simulator uses Gillespie's continuous-time direct method. Infection
#' propensities reuse the package force-of-infection convention:
#' rows of `contact_matrix` are recipient age groups and columns are infectious
#' source age groups. SIR infection events move `S -> I` with propensity
#' `lambda * S`, and recovery events move `I -> R` with propensity `gamma * I`.
#' SEIR infection events move `S -> E` with propensity `lambda * S`,
#' progression events move `E -> I` with propensity `sigma * E`, and recovery
#' events move `I -> R` with propensity `gamma * I`.
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
#' @param model Disease model. `SIRModel()` and `SEIRModel()` outputs are
#'   currently supported.
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
#'
#' @return If `return_events = FALSE`, a data frame with columns `time`,
#'   `compartment`, `age_group`, and `value`, ordered by time outermost,
#'   compartment next, and age group innermost. If `return_events = TRUE`, a
#'   list with `trajectory` and an event log containing `time`, `event`,
#'   `age_group`, `from`, `to`, and `rate`.
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
  demographic_process = NULL
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
      return_events = return_events
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
  return_events
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

    if (return_events) {
      event_count <- event_count + 1L
      event_rows[[event_count]] <- data.frame(
        time = current_time,
        event = event$event,
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

  if (!return_events) {
    return(trajectory)
  }

  events <- stochastic_event_log(event_rows)
  list(trajectory = trajectory, events = events)
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
  if (identical(model$model_type, "SEIR")) {
    return(stochastic_seir_propensities(
      state_vector = state_vector,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      susceptibility = susceptibility,
      infectiousness = infectiousness,
      population = population
    ))
  }

  stochastic_sir_propensities(
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

stochastic_sir_propensities <- function(
  state_vector,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  population
) {
  n_age_groups <- age_structure$n_age_groups
  S <- state_vector[seq_len(n_age_groups)]
  I <- state_vector[n_age_groups + seq_len(n_age_groups)]

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
    event = rep(c("infection", "recovery"), each = n_age_groups),
    age_group = rep(age_structure$age_groups, times = 2),
    age_index = rep(seq_len(n_age_groups), times = 2),
    from = rep(c("S", "I"), each = n_age_groups),
    to = rep(c("I", "R"), each = n_age_groups),
    rate = c(infection_rates, recovery_rates),
    stringsAsFactors = FALSE
  )
}

stochastic_seir_propensities <- function(
  state_vector,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  population
) {
  n_age_groups <- age_structure$n_age_groups
  S <- state_vector[seq_len(n_age_groups)]
  E <- state_vector[n_age_groups + seq_len(n_age_groups)]
  I <- state_vector[(2L * n_age_groups) + seq_len(n_age_groups)]

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
  progression_rates <- model$sigma * E
  recovery_rates <- model$gamma * I

  data.frame(
    event = rep(c("infection", "progression", "recovery"), each = n_age_groups),
    age_group = rep(age_structure$age_groups, times = 3),
    age_index = rep(seq_len(n_age_groups), times = 3),
    from = rep(c("S", "E", "I"), each = n_age_groups),
    to = rep(c("E", "I", "R"), each = n_age_groups),
    rate = c(infection_rates, progression_rates, recovery_rates),
    stringsAsFactors = FALSE
  )
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
  if (!identical(model$model_type, "SIR") && !identical(model$model_type, "SEIR")) {
    stop("simulate_stochastic() currently supports only SIRModel() and SEIRModel() models.", call. = FALSE)
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
