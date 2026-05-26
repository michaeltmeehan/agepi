#' Simulate a deterministic SIR or SEIR model
#'
#' Runs a simple deterministic prototype simulation. For infection-only runs,
#' deSolve is used by default when the suggested `deSolve` package is
#' installed; otherwise the simulation falls back to explicit Euler time steps.
#' Coupled demographic runs keep the existing Euler default because exact
#' schedule lookup is the default demographic policy. Set `method = "deSolve"`
#' or `method = "euler"` explicitly to choose a backend.
#'
#' The initial state may be supplied either as long-form state data or as a
#' numeric vector in the existing compartment-major, age-group-minor ordering.
#' Names on numeric state vectors are ignored, matching the state-mapping
#' helpers. The initial state is included in the returned output.
#'
#' Supported combinations are:
#'
#' - SIR + Euler, infection-only;
#' - SEIR + Euler, infection-only;
#' - SIR + deSolve, infection-only;
#' - SEIR + deSolve, infection-only;
#' - SIR + demography + Euler;
#' - SEIR + demography + Euler;
#' - SIR + demography + deSolve;
#' - SEIR + demography + deSolve.
#'
#' Both solver backends use the same transition-rate and derivative pathway:
#' static contact matrix, static `beta`, static susceptibility, and static
#' infectiousness. The deSolve backend also supports demographic coupling via
#' the existing demographic derivative path.
#'
#' A first-pass SIR/SEIR-demography coupling is available via
#' `demographic_process`. In this coupled mode, births enter the youngest
#' susceptible age group, background mortality applies independently to every
#' disease compartment, ageing moves each compartment independently using the
#' process [AgeingOperator()] convention, and net migration is allocated using
#' `migration_policy`. Fertility and migration-rate exposure use the current
#' total infection-state population by age group: `S + I + R` for SIR and
#' `S + E + I + R` for SEIR. This is not WPP projection matching, and the
#' default `S`-only migration rule is an allocation convention for age-total
#' net migration inputs rather than a mechanistic movement model. For adaptive
#' deSolve runs, `time_policy = "linear"` is generally recommended for
#' demographic schedules; `time_policy = "step"` gives piecewise-constant rates,
#' and `time_policy = "exact"` will usually fail unless every solver evaluation
#' time is present in the schedules.
#'
#' Euler updates are intentionally not truncated: if a step would produce
#' negative compartment values, simulation stops with an error.
#'
#' @param initial_state Long-form state data frame with columns `compartment`,
#'   `age_group`, and `value`, or numeric state vector.
#' @param times Numeric vector of finite, non-missing, strictly increasing time
#'   points. Must have length at least two.
#' @param model Disease model. `SIRModel()` and `SEIRModel()` output are
#'   supported for infection-only deterministic simulation.
#' @param age_structure Valid age structure.
#' @param contact_matrix Numeric contact matrix with rows as recipient age
#'   groups and columns as source age groups.
#' @param beta Non-negative finite transmission scaling parameter.
#' @param susceptibility Optional non-negative numeric vector by recipient age
#'   group.
#' @param infectiousness Optional non-negative numeric vector by source age
#'   group.
#' @param method Simulation method. `NULL` selects `"deSolve"` for
#'   infection-only runs when available and `"euler"` otherwise; coupled
#'   demographic runs preserve the existing Euler default. `"deSolve"` and
#'   `"ode"` request the optional `deSolve::ode()` backend. `"euler"` requests
#'   explicit Euler time steps.
#' @param demographic_process Optional [DemographicProcess()] object for
#'   first-pass deterministic SIR/SEIR/generic-demography coupling. Defaults to
#'   `NULL`, which preserves infection-only simulation.
#' @param time_policy Demographic schedule lookup policy used only when
#'   `demographic_process` is supplied. `"exact"` requires exact schedule times;
#'   `"step"` uses left-continuous interval-start lookup; `"linear"`
#'   interpolates rate-like demographic schedules through the same demographic
#'   derivative path.
#' @param migration_policy Net migration allocation policy used only when
#'   `demographic_process` is supplied. `"susceptible"` preserves the default
#'   behaviour by applying all net migration to `S`; `"proportional"` allocates
#'   net migration across compartments by current age-specific compartment
#'   shares; `"error"` allows zero migration but errors if any net migration is
#'   non-zero because age-total migration allocation is ambiguous.
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
  method = NULL,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error")
) {
  method_was_null <- is.null(method)
  method <- validate_simulation_method(method)
  migration_policy <- validate_migration_policy(migration_policy)
  validate_simulation_times(times)
  validate_disease_model(model)
  validate_age_structure(age_structure)
  time_policy <- validate_simulation_demography_inputs(
    demographic_process = demographic_process,
    time_policy = time_policy,
    method = method,
    model = model,
    age_structure = age_structure,
    times = times
  )
  if (method_was_null && !is.null(demographic_process)) {
    method <- "euler"
  }

  state_vector <- simulation_state_to_vector(
    initial_state,
    age_structure,
    model$compartments
  )
  validate_non_negative_simulation_state(state_vector, "initial_state")

  simulate_deterministic_integrated(
    state_vector = state_vector,
    times = times,
    method = method,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    demographic_process = demographic_process,
    time_policy = time_policy,
    migration_policy = migration_policy
  )
}

simulate_deterministic_integrated <- function(
  state_vector,
  times,
  method,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error")
) {
  integrate_state_trajectory(
    initial_state = state_vector,
    times = times,
    method = method,
    derivative = function(time, state) {
      deterministic_derivative(
        state_vector = state,
        time = time,
        model = model,
        age_structure = age_structure,
        contact_matrix = contact_matrix,
        beta = beta,
        susceptibility = susceptibility,
        infectiousness = infectiousness,
        demographic_process = demographic_process,
        time_policy = time_policy,
        migration_policy = migration_policy
      )
    },
    output = function(state, time) {
      simulation_state_output(
        state,
        time = time,
        age_structure = age_structure,
        compartments = model$compartments
      )
    },
    non_negative = validate_non_negative_euler_state,
    tcrit = if (is.null(demographic_process)) NULL else desolve_schedule_tcrit(demographic_process, times),
    desolve_error = paste(
      "method = \"deSolve\" requires the optional deSolve package.",
      "Install deSolve or use method = \"euler\"."
    )
  )
}

deterministic_derivative <- function(
  state_vector,
  time,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error")
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
  derivative <- rates_to_derivative(
    transition_rate_table = rates,
    compartments = model$compartments,
    age_structure = age_structure
  )
  infection_derivative <- derivative$derivative

  if (is.null(demographic_process)) {
    return(infection_derivative)
  }

  infection_derivative + compartment_demographic_derivative(
    state_vector = as.numeric(state_vector),
    time = time,
    model = model,
    age_structure = age_structure,
    demographic_process = demographic_process,
    time_policy = time_policy,
    migration_policy = migration_policy
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
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error")
) {
  simulate_deterministic_integrated(
    state_vector = state_vector,
    times = times,
    method = "euler",
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    demographic_process = demographic_process,
    time_policy = time_policy,
    migration_policy = migration_policy
  )
}

deterministic_euler_derivative <- function(
  state_vector,
  time,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error")
) {
  deterministic_derivative(
    state_vector = state_vector,
    time = time,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    demographic_process = demographic_process,
    time_policy = time_policy,
    migration_policy = migration_policy
  )
}

compartment_demographic_derivative <- function(
  state_vector,
  time,
  model,
  age_structure,
  demographic_process,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error")
) {
  time_policy <- validate_demographic_time_policy(time_policy)
  migration_policy <- validate_migration_policy(migration_policy)
  validate_disease_model(model)
  validate_age_structure(age_structure)
  validate_demographic_process(demographic_process)

  birth_compartment <- demographic_birth_compartment(model)
  migration_compartment <- demographic_migration_compartment(model)

  if (is.null(birth_compartment)) {
    stop("demographic_process coupling requires a birth_compartment or an S compartment.", call. = FALSE)
  }
  if (migration_policy == "susceptible" && is.null(migration_compartment)) {
    stop("migration_policy = \"susceptible\" requires a migration_compartment or an S compartment.", call. = FALSE)
  }

  state_long <- state_vector_to_long(state_vector, age_structure, model$compartments)
  age_groups <- age_structure$age_groups
  n_age_groups <- age_structure$n_age_groups

  compartment_indices <- lapply(seq_along(model$compartments), function(compartment_position) {
    ((compartment_position - 1) * n_age_groups) + seq_len(n_age_groups)
  })
  names(compartment_indices) <- model$compartments

  population <- numeric(n_age_groups)
  for (compartment in model$compartments) {
    population <- population +
      transition_compartment_values(state_long, age_structure, compartment)
  }

  derivative <- numeric(length(state_vector))
  birth_index <- compartment_indices[[birth_compartment]]

  ageing <- demographic_process$ageing_operator
  for (index in compartment_indices) {
    compartment_state <- as.numeric(state_vector[index])
    ageing_out <- ageing$departure_rate * compartment_state
    derivative[index] <- derivative[index] - ageing_out

    has_destination <- !is.na(ageing$destination_index)
    derivative[index[ageing$destination_index[has_destination]]] <-
      derivative[index[ageing$destination_index[has_destination]]] +
      ageing_out[has_destination]
  }

  fertility_rates <- fertility_rates_at(
    demographic_process$fertility_schedule,
    time,
    age_groups,
    time_policy
  )
  derivative[birth_index[1]] <- derivative[birth_index[1]] + sum(fertility_rates * population)

  mortality_rates <- mortality_rates_at(
    demographic_process$mortality_schedule,
    time,
    age_groups,
    time_policy
  )
  for (index in compartment_indices) {
    derivative[index] <- derivative[index] - mortality_rates * as.numeric(state_vector[index])
  }

  migration <- migration_values_at(
    demographic_process$migration_schedule,
    time,
    population,
    age_groups,
    time_policy
  )
  derivative <- derivative + compartment_migration_derivative(
    migration = migration,
    state_vector = state_vector,
    population = population,
    compartment_indices = compartment_indices,
    migration_compartment = migration_compartment,
    migration_policy = migration_policy
  )

  if (any(!is.finite(derivative))) {
    stop("compartment_demographic_derivative result must contain only finite values.", call. = FALSE)
  }

  derivative
}

compartment_migration_derivative <- function(
  migration,
  state_vector,
  population,
  compartment_indices,
  migration_compartment,
  migration_policy = c("susceptible", "proportional", "error")
) {
  migration_policy <- validate_migration_policy(migration_policy)
  migration_derivative <- numeric(length(state_vector))

  if (migration_policy == "susceptible") {
    migration_derivative[compartment_indices[[migration_compartment]]] <- migration
    return(migration_derivative)
  }

  non_zero_migration <- migration != 0
  if (migration_policy == "error") {
    if (any(non_zero_migration)) {
      stop(
        "migration_policy = \"error\" does not allow non-zero net migration; ",
        "age-total migration allocation across disease compartments is ambiguous.",
        call. = FALSE
      )
    }
    return(migration_derivative)
  }

  zero_population_non_zero_migration <- population == 0 & non_zero_migration
  if (any(zero_population_non_zero_migration)) {
    stop(
      "migration_policy = \"proportional\" cannot allocate non-zero net migration ",
      "when an age-specific total population is zero.",
      call. = FALSE
    )
  }

  positive_population <- population > 0
  for (index in compartment_indices) {
    shares <- numeric(length(population))
    shares[positive_population] <- as.numeric(state_vector[index][positive_population]) /
      population[positive_population]
    migration_derivative[index] <- migration * shares
  }

  migration_derivative
}

simulate_deterministic_desolve <- function(
  state_vector,
  times,
  model,
  age_structure,
  contact_matrix,
  beta,
  susceptibility,
  infectiousness,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error")
) {
  simulate_deterministic_integrated(
    state_vector = state_vector,
    times = times,
    method = "deSolve",
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    demographic_process = demographic_process,
    time_policy = time_policy,
    migration_policy = migration_policy
  )
}

validate_simulation_method <- function(method) {
  if (is.null(method)) {
    if (desolve_is_available()) {
      return("deSolve")
    }
    return("euler")
  }

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

validate_migration_policy <- function(migration_policy) {
  if (!is.character(migration_policy) || anyNA(migration_policy) || any(migration_policy == "")) {
    stop("migration_policy must be non-missing character value(s).", call. = FALSE)
  }

  tryCatch(
    match.arg(migration_policy, c("susceptible", "proportional", "error")),
    error = function(e) {
      stop("unsupported migration_policy: ", paste(migration_policy, collapse = ", "), call. = FALSE)
    }
  )
}

validate_simulation_demography_inputs <- function(
  demographic_process,
  time_policy,
  method,
  model,
  age_structure,
  times
) {
  if (is.null(demographic_process)) {
    return(validate_demographic_time_policy(time_policy))
  }

  time_policy <- validate_demographic_time_policy(time_policy)
  validate_demographic_process(demographic_process)

  validate_same_age_structure(
    age_structure,
    demographic_process$age_structure,
    "demographic_process"
  )
  validate_demography_schedule_coverage(
    demographic_process,
    times,
    time_policy,
    include_output_times = TRUE
  )

  time_policy
}

demographic_birth_compartment <- function(model) {
  if (model$model_type == "CompartmentModel") {
    return(model$birth_compartment)
  }

  if ("S" %in% model$compartments) {
    return("S")
  }

  NULL
}

demographic_migration_compartment <- function(model) {
  if (model$model_type == "CompartmentModel") {
    return(model$migration_compartment)
  }

  if ("S" %in% model$compartments) {
    return("S")
  }

  NULL
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
