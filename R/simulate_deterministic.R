#' Simulate a deterministic SIR or SEIR model
#'
#' Runs a simple deterministic prototype simulation. The package uses the
#' `deSolve` backend for ODE integration, with explicit Euler still available as
#' an alternative solver choice. Coupled demographic runs with the default
#' derivative ageing policy keep the existing Euler default because exact
#' schedule lookup is the default demographic policy. Annual-cohort operator
#' splitting follows the ordinary method default unless `method` is set
#' explicitly. Set `method = "deSolve"` or `method = "euler"` explicitly to
#' choose a backend.
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
#' - SEIR + demography + deSolve;
#' - deterministic annual-cohort operator splitting for SIR, SEIR, and generic
#'   deterministic compartment models.
#'
#' Both solver backends use the same transition-rate and derivative pathway:
#' static contact matrix and a resolved transmission `beta` drawn from the
#' simulation override or model default. Model-level infection-transition
#' susceptibility and infectiousness weights are used automatically. The
#' deSolve backend also supports demographic coupling via the existing
#' demographic derivative path.
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
#' net migration inputs rather than a mechanistic movement model. The
#' `proportional` policy allocates net migration across compartments by current
#' age-specific shares, and the `error` policy rejects non-zero net migration.
#' For adaptive deSolve runs, `time_policy = "linear"` is generally recommended for
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
#' @param beta Optional non-negative finite transmission scaling parameter.
#'   When omitted, infectious models use `model$beta` if present and error if
#'   no default is available. Models without infection transitions do not
#'   require beta; an explicit beta is validated and otherwise ignored.
#' @param method Simulation method. `NULL` selects `"deSolve"` for infection-
#'   only runs and `"euler"` for coupled demographic runs with
#'   `ageing_policy = "exponential"`. `"deSolve"` and `"ode"` request the
#'   `deSolve::ode()` backend. `"euler"` requests explicit Euler time steps.
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
#' @param ageing_policy Demographic ageing implementation. `"exponential"`
#'   preserves the existing derivative-based coupling. `"annual_cohort"` runs
#'   deterministic infection dynamics over each 1-year interval, then applies
#'   `annual_cohort_demographic_step()` compartment-wise at the annual boundary.
#'   Annual-cohort coupling is deterministic-only, opt-in, and requires annual
#'   output times. If `demographic_process` is supplied, it must use a complete
#'   1-year age grid ending in an open-ended age group.
#' @param output_age_structure Optional reporting age structure. When supplied,
#'   deterministic trajectories and cumulative-flow outputs are aggregated from
#'   the simulation age grid to this reporting grid using `AgeGridMapping()` and
#'   the existing count-preserving output aggregation helpers.
#' @param cumulative_flows Optional named list or data frame specifying disease
#'   transition flows to track as auxiliary cumulative state variables. Named
#'   list entries must contain `from` and `to` fields, each either a character
#'   scalar or equal-length character vector; vector entries are summed into one
#'   cumulative counter per name and age group. Data frames must contain `name`,
#'   `from`, and `to` columns. Cumulative flows track disease transition flows
#'   only, and are supported with or without `demographic_process`.
#'
#' @return Data frame with columns `time`, `compartment`, `age_group`, and
#'   `value`, ordered by time outermost, compartment next, and age group
#'   innermost when `cumulative_flows = NULL`. When `cumulative_flows` is
#'   supplied, a list is returned with `trajectory` containing that ordinary
#'   compartment output and `cumulative` containing columns `time`,
#'   `cumulative_name`, `transition_id`, `from`, `to`, `age_group`, and `value`.
#' @export
simulate_deterministic <- function(
  initial_state,
  times,
  model,
  age_structure,
  contact_matrix,
  beta = NULL,
  method = NULL,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error"),
  ageing_policy = c("exponential", "annual_cohort"),
  output_age_structure = NULL,
  cumulative_flows = NULL
) {
  method_was_null <- is.null(method)
  method <- validate_simulation_method(method)
  migration_policy <- validate_migration_policy(migration_policy)
  ageing_policy <- validate_demography_ageing_policy(ageing_policy)
  validate_simulation_times(times)
  validate_disease_model(model)
  validate_age_structure(age_structure)
  if (!is.null(output_age_structure)) {
    validate_age_structure(output_age_structure)
  }
  if (model_has_infection_process(model)) {
    beta <- resolve_transmission_beta(model, beta)
  } else if (!is.null(beta)) {
    beta <- validate_force_beta(beta)
  }
  time_policy <- validate_simulation_demography_inputs(
    demographic_process = demographic_process,
    time_policy = time_policy,
    method = method,
    model = model,
    age_structure = age_structure,
    times = times,
    ageing_policy = ageing_policy
  )
  if (method_was_null && !is.null(demographic_process) && ageing_policy == "exponential") {
    method <- "euler"
  }

  transition_context <- prepare_transition_rate_context_validated(
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta
  )

  state_vector <- simulation_state_to_vector(
    initial_state,
    age_structure,
    model$compartments
  )
  validate_non_negative_simulation_state(state_vector, "initial_state")

  cumulative_spec <- NULL
  if (!is.null(cumulative_flows)) {
    cumulative_spec <- prepare_deterministic_cumulative_flows(
      cumulative_flows = cumulative_flows,
      state_vector = state_vector,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta
    )
  }

  output <- if (ageing_policy == "annual_cohort") {
    simulate_deterministic_annual_cohort_split(
      state_vector = state_vector,
      times = times,
      method = method,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      demographic_process = demographic_process,
      time_policy = time_policy,
      migration_policy = migration_policy,
      cumulative_spec = cumulative_spec,
      transition_context = transition_context
    )
  } else {
    simulate_deterministic_integrated(
      state_vector = state_vector,
      times = times,
      method = method,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      demographic_process = demographic_process,
      time_policy = time_policy,
      migration_policy = migration_policy,
      cumulative_spec = cumulative_spec,
      transition_context = transition_context
    )
  }

  aggregate_deterministic_output_if_requested(output, age_structure, output_age_structure)
}

simulate_deterministic_annual_cohort_split <- function(
  state_vector,
  times,
  method,
  model,
  age_structure,
  contact_matrix,
  beta,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error"),
  cumulative_spec = NULL,
  transition_context = NULL
) {
  time_policy <- validate_demographic_time_policy(time_policy)
  migration_policy <- validate_migration_policy(migration_policy)
  validate_annual_cohort_simulation_times(times)

  ordinary_state_length <- length(state_vector)
  current_state <- state_vector
  if (!is.null(cumulative_spec)) {
    current_state <- c(current_state, numeric(nrow(cumulative_spec$state_order)))
  }

  output <- vector("list", length(times))
  output[[1]] <- deterministic_augmented_state_output(
    state = current_state,
    ordinary_state_length = ordinary_state_length,
    time = times[1],
    age_structure = age_structure,
    compartments = model$compartments,
    cumulative_spec = cumulative_spec
  )

  for (time_index in seq_len(length(times) - 1)) {
    interval_times <- times[time_index + 0:1]
    current_state <- deterministic_annual_interval_epidemic_step(
      state = current_state,
      times = interval_times,
      method = method,
      ordinary_state_length = ordinary_state_length,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      cumulative_spec = cumulative_spec,
      transition_context = transition_context
    )

    current_state <- truncate_near_zero_state(current_state)
    if (!is.null(demographic_process)) {
      current_state[seq_len(ordinary_state_length)] <- compartment_annual_cohort_demographic_step(
        state_vector = current_state[seq_len(ordinary_state_length)],
        time = times[time_index],
        model = model,
        age_structure = age_structure,
        demographic_process = demographic_process,
        time_policy = time_policy,
        migration_policy = migration_policy
      )
    }

    validate_non_negative_simulation_state(current_state, "annual_cohort state")
    output[[time_index + 1]] <- deterministic_augmented_state_output(
      state = current_state,
      ordinary_state_length = ordinary_state_length,
      time = times[time_index + 1],
      age_structure = age_structure,
      compartments = model$compartments,
      cumulative_spec = cumulative_spec
    )
  }

  if (is.null(cumulative_spec)) {
    return(bind_integrated_outputs(output))
  }

  bind_integrated_outputs(output)
}

deterministic_annual_interval_epidemic_step <- function(
  state,
  times,
  method,
  ordinary_state_length,
  model,
  age_structure,
  contact_matrix,
  beta,
  cumulative_spec = NULL,
  transition_context = NULL
) {
  if (is.null(cumulative_spec) && deterministic_epidemic_dynamics_disabled(model, beta)) {
    return(as.numeric(state))
  }

  integrated <- integrate_state_trajectory(
    initial_state = state,
    times = times,
    method = method,
    derivative = function(time, state) {
      deterministic_derivative_augmented(
        state_vector = state,
        ordinary_state_length = ordinary_state_length,
        time = time,
        model = model,
        age_structure = age_structure,
        contact_matrix = contact_matrix,
        beta = beta,
        demographic_process = NULL,
        cumulative_spec = cumulative_spec,
        transition_context = transition_context
      )
    },
    output = function(state, time) {
      data.frame(
        time = unname(time),
        state_index = seq_len(length(state)),
        value = as.numeric(state),
        stringsAsFactors = FALSE
      )
    },
    non_negative = validate_non_negative_euler_state,
    desolve_error = paste(
      "ageing_policy = \"annual_cohort\" requires the deSolve package",
      "for method = \"deSolve\". Install deSolve or use method = \"euler\"."
    )
  )

  as.numeric(integrated$value[integrated$time == times[2]])
}

deterministic_epidemic_dynamics_disabled <- function(model, beta) {
  validate_disease_model(model)
  if (!model_has_infection_process(model)) {
    return(nrow(model$transitions) == 0)
  }

  beta_value <- resolve_transmission_beta(model, beta)
  beta_is_zero <- is.numeric(beta_value) && length(beta_value) == 1 && !is.na(beta_value) && is.finite(beta_value) && beta_value == 0

  if (model$model_type == "SIR") {
    return(beta_is_zero && model$gamma == 0)
  }

  if (model$model_type == "SEIR") {
    return(beta_is_zero && model$sigma == 0 && model$gamma == 0)
  }

  if (model$model_type == "CompartmentModel") {
    return(beta_is_zero && nrow(model$infection_transitions) == 0 && nrow(model$transitions) == 0)
  }

  FALSE
}

compartment_annual_cohort_demographic_step <- function(
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
  validate_same_age_structure(age_structure, demographic_process$age_structure, "demographic_process")
  validate_simulate_demography_annual_cohort_age_structure(age_structure)

  birth_compartment <- demographic_birth_compartment(model)
  migration_compartment <- demographic_migration_compartment(model)
  if (is.null(birth_compartment)) {
    stop("annual cohort demographic coupling requires a birth_compartment or an S compartment.", call. = FALSE)
  }
  if (migration_policy == "susceptible" && is.null(migration_compartment)) {
    stop("migration_policy = \"susceptible\" requires a migration_compartment or an S compartment.", call. = FALSE)
  }

  age_groups <- age_structure$age_groups
  n_age_groups <- age_structure$n_age_groups
  compartments <- model$compartments
  compartment_indices <- lapply(seq_along(compartments), function(compartment_position) {
    ((compartment_position - 1) * n_age_groups) + seq_len(n_age_groups)
  })
  names(compartment_indices) <- compartments

  state_by_compartment <- lapply(compartment_indices, function(index) {
    stats::setNames(as.numeric(state_vector[index]), age_groups)
  })
  population <- Reduce(`+`, state_by_compartment)

  fertility <- stats::setNames(
    fertility_rates_at(demographic_process$fertility_schedule, time, age_groups, time_policy),
    age_groups
  )
  mortality <- stats::setNames(
    mortality_rates_at(demographic_process$mortality_schedule, time, age_groups, time_policy),
    age_groups
  )
  migration <- annual_cohort_migration_values_at(
    demographic_process$migration_schedule,
    time = time,
    age_groups = age_groups,
    time_policy = time_policy
  )
  zero <- stats::setNames(numeric(n_age_groups), age_groups)

  stepped_by_compartment <- lapply(state_by_compartment, function(compartment_state) {
    step <- annual_cohort_demographic_step(
      population = compartment_state,
      age_structure = age_structure,
      fertility = zero,
      mortality = mortality,
      migration = zero,
      migration_type = "count"
    )
    stats::setNames(step$population, age_groups)
  })

  births <- sum(fertility * demographic_process$fertility_exposure_fraction * population)
  stepped_by_compartment[[birth_compartment]][1] <-
    stepped_by_compartment[[birth_compartment]][1] + births

  total_after_births <- Reduce(`+`, stepped_by_compartment)
  total_step <- annual_cohort_demographic_step(
    population = population,
    age_structure = age_structure,
    fertility = fertility,
    mortality = mortality,
    migration = migration$values,
    fertility_exposure_fraction = demographic_process$fertility_exposure_fraction,
    migration_type = migration$type
  )
  net_migration <- stats::setNames(total_step$population - total_after_births, age_groups)

  stepped_by_compartment <- allocate_annual_cohort_migration(
    stepped_by_compartment = stepped_by_compartment,
    net_migration = net_migration,
    total_after_births = total_after_births,
    migration_compartment = migration_compartment,
    migration_policy = migration_policy
  )

  result <- numeric(length(state_vector))
  for (compartment in compartments) {
    result[compartment_indices[[compartment]]] <- stepped_by_compartment[[compartment]]
  }
  validate_non_negative_simulation_state(result, "annual cohort demographic state")
  result
}

allocate_annual_cohort_migration <- function(
  stepped_by_compartment,
  net_migration,
  total_after_births,
  migration_compartment,
  migration_policy = c("susceptible", "proportional", "error")
) {
  migration_policy <- validate_migration_policy(migration_policy)

  if (migration_policy == "susceptible") {
    stepped_by_compartment[[migration_compartment]] <-
      stepped_by_compartment[[migration_compartment]] + net_migration
    return(stepped_by_compartment)
  }

  non_zero_migration <- net_migration != 0
  if (migration_policy == "error") {
    if (any(non_zero_migration)) {
      stop(
        "migration_policy = \"error\" does not allow non-zero net migration; ",
        "age-total migration allocation across disease compartments is ambiguous.",
        call. = FALSE
      )
    }
    return(stepped_by_compartment)
  }

  zero_population_non_zero_migration <- total_after_births == 0 & non_zero_migration
  if (any(zero_population_non_zero_migration)) {
    stop(
      "migration_policy = \"proportional\" cannot allocate non-zero net migration ",
      "when an age-specific total population is zero.",
      call. = FALSE
    )
  }

  positive_population <- total_after_births > 0
  for (compartment in names(stepped_by_compartment)) {
    shares <- numeric(length(total_after_births))
    shares[positive_population] <- stepped_by_compartment[[compartment]][positive_population] /
      total_after_births[positive_population]
    stepped_by_compartment[[compartment]] <- stepped_by_compartment[[compartment]] +
      net_migration * shares
  }

  stepped_by_compartment
}

truncate_near_zero_state <- function(state, tolerance = sqrt(.Machine$double.eps)) {
  state <- as.numeric(state)
  state[state < 0 & state >= -tolerance] <- 0
  state
}

aggregate_deterministic_output_if_requested <- function(output, age_structure, output_age_structure = NULL) {
  if (is.null(output_age_structure)) {
    return(output)
  }
  mapping <- AgeGridMapping(age_structure, output_age_structure, open_ended = "include")

  if (is.list(output) && !is.data.frame(output)) {
    output$trajectory <- aggregate_epidemic_trajectory_age_grid(output$trajectory, mapping)
    output$cumulative <- aggregate_cumulative_flows_age_grid(output$cumulative, mapping)
    return(output)
  }

  aggregate_epidemic_trajectory_age_grid(output, mapping)
}

simulate_deterministic_integrated <- function(
  state_vector,
  times,
  method,
  model,
  age_structure,
  contact_matrix,
  beta,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error"),
  cumulative_spec = NULL,
  transition_context = NULL
) {
  ordinary_state_length <- length(state_vector)
  if (!is.null(cumulative_spec)) {
    state_vector <- c(state_vector, numeric(nrow(cumulative_spec$state_order)))
  }

  integrated <- integrate_state_trajectory(
    initial_state = state_vector,
    times = times,
    method = method,
    derivative = function(time, state) {
      deterministic_derivative_augmented(
        state_vector = state,
        ordinary_state_length = ordinary_state_length,
        time = time,
        model = model,
        age_structure = age_structure,
        contact_matrix = contact_matrix,
        beta = beta,
        demographic_process = demographic_process,
        time_policy = time_policy,
        migration_policy = migration_policy,
        cumulative_spec = cumulative_spec,
        transition_context = transition_context
      )
    },
    output = function(state, time) {
      deterministic_augmented_state_output(
        state = state,
        ordinary_state_length = ordinary_state_length,
        time = time,
        age_structure = age_structure,
        compartments = model$compartments,
        cumulative_spec = cumulative_spec
      )
    },
    non_negative = validate_non_negative_euler_state,
    tcrit = if (is.null(demographic_process)) NULL else desolve_schedule_tcrit(demographic_process, times),
    desolve_error = paste(
      "method = \"deSolve\" requires the deSolve package.",
      "Install deSolve or use method = \"euler\"."
    )
  )

  if (is.null(cumulative_spec)) {
    return(integrated)
  }

  list(
    trajectory = integrated$trajectory,
    cumulative = integrated$cumulative
  )
}

deterministic_derivative_augmented <- function(
  state_vector,
  ordinary_state_length,
  time,
  model,
  age_structure,
  contact_matrix,
  beta,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error"),
  cumulative_spec = NULL,
  transition_context = NULL
) {
  ordinary_state <- as.numeric(state_vector)[seq_len(ordinary_state_length)]
  if (is.null(cumulative_spec)) {
    return(deterministic_derivative(
      state_vector = ordinary_state,
      time = time,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta,
      demographic_process = demographic_process,
      time_policy = time_policy,
      migration_policy = migration_policy,
      transition_context = transition_context
    ))
  }

  if (is.null(transition_context)) {
    rates <- transition_rates(
      state = ordinary_state,
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta
    )
    derivative <- rates_to_derivative(
      transition_rate_table = rates,
      compartments = model$compartments,
      age_structure = age_structure
    )
  } else {
    rates <- transition_rates_from_state_vector(ordinary_state, transition_context)
    derivative <- rates_to_derivative_from_rates(rates, transition_context)
  }
  compartment_derivative <- derivative$derivative
  if (!is.null(demographic_process)) {
    compartment_derivative <- compartment_derivative +
      compartment_demographic_derivative(
        state_vector = ordinary_state,
        time = time,
        model = model,
        age_structure = age_structure,
        demographic_process = demographic_process,
        time_policy = time_policy,
        migration_policy = migration_policy
      )
  }

  c(
    compartment_derivative,
    cumulative_flow_derivative(
      transition_rate_table = rates,
      cumulative_spec = cumulative_spec
    )
  )
}

prepare_deterministic_cumulative_flows <- function(
  cumulative_flows,
  state_vector,
  model,
  age_structure,
  contact_matrix,
  beta
) {
  rates <- transition_rates(
    state = state_vector,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta
  )
  flows <- validate_cumulative_flows(cumulative_flows, rates)
  state_spec <- cumulative_flow_state_spec(flows, age_structure$age_groups)

  list(
    flows = flows,
    state_order = state_spec$state_order,
    state_transitions = state_spec$transition_sets
  )
}

cumulative_flow_derivative <- function(transition_rate_table, cumulative_spec) {
  state_order <- cumulative_spec$state_order
  derivative <- numeric(nrow(state_order))

  for (i in seq_len(nrow(state_order))) {
    transition_ids <- if (is.null(cumulative_spec$state_transitions)) {
      state_order$transition_id[i]
    } else {
      cumulative_spec$state_transitions[[i]]
    }
    matched <- transition_rate_table$transition_id %in% transition_ids &
      transition_rate_table$age_group == state_order$age_group[i]
    if (sum(matched) != length(transition_ids)) {
      stop(
        "cumulative flow '",
        state_order$cumulative_name[i],
        "' did not match the requested transition-rate row(s) for age group ",
        state_order$age_group[i],
        ".",
        call. = FALSE
      )
    }
    derivative[i] <- sum(transition_rate_table$rate[matched])
  }

  derivative
}

deterministic_augmented_state_output <- function(
  state,
  ordinary_state_length,
  time,
  age_structure,
  compartments,
  cumulative_spec = NULL
) {
  ordinary_state <- as.numeric(state)[seq_len(ordinary_state_length)]
  trajectory <- simulation_state_output(
    ordinary_state,
    time = time,
    age_structure = age_structure,
    compartments = compartments
  )

  if (is.null(cumulative_spec)) {
    return(trajectory)
  }

  cumulative <- cumulative_spec$state_order
  cumulative$time <- time
  cumulative$value <- as.numeric(state)[ordinary_state_length + seq_len(nrow(cumulative))]
  cumulative <- cumulative[, c(
    "time",
    "cumulative_name",
    "transition_id",
    "from",
    "to",
    "age_group",
    "value"
  )]

  list(
    trajectory = trajectory,
    cumulative = cumulative
  )
}

deterministic_derivative <- function(
  state_vector,
  time,
  model,
  age_structure,
  contact_matrix,
  beta,
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error"),
  transition_context = NULL
) {
  if (is.null(transition_context)) {
    rates <- transition_rates(
      state = as.numeric(state_vector),
      model = model,
      age_structure = age_structure,
      contact_matrix = contact_matrix,
      beta = beta
    )
    derivative <- rates_to_derivative(
      transition_rate_table = rates,
      compartments = model$compartments,
      age_structure = age_structure
    )
  } else {
    rates <- transition_rates_from_state_vector(as.numeric(state_vector), transition_context)
    derivative <- rates_to_derivative_from_rates(rates, transition_context)
  }
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
  demographic_process = NULL,
  time_policy = c("exact", "step", "linear"),
  migration_policy = c("susceptible", "proportional", "error"),
  transition_context = NULL
) {
  deterministic_derivative(
    state_vector = state_vector,
    time = time,
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = beta,
    demographic_process = demographic_process,
    time_policy = time_policy,
    migration_policy = migration_policy,
    transition_context = transition_context
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
  derivative[birth_index[1]] <- derivative[birth_index[1]] +
    sum(fertility_rates * demographic_process$fertility_exposure_fraction * population)

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
  times,
  ageing_policy = c("exponential", "annual_cohort")
) {
  ageing_policy <- validate_demography_ageing_policy(ageing_policy)
  if (is.null(demographic_process)) {
    if (ageing_policy == "annual_cohort") {
      validate_annual_cohort_simulation_times(times)
    }
    return(validate_demographic_time_policy(time_policy))
  }

  time_policy <- validate_demographic_time_policy(time_policy)
  validate_demographic_process(demographic_process)

  validate_same_age_structure(
    age_structure,
    demographic_process$age_structure,
    "demographic_process"
  )
  if (ageing_policy == "annual_cohort") {
    validate_annual_cohort_simulation_times(times)
    validate_simulate_demography_annual_cohort_age_structure(age_structure)
  }
  validate_demography_schedule_coverage(
    demographic_process,
    times,
    time_policy,
    include_output_times = ageing_policy == "exponential"
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
