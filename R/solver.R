integrate_state_trajectory <- function(initial_state,
                                       times,
                                       method,
                                       derivative,
                                       output,
                                       non_negative = NULL,
                                       tcrit = NULL,
                                       desolve_error = NULL) {
  method <- validate_simulation_method(method)

  if (method == "euler") {
    return(integrate_state_trajectory_euler(
      initial_state = initial_state,
      times = times,
      derivative = derivative,
      output = output,
      non_negative = non_negative
    ))
  }

  integrate_state_trajectory_desolve(
    initial_state = initial_state,
    times = times,
    derivative = derivative,
    output = output,
    tcrit = tcrit,
    desolve_error = desolve_error
  )
}

integrate_state_trajectory_values <- function(initial_state,
                                              times,
                                              method,
                                              derivative,
                                              non_negative = NULL,
                                              tcrit = NULL,
                                              desolve_error = NULL) {
  method <- validate_simulation_method(method)

  if (method == "euler") {
    return(integrate_state_trajectory_values_euler(
      initial_state = initial_state,
      times = times,
      derivative = derivative,
      non_negative = non_negative
    ))
  }

  integrate_state_trajectory_values_desolve(
    initial_state = initial_state,
    times = times,
    derivative = derivative,
    tcrit = tcrit,
    desolve_error = desolve_error
  )
}

integrate_state_trajectory_euler <- function(initial_state,
                                             times,
                                             derivative,
                                             output,
                                             non_negative = NULL) {
  values <- vector("list", length(times))
  values[[1]] <- output(initial_state, times[1])

  current_state <- initial_state
  for (i in seq_len(length(times) - 1)) {
    dt <- times[i + 1] - times[i]
    next_state <- as.numeric(current_state) + dt * as.numeric(derivative(times[i], current_state))

    if (!is.null(non_negative)) {
      non_negative(next_state, times[i + 1])
    }

    current_state <- next_state
    values[[i + 1]] <- output(current_state, times[i + 1])
  }

  bind_integrated_outputs(values)
}

integrate_state_trajectory_values_euler <- function(initial_state,
                                                    times,
                                                    derivative,
                                                    non_negative = NULL) {
  trajectory <- matrix(NA_real_, nrow = length(times), ncol = length(initial_state))
  trajectory[1, ] <- as.numeric(initial_state)

  current_state <- as.numeric(initial_state)
  for (i in seq_len(length(times) - 1)) {
    dt <- times[i + 1] - times[i]
    next_state <- current_state + dt * as.numeric(derivative(times[i], current_state))

    if (!is.null(non_negative)) {
      non_negative(next_state, times[i + 1])
    }

    current_state <- next_state
    trajectory[i + 1, ] <- current_state
  }

  trajectory
}

integrate_state_trajectory_desolve <- function(initial_state,
                                               times,
                                               derivative,
                                               output,
                                               tcrit = NULL,
                                               desolve_error = NULL) {
  if (!desolve_is_available()) {
    if (is.null(desolve_error)) {
      desolve_error <- paste(
        "method = \"deSolve\" requires the deSolve package.",
        "Install deSolve or use method = \"euler\"."
      )
    }
    stop(desolve_error, call. = FALSE)
  }

  critical_times <- desolve_tcrit(times, tcrit)
  solved <- if (length(critical_times) == 0) {
    integrate_state_trajectory_desolve_interval(
      initial_state = initial_state,
      times = times,
      derivative = derivative
    )
  } else {
    integrate_state_trajectory_desolve_segmented(
      initial_state = initial_state,
      times = times,
      derivative = derivative,
      critical_times = critical_times
    )
  }

  values <- vector("list", length(times))
  state_columns <- seq_len(length(initial_state)) + 1
  for (i in seq_along(times)) {
    values[[i]] <- output(as.numeric(solved[i, state_columns]), solved[i, "time"])
  }

  bind_integrated_outputs(values)
}

integrate_state_trajectory_values_desolve <- function(initial_state,
                                                      times,
                                                      derivative,
                                                      tcrit = NULL,
                                                      desolve_error = NULL) {
  if (!desolve_is_available()) {
    if (is.null(desolve_error)) {
      desolve_error <- paste(
        "method = \"deSolve\" requires the deSolve package.",
        "Install deSolve or use method = \"euler\"."
      )
    }
    stop(desolve_error, call. = FALSE)
  }

  critical_times <- desolve_tcrit(times, tcrit)
  solved <- if (length(critical_times) == 0) {
    integrate_state_trajectory_desolve_interval(
      initial_state = initial_state,
      times = times,
      derivative = derivative
    )
  } else {
    integrate_state_trajectory_desolve_segmented(
      initial_state = initial_state,
      times = times,
      derivative = derivative,
      critical_times = critical_times
    )
  }

  solved[, -1, drop = FALSE]
}

integrate_state_trajectory_desolve_interval <- function(initial_state,
                                                        times,
                                                        derivative,
                                                        tcrit = NULL) {
  ode_args <- list(
    y = as.numeric(initial_state),
    times = times,
    func = function(time, state, parms) {
      list(as.numeric(derivative(time, state)))
    },
    parms = NULL
  )
  if (!is.null(tcrit)) {
    ode_args$tcrit <- tcrit
  }

  do.call(deSolve::ode, ode_args)
}

integrate_state_trajectory_desolve_segmented <- function(initial_state,
                                                         times,
                                                         derivative,
                                                         critical_times) {
  segment_breaks <- unique(c(min(times), critical_times))
  segment_breaks <- sort(segment_breaks)
  current_state <- as.numeric(initial_state)
  solved_segments <- list()

  for (segment_index in seq_len(length(segment_breaks) - 1)) {
    segment_start <- segment_breaks[segment_index]
    segment_end <- segment_breaks[segment_index + 1]
    segment_times <- times[times >= segment_start & times <= segment_end]
    segment_times <- sort(unique(c(segment_start, segment_times, segment_end)))

    solved <- integrate_state_trajectory_desolve_interval(
      initial_state = current_state,
      times = segment_times,
      derivative = derivative,
      tcrit = segment_end
    )

    if (segment_index > 1) {
      solved <- solved[solved[, "time"] != segment_start, , drop = FALSE]
    }
    solved_segments[[segment_index]] <- solved
    current_state <- as.numeric(solved[nrow(solved), -1])
  }

  solved <- do.call(rbind, solved_segments)
  requested_rows <- match(times, solved[, "time"])
  if (anyNA(requested_rows)) {
    stop("deSolve segmented integration did not return all requested output times.", call. = FALSE)
  }
  solved[requested_rows, , drop = FALSE]
}

state_trajectory_to_data_frame <- function(state_matrix, times, state_template) {
  if (nrow(state_matrix) != length(times)) {
    stop("state_matrix row count must match length(times).", call. = FALSE)
  }

  state_count <- ncol(state_matrix)
  data.frame(
    time = rep(times, each = state_count),
    compartment = rep(state_template$compartment, times = length(times)),
    age_group = rep(state_template$age_group, times = length(times)),
    value = as.numeric(t(state_matrix)),
    stringsAsFactors = FALSE
  )
}

bind_integrated_outputs <- function(values) {
  if (length(values) > 0 && is.list(values[[1]]) && !is.data.frame(values[[1]])) {
    result <- lapply(names(values[[1]]), function(component) {
      component_values <- lapply(values, `[[`, component)
      bound <- do.call(rbind, component_values)
      row.names(bound) <- NULL
      bound
    })
    names(result) <- names(values[[1]])
    return(result)
  }

  result <- do.call(rbind, values)
  row.names(result) <- NULL
  result
}

desolve_tcrit <- function(times, tcrit = NULL) {
  critical_times <- unique(c(tcrit, max(times)))
  critical_times <- critical_times[is.finite(critical_times)]
  critical_times <- critical_times[critical_times > min(times) & critical_times <= max(times)]
  sort(critical_times)
}

desolve_schedule_tcrit <- function(process, times) {
  schedule_times <- unique(c(
    process$fertility_schedule$times,
    process$mortality_schedule$times,
    process$migration_schedule$times
  ))
  desolve_tcrit(times, schedule_times)
}

desolve_is_available <- function() {
  requireNamespace("deSolve", quietly = TRUE)
}
