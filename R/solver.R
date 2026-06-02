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

integrate_state_trajectory_desolve <- function(initial_state,
                                               times,
                                               derivative,
                                               output,
                                               tcrit = NULL,
                                               desolve_error = NULL) {
  if (!desolve_is_available()) {
    if (is.null(desolve_error)) {
      desolve_error <- paste(
        "method = \"deSolve\" requires the optional deSolve package.",
        "Install deSolve or use method = \"euler\"."
      )
    }
    stop(desolve_error, call. = FALSE)
  }

  ode_args <- list(
    y = as.numeric(initial_state),
    times = times,
    func = function(time, state, parms) {
      list(as.numeric(derivative(time, state)))
    },
    parms = NULL
  )
  if (!is.null(tcrit)) {
    ode_args$tcrit <- desolve_tcrit(times, tcrit)
  }

  solved <- do.call(deSolve::ode, ode_args)

  values <- vector("list", length(times))
  state_columns <- seq_len(length(initial_state)) + 1
  for (i in seq_along(times)) {
    values[[i]] <- output(as.numeric(solved[i, state_columns]), solved[i, "time"])
  }

  bind_integrated_outputs(values)
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
