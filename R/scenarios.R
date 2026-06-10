#' Extract an initial state from simulation output
#'
#' Extracts ordinary compartment rows at an exact simulation time and returns a
#' long-form initial-state data frame that can be passed to
#' [simulate_deterministic()] or [simulate_stochastic()]. The input may be a
#' trajectory data frame or a simulation result list containing a `trajectory`
#' element. Cumulative-flow tables are ignored.
#'
#' @param simulation Simulation trajectory data frame, or a list containing a
#'   `trajectory` data frame.
#' @param time Single numeric branch time. The time must exist exactly in the
#'   trajectory.
#'
#' @return Data frame with columns `compartment`, `age_group`, and `value`.
#' @export
initial_state_from_simulation <- function(simulation, time) {
  validate_initial_state_time(time)
  trajectory <- simulation_trajectory(simulation)
  validate_simulation_trajectory_for_state(trajectory)

  rows <- trajectory[trajectory$time == time, , drop = FALSE]
  if (nrow(rows) == 0) {
    stop("simulation trajectory does not contain requested time: ", time, call. = FALSE)
  }

  keys <- rows[, c("compartment", "age_group"), drop = FALSE]
  duplicate_rows <- duplicated(keys)
  if (any(duplicate_rows)) {
    duplicated_key <- keys[duplicate_rows, , drop = FALSE]
    stop(
      "simulation trajectory contains multiple rows for compartment-age pair at requested time: ",
      duplicated_key$compartment[1],
      "/",
      duplicated_key$age_group[1],
      call. = FALSE
    )
  }

  state <- rows[, c("compartment", "age_group", "value"), drop = FALSE]
  row.names(state) <- NULL
  state
}

#' Create a scenario
#'
#' Constructs a simple scenario object describing named argument overrides and
#' optional metadata for future scenario execution helpers.
#'
#' @param name Non-empty scenario name.
#' @param description Optional character description.
#' @param overrides Named list of simulation argument replacements.
#' @param modifier `NULL` or a function used by future scenario runners to
#'   modify a resolved simulation argument list.
#' @param metadata Named list of scenario metadata.
#'
#' @return A named list with class `"agepi_scenario"`.
#' @export
Scenario <- function(
  name,
  description = NULL,
  overrides = list(),
  modifier = NULL,
  metadata = list()
) {
  name <- validate_scenario_name(name)
  description <- validate_scenario_description(description)
  overrides <- validate_named_list(overrides, "overrides")
  validate_scenario_modifier(modifier)
  metadata <- validate_named_list(metadata, "metadata")

  structure(
    list(
      name = name,
      description = description,
      overrides = overrides,
      modifier = modifier,
      metadata = metadata
    ),
    class = "agepi_scenario"
  )
}

#' Create a scenario set
#'
#' Constructs a simple collection of named [Scenario()] objects for future
#' scenario execution helpers.
#'
#' @param ... One or more [Scenario()] objects.
#' @param baseline_name Optional name of the baseline scenario.
#' @param branch_time `NULL` or a single numeric branch time.
#' @param cumulative_policy How cumulative counters should be handled after a
#'   branch. Must be `"reset"` or `"continue"`.
#'
#' @return A named list with class `"agepi_scenario_set"`.
#' @export
ScenarioSet <- function(
  ...,
  baseline_name = NULL,
  branch_time = NULL,
  cumulative_policy = c("reset", "continue")
) {
  scenarios <- list(...)
  if (length(scenarios) == 1 && is.list(scenarios[[1]]) && !inherits(scenarios[[1]], "agepi_scenario")) {
    scenarios <- scenarios[[1]]
  }
  validate_scenarios(scenarios)

  scenario_names <- vapply(scenarios, function(scenario) scenario$name, character(1))
  duplicated_names <- unique(scenario_names[duplicated(scenario_names)])
  if (length(duplicated_names) > 0) {
    stop(
      "scenario names must be unique; duplicate name(s): ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }
  names(scenarios) <- scenario_names

  baseline_name <- validate_baseline_name(baseline_name, scenario_names)
  branch_time <- validate_branch_time(branch_time)
  cumulative_policy <- validate_cumulative_policy(cumulative_policy)

  structure(
    list(
      scenarios = scenarios,
      baseline_name = baseline_name,
      branch_time = branch_time,
      cumulative_policy = cumulative_policy
    ),
    class = "agepi_scenario_set"
  )
}

simulation_trajectory <- function(simulation) {
  if (is.data.frame(simulation)) {
    return(simulation)
  }

  if (is.list(simulation) && is.data.frame(simulation$trajectory)) {
    return(simulation$trajectory)
  }

  stop("simulation must be a trajectory data frame or a list with a trajectory data frame.", call. = FALSE)
}

validate_initial_state_time <- function(time) {
  if (!is.numeric(time) || length(time) != 1 || is.na(time) || !is.finite(time)) {
    stop("time must be a single finite numeric value.", call. = FALSE)
  }

  invisible(time)
}

validate_simulation_trajectory_for_state <- function(trajectory) {
  required_columns <- c("time", "compartment", "age_group", "value")
  missing_columns <- setdiff(required_columns, names(trajectory))
  if (length(missing_columns) > 0) {
    stop(
      "simulation trajectory is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(trajectory$time)) {
    stop("simulation trajectory time must be numeric.", call. = FALSE)
  }

  if (!is.numeric(trajectory$value)) {
    stop("simulation trajectory value must be numeric.", call. = FALSE)
  }

  if (anyNA(trajectory$time) || anyNA(trajectory$compartment) || anyNA(trajectory$age_group) || anyNA(trajectory$value)) {
    stop("simulation trajectory required columns cannot contain missing values.", call. = FALSE)
  }

  invisible(trajectory)
}

validate_scenario_name <- function(name) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || name == "") {
    stop("name must be a single non-empty string.", call. = FALSE)
  }

  name
}

validate_scenario_description <- function(description) {
  if (is.null(description)) {
    return(NULL)
  }

  if (!is.character(description) || length(description) != 1 || is.na(description)) {
    stop("description must be NULL or a single string.", call. = FALSE)
  }

  description
}

validate_named_list <- function(x, label) {
  if (!is.list(x) || is.data.frame(x)) {
    stop(label, " must be a named list.", call. = FALSE)
  }

  if (length(x) == 0) {
    return(x)
  }

  x_names <- names(x)
  if (is.null(x_names) || length(x_names) != length(x) || anyNA(x_names) || any(x_names == "")) {
    stop(label, " must be a named list.", call. = FALSE)
  }

  x
}

validate_scenario_modifier <- function(modifier) {
  if (!is.null(modifier) && !is.function(modifier)) {
    stop("modifier must be NULL or a function.", call. = FALSE)
  }

  invisible(modifier)
}

validate_scenarios <- function(scenarios) {
  if (length(scenarios) == 0) {
    stop("ScenarioSet requires at least one scenario.", call. = FALSE)
  }

  invalid <- !vapply(scenarios, inherits, logical(1), what = "agepi_scenario")
  if (any(invalid)) {
    stop("ScenarioSet entries must be Scenario objects.", call. = FALSE)
  }

  invisible(scenarios)
}

validate_baseline_name <- function(baseline_name, scenario_names) {
  if (is.null(baseline_name)) {
    return(NULL)
  }

  if (!is.character(baseline_name) || length(baseline_name) != 1 || is.na(baseline_name) || baseline_name == "") {
    stop("baseline_name must be NULL or a single non-empty string.", call. = FALSE)
  }

  if (!baseline_name %in% scenario_names) {
    stop("baseline_name must refer to one of the supplied scenarios.", call. = FALSE)
  }

  baseline_name
}

validate_branch_time <- function(branch_time) {
  if (is.null(branch_time)) {
    return(NULL)
  }

  if (!is.numeric(branch_time) || length(branch_time) != 1 || is.na(branch_time) || !is.finite(branch_time)) {
    stop("branch_time must be NULL or a single finite numeric value.", call. = FALSE)
  }

  branch_time
}

validate_cumulative_policy <- function(cumulative_policy) {
  allowed <- c("reset", "continue")
  if (!is.character(cumulative_policy) ||
      length(cumulative_policy) != 1 ||
      is.na(cumulative_policy) ||
      !cumulative_policy %in% allowed) {
    stop(
      "cumulative_policy must be one of: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }

  cumulative_policy
}
