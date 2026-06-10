#' Calibration target for already-run simulation outputs
#'
#' Defines one observed data target to compare against an existing deterministic
#' trajectory or cumulative-flow simulation output. This first calibration
#' scaffold evaluates simulation outputs only; it does not run simulations,
#' inject parameters, optimise, or sample.
#'
#' @param name Character scalar target name.
#' @param data Observed data frame. Must contain `time` and `value`; may contain
#'   `age_group`. Gaussian targets need an uncertainty from `sd`, `se`, or
#'   `sigma`. Lognormal targets need `sdlog` or `sigma`.
#' @param type Output source, either `"trajectory"` or `"cumulative"`.
#' @param observation_model Observation family: `"gaussian"`, `"lognormal"`, or
#'   `"poisson"`.
#' @param compartment Compartment to select for trajectory targets.
#' @param cumulative_name,transition_id,from,to Optional cumulative-flow
#'   selectors for cumulative targets.
#' @param sigma,sd,se,sdlog Optional scalar observation uncertainty. These are
#'   used when the corresponding column is absent from `data`.
#' @param weight Numeric scalar multiplier for this target contribution.
#' @param value_column Name of the observed value column in `data`.
#'
#' @return An S3 list with class `"agepi_calibration_target"`.
#' @export
#'
#' @examples
#' target <- CalibrationTarget(
#'   name = "infectious",
#'   data = data.frame(time = c(0, 1), value = c(10, 12), sd = c(1, 1)),
#'   type = "trajectory",
#'   observation_model = "gaussian",
#'   compartment = "I"
#' )
CalibrationTarget <- function(
  name,
  data,
  type = c("trajectory", "cumulative"),
  observation_model = c("gaussian", "lognormal", "poisson"),
  compartment = NULL,
  cumulative_name = NULL,
  transition_id = NULL,
  from = NULL,
  to = NULL,
  sigma = NULL,
  sd = NULL,
  se = NULL,
  sdlog = NULL,
  weight = 1,
  value_column = "value"
) {
  type <- match.arg(type)
  observation_model <- validate_calibration_observation_model(observation_model)
  validate_calibration_target_name(name)
  validate_calibration_target_data(data, value_column)
  validate_calibration_target_weight(weight)

  if (type == "trajectory") {
    validate_calibration_selector(compartment, "compartment", required = TRUE)
  } else {
    validate_calibration_selector(cumulative_name, "cumulative_name")
    validate_calibration_selector(transition_id, "transition_id")
    validate_calibration_selector(from, "from")
    validate_calibration_selector(to, "to")
    if (all(vapply(
      list(cumulative_name, transition_id, from, to),
      is.null,
      logical(1)
    ))) {
      stop(
        "cumulative calibration targets must specify at least one cumulative-flow selector.",
        call. = FALSE
      )
    }
  }

  observation <- list(
    family = observation_model,
    sigma = validate_optional_calibration_scalar(sigma, "sigma"),
    sd = validate_optional_calibration_scalar(sd, "sd"),
    se = validate_optional_calibration_scalar(se, "se"),
    sdlog = validate_optional_calibration_scalar(sdlog, "sdlog")
  )

  target <- list(
    name = name,
    data = data,
    type = type,
    observation_model = observation,
    selectors = list(
      compartment = compartment,
      cumulative_name = cumulative_name,
      transition_id = transition_id,
      from = from,
      to = to
    ),
    weight = weight,
    value_column = value_column
  )
  class(target) <- "agepi_calibration_target"
  target
}

#' Evaluate one calibration target
#'
#' @param target A target created by [CalibrationTarget()].
#' @param simulation_output A deterministic trajectory data frame, or a result
#'   list containing `trajectory` and/or `cumulative`.
#'
#' @return A list containing `target`, `log_likelihood`, `objective`,
#'   `n_observations`, `n_used`, `n_missing`, `missing`, and pointwise `details`.
#' @export
evaluate_calibration_target <- function(target, simulation_output) {
  validate_calibration_target_object(target)
  output <- calibration_extract_output(simulation_output, target$type)
  output <- calibration_filter_output(output, target)
  predicted <- calibration_prediction_table(output, target$data)
  mapped <- calibration_join_observed_predicted(target, predicted)

  usable <- !is.na(mapped$predicted)
  missing <- mapped[!usable, , drop = FALSE]
  evaluated <- mapped[usable, , drop = FALSE]

  if (nrow(evaluated) > 0) {
    contributions <- calibration_log_likelihood(target, evaluated)
    evaluated$log_likelihood <- contributions
  } else {
    contributions <- numeric()
  }

  log_likelihood <- sum(contributions)
  objective <- -target$weight * log_likelihood

  list(
    target = target$name,
    log_likelihood = log_likelihood,
    objective = objective,
    weight = target$weight,
    n_observations = nrow(mapped),
    n_used = nrow(evaluated),
    n_missing = nrow(missing),
    missing = missing,
    details = evaluated
  )
}

#' Evaluate a summed calibration objective
#'
#' @param targets A single [CalibrationTarget()] or a list of targets.
#' @param simulation_output A deterministic trajectory data frame, or a result
#'   list containing `trajectory` and/or `cumulative`.
#'
#' @return A list containing `log_likelihood`, `objective`, `n_observations`,
#'   `n_used`, and per-target `targets` results.
#' @export
evaluate_calibration_objective <- function(targets, simulation_output) {
  if (inherits(targets, "agepi_calibration_target")) {
    targets <- list(targets)
  }
  if (!is.list(targets) || length(targets) == 0) {
    stop("targets must be a CalibrationTarget or a non-empty list of targets.", call. = FALSE)
  }

  results <- lapply(targets, evaluate_calibration_target, simulation_output = simulation_output)
  list(
    log_likelihood = sum(vapply(results, `[[`, numeric(1), "log_likelihood")),
    objective = sum(vapply(results, `[[`, numeric(1), "objective")),
    n_observations = sum(vapply(results, `[[`, integer(1), "n_observations")),
    n_used = sum(vapply(results, `[[`, integer(1), "n_used")),
    targets = results
  )
}

validate_calibration_target_name <- function(name) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || name == "") {
    stop("name must be a non-empty character scalar.", call. = FALSE)
  }
}

validate_calibration_target_data <- function(data, value_column) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  if (!is.character(value_column) || length(value_column) != 1 ||
      is.na(value_column) || value_column == "") {
    stop("value_column must be a non-empty character scalar.", call. = FALSE)
  }

  required_columns <- c("time", value_column)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "calibration target data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(data[[value_column]])) {
    stop("calibration target value column must be numeric.", call. = FALSE)
  }
  if (anyNA(data[[value_column]]) || any(!is.finite(data[[value_column]]))) {
    stop("calibration target values must be finite and non-missing.", call. = FALSE)
  }

  invisible(data)
}

validate_calibration_target_weight <- function(weight) {
  if (!is.numeric(weight) || length(weight) != 1 || is.na(weight) ||
      !is.finite(weight) || weight < 0) {
    stop("weight must be a finite non-negative numeric scalar.", call. = FALSE)
  }
}

validate_calibration_observation_model <- function(observation_model) {
  if (!is.character(observation_model) || length(observation_model) != 1 ||
      is.na(observation_model) || observation_model == "") {
    stop("observation_model must be a non-empty character scalar.", call. = FALSE)
  }

  supported <- c("gaussian", "lognormal", "poisson")
  if (!observation_model %in% supported) {
    stop("unsupported observation model: ", observation_model, call. = FALSE)
  }

  observation_model
}

validate_calibration_selector <- function(x, name, required = FALSE) {
  if (is.null(x)) {
    if (required) {
      stop(name, " must be specified.", call. = FALSE)
    }
    return(invisible(x))
  }
  if (!is.character(x) || length(x) != 1 || is.na(x) || x == "") {
    stop(name, " must be a non-empty character scalar.", call. = FALSE)
  }
  invisible(x)
}

validate_optional_calibration_scalar <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x) || x <= 0) {
    stop(name, " must be a positive finite numeric scalar.", call. = FALSE)
  }
  as.numeric(x)
}

validate_calibration_target_object <- function(target) {
  if (!inherits(target, "agepi_calibration_target")) {
    stop("target must be created by CalibrationTarget().", call. = FALSE)
  }
  if (!is.character(target$type) || length(target$type) != 1 ||
      !target$type %in% c("trajectory", "cumulative")) {
    stop("unsupported calibration output type: ", target$type, call. = FALSE)
  }
  invisible(target)
}

calibration_extract_output <- function(simulation_output, type) {
  if (type == "trajectory") {
    if (is.data.frame(simulation_output)) {
      output <- simulation_output
    } else if (is.list(simulation_output) && !is.null(simulation_output$trajectory)) {
      output <- simulation_output$trajectory
    } else {
      stop(
        "simulation_output must be a trajectory data frame or a list with $trajectory.",
        call. = FALSE
      )
    }
    validate_calibration_output(output, c("time", "compartment", "age_group", "value"), "trajectory output")
    return(output)
  }

  if (is.data.frame(simulation_output)) {
    output <- simulation_output
  } else if (is.list(simulation_output) && !is.null(simulation_output$cumulative)) {
    output <- simulation_output$cumulative
  } else {
    stop(
      "simulation_output must be a cumulative data frame or a list with $cumulative.",
      call. = FALSE
    )
  }
  validate_calibration_output(
    output,
    c("time", "cumulative_name", "transition_id", "from", "to", "age_group", "value"),
    "cumulative output"
  )
  output
}

validate_calibration_output <- function(output, required_columns, label) {
  if (!is.data.frame(output)) {
    stop(label, " must be a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(required_columns, names(output))
  if (length(missing_columns) > 0) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(output$value)) {
    stop(label, " value must be numeric.", call. = FALSE)
  }
  if (anyNA(output$value) || any(!is.finite(output$value))) {
    stop(label, " value must be finite and non-missing.", call. = FALSE)
  }

  invisible(output)
}

calibration_filter_output <- function(output, target) {
  selectors <- target$selectors
  matches <- rep(TRUE, nrow(output))
  if (target$type == "trajectory") {
    matches <- matches & output$compartment == selectors$compartment
  } else {
    for (column in c("cumulative_name", "transition_id", "from", "to")) {
      if (!is.null(selectors[[column]])) {
        matches <- matches & output[[column]] == selectors[[column]]
      }
    }
  }

  selected <- output[matches, , drop = FALSE]
  if (nrow(selected) == 0) {
    stop("no model outputs matched calibration target '", target$name, "'.", call. = FALSE)
  }
  selected
}

calibration_prediction_table <- function(output, observed) {
  join_columns <- intersect(c("time", "age_group"), names(observed))
  grouping_key <- do.call(paste, c(output[join_columns], list(sep = "\r")))
  group_levels <- unique(grouping_key)
  group_index <- match(grouping_key, group_levels)
  first_rows <- match(group_levels, grouping_key)

  predicted <- output[first_rows, join_columns, drop = FALSE]
  predicted$predicted <- as.numeric(tapply(
    output$value,
    factor(group_index, levels = seq_along(group_levels)),
    sum
  ))
  row.names(predicted) <- NULL
  predicted
}

calibration_join_observed_predicted <- function(target, predicted) {
  observed <- target$data
  observed$observed <- observed[[target$value_column]]
  join_columns <- intersect(c("time", "age_group"), names(observed))

  mapped <- merge(observed, predicted, by = join_columns, all.x = TRUE, sort = FALSE)
  mapped <- mapped[order(match(do.call(paste, c(mapped[join_columns], list(sep = "\r"))),
                              do.call(paste, c(observed[join_columns], list(sep = "\r"))))), ,
                   drop = FALSE]
  row.names(mapped) <- NULL
  mapped
}

calibration_log_likelihood <- function(target, mapped) {
  family <- target$observation_model$family
  observed <- mapped$observed
  predicted <- mapped$predicted

  if (any(!is.finite(predicted))) {
    stop("model predictions must be finite and non-missing.", call. = FALSE)
  }

  if (family == "gaussian") {
    sigma <- calibration_observation_scale(target, mapped, c("sd", "se", "sigma"))
    return(stats::dnorm(observed, mean = predicted, sd = sigma, log = TRUE))
  }

  if (family == "lognormal") {
    if (any(observed <= 0) || any(predicted <= 0)) {
      stop("lognormal calibration targets require positive observed values and predictions.", call. = FALSE)
    }
    sigma <- calibration_observation_scale(target, mapped, c("sdlog", "sigma"))
    return(stats::dlnorm(observed, meanlog = log(predicted), sdlog = sigma, log = TRUE))
  }

  if (family == "poisson") {
    if (any(observed < 0)) {
      stop("poisson calibration targets cannot contain negative observed counts.", call. = FALSE)
    }
    if (any(abs(observed - round(observed)) > sqrt(.Machine$double.eps))) {
      stop("poisson calibration targets require integer observed counts.", call. = FALSE)
    }
    if (any(predicted < 0)) {
      stop("poisson calibration targets require non-negative model predictions.", call. = FALSE)
    }
    return(stats::dpois(observed, lambda = predicted, log = TRUE))
  }

  stop("unsupported observation model: ", family, call. = FALSE)
}

calibration_observation_scale <- function(target, mapped, candidates) {
  for (column in candidates) {
    if (column %in% names(mapped)) {
      scale <- mapped[[column]]
      validate_calibration_scale_vector(scale, column)
      return(scale)
    }
  }

  for (field in candidates) {
    scale <- target$observation_model[[field]]
    if (!is.null(scale)) {
      return(rep(scale, nrow(mapped)))
    }
  }

  stop(
    target$observation_model$family,
    " calibration targets require ",
    paste(candidates, collapse = " or "),
    ".",
    call. = FALSE
  )
}

validate_calibration_scale_vector <- function(x, name) {
  if (!is.numeric(x)) {
    stop(name, " must be numeric.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
    stop(name, " must contain positive finite non-missing values.", call. = FALSE)
  }
}
