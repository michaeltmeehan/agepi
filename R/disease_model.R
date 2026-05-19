#' Construct a minimal SIR disease model
#'
#' Creates a simple, inspectable disease-model object with compartments
#' `S`, `I`, and `R`; transitions `S -> I` and `I -> R`; and recovery
#' rate `gamma`.
#'
#' @param gamma Non-negative finite numeric scalar recovery rate.
#'
#' @return A `DiseaseModel` list describing an SIR model.
#' @examples
#' model <- SIRModel(gamma = 0.25)
#' model$compartments
#' @export
SIRModel <- function(gamma) {
  if (missing(gamma)) {
    stop("gamma is required.", call. = FALSE)
  }

  gamma <- validate_sir_gamma(gamma)

  model <- list(
    model_type = "SIR",
    compartments = c("S", "I", "R"),
    transitions = data.frame(
      from = c("S", "I"),
      to = c("I", "R"),
      stringsAsFactors = FALSE
    ),
    gamma = gamma
  )

  class(model) <- "DiseaseModel"
  validate_disease_model(model)
  model
}

#' Construct a minimal SEIR disease model
#'
#' Creates a simple, inspectable disease-model object with compartments
#' `S`, `E`, `I`, and `R`; transitions `S -> E`, `E -> I`, and `I -> R`;
#' progression rate `sigma`; and recovery rate `gamma`.
#'
#' @param sigma Non-negative finite numeric scalar progression rate.
#' @param gamma Non-negative finite numeric scalar recovery rate.
#'
#' @return A `DiseaseModel` list describing an SEIR model.
#' @examples
#' model <- SEIRModel(sigma = 0.2, gamma = 0.25)
#' model$transitions
#' @export
SEIRModel <- function(sigma, gamma) {
  if (missing(sigma)) {
    stop("sigma is required.", call. = FALSE)
  }
  if (missing(gamma)) {
    stop("gamma is required.", call. = FALSE)
  }

  sigma <- validate_seir_sigma(sigma)
  gamma <- validate_sir_gamma(gamma)

  model <- list(
    model_type = "SEIR",
    compartments = c("S", "E", "I", "R"),
    transitions = data.frame(
      from = c("S", "E", "I"),
      to = c("E", "I", "R"),
      stringsAsFactors = FALSE
    ),
    sigma = sigma,
    gamma = gamma
  )

  class(model) <- "DiseaseModel"
  validate_disease_model(model)
  model
}

#' Validate a disease model
#'
#' Checks the minimal disease-model fields currently required by `agepi`.
#' At present, SIR and SEIR models are supported.
#'
#' @param model Disease-model object to validate.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_disease_model <- function(model) {
  if (!is.list(model)) {
    stop("model must be a list.", call. = FALSE)
  }

  base_required_fields <- c("model_type", "compartments", "transitions")
  missing_base_fields <- setdiff(base_required_fields, names(model))
  if (length(missing_base_fields) > 0) {
    stop(
      "model is missing required field(s): ",
      paste(missing_base_fields, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.character(model$model_type) || length(model$model_type) != 1 ||
      anyNA(model$model_type) || model$model_type == "") {
    stop("model_type must be a non-missing character scalar.", call. = FALSE)
  }

  if (!model$model_type %in% c("SIR", "SEIR")) {
    stop("unsupported disease model type: ", model$model_type, call. = FALSE)
  }

  required_fields <- disease_model_required_fields(model$model_type)
  missing_fields <- setdiff(required_fields, names(model))
  if (length(missing_fields) > 0) {
    stop(
      "model is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  validate_compartments(model$compartments)

  if (model$model_type == "SIR") {
    if (!identical(model$compartments, c("S", "I", "R"))) {
      stop("SIR model compartments must be S, I, R.", call. = FALSE)
    }

    validate_sir_transitions(model$transitions)
    validate_sir_gamma(model$gamma)

    return(invisible(model))
  }

  if (!identical(model$compartments, c("S", "E", "I", "R"))) {
    stop("SEIR model compartments must be S, E, I, R.", call. = FALSE)
  }

  validate_seir_transitions(model$transitions)
  validate_seir_sigma(model$sigma)
  validate_sir_gamma(model$gamma)

  invisible(model)
}

disease_model_required_fields <- function(model_type) {
  if (model_type == "SIR") {
    return(c("model_type", "compartments", "transitions", "gamma"))
  }

  c("model_type", "compartments", "transitions", "sigma", "gamma")
}

validate_sir_gamma <- function(gamma) {
  if (!is.numeric(gamma) || length(gamma) != 1 || anyNA(gamma) || !is.finite(gamma)) {
    stop("gamma must be a finite numeric scalar.", call. = FALSE)
  }

  if (gamma < 0) {
    stop("gamma cannot be negative.", call. = FALSE)
  }

  as.numeric(gamma)
}

validate_seir_sigma <- function(sigma) {
  if (!is.numeric(sigma) || length(sigma) != 1 || anyNA(sigma) || !is.finite(sigma)) {
    stop("sigma must be a finite numeric scalar.", call. = FALSE)
  }

  if (sigma < 0) {
    stop("sigma cannot be negative.", call. = FALSE)
  }

  as.numeric(sigma)
}

validate_sir_transitions <- function(transitions) {
  validate_transition_schema(transitions)

  expected <- data.frame(
    from = c("S", "I"),
    to = c("I", "R"),
    stringsAsFactors = FALSE
  )

  observed <- transitions[, c("from", "to")]
  row.names(observed) <- NULL

  if (!identical(observed, expected)) {
    stop("SIR model transitions must be S -> I and I -> R.", call. = FALSE)
  }

  invisible(transitions)
}

validate_seir_transitions <- function(transitions) {
  validate_transition_schema(transitions)

  expected <- data.frame(
    from = c("S", "E", "I"),
    to = c("E", "I", "R"),
    stringsAsFactors = FALSE
  )

  observed <- transitions[, c("from", "to")]
  row.names(observed) <- NULL

  if (!identical(observed, expected)) {
    stop("SEIR model transitions must be S -> E, E -> I, and I -> R.", call. = FALSE)
  }

  invisible(transitions)
}

validate_transition_schema <- function(transitions) {
  if (!is.data.frame(transitions)) {
    stop("model transitions must be a data frame.", call. = FALSE)
  }

  required_columns <- c("from", "to")
  missing_columns <- setdiff(required_columns, names(transitions))
  if (length(missing_columns) > 0) {
    stop(
      "model transitions are missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(transitions)
}
