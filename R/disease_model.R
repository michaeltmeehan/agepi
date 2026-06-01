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

#' Construct a generic age-structured compartment model
#'
#' Creates a simple generic disease-model object for deterministic
#' age-structured compartmental models. Infection transitions use the existing
#' [force_of_infection()] convention. Other transitions are per-capita flows
#' from `from` to `to`, with rates supplied either as finite non-negative
#' scalars or as age-specific vectors when transition rates are evaluated.
#'
#' `SIRModel()` and `SEIRModel()` remain supported convenience constructors for
#' the built-in SIR and SEIR workflows. `CompartmentModel()` is intended for
#' explicit custom models with the same compartment-major state ordering.
#'
#' @param compartments Unique non-empty character vector of compartment names.
#' @param infection_transitions Data frame with columns `from` and `to`.
#'   Infection flow is `lambda * state[from]` within each recipient age group.
#' @param transitions Optional data frame with columns `from`, `to`, and
#'   `rate`, where `rate` is a non-negative per-capita transition rate.
#' @param infectious_compartments Character vector naming compartment(s) that
#'   contribute to infectious pressure. Defaults to `"I"` when present.
#' @param infectiousness_weights Optional named non-negative numeric weights for
#'   `infectious_compartments`. Names must match infectious compartments.
#'   Defaults to one for each infectious compartment.
#' @param birth_compartment Optional compartment receiving demographic births.
#'   Defaults to `"S"` when present.
#' @param migration_compartment Optional compartment receiving net migration
#'   under `migration_policy = "susceptible"`. Defaults to `birth_compartment`.
#'
#' @return A `DiseaseModel` list describing a generic compartment model.
#' @examples
#' model <- CompartmentModel(
#'   compartments = c("S", "E", "I", "R"),
#'   infection_transitions = data.frame(from = "S", to = "E"),
#'   transitions = data.frame(
#'     from = c("E", "I"),
#'     to = c("I", "R"),
#'     rate = c(0.2, 0.25)
#'   ),
#'   infectious_compartments = "I"
#' )
#' @export
CompartmentModel <- function(
  compartments,
  infection_transitions = NULL,
  transitions = NULL,
  infectious_compartments = NULL,
  infectiousness_weights = NULL,
  birth_compartment = NULL,
  migration_compartment = NULL
) {
  validate_compartments(compartments)

  if (is.null(infection_transitions)) {
    infection_transitions <- data.frame(
      from = character(),
      to = character(),
      stringsAsFactors = FALSE
    )
  }
  infection_transitions <- validate_generic_transitions(
    infection_transitions,
    compartments = compartments,
    name = "infection_transitions",
    require_rate = FALSE
  )

  if (is.null(transitions)) {
    transitions <- data.frame(
      from = character(),
      to = character(),
      rate = numeric(),
      stringsAsFactors = FALSE
    )
  }
  transitions <- validate_generic_transitions(
    transitions,
    compartments = compartments,
    name = "transitions",
    require_rate = TRUE
  )

  if (is.null(infectious_compartments)) {
    infectious_compartments <- if ("I" %in% compartments) "I" else character()
  }
  infectious_compartments <- validate_generic_compartment_subset(
    infectious_compartments,
    compartments,
    "infectious_compartments",
    allow_empty = nrow(infection_transitions) == 0
  )

  if (nrow(infection_transitions) > 0 && length(infectious_compartments) == 0) {
    stop("infectious_compartments must name at least one compartment when infection_transitions are supplied.", call. = FALSE)
  }

  infectiousness_weights <- validate_generic_infectiousness_weights(
    infectiousness_weights,
    infectious_compartments,
    require_positive = nrow(infection_transitions) > 0
  )

  birth_compartment <- validate_optional_generic_compartment(
    birth_compartment,
    compartments,
    "birth_compartment"
  )
  if (is.null(birth_compartment) && "S" %in% compartments) {
    birth_compartment <- "S"
  }

  migration_compartment <- validate_optional_generic_compartment(
    migration_compartment,
    compartments,
    "migration_compartment"
  )
  if (is.null(migration_compartment)) {
    migration_compartment <- birth_compartment
  }

  model <- list(
    model_type = "CompartmentModel",
    compartments = compartments,
    transitions = transitions,
    infection_transitions = infection_transitions,
    infectious_compartments = infectious_compartments,
    infectiousness_weights = infectiousness_weights,
    birth_compartment = birth_compartment,
    migration_compartment = migration_compartment
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
#' SIR, SEIR, and generic `CompartmentModel()` models are supported.
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

  if (!model$model_type %in% c("SIR", "SEIR", "CompartmentModel")) {
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

  if (model$model_type == "CompartmentModel") {
    validate_generic_transitions(
      model$infection_transitions,
      compartments = model$compartments,
      name = "infection_transitions",
      require_rate = FALSE
    )
    validate_generic_transitions(
      model$transitions,
      compartments = model$compartments,
      name = "transitions",
      require_rate = TRUE
    )
    validate_generic_compartment_subset(
      model$infectious_compartments,
      model$compartments,
      "infectious_compartments",
      allow_empty = nrow(model$infection_transitions) == 0
    )
    validate_generic_infectiousness_weights(
      model$infectiousness_weights,
      model$infectious_compartments,
      require_positive = nrow(model$infection_transitions) > 0
    )
    validate_optional_generic_compartment(
      model$birth_compartment,
      model$compartments,
      "birth_compartment"
    )
    validate_optional_generic_compartment(
      model$migration_compartment,
      model$compartments,
      "migration_compartment"
    )

    return(invisible(model))
  }

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

  if (model_type == "CompartmentModel") {
    return(c(
      "model_type", "compartments", "transitions", "infection_transitions",
      "infectious_compartments", "infectiousness_weights",
      "birth_compartment", "migration_compartment"
    ))
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

validate_generic_transitions <- function(transitions, compartments, name, require_rate) {
  if (!is.data.frame(transitions)) {
    stop(name, " must be a data frame.", call. = FALSE)
  }

  required_columns <- c("from", "to")
  if (require_rate) {
    required_columns <- c(required_columns, "rate")
  }
  missing_columns <- setdiff(required_columns, names(transitions))
  if (length(missing_columns) > 0) {
    stop(
      name,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(transitions$from) || anyNA(transitions$to) ||
      any(transitions$from == "") || any(transitions$to == "")) {
    stop(name, " from and to cannot contain missing or empty values.", call. = FALSE)
  }

  unknown_from <- setdiff(unique(transitions$from), compartments)
  if (length(unknown_from) > 0) {
    stop(
      name,
      " contains unknown source compartment value(s): ",
      paste(unknown_from, collapse = ", "),
      call. = FALSE
    )
  }

  unknown_to <- setdiff(unique(transitions$to), compartments)
  if (length(unknown_to) > 0) {
    stop(
      name,
      " contains unknown destination compartment value(s): ",
      paste(unknown_to, collapse = ", "),
      call. = FALSE
    )
  }

  same_compartment <- transitions$from == transitions$to
  if (any(same_compartment)) {
    stop(name, " cannot contain self-transitions.", call. = FALSE)
  }

  transition_keys <- transitions[, c("from", "to"), drop = FALSE]
  duplicate_rows <- duplicated(transition_keys)
  if (any(duplicate_rows)) {
    duplicated_key <- transition_keys[which(duplicate_rows)[1], , drop = FALSE]
    stop(
      name,
      " contains duplicate transition: ",
      duplicated_key$from,
      "->",
      duplicated_key$to,
      call. = FALSE
    )
  }

  if (require_rate) {
    validate_generic_rate_column(transitions$rate, paste0(name, " rate"))
  }

  transitions[, required_columns, drop = FALSE]
}

validate_generic_rate_column <- function(rate, name) {
  if (!is.numeric(rate) && !is.list(rate)) {
    stop(name, " must be numeric.", call. = FALSE)
  }

  rate_list <- if (is.list(rate)) rate else as.list(rate)
  for (value in rate_list) {
    if (!is.numeric(value) || is.matrix(value) || is.data.frame(value) ||
        length(value) == 0 || anyNA(value) || any(!is.finite(value))) {
      stop(name, " must contain finite non-missing numeric value(s).", call. = FALSE)
    }
    if (any(value < 0)) {
      stop(name, " cannot contain negative values.", call. = FALSE)
    }
  }

  invisible(rate)
}

validate_generic_compartment_subset <- function(x, compartments, name, allow_empty) {
  if (!is.character(x)) {
    stop(name, " must be a character vector.", call. = FALSE)
  }

  if (!allow_empty && length(x) == 0) {
    stop(name, " must contain at least one compartment.", call. = FALSE)
  }

  if (anyNA(x) || any(x == "")) {
    stop(name, " cannot contain missing or empty values.", call. = FALSE)
  }

  duplicated_values <- unique(x[duplicated(x)])
  if (length(duplicated_values) > 0) {
    stop(
      name,
      " must be unique; duplicate compartment(s): ",
      paste(duplicated_values, collapse = ", "),
      call. = FALSE
    )
  }

  unknown <- setdiff(x, compartments)
  if (length(unknown) > 0) {
    stop(
      name,
      " contains unknown compartment value(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  x
}

validate_generic_infectiousness_weights <- function(
  weights,
  infectious_compartments,
  require_positive = FALSE
) {
  if (is.null(weights)) {
    weights <- rep(1, length(infectious_compartments))
    names(weights) <- infectious_compartments
  }

  if (!is.numeric(weights) || is.matrix(weights) || is.data.frame(weights) ||
      length(weights) != length(infectious_compartments) ||
      anyNA(weights) || any(!is.finite(weights))) {
    stop(
      "infectiousness_weights must be a finite numeric vector with one value per infectious_compartment.",
      call. = FALSE
    )
  }

  if (length(weights) > 0) {
    if (is.null(names(weights)) || anyNA(names(weights)) || any(names(weights) == "")) {
      stop("infectiousness_weights must be named by infectious compartment.", call. = FALSE)
    }

    duplicated_names <- unique(names(weights)[duplicated(names(weights))])
    if (length(duplicated_names) > 0) {
      stop(
        "infectiousness_weights names must be unique; duplicate compartment(s): ",
        paste(duplicated_names, collapse = ", "),
        call. = FALSE
      )
    }

    unknown_names <- setdiff(names(weights), infectious_compartments)
    if (length(unknown_names) > 0) {
      stop(
        "infectiousness_weights contains unknown infectious compartment value(s): ",
        paste(unknown_names, collapse = ", "),
        call. = FALSE
      )
    }

    missing_names <- setdiff(infectious_compartments, names(weights))
    if (length(missing_names) > 0) {
      stop(
        "infectiousness_weights is missing infectious compartment value(s): ",
        paste(missing_names, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (any(weights < 0)) {
    stop("infectiousness_weights cannot contain negative values.", call. = FALSE)
  }

  if (require_positive && !any(weights > 0)) {
    stop(
      "infectiousness_weights must contain at least one positive value when infection_transitions are supplied.",
      call. = FALSE
    )
  }

  weights <- as.numeric(weights[infectious_compartments])
  names(weights) <- infectious_compartments
  weights
}

validate_optional_generic_compartment <- function(x, compartments, name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (!is.character(x) || length(x) != 1 || anyNA(x) || x == "") {
    stop(name, " must be a non-missing character scalar.", call. = FALSE)
  }

  if (!x %in% compartments) {
    stop(name, " must name a model compartment.", call. = FALSE)
  }

  x
}
