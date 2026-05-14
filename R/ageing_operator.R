#' Construct an ODE-compatible ageing operator
#'
#' Creates an explicit ageing operator from an age structure. Finite age bins
#' have departure rates equal to `1 / width`, where width is calculated from
#' inclusive age bounds as `upper_bound - lower_bound + 1`. The final
#' open-ended bin is terminal: it has zero departure and no destination.
#'
#' @param age_structure Age structure validated by [validate_age_structure()].
#'
#' @return An `agepi_ageing_operator` list.
#' @examples
#' ageing <- AgeingOperator(wpp_age_structure_5year())
#' head(ageing$departure_rate)
#' @export
AgeingOperator <- function(age_structure) {
  validate_age_structure(age_structure)
  validate_contiguous_age_structure(age_structure)

  n_age_groups <- age_structure$n_age_groups
  widths <- age_structure$upper_bounds - age_structure$lower_bounds + 1
  widths[n_age_groups] <- Inf

  departure_rate <- 1 / widths
  departure_rate[n_age_groups] <- 0

  destination_index <- ageing_destination_index(n_age_groups)

  ageing_operator <- list(
    age_structure = age_structure,
    age_groups = age_structure$age_groups,
    n_age_groups = n_age_groups,
    lower_bounds = age_structure$lower_bounds,
    upper_bounds = age_structure$upper_bounds,
    widths = widths,
    departure_rate = departure_rate,
    destination_index = destination_index
  )

  class(ageing_operator) <- c("agepi_ageing_operator", "list")
  validate_ageing_operator(ageing_operator)
  ageing_operator
}

#' Validate an ageing operator
#'
#' Checks that an ageing operator has internally consistent dimensions,
#' non-negative finite departure rates, next-bin destinations for finite bins,
#' and terminal final-bin behaviour.
#'
#' @param x Object to validate.
#'
#' @return The input invisibly if validation succeeds.
#' @examples
#' ageing <- AgeingOperator(wpp_age_structure_1year())
#' validate_ageing_operator(ageing)
#' @export
validate_ageing_operator <- function(x) {
  if (!inherits(x, "agepi_ageing_operator")) {
    stop("x must be an agepi_ageing_operator object.", call. = FALSE)
  }

  if (!is.list(x)) {
    stop("x must be a list.", call. = FALSE)
  }

  required_fields <- c(
    "age_structure",
    "age_groups",
    "n_age_groups",
    "lower_bounds",
    "upper_bounds",
    "widths",
    "departure_rate",
    "destination_index"
  )
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop(
      "ageing_operator is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  validate_age_structure(x$age_structure)

  n_age_groups <- x$n_age_groups
  if (!is.numeric(n_age_groups) || length(n_age_groups) != 1 || anyNA(n_age_groups) || !is.finite(n_age_groups)) {
    stop("ageing_operator n_age_groups must be a finite numeric scalar.", call. = FALSE)
  }

  if (n_age_groups != as.integer(n_age_groups)) {
    stop("ageing_operator n_age_groups must be a whole number.", call. = FALSE)
  }

  if (n_age_groups != x$age_structure$n_age_groups) {
    stop("ageing_operator n_age_groups must match age_structure$n_age_groups.", call. = FALSE)
  }

  validate_ageing_operator_lengths(x, n_age_groups)

  if (!identical(x$age_groups, x$age_structure$age_groups)) {
    stop("ageing_operator age_groups must match age_structure$age_groups.", call. = FALSE)
  }

  if (!identical(x$lower_bounds, x$age_structure$lower_bounds) ||
      !identical(x$upper_bounds, x$age_structure$upper_bounds)) {
    stop("ageing_operator bounds must match age_structure bounds.", call. = FALSE)
  }

  if (!is.numeric(x$widths) || anyNA(x$widths) || any(x$widths <= 0)) {
    stop("ageing_operator widths must be positive numeric values.", call. = FALSE)
  }

  if (!is.infinite(x$widths[n_age_groups])) {
    stop("ageing_operator final width must be infinite.", call. = FALSE)
  }

  if (!is.numeric(x$departure_rate) || anyNA(x$departure_rate) || any(!is.finite(x$departure_rate))) {
    stop("ageing_operator departure_rate must contain finite numeric values.", call. = FALSE)
  }

  if (any(x$departure_rate < 0)) {
    stop("ageing_operator departure_rate cannot contain negative values.", call. = FALSE)
  }

  if (x$departure_rate[n_age_groups] != 0) {
    stop("ageing_operator final bin must have zero departure_rate.", call. = FALSE)
  }

  expected_destination <- ageing_destination_index(n_age_groups)
  if (!is.integer(x$destination_index) && !is.numeric(x$destination_index)) {
    stop("ageing_operator destination_index must be numeric or integer.", call. = FALSE)
  }

  if (!identical(as.integer(x$destination_index), expected_destination)) {
    stop("ageing_operator destination_index must point to the next age group, with NA for the final bin.", call. = FALSE)
  }

  invisible(x)
}

validate_contiguous_age_structure <- function(age_structure) {
  if (age_structure$n_age_groups == 1) {
    if (!is.infinite(age_structure$upper_bounds[1])) {
      stop("age_structure must end with an open-ended final age bin.", call. = FALSE)
    }
    return(invisible(age_structure))
  }

  if (!is.infinite(age_structure$upper_bounds[age_structure$n_age_groups])) {
    stop("age_structure must end with an open-ended final age bin.", call. = FALSE)
  }

  expected_next_lower <- age_structure$upper_bounds[-age_structure$n_age_groups] + 1
  if (any(age_structure$lower_bounds[-1] != expected_next_lower)) {
    stop("age_structure age bins must be contiguous for ageing.", call. = FALSE)
  }

  invisible(age_structure)
}

validate_ageing_operator_lengths <- function(x, n_age_groups) {
  length_fields <- c(
    "age_groups",
    "lower_bounds",
    "upper_bounds",
    "widths",
    "departure_rate",
    "destination_index"
  )

  for (field in length_fields) {
    if (length(x[[field]]) != n_age_groups) {
      stop(
        "ageing_operator ",
        field,
        " length must equal n_age_groups.",
        call. = FALSE
      )
    }
  }

  invisible(x)
}

ageing_destination_index <- function(n_age_groups) {
  if (n_age_groups == 1) {
    return(NA_integer_)
  }

  c(seq.int(2, n_age_groups), NA_integer_)
}
