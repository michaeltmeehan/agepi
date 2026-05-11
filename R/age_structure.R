#' Construct an age structure
#'
#' Creates a small list object defining explicit model age groups.
#' Age groups must be sorted, non-overlapping, and use unique labels.
#' The final age group may be open-ended with `upper_bounds = Inf`.
#'
#' @param age_groups Character vector of age-group labels.
#' @param lower_bounds Numeric vector of inclusive lower age bounds.
#' @param upper_bounds Numeric vector of inclusive upper age bounds.
#'
#' @return An `AgeStructure` list with `lower`, `upper`, and `labels` fields.
#' @export
AgeStructure <- function(age_groups, lower_bounds = NULL, upper_bounds = NULL) {
  if (is.null(lower_bounds) || is.null(upper_bounds)) {
    stop(
      "AgeStructure() requires explicit lower_bounds and upper_bounds.",
      call. = FALSE
    )
  }

  age_structure <- list(
    lower = lower_bounds,
    upper = upper_bounds,
    labels = age_groups
  )

  class(age_structure) <- "AgeStructure"
  validate_age_structure(age_structure)
  age_structure
}

#' Validate an age structure
#'
#' Checks that an age structure has explicit numeric bounds and unique,
#' non-missing character labels. Age bins must be sorted and non-overlapping.
#'
#' @param age_structure An object with `lower`, `upper`, and `labels` fields.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_age_structure <- function(age_structure) {
  if (!is.list(age_structure)) {
    stop("age_structure must be a list.", call. = FALSE)
  }

  required_fields <- c("lower", "upper", "labels")
  missing_fields <- setdiff(required_fields, names(age_structure))
  if (length(missing_fields) > 0) {
    stop(
      "age_structure is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  lower <- age_structure$lower
  upper <- age_structure$upper
  labels <- age_structure$labels

  if (!is.numeric(lower) || !is.numeric(upper)) {
    stop("age_structure lower and upper bounds must be numeric.", call. = FALSE)
  }

  if (!is.character(labels)) {
    stop("age_structure labels must be a character vector.", call. = FALSE)
  }

  lengths <- c(length(lower), length(upper), length(labels))
  if (length(unique(lengths)) != 1) {
    stop(
      "age_structure lower, upper, and labels must have the same length.",
      call. = FALSE
    )
  }

  if (length(labels) == 0) {
    stop("age_structure must contain at least one age group.", call. = FALSE)
  }

  if (anyNA(lower) || anyNA(upper) || anyNA(labels)) {
    stop("age_structure lower, upper, and labels cannot contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(lower))) {
    stop("age_structure lower bounds must be finite.", call. = FALSE)
  }

  if (any(!is.finite(upper[-length(upper)]))) {
    stop(
      "age_structure upper bounds must be finite except for the final open-ended bin.",
      call. = FALSE
    )
  }

  if (length(upper) > 1 && is.infinite(upper[-length(upper)][1])) {
    stop("Only the final age bin may be open-ended.", call. = FALSE)
  }

  if (any(upper <= lower)) {
    stop("Each age_structure upper bound must be greater than its lower bound.", call. = FALSE)
  }

  if (any(diff(lower) <= 0)) {
    stop("age_structure lower bounds must be strictly increasing.", call. = FALSE)
  }

  if (length(lower) > 1 && any(lower[-1] <= upper[-length(upper)])) {
    stop("age_structure age bins must be sorted and non-overlapping.", call. = FALSE)
  }

  if (any(labels == "")) {
    stop("age_structure labels cannot be empty strings.", call. = FALSE)
  }

  duplicated_labels <- unique(labels[duplicated(labels)])
  if (length(duplicated_labels) > 0) {
    stop(
      "age_structure labels must be unique; duplicate label(s): ",
      paste(duplicated_labels, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(age_structure)
}
