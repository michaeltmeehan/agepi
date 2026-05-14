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
#' @return An `AgeStructure` list with `age_groups`, `n_age_groups`,
#'   `lower_bounds`, and `upper_bounds` fields.
#' @export
AgeStructure <- function(age_groups, lower_bounds = NULL, upper_bounds = NULL) {
  if (is.null(lower_bounds) || is.null(upper_bounds)) {
    stop(
      "AgeStructure() requires explicit lower_bounds and upper_bounds.",
      call. = FALSE
    )
  }

  age_structure <- list(
    age_groups = age_groups,
    n_age_groups = length(age_groups),
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds
  )

  class(age_structure) <- "AgeStructure"
  validate_age_structure(age_structure)
  age_structure
}

#' Validate an age structure
#'
#' Checks that an age structure has explicit numeric bounds and unique,
#' non-missing character age groups. Age bins must be sorted and
#' non-overlapping.
#'
#' @param age_structure An object with `age_groups`, `n_age_groups`,
#'   `lower_bounds`, and `upper_bounds` fields.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_age_structure <- function(age_structure) {
  if (!is.list(age_structure)) {
    stop("age_structure must be a list.", call. = FALSE)
  }

  required_fields <- c("age_groups", "n_age_groups", "lower_bounds", "upper_bounds")
  missing_fields <- setdiff(required_fields, names(age_structure))
  if (length(missing_fields) > 0) {
    stop(
      "age_structure is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  lower_bounds <- age_structure$lower_bounds
  upper_bounds <- age_structure$upper_bounds
  age_groups <- age_structure$age_groups
  n_age_groups <- age_structure$n_age_groups

  if (!is.numeric(lower_bounds) || !is.numeric(upper_bounds)) {
    stop("age_structure lower_bounds and upper_bounds must be numeric.", call. = FALSE)
  }

  if (!is.character(age_groups)) {
    stop("age_structure age_groups must be a character vector.", call. = FALSE)
  }

  if (!is.integer(n_age_groups) && !is.numeric(n_age_groups)) {
    stop("age_structure n_age_groups must be a numeric scalar.", call. = FALSE)
  }

  if (length(n_age_groups) != 1 || anyNA(n_age_groups) || !is.finite(n_age_groups)) {
    stop("age_structure n_age_groups must be a finite numeric scalar.", call. = FALSE)
  }

  if (n_age_groups != as.integer(n_age_groups)) {
    stop("age_structure n_age_groups must be a whole number.", call. = FALSE)
  }

  if (n_age_groups != length(age_groups)) {
    stop("age_structure n_age_groups must equal length(age_groups).", call. = FALSE)
  }

  lengths <- c(length(lower_bounds), length(upper_bounds), length(age_groups))
  if (length(unique(lengths)) != 1) {
    stop(
      "age_structure lower_bounds, upper_bounds, and age_groups must have the same length.",
      call. = FALSE
    )
  }

  if (length(age_groups) == 0) {
    stop("age_structure must contain at least one age group.", call. = FALSE)
  }

  if (anyNA(lower_bounds) || anyNA(upper_bounds) || anyNA(age_groups)) {
    stop("age_structure lower_bounds, upper_bounds, and age_groups cannot contain missing values.", call. = FALSE)
  }

  if (any(!is.finite(lower_bounds))) {
    stop("age_structure lower_bounds must be finite.", call. = FALSE)
  }

  if (any(!is.finite(upper_bounds[-length(upper_bounds)]))) {
    stop(
      "age_structure upper_bounds must be finite except for the final open-ended bin.",
      call. = FALSE
    )
  }

  if (length(upper_bounds) > 1 && is.infinite(upper_bounds[-length(upper_bounds)][1])) {
    stop("Only the final age bin may be open-ended.", call. = FALSE)
  }

  if (any(upper_bounds < lower_bounds)) {
    stop("Each age_structure upper_bound must be greater than or equal to its lower_bound.", call. = FALSE)
  }

  if (any(diff(lower_bounds) <= 0)) {
    stop("age_structure lower_bounds must be strictly increasing.", call. = FALSE)
  }

  if (length(lower_bounds) > 1 && any(lower_bounds[-1] <= upper_bounds[-length(upper_bounds)])) {
    stop("age_structure age bins must be sorted and non-overlapping.", call. = FALSE)
  }

  if (any(age_groups == "")) {
    stop("age_structure age_groups cannot be empty strings.", call. = FALSE)
  }

  duplicated_age_groups <- unique(age_groups[duplicated(age_groups)])
  if (length(duplicated_age_groups) > 0) {
    stop(
      "age_structure age_groups must be unique; duplicate age_group(s): ",
      paste(duplicated_age_groups, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(age_structure)
}
