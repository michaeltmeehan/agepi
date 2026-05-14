#' Construct a WPP-style single-year age structure
#'
#' Creates a single-year demographic age structure with a final open-ended
#' terminal age group. For `max_age = 100`, age groups are `"0"`, `"1"`,
#' ..., `"99"`, and `"100+"`.
#'
#' @param max_age Positive whole-number lower bound for the final open-ended
#'   age group.
#'
#' @return An [AgeStructure()] object.
#' @examples
#' ages <- wpp_age_structure_1year()
#' tail(ages$age_groups)
#' @export
wpp_age_structure_1year <- function(max_age = 100) {
  validate_wpp_max_age(max_age)

  finite_ages <- seq.int(0, max_age - 1)
  AgeStructure(
    age_groups = c(as.character(finite_ages), paste0(max_age, "+")),
    lower_bounds = c(finite_ages, max_age),
    upper_bounds = c(finite_ages, Inf)
  )
}

#' Construct a WPP-style five-year age structure
#'
#' Creates a five-year demographic age structure with a final open-ended
#' terminal age group. For `max_age = 100`, age groups are `"0-4"`, `"5-9"`,
#' ..., `"95-99"`, and `"100+"`.
#'
#' @param max_age Positive whole-number lower bound for the final open-ended
#'   age group. Must be a multiple of 5.
#'
#' @return An [AgeStructure()] object.
#' @examples
#' ages <- wpp_age_structure_5year()
#' tail(ages$age_groups)
#' @export
wpp_age_structure_5year <- function(max_age = 100) {
  validate_wpp_max_age(max_age)
  if (max_age %% 5 != 0) {
    stop("max_age must be a positive multiple of 5.", call. = FALSE)
  }

  finite_lower_bounds <- seq.int(0, max_age - 5, by = 5)
  finite_upper_bounds <- finite_lower_bounds + 4
  AgeStructure(
    age_groups = c(
      paste0(finite_lower_bounds, "-", finite_upper_bounds),
      paste0(max_age, "+")
    ),
    lower_bounds = c(finite_lower_bounds, max_age),
    upper_bounds = c(finite_upper_bounds, Inf)
  )
}

validate_wpp_max_age <- function(max_age) {
  if (!is.numeric(max_age) || length(max_age) != 1 || anyNA(max_age) || !is.finite(max_age)) {
    stop("max_age must be a finite numeric scalar.", call. = FALSE)
  }

  if (max_age != as.integer(max_age)) {
    stop("max_age must be a whole number.", call. = FALSE)
  }

  if (max_age <= 0) {
    stop("max_age must be positive.", call. = FALSE)
  }

  invisible(max_age)
}
