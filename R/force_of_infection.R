#' Compute age-specific force of infection
#'
#' Computes a reusable age-structured force of infection without depending on a
#' disease-model object. Rows of `contact_matrix` are recipient age groups and
#' columns are source age groups.
#'
#' The calculation is:
#' `lambda = beta * susceptibility * (contact_matrix %*% (infectiousness * infectious / population))`.
#'
#' @param infectious Numeric vector of infectious counts by source age group.
#' @param population Numeric vector of total population by source age group.
#' @param contact_matrix Numeric square matrix. Rows are recipient age groups;
#'   columns are source age groups.
#' @param beta Non-negative finite transmission scaling parameter.
#' @param susceptibility Optional non-negative numeric vector by recipient age
#'   group. Defaults to one for all age groups.
#' @param infectiousness Optional non-negative numeric vector by source age
#'   group. Defaults to one for all age groups.
#' @param age_structure Optional valid age structure used to validate lengths
#'   and name the output by age group.
#'
#' @return Numeric vector of force-of-infection values by recipient age group.
#' @export
force_of_infection <- function(
  infectious,
  population,
  contact_matrix,
  beta = 1,
  susceptibility = NULL,
  infectiousness = NULL,
  age_structure = NULL
) {
  n_age_groups <- force_of_infection_n_age_groups(
    infectious = infectious,
    population = population,
    contact_matrix = contact_matrix,
    age_structure = age_structure
  )

  infectious <- validate_force_vector(
    infectious,
    n_age_groups,
    "infectious",
    allow_zero = TRUE
  )
  population <- validate_force_vector(
    population,
    n_age_groups,
    "population",
    allow_zero = FALSE
  )
  contact_matrix <- validate_contact_matrix(contact_matrix, age_structure)
  validate_contact_matrix_length(contact_matrix, n_age_groups)
  beta <- validate_force_beta(beta)
  susceptibility <- validate_optional_force_vector(
    susceptibility,
    n_age_groups,
    "susceptibility"
  )
  infectiousness <- validate_optional_force_vector(
    infectiousness,
    n_age_groups,
    "infectiousness"
  )

  infectious_fraction <- infectiousness * infectious / population
  lambda <- as.numeric(beta * susceptibility * (contact_matrix %*% infectious_fraction))

  if (!is.null(age_structure)) {
    names(lambda) <- age_structure$age_groups
  }

  lambda
}

force_of_infection_n_age_groups <- function(
  infectious,
  population,
  contact_matrix,
  age_structure
) {
  if (!is.null(age_structure)) {
    validate_age_structure(age_structure)
    return(age_structure$n_age_groups)
  }

  if (!is.numeric(infectious)) {
    stop("infectious must be a numeric vector.", call. = FALSE)
  }

  length(infectious)
}

validate_force_vector <- function(x, expected_length, name, allow_zero) {
  if (!is.numeric(x) || is.matrix(x) || is.data.frame(x)) {
    stop(name, " must be a numeric vector.", call. = FALSE)
  }

  if (length(x) != expected_length) {
    stop(
      name,
      " length must match the number of age groups: ",
      expected_length,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop(name, " cannot contain missing or non-finite values.", call. = FALSE)
  }

  if (any(x < 0)) {
    stop(name, " cannot contain negative values.", call. = FALSE)
  }

  if (!allow_zero && any(x == 0)) {
    stop(name, " must contain positive values.", call. = FALSE)
  }

  as.numeric(x)
}

validate_optional_force_vector <- function(x, expected_length, name) {
  if (is.null(x)) {
    return(rep(1, expected_length))
  }

  validate_force_vector(
    x,
    expected_length = expected_length,
    name = name,
    allow_zero = TRUE
  )
}

#' Validate a contact matrix
#'
#' Checks that a contact matrix is numeric, finite, non-missing,
#' non-negative, square, and optionally consistent with an age structure.
#'
#' @param contact_matrix Numeric square matrix. Rows are recipient age groups;
#'   columns are source age groups.
#' @param age_structure Optional valid age structure used to validate matrix
#'   dimensions.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_contact_matrix <- function(contact_matrix, age_structure = NULL) {
  if (!is.numeric(contact_matrix) || !is.matrix(contact_matrix)) {
    stop("contact_matrix must be a numeric matrix.", call. = FALSE)
  }

  if (anyNA(contact_matrix) || any(!is.finite(contact_matrix))) {
    stop("contact_matrix cannot contain missing or non-finite values.", call. = FALSE)
  }

  if (any(contact_matrix < 0)) {
    stop("contact_matrix cannot contain negative values.", call. = FALSE)
  }

  matrix_dimensions <- dim(contact_matrix)
  if (matrix_dimensions[1] != matrix_dimensions[2]) {
    stop("contact_matrix must be square.", call. = FALSE)
  }

  if (!is.null(age_structure)) {
    validate_age_structure(age_structure)
    validate_contact_matrix_length(contact_matrix, age_structure$n_age_groups)
  }

  invisible(contact_matrix)
}

validate_contact_matrix_length <- function(contact_matrix, expected_length) {
  matrix_dimensions <- dim(contact_matrix)
  if (matrix_dimensions[1] != expected_length) {
    stop(
      "contact_matrix dimensions must match the number of age groups: ",
      expected_length,
      ".",
      call. = FALSE
    )
  }

  invisible(contact_matrix)
}

validate_force_beta <- function(beta) {
  if (!is.numeric(beta) || length(beta) != 1 || anyNA(beta) || !is.finite(beta)) {
    stop("beta must be a finite numeric scalar.", call. = FALSE)
  }

  if (beta < 0) {
    stop("beta cannot be negative.", call. = FALSE)
  }

  as.numeric(beta)
}
