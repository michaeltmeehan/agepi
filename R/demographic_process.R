#' Construct a demographic process configuration
#'
#' Combines an age structure, ageing operator, and optional fertility,
#' mortality, and migration schedules. This object only stores and validates
#' demographic process inputs; it does not implement derivatives, simulation,
#' interpolation, or infection coupling.
#'
#' @param age_structure Age structure validated by [validate_age_structure()].
#' @param ageing_operator Optional ageing operator. Defaults to
#'   [AgeingOperator()] for `age_structure`.
#' @param fertility_schedule Optional [FertilitySchedule()] object.
#' @param fertility_exposure_fraction Single non-negative finite multiplier
#'   applied to fertility exposure when computing births. The default, `1`,
#'   preserves the historical total-population exposure convention. Use values
#'   such as `0.5` only as an approximate female-exposure correction for
#'   total-population simulations with female ASFR; this is not a substitute for
#'   sex-structured demography.
#' @param mortality_schedule Optional [MortalitySchedule()] object.
#' @param migration_schedule Optional [MigrationSchedule()] object.
#' @param mode Process mode, either `"closed"` or `"migration"`.
#'
#' @return An `agepi_demographic_process` list.
#' @examples
#' ages <- wpp_age_structure_5year()
#' process <- DemographicProcess(age_structure = ages, mode = "closed")
#' process$mode
#' @export
DemographicProcess <- function(
  age_structure,
  ageing_operator = NULL,
  fertility_schedule = NULL,
  fertility_exposure_fraction = 1,
  mortality_schedule = NULL,
  migration_schedule = NULL,
  mode = c("closed", "migration")
) {
  mode <- match.arg(mode)
  validate_age_structure(age_structure)
  fertility_exposure_fraction <- validate_fertility_exposure_fraction(fertility_exposure_fraction)

  if (is.null(ageing_operator)) {
    ageing_operator <- AgeingOperator(age_structure)
  }

  process <- list(
    age_structure = age_structure,
    ageing_operator = ageing_operator,
    fertility_schedule = fertility_schedule,
    fertility_exposure_fraction = fertility_exposure_fraction,
    mortality_schedule = mortality_schedule,
    migration_schedule = migration_schedule,
    mode = mode,
    times = infer_demographic_process_common_times(
      fertility_schedule,
      mortality_schedule,
      migration_schedule
    )
  )
  class(process) <- c("agepi_demographic_process", "list")
  validate_demographic_process(process)
  process
}

#' Validate a demographic process configuration
#'
#' @param x Object to validate.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_demographic_process <- function(x) {
  if (!inherits(x, "agepi_demographic_process")) {
    stop("x must be an agepi_demographic_process object.", call. = FALSE)
  }

  if (!is.list(x)) {
    stop("demographic_process must be a list.", call. = FALSE)
  }

  required_fields <- c(
    "age_structure",
    "ageing_operator",
    "fertility_schedule",
    "fertility_exposure_fraction",
    "mortality_schedule",
    "migration_schedule",
    "mode",
    "times"
  )
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop(
      "demographic_process is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  validate_age_structure(x$age_structure)
  validate_ageing_operator(x$ageing_operator)
  validate_same_age_structure(
    x$age_structure,
    x$ageing_operator$age_structure,
    "ageing_operator"
  )

  if (!is.character(x$mode) || length(x$mode) != 1 || anyNA(x$mode)) {
    stop("demographic_process mode must be a single character value.", call. = FALSE)
  }

  if (!x$mode %in% c("closed", "migration")) {
    stop("demographic_process mode must be closed or migration.", call. = FALSE)
  }

  if (x$mode == "closed" && !is.null(x$migration_schedule)) {
    stop("closed demographic_process mode requires migration_schedule to be NULL.", call. = FALSE)
  }

  validate_fertility_exposure_fraction(x$fertility_exposure_fraction)

  validate_optional_process_schedule(
    x$fertility_schedule,
    "fertility_schedule",
    "agepi_fertility_schedule",
    validate_fertility_schedule,
    x$age_structure
  )
  validate_optional_process_schedule(
    x$mortality_schedule,
    "mortality_schedule",
    "agepi_mortality_schedule",
    validate_mortality_schedule,
    x$age_structure
  )
  validate_optional_process_schedule(
    x$migration_schedule,
    "migration_schedule",
    "agepi_migration_schedule",
    validate_migration_schedule,
    x$age_structure
  )

  expected_times <- infer_demographic_process_common_times(
    x$fertility_schedule,
    x$mortality_schedule,
    x$migration_schedule
  )
  if (!identical(x$times, expected_times)) {
    stop("demographic_process times must be the common schedule time grid or NULL.", call. = FALSE)
  }

  invisible(x)
}

validate_fertility_exposure_fraction <- function(fertility_exposure_fraction) {
  if (!is.numeric(fertility_exposure_fraction) || length(fertility_exposure_fraction) != 1) {
    stop("fertility_exposure_fraction must be a single numeric value.", call. = FALSE)
  }

  if (anyNA(fertility_exposure_fraction) || !is.finite(fertility_exposure_fraction)) {
    stop("fertility_exposure_fraction must be finite and non-missing.", call. = FALSE)
  }

  if (fertility_exposure_fraction < 0) {
    stop("fertility_exposure_fraction must be non-negative.", call. = FALSE)
  }

  as.numeric(fertility_exposure_fraction)
}

validate_optional_process_schedule <- function(
  schedule,
  field_name,
  class_name,
  validator,
  age_structure
) {
  if (is.null(schedule)) {
    return(invisible(NULL))
  }

  if (!inherits(schedule, class_name)) {
    stop(field_name, " must be a ", class_name, " object.", call. = FALSE)
  }

  validator(schedule)
  validate_same_age_structure(age_structure, schedule$age_structure, field_name)
  invisible(schedule)
}

validate_same_age_structure <- function(expected, observed, object_name) {
  validate_age_structure(expected)
  validate_age_structure(observed)

  if (!identical(expected$age_groups, observed$age_groups) ||
      !identical(expected$lower_bounds, observed$lower_bounds) ||
      !identical(expected$upper_bounds, observed$upper_bounds)) {
    stop(object_name, " must use the same age_structure.", call. = FALSE)
  }

  invisible(observed)
}

infer_demographic_process_common_times <- function(...) {
  schedules <- list(...)
  schedules <- schedules[!vapply(schedules, is.null, logical(1))]

  if (length(schedules) == 0) {
    return(NULL)
  }

  first_times <- schedules[[1]]$times
  all_same <- all(vapply(
    schedules,
    function(schedule) identical(schedule$times, first_times),
    logical(1)
  ))

  if (all_same) {
    return(first_times)
  }

  NULL
}
