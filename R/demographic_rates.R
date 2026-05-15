#' Construct a fertility schedule
#'
#' Stores age-specific fertility rates for later demographic ODE use.
#' `fertility_rate` is interpreted as annual births per female person-year.
#' Fertility schedules may contain only reproductive age groups; full age
#' coverage is not required.
#'
#' @param fertility Data frame with `time`, `age_group`, and
#'   `fertility_rate` columns. Optional columns such as `sex` and
#'   `rate_source` are retained.
#' @param age_structure Age structure validated by [validate_age_structure()].
#'
#' @return An `agepi_fertility_schedule` list.
#' @examples
#' ages <- wpp_age_structure_5year()
#' fertility <- FertilitySchedule(
#'   data.frame(
#'     time = c(2020, 2020),
#'     age_group = c("20-24", "25-29"),
#'     fertility_rate = c(0.08, 0.1)
#'   ),
#'   ages
#' )
#' fertility$rate_convention
#' @export
FertilitySchedule <- function(fertility, age_structure) {
  validate_fertility_schedule_table(fertility, age_structure)

  fertility <- sort_demographic_rate_table(
    fertility,
    age_structure,
    value_columns = "fertility_rate"
  )
  times <- sort(unique(fertility$time))

  schedule <- list(
    data = fertility,
    age_structure = age_structure,
    times = times,
    age_groups = age_structure$age_groups,
    n_times = length(times),
    n_age_groups = age_structure$n_age_groups,
    rate_convention = "births_per_female_person_year"
  )
  class(schedule) <- c("agepi_fertility_schedule", "list")
  validate_fertility_schedule(schedule)
  schedule
}

#' Validate a fertility schedule
#'
#' @param x Object to validate.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_fertility_schedule <- function(x) {
  validate_demographic_schedule_object(
    x,
    class_name = "agepi_fertility_schedule",
    value_column = "fertility_rate",
    value_label = "fertility_rate",
    require_full_coverage = FALSE,
    allow_negative = FALSE
  )

  if (!identical(x$rate_convention, "births_per_female_person_year")) {
    stop("fertility_schedule rate_convention must be births_per_female_person_year.", call. = FALSE)
  }

  invisible(x)
}

#' Construct a mortality schedule
#'
#' Stores age-specific mortality rates for later demographic ODE use.
#' `mortality_rate` is interpreted as an annual continuous-time hazard. It may
#' also be used as a first-pass approximation to WPP-style central death rates.
#'
#' @param mortality Data frame with `time`, `age_group`, and
#'   `mortality_rate` columns.
#' @param age_structure Age structure validated by [validate_age_structure()].
#'
#' @return An `agepi_mortality_schedule` list.
#' @examples
#' ages <- wpp_age_structure_5year()
#' mortality <- MortalitySchedule(
#'   data.frame(
#'     time = 2020,
#'     age_group = ages$age_groups,
#'     mortality_rate = rep(0.01, ages$n_age_groups)
#'   ),
#'   ages
#' )
#' mortality$rate_convention
#' @export
MortalitySchedule <- function(mortality, age_structure) {
  validate_mortality_schedule_table(mortality, age_structure)

  mortality <- sort_demographic_rate_table(
    mortality,
    age_structure,
    value_columns = "mortality_rate"
  )
  times <- sort(unique(mortality$time))

  schedule <- list(
    data = mortality,
    age_structure = age_structure,
    times = times,
    age_groups = age_structure$age_groups,
    n_times = length(times),
    n_age_groups = age_structure$n_age_groups,
    rate_convention = "annual_hazard"
  )
  class(schedule) <- c("agepi_mortality_schedule", "list")
  validate_mortality_schedule(schedule)
  schedule
}

#' Validate a mortality schedule
#'
#' @param x Object to validate.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_mortality_schedule <- function(x) {
  validate_demographic_schedule_object(
    x,
    class_name = "agepi_mortality_schedule",
    value_column = "mortality_rate",
    value_label = "mortality_rate",
    require_full_coverage = TRUE,
    allow_negative = FALSE
  )

  if (!identical(x$rate_convention, "annual_hazard")) {
    stop("mortality_schedule rate_convention must be annual_hazard.", call. = FALSE)
  }

  invisible(x)
}

#' Construct a migration schedule
#'
#' Stores age-specific net migration inputs for later demographic ODE use.
#' Supply exactly one of `migration_rate` or `migration_count`. Migration values
#' may be negative, but must be finite.
#'
#' @param migration Data frame with `time`, `age_group`, and exactly one of
#'   `migration_rate` or `migration_count`.
#' @param age_structure Age structure validated by [validate_age_structure()].
#'
#' @return An `agepi_migration_schedule` list.
#' @examples
#' ages <- wpp_age_structure_5year()
#' migration <- MigrationSchedule(
#'   data.frame(
#'     time = 2020,
#'     age_group = ages$age_groups,
#'     migration_count = rep(0, ages$n_age_groups)
#'   ),
#'   ages
#' )
#' migration$migration_type
#' @export
MigrationSchedule <- function(migration, age_structure) {
  migration_type <- validate_migration_schedule_table(migration, age_structure)
  value_column <- paste0("migration_", migration_type)

  migration <- sort_demographic_rate_table(
    migration,
    age_structure,
    value_columns = value_column
  )
  times <- sort(unique(migration$time))

  schedule <- list(
    data = migration,
    age_structure = age_structure,
    times = times,
    age_groups = age_structure$age_groups,
    n_times = length(times),
    n_age_groups = age_structure$n_age_groups,
    migration_type = migration_type
  )
  class(schedule) <- c("agepi_migration_schedule", "list")
  validate_migration_schedule(schedule)
  schedule
}

#' Validate a migration schedule
#'
#' @param x Object to validate.
#'
#' @return The input invisibly if validation succeeds.
#' @export
validate_migration_schedule <- function(x) {
  if (!inherits(x, "agepi_migration_schedule")) {
    stop("x must be an agepi_migration_schedule object.", call. = FALSE)
  }

  if (!is.list(x) || !is.character(x$migration_type) || length(x$migration_type) != 1) {
    stop("migration_schedule must contain a single migration_type.", call. = FALSE)
  }

  if (!x$migration_type %in% c("rate", "count")) {
    stop("migration_schedule migration_type must be rate or count.", call. = FALSE)
  }

  validate_demographic_schedule_object(
    x,
    class_name = "agepi_migration_schedule",
    value_column = paste0("migration_", x$migration_type),
    value_label = paste0("migration_", x$migration_type),
    require_full_coverage = TRUE,
    allow_negative = TRUE
  )

  invisible(x)
}

validate_fertility_schedule_table <- function(fertility, age_structure) {
  validate_demographic_rate_table(
    fertility,
    age_structure,
    required_columns = c("time", "age_group", "fertility_rate"),
    value_column = "fertility_rate",
    value_label = "fertility_rate",
    require_full_coverage = FALSE,
    allow_negative = FALSE
  )
}

validate_mortality_schedule_table <- function(mortality, age_structure) {
  validate_demographic_rate_table(
    mortality,
    age_structure,
    required_columns = c("time", "age_group", "mortality_rate"),
    value_column = "mortality_rate",
    value_label = "mortality_rate",
    require_full_coverage = TRUE,
    allow_negative = FALSE
  )
}

validate_migration_schedule_table <- function(migration, age_structure) {
  if (!is.data.frame(migration)) {
    stop("migration must be a data frame.", call. = FALSE)
  }

  has_rate <- "migration_rate" %in% names(migration)
  has_count <- "migration_count" %in% names(migration)
  if (has_rate && has_count) {
    stop("migration must supply exactly one of migration_rate or migration_count, not both.", call. = FALSE)
  }
  if (!has_rate && !has_count) {
    stop("migration must supply exactly one of migration_rate or migration_count.", call. = FALSE)
  }

  migration_type <- if (has_rate) "rate" else "count"
  validate_demographic_rate_table(
    migration,
    age_structure,
    required_columns = c("time", "age_group", paste0("migration_", migration_type)),
    value_column = paste0("migration_", migration_type),
    value_label = paste0("migration_", migration_type),
    require_full_coverage = TRUE,
    allow_negative = TRUE
  )
  migration_type
}

validate_demographic_rate_table <- function(
  data,
  age_structure,
  required_columns,
  value_column,
  value_label,
  require_full_coverage,
  allow_negative
) {
  validate_age_structure(age_structure)

  if (!is.data.frame(data)) {
    stop("schedule data must be a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "schedule data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  time <- data$time
  age_group <- data$age_group
  value <- data[[value_column]]

  if (!is.numeric(time) || anyNA(time) || any(!is.finite(time))) {
    stop("schedule time must be finite and non-missing numeric values.", call. = FALSE)
  }

  if (anyNA(age_group)) {
    stop("schedule age_group cannot contain missing values.", call. = FALSE)
  }

  age_group <- as.character(age_group)
  unknown_ages <- setdiff(unique(age_group), age_structure$age_groups)
  if (length(unknown_ages) > 0) {
    stop(
      "schedule contains age_group value(s) not in age_structure: ",
      paste(unknown_ages, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
    stop(value_label, " must contain finite non-missing numeric values.", call. = FALSE)
  }

  if (!allow_negative && any(value < 0)) {
    stop(value_label, " cannot contain negative values.", call. = FALSE)
  }

  duplicate_rows <- duplicated(data.frame(time = time, age_group = age_group))
  if (any(duplicate_rows)) {
    stop("schedule contains duplicate time-age_group rows.", call. = FALSE)
  }

  if (require_full_coverage) {
    for (this_time in unique(time)) {
      observed_age_groups <- age_group[time == this_time]
      missing_age_groups <- setdiff(age_structure$age_groups, observed_age_groups)
      if (length(missing_age_groups) > 0) {
        stop(
          "schedule is missing age_group value(s) at time ",
          this_time,
          ": ",
          paste(missing_age_groups, collapse = ", "),
          call. = FALSE
        )
      }

      if (length(observed_age_groups) != age_structure$n_age_groups) {
        stop("schedule must contain exactly one row per age group at each time point.", call. = FALSE)
      }
    }
  }

  invisible(data)
}

sort_demographic_rate_table <- function(data, age_structure, value_columns) {
  retained_columns <- unique(c(
    "time",
    "age_group",
    value_columns,
    setdiff(names(data), c("time", "age_group", value_columns))
  ))
  data <- data[, retained_columns, drop = FALSE]
  data$age_group <- as.character(data$age_group)
  row_order <- order(data$time, match(data$age_group, age_structure$age_groups))
  data <- data[row_order, , drop = FALSE]
  row.names(data) <- NULL
  data
}

validate_demographic_schedule_object <- function(
  x,
  class_name,
  value_column,
  value_label,
  require_full_coverage,
  allow_negative
) {
  if (!inherits(x, class_name)) {
    stop("x must be a ", class_name, " object.", call. = FALSE)
  }

  if (!is.list(x)) {
    stop("schedule must be a list.", call. = FALSE)
  }

  required_fields <- c(
    "data",
    "age_structure",
    "times",
    "age_groups",
    "n_times",
    "n_age_groups"
  )
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop(
      "schedule is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  validate_demographic_rate_table(
    x$data,
    x$age_structure,
    required_columns = c("time", "age_group", value_column),
    value_column = value_column,
    value_label = value_label,
    require_full_coverage = require_full_coverage,
    allow_negative = allow_negative
  )

  expected_times <- sort(unique(x$data$time))
  if (!identical(x$times, expected_times)) {
    stop("schedule times must match sorted unique data times.", call. = FALSE)
  }

  if (!identical(x$age_groups, x$age_structure$age_groups)) {
    stop("schedule age_groups must match age_structure$age_groups.", call. = FALSE)
  }

  if (!identical(x$n_times, length(expected_times))) {
    stop("schedule n_times must equal length(times).", call. = FALSE)
  }

  if (!identical(x$n_age_groups, x$age_structure$n_age_groups)) {
    stop("schedule n_age_groups must match age_structure$n_age_groups.", call. = FALSE)
  }

  invisible(x)
}
