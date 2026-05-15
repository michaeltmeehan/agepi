#' Standardise WPP-like mortality inputs
#'
#' Converts a WPP-like age-specific mortality table into a
#' [MortalitySchedule()]. This helper is an adapter only: it standardises
#' columns and age labels, then delegates validation and storage to
#' [MortalitySchedule()]. It does not implement projection dynamics, simulate
#' demography, or guarantee reproduction of WPP projections. The optional
#' `wpp2024` package is not required when supplying data directly.
#'
#' The mortality input is interpreted as an annual hazard, compatible with the
#' [MortalitySchedule()] convention. WPP central death rates (`mx`) may be used
#' as a first-pass annual hazard approximation.
#'
#' @param data Data frame containing WPP-like mortality data.
#' @param age_structure Target age-structure object validated by
#'   [validate_age_structure()].
#' @param time_col,age_col,mortality_col Column names containing time/year, age
#'   labels, and mortality rates.
#'
#' @return An `agepi_mortality_schedule` object.
#' @examples
#' ages <- wpp_age_structure_5year(max_age = 10)
#' mortality_like <- data.frame(
#'   year = rep(2020, ages$n_age_groups),
#'   age = ages$age_groups,
#'   mx = c(0.01, 0.002, 0.03)
#' )
#' standardise_wpp_mortality(
#'   mortality_like,
#'   age_structure = ages,
#'   time_col = "year",
#'   age_col = "age",
#'   mortality_col = "mx"
#' )
#' @export
standardise_wpp_mortality <- function(data,
                                      age_structure,
                                      time_col,
                                      age_col,
                                      mortality_col) {
  validate_wpp_standardisation_inputs(
    data,
    age_structure,
    required_columns = c(time_col, age_col, mortality_col)
  )

  mortality <- data.frame(
    time = data[[time_col]],
    age_group = parse_wpp_age_labels(data[[age_col]], age_structure),
    mortality_rate = data[[mortality_col]],
    stringsAsFactors = FALSE
  )

  MortalitySchedule(mortality, age_structure)
}

#' Standardise WPP-like fertility inputs
#'
#' Converts a WPP-like age-specific fertility table into a
#' [FertilitySchedule()]. This helper accepts already-computed age-specific
#' fertility rates. It does not derive rates from total fertility rates or WPP
#' percentage age schedules.
#'
#' This helper is an adapter only: it standardises columns and age labels, then
#' delegates validation and storage to [FertilitySchedule()]. It does not
#' implement projection dynamics, simulate demography, or guarantee reproduction
#' of WPP projections. The optional `wpp2024` package is not required when
#' supplying data directly.
#'
#' Fertility rates are interpreted as annual births per female person-year.
#'
#' @param data Data frame containing WPP-like fertility data.
#' @param age_structure Target age-structure object validated by
#'   [validate_age_structure()].
#' @param time_col,age_col Column names containing time/year and age labels.
#' @param fertility_col Column name containing age-specific fertility rates.
#'
#' @return An `agepi_fertility_schedule` object.
#' @examples
#' ages <- wpp_age_structure_5year()
#' fertility_like <- data.frame(
#'   year = c(2020, 2020),
#'   age = c("20-24", "25-29"),
#'   asfr = c(0.08, 0.1)
#' )
#' standardise_wpp_fertility(
#'   fertility_like,
#'   age_structure = ages,
#'   time_col = "year",
#'   age_col = "age",
#'   fertility_col = "asfr"
#' )
#' @export
standardise_wpp_fertility <- function(data,
                                      age_structure,
                                      time_col,
                                      age_col,
                                      fertility_col) {
  validate_wpp_standardisation_inputs(
    data,
    age_structure,
    required_columns = c(time_col, age_col, fertility_col)
  )

  fertility <- data.frame(
    time = data[[time_col]],
    age_group = parse_wpp_age_labels(data[[age_col]], age_structure),
    fertility_rate = data[[fertility_col]],
    stringsAsFactors = FALSE
  )

  FertilitySchedule(fertility, age_structure)
}

#' Standardise WPP-like migration inputs
#'
#' Converts a WPP-like age-specific migration table into a
#' [MigrationSchedule()]. Migration values may be negative. Use
#' `migration_type = "count"` for net migrant counts and
#' `migration_type = "rate"` for annual per-person migration rates.
#'
#' This helper is an adapter only: it standardises columns and age labels, then
#' delegates validation and storage to [MigrationSchedule()]. It does not
#' implement projection dynamics, simulate demography, or guarantee reproduction
#' of WPP projections. The optional `wpp2024` package is not required when
#' supplying data directly.
#'
#' @param data Data frame containing WPP-like migration data.
#' @param age_structure Target age-structure object validated by
#'   [validate_age_structure()].
#' @param time_col,age_col Column names containing time/year and age labels.
#' @param migration_col Column name containing age-specific migration values.
#' @param migration_type Whether `migration_col` contains `"count"` or `"rate"`
#'   values.
#'
#' @return An `agepi_migration_schedule` object.
#' @examples
#' ages <- wpp_age_structure_5year(max_age = 10)
#' migration_like <- data.frame(
#'   year = rep(2020, ages$n_age_groups),
#'   age = ages$age_groups,
#'   net_migration = c(-10, 5, 2)
#' )
#' standardise_wpp_migration(
#'   migration_like,
#'   age_structure = ages,
#'   time_col = "year",
#'   age_col = "age",
#'   migration_col = "net_migration",
#'   migration_type = "count"
#' )
#' @export
standardise_wpp_migration <- function(data,
                                      age_structure,
                                      time_col,
                                      age_col,
                                      migration_col,
                                      migration_type = c("count", "rate")) {
  migration_type <- match.arg(migration_type)
  validate_wpp_standardisation_inputs(
    data,
    age_structure,
    required_columns = c(time_col, age_col, migration_col)
  )

  value_column <- paste0("migration_", migration_type)
  migration <- data.frame(
    time = data[[time_col]],
    age_group = parse_wpp_age_labels(data[[age_col]], age_structure),
    value = data[[migration_col]],
    stringsAsFactors = FALSE
  )
  names(migration)[names(migration) == "value"] <- value_column

  MigrationSchedule(migration, age_structure)
}

validate_wpp_standardisation_inputs <- function(data, age_structure, required_columns) {
  validate_age_structure(age_structure)

  if (!is.data.frame(data)) {
    stop("WPP data must be a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "WPP data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(data)
}

parse_wpp_age_labels <- function(age, age_structure) {
  validate_age_structure(age_structure)

  if (anyNA(age)) {
    stop("WPP age labels cannot contain missing values.", call. = FALSE)
  }

  labels <- trimws(as.character(age))
  parsed <- vapply(
    labels,
    parse_one_wpp_age_label,
    character(1),
    age_structure = age_structure,
    USE.NAMES = FALSE
  )

  parsed
}

parse_one_wpp_age_label <- function(label, age_structure) {
  if (label %in% age_structure$age_groups) {
    return(label)
  }

  lower_bound <- parse_wpp_age_lower_bound(label)
  if (is.na(lower_bound)) {
    stop("Unsupported WPP age label: ", label, call. = FALSE)
  }

  matches <- which(age_structure$lower_bounds == lower_bound)
  if (length(matches) != 1) {
    stop("WPP age label is not in age_structure: ", label, call. = FALSE)
  }

  age_structure$age_groups[matches]
}

parse_wpp_age_lower_bound <- function(label) {
  if (grepl("^[0-9]+$", label)) {
    return(as.numeric(label))
  }

  if (grepl("^[0-9]+\\+$", label)) {
    return(as.numeric(sub("\\+$", "", label)))
  }

  if (grepl("^[0-9]+-[0-9]+$", label)) {
    return(as.numeric(sub("-[0-9]+$", "", label)))
  }

  NA_real_
}
