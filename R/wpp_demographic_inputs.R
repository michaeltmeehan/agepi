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
#' Inputs should be pre-filtered to one country, location, or entity before
#' calling this adapter. Multi-location-shaped inputs are not selected or
#' grouped here; duplicate time-age rows are rejected by [MortalitySchedule()].
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

#' Convert WPP-style central death rates to a mortality schedule
#'
#' Narrow adapter for WPP-style age-specific central death rates (`mx`). Values
#' are interpreted as annual per-capita central death rates and stored using the
#' package's [MortalitySchedule()] `annual_hazard` convention. This is an
#' approximation suitable for agepi's current continuous-time demographic
#' derivative; it does not convert survival probabilities or age-interval death
#' probabilities.
#'
#' Inputs should be pre-filtered to one country, location, or entity before
#' calling this adapter. Multi-location-shaped inputs are not selected or
#' grouped here; duplicate time-age rows are rejected.
#'
#' @param data Data frame containing WPP-like mortality data.
#' @param age_structure Target age-structure object validated by
#'   [validate_age_structure()].
#' @param time_col,age_col,mx_col Column names containing time/year, age labels,
#'   and WPP-style central death rates.
#'
#' @return An `agepi_mortality_schedule` object.
#' @export
mortality_from_wpp_mx <- function(data,
                                  age_structure,
                                  time_col,
                                  age_col,
                                  mx_col) {
  standardise_wpp_mortality(
    data = data,
    age_structure = age_structure,
    time_col = time_col,
    age_col = age_col,
    mortality_col = mx_col
  )
}

#' Convert WPP-style mortality data to a mortality schedule
#'
#' Dispatches supported WPP-style mortality quantities to the appropriate
#' narrow adapter. Currently only `quantity = "mx"` is supported. Age-specific
#' death probabilities (`qx`) and survival probabilities are intentionally
#' rejected because their period and interval conventions are ambiguous without
#' additional metadata.
#'
#' @inheritParams mortality_from_wpp_mx
#' @param time_col,age_col Column names containing time/year and age labels.
#' @param mortality_col Column containing the selected mortality quantity.
#' @param quantity Mortality quantity convention. Currently only `"mx"` is
#'   supported.
#'
#' @return An `agepi_mortality_schedule` object.
#' @export
mortality_from_wpp <- function(data,
                               age_structure,
                               time_col,
                               age_col,
                               mortality_col,
                               quantity = c("mx")) {
  quantity <- match.arg(quantity)

  mortality_from_wpp_mx(
    data = data,
    age_structure = age_structure,
    time_col = time_col,
    age_col = age_col,
    mx_col = mortality_col
  )
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

#' Convert WPP-style percent ASFR weights to fertility rates
#'
#' Converts WPP 2024 `percentASFR1dt`-style fertility schedules into an
#' [FertilitySchedule()]. The supplied age-specific values are interpreted as a
#' distribution of total fertility across maternal age groups, not as
#' per-capita fertility rates. For each time and maternal age bin, the returned
#' annual fertility rate is `TFR * fraction / age_bin_width`, where percentages
#' are first divided by 100. Maternal age bins must have finite widths;
#' open-ended maternal bins are rejected.
#'
#' This is an adapter only. It does not implement WPP projection dynamics,
#' interpolation, cohort fertility, sex-specific exposure correction, or any
#' complete WPP projection system.
#'
#' @param data Data frame containing WPP-like percent-ASFR or fraction-ASFR
#'   weights.
#' @param age_structure Target age-structure object validated by
#'   [validate_age_structure()].
#' @param time_col,age_col,weight_col Column names containing time/year, maternal
#'   age labels, and age-specific fertility weights.
#' @param tfr_data Data frame containing total fertility rates by time. Inputs
#'   should be pre-filtered to one country, location, or other entity before
#'   calling this adapter.
#' @param tfr_time_col,tfr_col Column names in `tfr_data` containing time/year
#'   and total fertility rate values.
#' @param weight_type Whether `weight_col` contains `"percent"` values summing
#'   to 100 or `"fraction"` values summing to 1 within each time point.
#' @param maternal_age_groups Character vector of age groups that must be present
#'   at every time point. Defaults to WPP's usual five-year reproductive ages
#'   from 15-19 through 45-49 when those labels are present in `age_structure`,
#'   otherwise all age groups appearing in `data`.
#' @param tolerance Numeric tolerance for checking weight sums.
#'
#' @return An `agepi_fertility_schedule` object with annual births per female
#'   person-year.
#' @examples
#' ages <- wpp_age_structure_5year()
#' weights <- data.frame(
#'   year = rep(2020, 7),
#'   age = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
#'   percent_asfr = c(5, 20, 30, 25, 15, 4, 1)
#' )
#' tfr <- data.frame(year = 2020, tfr = 2.1)
#' fertility_from_wpp_percent_asfr(
#'   weights,
#'   age_structure = ages,
#'   time_col = "year",
#'   age_col = "age",
#'   weight_col = "percent_asfr",
#'   tfr_data = tfr,
#'   tfr_time_col = "year",
#'   tfr_col = "tfr"
#' )
#' @export
fertility_from_wpp_percent_asfr <- function(data,
                                            age_structure,
                                            time_col,
                                            age_col,
                                            weight_col,
                                            tfr_data,
                                            tfr_time_col,
                                            tfr_col,
                                            weight_type = c("percent", "fraction"),
                                            maternal_age_groups = NULL,
                                            tolerance = 1e-6) {
  weight_type <- match.arg(weight_type)
  validate_wpp_standardisation_inputs(
    data,
    age_structure,
    required_columns = c(time_col, age_col, weight_col)
  )
  validate_wpp_tfr_inputs(tfr_data, tfr_time_col, tfr_col)
  validate_wpp_weight_tolerance(tolerance)

  fertility_weights <- data.frame(
    time = data[[time_col]],
    age_group = parse_wpp_age_labels(data[[age_col]], age_structure),
    weight = data[[weight_col]],
    stringsAsFactors = FALSE
  )

  maternal_age_groups <- infer_wpp_maternal_age_groups(
    maternal_age_groups,
    fertility_weights$age_group,
    age_structure
  )
  validate_wpp_asfr_weight_table(
    fertility_weights,
    age_structure,
    maternal_age_groups,
    weight_type,
    tolerance
  )

  tfr <- data.frame(
    time = tfr_data[[tfr_time_col]],
    tfr = tfr_data[[tfr_col]],
    stringsAsFactors = FALSE
  )
  validate_wpp_tfr_table(tfr)

  fertility_rates <- merge(fertility_weights, tfr, by = "time", all.x = TRUE, sort = FALSE)
  if (anyNA(fertility_rates$tfr)) {
    missing_times <- unique(fertility_rates$time[is.na(fertility_rates$tfr)])
    stop(
      "tfr_data is missing TFR value(s) for time(s): ",
      paste(missing_times, collapse = ", "),
      call. = FALSE
    )
  }

  age_widths <- age_bin_widths(age_structure, fertility_rates$age_group)
  fraction <- fertility_rates$weight
  if (identical(weight_type, "percent")) {
    fraction <- fraction / 100
  }

  fertility <- data.frame(
    time = fertility_rates$time,
    age_group = fertility_rates$age_group,
    fertility_rate = fertility_rates$tfr * fraction / age_widths,
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

#' Load a dataset from the optional wpp2024 package
#'
#' @param dataset Single dataset name to load from `wpp2024`.
#'
#' @return A data frame.
#' @export
wpp_dataset <- function(dataset) {
  if (!is.character(dataset) || length(dataset) != 1 || is.na(dataset) || !nzchar(dataset)) {
    stop("dataset must be a single non-empty string.", call. = FALSE)
  }

  if (!requireNamespace("wpp2024", quietly = TRUE)) {
    stop(
      "Package wpp2024 is required to load dataset '",
      dataset,
      "'. Install it with pak::pkg_install(\"PPgp/wpp2024\") or remotes::install_github(\"PPgp/wpp2024\"), or supply data directly.",
      call. = FALSE
    )
  }

  dataset_environment <- new.env(parent = baseenv())
  utils::data(list = dataset, package = "wpp2024", envir = dataset_environment)
  if (!exists(dataset, envir = dataset_environment, inherits = FALSE)) {
    stop("Dataset not found in wpp2024: ", dataset, call. = FALSE)
  }

  as.data.frame(get(dataset, envir = dataset_environment, inherits = FALSE))
}

#' Select one location from WPP-like data
#'
#' @param data Data frame containing WPP-like rows.
#' @param location Single location value to select.
#' @param location_col Column containing location names. Defaults to `"name"`.
#' @param columns Optional columns to retain after filtering.
#'
#' @return A filtered data frame.
#' @export
wpp_location_rows <- function(data,
                              location,
                              location_col = "name",
                              columns = names(data)) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }

  if (length(location) != 1 || anyNA(location)) {
    stop("location must be a single non-missing value.", call. = FALSE)
  }

  if (!is.character(location_col) || length(location_col) != 1 || is.na(location_col) || !nzchar(location_col)) {
    stop("location_col must be a single non-empty string.", call. = FALSE)
  }

  missing_columns <- setdiff(c(location_col, columns), names(data))
  if (length(missing_columns) > 0) {
    stop(
      "data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  data[data[[location_col]] == location, columns, drop = FALSE]
}

#' Collapse WPP-style age-specific counts into an open-ended age group
#'
#' Sums all rows with age lower bound greater than or equal to `open_age` into
#' a single `open_age+` row within each non-age grouping.
#'
#' @param data Data frame containing WPP-like age-specific counts.
#' @param value_col Column containing count values.
#' @param age_col Column containing WPP-style age labels. Defaults to `"age"`.
#' @param open_age Lower bound of the open-ended terminal age group.
#'
#' @return A data frame with collapsed count rows.
#' @export
collapse_wpp_open_age_counts <- function(data,
                                         value_col,
                                         age_col = "age",
                                         open_age) {
  collapse_wpp_open_age_values(
    data = data,
    value_col = value_col,
    age_col = age_col,
    open_age = open_age,
    mode = "count"
  )
}

#' Collapse WPP-style age-specific rates into an open-ended age group
#'
#' Keeps ages up to `open_age`, relabels the row whose lower bound is
#' `open_age` as `open_age+`, and drops finer ages above `open_age`.
#'
#' @inheritParams collapse_wpp_open_age_counts
#'
#' @return A data frame with one terminal open-ended rate row.
#' @export
collapse_wpp_open_age_rates <- function(data,
                                        value_col,
                                        age_col = "age",
                                        open_age) {
  collapse_wpp_open_age_values(
    data = data,
    value_col = value_col,
    age_col = age_col,
    open_age = open_age,
    mode = "rate"
  )
}

#' Build an agepi demographic process from WPP 2024 inputs
#'
#' Loads WPP 2024 national population, percent-ASFR fertility, total fertility,
#' mortality, and optionally migration inputs for one country. Annual WPP inputs
#' are used directly as stepwise schedules; no interpolation is introduced.
#'
#' @param country Single country/location name in WPP's `name` column.
#' @param years Numeric vector of annual years to use.
#' @param age_structure Target age structure.
#' @param migration Whether to include WPP net migration counts.
#' @param fertility_exposure_fraction Female exposure fraction passed to
#'   [build_demographic_process()].
#'
#' @return A list with top-level fields `country`, `years`, `age_structure`,
#'   `population`, `demographic_process`, `schedules`, and `inputs`.
#'   `schedules` contains the WPP-derived fertility, mortality, and optional
#'   migration schedules. `inputs` contains the selected and standardised WPP
#'   source tables used to build those objects.
#' @examples
#' if (requireNamespace("wpp2024", quietly = TRUE)) {
#'   ages <- wpp_age_structure_1year(max_age = 95)
#'   wpp_inputs <- demographic_process_from_wpp(
#'     country = "Kiribati",
#'     years = 2025:2026,
#'     age_structure = ages
#'   )
#'   demography_population_vector(wpp_inputs$population, time = 2025)
#' }
#' @export
demographic_process_from_wpp <- function(country,
                                         years,
                                         age_structure,
                                         migration = TRUE,
                                         fertility_exposure_fraction = 0.5) {
  validate_age_structure(age_structure)
  validate_wpp_workflow_arguments(country, years, migration, fertility_exposure_fraction)

  open_age <- wpp_workflow_open_age(age_structure)

  population_data <- wpp_dataset("popprojAge1dt")
  fertility_weight_data <- wpp_dataset("percentASFR1dt")
  tfr_data <- wpp_dataset("tfrproj1dt")
  mortality_data <- wpp_dataset("mx1dt")

  population_input <- wpp_location_rows(
    population_data,
    location = country,
    columns = c("name", "year", "age", "pop")
  )
  population_input$year <- as.numeric(population_input$year)
  population_input <- population_input[population_input$year %in% years, ]
  population_input$pop <- 1000 * as.numeric(population_input$pop)
  population_input <- collapse_wpp_open_age_counts(
    population_input,
    "pop",
    open_age = open_age
  )

  population <- population_from_wpp(
    data = population_input,
    age_structure = age_structure,
    time_col = "year",
    age_group_col = "age",
    population_col = "pop",
    location = country,
    location_col = "name"
  )

  fertility_weights <- wpp_location_rows(
    fertility_weight_data,
    location = country,
    columns = c("year", "age", "pasfr")
  )
  fertility_weights$year <- as.numeric(fertility_weights$year)
  fertility_weights <- fertility_weights[fertility_weights$year %in% years, ]

  tfr_source <- wpp_location_rows(
    tfr_data,
    location = country,
    columns = c("year", "tfr")
  )
  tfr_source$year <- as.numeric(tfr_source$year)
  tfr_source$tfr <- as.numeric(tfr_source$tfr)
  tfr_years <- sort(unique(fertility_weights$year))
  tfr_input <- tfr_source[tfr_source$year %in% tfr_years, , drop = FALSE]
  missing_tfr_years <- setdiff(tfr_years, tfr_input$year)
  if (length(missing_tfr_years) > 0) {
    stop(
      "Annual WPP TFR input is missing year(s): ",
      paste(missing_tfr_years, collapse = ", "),
      call. = FALSE
    )
  }

  fertility_schedule <- fertility_from_wpp_percent_asfr(
    data = fertility_weights,
    age_structure = age_structure,
    time_col = "year",
    age_col = "age",
    weight_col = "pasfr",
    tfr_data = tfr_input,
    tfr_time_col = "year",
    tfr_col = "tfr",
    weight_type = "percent",
    maternal_age_groups = wpp_workflow_maternal_age_groups(age_structure),
    tolerance = 1e-5
  )

  mortality_input <- wpp_location_rows(
    mortality_data,
    location = country,
    columns = c("year", "age", "mxB")
  )
  mortality_input$year <- as.numeric(mortality_input$year)
  mortality_input <- mortality_input[mortality_input$year %in% years, ]
  mortality_input <- collapse_wpp_open_age_rates(
    mortality_input,
    "mxB",
    open_age = open_age
  )

  mortality_schedule <- mortality_from_wpp(
    data = mortality_input,
    age_structure = age_structure,
    time_col = "year",
    age_col = "age",
    mortality_col = "mxB",
    quantity = "mx"
  )

  migration_schedule <- NULL
  migration_input <- NULL
  process_mode <- "closed"
  if (isTRUE(migration)) {
    migration_data <- wpp_dataset("migprojAge1dt")
    migration_input <- wpp_location_rows(
      migration_data,
      location = country,
      columns = c("year", "age", "mig")
    )
    migration_input$year <- as.numeric(migration_input$year)
    migration_input <- migration_input[migration_input$year %in% years, ]
    migration_input$mig <- 1000 * as.numeric(migration_input$mig)
    migration_input <- collapse_wpp_open_age_counts(
      migration_input,
      "mig",
      open_age = open_age
    )

    migration_schedule <- standardise_wpp_migration(
      data = migration_input,
      age_structure = age_structure,
      time_col = "year",
      age_col = "age",
      migration_col = "mig",
      migration_type = "count"
    )
    process_mode <- "migration"
  }

  demographic_process <- build_demographic_process(
    age_structure = age_structure,
    fertility_schedule = fertility_schedule,
    fertility_exposure_fraction = fertility_exposure_fraction,
    mortality_schedule = mortality_schedule,
    migration_schedule = migration_schedule,
    mode = process_mode
  )

  result <- list(
    country = country,
    years = years,
    age_structure = age_structure,
    population = population,
    demographic_process = demographic_process,
    schedules = list(
      fertility = fertility_schedule,
      mortality = mortality_schedule,
      migration = migration_schedule
    ),
    inputs = list(
      population = population_input,
      fertility_weights = fertility_weights,
      tfr = tfr_input,
      mortality = mortality_input,
      migration = migration_input
    )
  )
  class(result) <- c("agepi_wpp_demographic_process_inputs", "list")
  result
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

validate_wpp_tfr_inputs <- function(tfr_data, tfr_time_col, tfr_col) {
  if (!is.data.frame(tfr_data)) {
    stop("tfr_data must be a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(c(tfr_time_col, tfr_col), names(tfr_data))
  if (length(missing_columns) > 0) {
    stop(
      "tfr_data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(tfr_data)
}

validate_wpp_weight_tolerance <- function(tolerance) {
  if (!is.numeric(tolerance) || length(tolerance) != 1 || anyNA(tolerance) || !is.finite(tolerance)) {
    stop("tolerance must be a single finite numeric value.", call. = FALSE)
  }

  if (tolerance < 0) {
    stop("tolerance cannot be negative.", call. = FALSE)
  }

  invisible(tolerance)
}

infer_wpp_maternal_age_groups <- function(maternal_age_groups, observed_age_groups, age_structure) {
  if (!is.null(maternal_age_groups)) {
    maternal_age_groups <- as.character(maternal_age_groups)
    unknown_age_groups <- setdiff(maternal_age_groups, age_structure$age_groups)
    if (length(unknown_age_groups) > 0) {
      stop(
        "maternal_age_groups contains value(s) not in age_structure: ",
        paste(unknown_age_groups, collapse = ", "),
        call. = FALSE
      )
    }
    return(maternal_age_groups)
  }

  wpp_reproductive_ages <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")
  default_age_groups <- intersect(wpp_reproductive_ages, age_structure$age_groups)
  if (length(default_age_groups) > 0) {
    return(default_age_groups)
  }

  age_structure$age_groups[age_structure$age_groups %in% unique(observed_age_groups)]
}

validate_wpp_asfr_weight_table <- function(weights,
                                           age_structure,
                                           maternal_age_groups,
                                           weight_type,
                                           tolerance) {
  validate_demographic_rate_table(
    weights,
    age_structure,
    required_columns = c("time", "age_group", "weight"),
    value_column = "weight",
    value_label = "fertility weight",
    require_full_coverage = FALSE,
    allow_negative = FALSE
  )

  if (length(maternal_age_groups) == 0) {
    stop("maternal_age_groups must contain at least one age group.", call. = FALSE)
  }

  extra_age_groups <- setdiff(unique(weights$age_group), maternal_age_groups)
  if (length(extra_age_groups) > 0) {
    stop(
      "fertility weights contain age_group value(s) outside maternal_age_groups: ",
      paste(extra_age_groups, collapse = ", "),
      call. = FALSE
    )
  }

  expected_sum <- if (identical(weight_type, "percent")) 100 else 1
  for (this_time in sort(unique(weights$time))) {
    observed_age_groups <- weights$age_group[weights$time == this_time]
    missing_age_groups <- setdiff(maternal_age_groups, observed_age_groups)
    if (length(missing_age_groups) > 0) {
      stop(
        "fertility weights are missing age_group value(s) at time ",
        this_time,
        ": ",
        paste(missing_age_groups, collapse = ", "),
        call. = FALSE
      )
    }

    weight_sum <- sum(weights$weight[weights$time == this_time])
    if (abs(weight_sum - expected_sum) > tolerance) {
      stop(
        "fertility weights at time ",
        this_time,
        " must sum to ",
        expected_sum,
        ".",
        call. = FALSE
      )
    }
  }

  age_bin_widths(age_structure, maternal_age_groups)
  invisible(weights)
}

validate_wpp_tfr_table <- function(tfr) {
  if (!is.numeric(tfr$time) || anyNA(tfr$time) || any(!is.finite(tfr$time))) {
    stop("TFR time must contain finite non-missing numeric values.", call. = FALSE)
  }

  if (!is.numeric(tfr$tfr) || anyNA(tfr$tfr) || any(!is.finite(tfr$tfr))) {
    stop("TFR values must contain finite non-missing numeric values.", call. = FALSE)
  }

  if (any(tfr$tfr < 0)) {
    stop("TFR values cannot be negative.", call. = FALSE)
  }

  if (any(duplicated(tfr$time))) {
    stop("tfr_data contains duplicate time rows.", call. = FALSE)
  }

  invisible(tfr)
}

age_bin_widths <- function(age_structure, age_groups) {
  age_index <- match(age_groups, age_structure$age_groups)
  widths <- age_structure$upper_bounds[age_index] - age_structure$lower_bounds[age_index] + 1

  if (any(!is.finite(widths)) || any(widths <= 0)) {
    stop("maternal age-bin widths must be finite and positive.", call. = FALSE)
  }

  widths
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

collapse_wpp_open_age_values <- function(data, value_col, age_col, open_age, mode) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }

  if (!is.character(value_col) || length(value_col) != 1 || is.na(value_col) || !nzchar(value_col)) {
    stop("value_col must be a single non-empty string.", call. = FALSE)
  }

  if (!is.character(age_col) || length(age_col) != 1 || is.na(age_col) || !nzchar(age_col)) {
    stop("age_col must be a single non-empty string.", call. = FALSE)
  }

  if (!is.numeric(open_age) || length(open_age) != 1 || anyNA(open_age) || !is.finite(open_age)) {
    stop("open_age must be a single finite numeric value.", call. = FALSE)
  }

  missing_columns <- setdiff(c(age_col, value_col), names(data))
  if (length(missing_columns) > 0) {
    stop(
      "data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  data[[value_col]] <- suppressWarnings(as.numeric(data[[value_col]]))
  if (anyNA(data[[value_col]]) || any(!is.finite(data[[value_col]]))) {
    stop("value_col must contain finite numeric values.", call. = FALSE)
  }

  age_lower <- vapply(
    as.character(data[[age_col]]),
    parse_wpp_age_lower_bound,
    numeric(1)
  )
  if (anyNA(age_lower)) {
    stop("Unsupported WPP age label in ", age_col, ".", call. = FALSE)
  }

  data$.agepi_age_lower <- age_lower

  if (identical(mode, "rate")) {
    data <- data[data$.agepi_age_lower <= open_age, , drop = FALSE]
  }

  data[[age_col]] <- ifelse(
    data$.agepi_age_lower >= open_age,
    paste0(open_age, "+"),
    as.character(data$.agepi_age_lower)
  )

  if (identical(mode, "rate")) {
    data$.agepi_age_lower <- NULL
    row.names(data) <- NULL
    return(data)
  }

  group_cols <- setdiff(names(data), c(value_col, ".agepi_age_lower"))
  collapsed <- stats::aggregate(
    data[[value_col]],
    data[group_cols],
    sum
  )
  names(collapsed) <- c(group_cols, value_col)
  row.names(collapsed) <- NULL
  collapsed
}

validate_wpp_workflow_arguments <- function(country, years, migration, fertility_exposure_fraction) {
  if (!is.character(country) || length(country) != 1 || is.na(country) || !nzchar(country)) {
    stop("country must be a single non-empty string.", call. = FALSE)
  }

  if (!is.numeric(years) || length(years) == 0 || anyNA(years) || any(!is.finite(years))) {
    stop("years must be a non-empty numeric vector of finite non-missing values.", call. = FALSE)
  }

  if (any(duplicated(years))) {
    stop("years cannot contain duplicate values.", call. = FALSE)
  }

  if (!is.logical(migration) || length(migration) != 1 || anyNA(migration)) {
    stop("migration must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(fertility_exposure_fraction) || length(fertility_exposure_fraction) != 1 ||
      anyNA(fertility_exposure_fraction) || !is.finite(fertility_exposure_fraction) ||
      fertility_exposure_fraction < 0 || fertility_exposure_fraction > 1) {
    stop("fertility_exposure_fraction must be a finite numeric scalar between 0 and 1.", call. = FALSE)
  }

  invisible(NULL)
}

wpp_workflow_open_age <- function(age_structure) {
  open_ages <- age_structure$lower_bounds[is.infinite(age_structure$upper_bounds)]
  if (length(open_ages) != 1) {
    stop("age_structure must contain exactly one open-ended terminal age group.", call. = FALSE)
  }
  open_ages
}

wpp_workflow_maternal_age_groups <- function(age_structure) {
  age_structure$age_groups[
    age_structure$lower_bounds >= 10 &
      age_structure$lower_bounds <= 54
  ]
}
