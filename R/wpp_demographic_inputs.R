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
