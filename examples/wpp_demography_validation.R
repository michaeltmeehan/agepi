# WPP demographic validation for agepi annual-cohort demography.
#
# This script benchmarks agepi's standalone annual-cohort demographic
# machinery against WPP-derived Kiribati inputs. It is not an exact attempt to
# reproduce WPP projections. The cleanest validation path starts from WPP
# single-year age cohorts, advances the population in annual time steps with
# simulate_demography(..., ageing_policy = "annual_cohort"), and compares back
# to WPP single-year age distributions.
#
# Mechanically, births are generated from WPP-derived total fertility and
# percent-ASFR schedules, deaths are generated from WPP single-year mx schedules,
# and annual cohort ageing moves survivors to the next 1-year age group at each
# year-end. Schedule values are linearly interpolated between WPP schedule years.
# Fertility, mortality, and benchmark population inputs all use WPP single-year
# age data. The maternal ages included in fertility schedules are inferred from
# WPP percent-ASFR rows rather than hard-coded.
#
# Limitations: no sex structure, no migration, simplified female exposure,
# interpolated WPP fertility and mortality schedules, and no full WPP
# projection machinery.

# Phase 1: package loading and required function checks ---------------------

if ("package:agepi" %in% search()) {
  # Already loaded by devtools::load_all() or library(agepi).
} else if (dir.exists("R")) {
  invisible(lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), source))
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else if (requireNamespace("agepi", quietly = TRUE)) {
  library(agepi)
} else {
  stop(
    "Package agepi is not installed. Run this script from the package root ",
    "or install agepi first.",
    call. = FALSE
  )
}

required_functions <- c(
  "wpp_age_structure_1year",
  "population_from_wpp",
  "fertility_from_wpp_percent_asfr",
  "mortality_from_wpp",
  "build_demographic_process",
  "simulate_demography",
  "demography_population_table",
  "demography_population_vector",
  "standardise_wpp_migration"
)

missing_functions <- required_functions[!vapply(
  required_functions,
  exists,
  logical(1),
  mode = "function"
)]
if (length(missing_functions) > 0) {
  stop(
    "Required agepi function(s) are unavailable: ",
    paste(missing_functions, collapse = ", "),
    call. = FALSE
  )
}

if (!requireNamespace("wpp2024", quietly = TRUE)) {
  message(
    "Skipping WPP-connected demography benchmark example: optional package wpp2024 ",
    "is not installed. Install wpp2024 separately to run this integration."
  )
} else {

# Phase 2: helper functions -------------------------------------------------

section <- function(title) {
  cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
}

load_wpp_dataset <- function(name) {
  dataset_environment <- new.env(parent = baseenv())
  tryCatch(
    {
      utils::data(list = name, package = "wpp2024", envir = dataset_environment)
      if (!exists(name, envir = dataset_environment, inherits = FALSE)) {
        stop("Dataset not found in wpp2024: ", name, call. = FALSE)
      }
      as.data.frame(get(name, envir = dataset_environment, inherits = FALSE))
    },
    error = function(error) {
      stop(
        "Could not load wpp2024 dataset '", name, "': ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
}

try_load_wpp_dataset <- function(name) {
  dataset_environment <- new.env(parent = baseenv())
  result <- tryCatch(
    {
      utils::data(list = name, package = "wpp2024", envir = dataset_environment)
      if (!exists(name, envir = dataset_environment, inherits = FALSE)) {
        stop("Dataset not found in wpp2024: ", name, call. = FALSE)
      }
      as.data.frame(get(name, envir = dataset_environment, inherits = FALSE))
    },
    error = function(error) error
  )

  if (inherits(result, "error")) {
    return(list(found = FALSE, reason = conditionMessage(result), data = NULL))
  }

  list(found = TRUE, reason = NULL, data = result)
}

wpp_population_input <- function(data, country, years) {
  population_input <- data[
    data$name == country & data$year %in% years,
    c("name", "year", "age", "pop"),
    drop = FALSE
  ]
  population_input$year <- as.numeric(population_input$year)
  population_input$pop <- as.numeric(population_input$pop)
  population_input
}

build_wpp_population <- function(data, age_structure, country, years) {
  population_input <- wpp_population_input(data, country, years)
  if (nrow(population_input) == 0) {
    stop(
      "No WPP population projection rows were found for ",
      country,
      " over ",
      paste(range(years), collapse = "-"),
      ".",
      call. = FALSE
    )
  }

  population_from_wpp(
    data = population_input,
    age_structure = age_structure,
    time_col = "year",
    age_group_col = "age",
    population_col = "pop",
    location = country,
    location_col = "name"
  )
}

try_build_wpp_population <- function(data, age_structure, country, years, dataset_name) {
  unavailable <- function(reason) {
    list(found = FALSE, reason = reason, demography = NULL, table = NULL)
  }

  if (is.null(data)) {
    return(unavailable(paste0("WPP dataset ", dataset_name, " is not available locally.")))
  }

  required_columns <- c("name", "year", "age", "pop")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    return(unavailable(paste(
      "WPP dataset",
      dataset_name,
      "is missing required column(s):",
      paste(missing_columns, collapse = ", ")
    )))
  }

  result <- tryCatch(
    build_wpp_population(data, age_structure, country, years),
    error = function(error) error
  )
  if (inherits(result, "error")) {
    return(unavailable(conditionMessage(result)))
  }

  list(
    found = TRUE,
    reason = NULL,
    demography = result,
    table = demography_population_table(result)
  )
}

try_build_wpp_migration_schedule <- function(data, age_structure, country, interval_years) {
  unavailable <- function(reason) {
    list(found = FALSE, reason = reason, schedule = NULL, input = NULL)
  }

  if (is.null(data)) {
    return(unavailable("WPP dataset migprojAge1dt is not available locally."))
  }

  required_columns <- c("name", "year", "age", "mig")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    return(unavailable(paste(
      "WPP dataset migprojAge1dt is missing required column(s):",
      paste(missing_columns, collapse = ", ")
    )))
  }

  migration_input <- data[
    data$name == country & data$year %in% interval_years,
    required_columns,
    drop = FALSE
  ]
  if (nrow(migration_input) == 0) {
    return(unavailable(paste(
      "No WPP migration projection rows were found for",
      country,
      "over",
      paste(range(interval_years), collapse = "-"),
      "."
    )))
  }

  migration_input$year <- as.numeric(migration_input$year)
  migration_input$mig <- as.numeric(migration_input$mig)
  result <- tryCatch(
    standardise_wpp_migration(
      data = migration_input,
      age_structure = age_structure,
      time_col = "year",
      age_col = "age",
      migration_col = "mig",
      migration_type = "count"
    ),
    error = function(error) error
  )
  if (inherits(result, "error")) {
    return(unavailable(conditionMessage(result)))
  }

  list(found = TRUE, reason = NULL, schedule = result, input = migration_input)
}

schedule_rates_at <- function(schedule, time, value_column, age_groups, fill_value = 0) {
  values <- rep(fill_value, length(age_groups))
  names(values) <- age_groups

  if (is.null(schedule)) {
    return(values)
  }

  for (age_group in unique(schedule$data$age_group)) {
    age_rows <- schedule$data$age_group == age_group
    age_times <- schedule$data$time[age_rows]
    age_values <- schedule$data[[value_column]][age_rows]
    age_order <- order(age_times)
    values[age_group] <- stats::approx(
      x = age_times[age_order],
      y = age_values[age_order],
      xout = time,
      ties = "ordered"
    )$y
  }

  values
}

compare_total_population <- function(simulated, wpp) {
  simulated_totals <- aggregate(population ~ time, simulated, sum)
  names(simulated_totals)[names(simulated_totals) == "population"] <- "simulated_population"

  wpp_totals <- aggregate(population ~ time, wpp, sum)
  names(wpp_totals)[names(wpp_totals) == "population"] <- "wpp_population"

  comparison <- merge(simulated_totals, wpp_totals, by = "time", all = FALSE)
  comparison$absolute_difference <- abs(
    comparison$simulated_population - comparison$wpp_population
  )
  comparison$relative_difference <- comparison$absolute_difference /
    pmax(abs(comparison$wpp_population), .Machine$double.eps)
  comparison
}

migration_adjusted_total_comparison <- function(closed_simulated,
                                                wpp,
                                                migration_schedule,
                                                migration_simulated = NULL,
                                                simulation_years) {
  closed_totals <- aggregate(population ~ time, closed_simulated, sum)
  names(closed_totals)[names(closed_totals) == "population"] <- "closed_simulated_total"

  wpp_totals <- aggregate(population ~ time, wpp, sum)
  names(wpp_totals)[names(wpp_totals) == "population"] <- "wpp_total_population"

  comparison <- merge(wpp_totals, closed_totals, by = "time", all = FALSE, sort = FALSE)

  interval_years <- simulation_years[-length(simulation_years)]
  annual_migration <- rep(0, length(interval_years))
  names(annual_migration) <- as.character(interval_years)
  if (!is.null(migration_schedule)) {
    annual_migration <- vapply(
      interval_years,
      function(year) {
        sum(schedule_rates_at(
          migration_schedule,
          year,
          "migration_count",
          migration_schedule$age_groups,
          fill_value = 0
        ))
      },
      numeric(1)
    )
    names(annual_migration) <- as.character(interval_years)
  }

  comparison$cumulative_wpp_net_migration <- vapply(
    comparison$time,
    function(year) {
      elapsed_years <- interval_years[interval_years < year]
      if (length(elapsed_years) == 0) {
        return(0)
      }
      sum(annual_migration[as.character(elapsed_years)])
    },
    numeric(1)
  )

  if (!is.null(migration_simulated)) {
    migration_totals <- aggregate(population ~ time, migration_simulated, sum)
    names(migration_totals)[names(migration_totals) == "population"] <- "migration_adjusted_simulated_total"
    comparison <- merge(comparison, migration_totals, by = "time", all.x = TRUE, sort = FALSE)
  } else {
    comparison$migration_adjusted_simulated_total <-
      comparison$closed_simulated_total + comparison$cumulative_wpp_net_migration
  }

  comparison$closed_signed_error <-
    comparison$closed_simulated_total - comparison$wpp_total_population
  comparison$closed_absolute_error <- abs(comparison$closed_signed_error)
  comparison$closed_relative_error <- comparison$closed_signed_error /
    pmax(abs(comparison$wpp_total_population), .Machine$double.eps)
  comparison$closed_absolute_relative_error <- comparison$closed_absolute_error /
    pmax(abs(comparison$wpp_total_population), .Machine$double.eps)
  comparison$migration_adjusted_signed_error <-
    comparison$migration_adjusted_simulated_total - comparison$wpp_total_population
  comparison$migration_adjusted_absolute_error <- abs(comparison$migration_adjusted_signed_error)
  comparison$migration_adjusted_relative_error <-
    comparison$migration_adjusted_signed_error /
    pmax(abs(comparison$wpp_total_population), .Machine$double.eps)
  comparison$migration_adjusted_absolute_relative_error <-
    comparison$migration_adjusted_absolute_error /
    pmax(abs(comparison$wpp_total_population), .Machine$double.eps)

  comparison
}

headline_total_comparison <- function(total_comparison, years) {
  result <- total_comparison[total_comparison$time %in% years, c(
    "time",
    "wpp_total_population",
    "closed_simulated_total",
    "closed_signed_error",
    "closed_relative_error",
    "migration_adjusted_simulated_total",
    "migration_adjusted_signed_error",
    "migration_adjusted_relative_error",
    "cumulative_wpp_net_migration"
  )]
  names(result) <- c(
    "year",
    "wpp_total_population",
    "closed_simulated_total",
    "closed_signed_error",
    "closed_relative_error",
    "migration_simulated_total",
    "migration_signed_error",
    "migration_relative_error",
    "cumulative_wpp_net_migration"
  )
  result
}

headline_improvement_metrics <- function(total_comparison, final_year) {
  final_row <- total_comparison[total_comparison$time == final_year, , drop = FALSE]
  percent_reduction <- 100 * (
    final_row$closed_absolute_error - final_row$migration_adjusted_absolute_error
  ) / pmax(final_row$closed_absolute_error, .Machine$double.eps)

  data.frame(
    metric = c(
      "closed_2050_signed_error",
      "closed_2050_relative_error",
      "migration_2050_signed_error",
      "migration_2050_relative_error",
      "percent_reduction_absolute_2050_total_error",
      "closed_max_absolute_relative_error_all_years",
      "migration_max_absolute_relative_error_all_years"
    ),
    value = c(
      final_row$closed_signed_error,
      final_row$closed_relative_error,
      final_row$migration_adjusted_signed_error,
      final_row$migration_adjusted_relative_error,
      percent_reduction,
      max(total_comparison$closed_absolute_relative_error, na.rm = TRUE),
      max(total_comparison$migration_adjusted_absolute_relative_error, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

age_error_table <- function(simulated, wpp, year, tiny_denominator) {
  simulated_year <- simulated[simulated$time == year, c("age_group", "population")]
  wpp_year <- wpp[wpp$time == year, c("age_group", "population")]
  comparison <- merge(
    simulated_year,
    wpp_year,
    by = "age_group",
    suffixes = c("_simulated", "_wpp"),
    sort = FALSE
  )
  comparison$signed_error <- comparison$population_simulated - comparison$population_wpp
  comparison$absolute_error <- abs(comparison$signed_error)
  relative_rows <- abs(comparison$population_wpp) > tiny_denominator
  comparison$relative_error <- NA_real_
  comparison$relative_error[relative_rows] <- comparison$signed_error[relative_rows] /
    comparison$population_wpp[relative_rows]
  comparison
}

age_residual_summary <- function(comparison, year) {
  over_row <- comparison[which.max(comparison$signed_error), , drop = FALSE]
  under_row <- comparison[which.min(comparison$signed_error), , drop = FALSE]

  data.frame(
    time = year,
    sum_positive_errors = sum(comparison$signed_error[comparison$signed_error > 0]),
    sum_negative_errors = sum(comparison$signed_error[comparison$signed_error < 0]),
    total_absolute_error = sum(comparison$absolute_error),
    age_largest_overprediction = over_row$age_group,
    largest_overprediction = over_row$signed_error,
    age_largest_underprediction = under_row$age_group,
    largest_underprediction = under_row$signed_error,
    stringsAsFactors = FALSE
  )
}

signed_age_residual_diagnostics <- function(simulated, wpp, years, tiny_denominator) {
  residuals <- lapply(
    years,
    function(year) {
      comparison <- age_error_table(simulated, wpp, year, tiny_denominator)
      comparison$time <- year
      comparison[, c(
        "time",
        "age_group",
        "population_simulated",
        "population_wpp",
        "signed_error",
        "absolute_error",
        "relative_error"
      )]
    }
  )
  do.call(rbind, residuals)
}

signed_age_residual_summaries <- function(simulated, wpp, years, tiny_denominator) {
  summaries <- lapply(
    years,
    function(year) {
      comparison <- age_error_table(simulated, wpp, year, tiny_denominator)
      age_residual_summary(comparison, year)
    }
  )
  do.call(rbind, summaries)
}

age_distribution_diagnostics <- function(simulated, wpp, years, tiny_denominator) {
  diagnostics <- lapply(
    years,
    function(year) {
      comparison <- age_error_table(
        simulated = simulated,
        wpp = wpp,
        year = year,
        tiny_denominator = tiny_denominator
      )
      data.frame(
        time = year,
        total_absolute_age_distribution_error = sum(comparison$absolute_error),
        mean_absolute_age_error = mean(comparison$absolute_error),
        maximum_absolute_age_error = max(comparison$absolute_error),
        age_largest_absolute_error = comparison$age_group[which.max(comparison$absolute_error)],
        maximum_relative_age_error = max(abs(comparison$relative_error), na.rm = TRUE),
        ages_excluded_from_relative_error = sum(is.na(comparison$relative_error)),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, diagnostics)
}

side_by_side_age_distribution_diagnostics <- function(closed_simulated,
                                                      migration_simulated,
                                                      wpp,
                                                      years,
                                                      tiny_denominator) {
  closed <- age_distribution_diagnostics(
    closed_simulated,
    wpp,
    years,
    tiny_denominator
  )
  names(closed)[names(closed) == "total_absolute_age_distribution_error"] <-
    "closed_total_absolute_age_error"
  names(closed)[names(closed) == "mean_absolute_age_error"] <-
    "closed_mean_absolute_age_error"
  names(closed)[names(closed) == "maximum_absolute_age_error"] <-
    "closed_maximum_absolute_age_error"
  names(closed)[names(closed) == "age_largest_absolute_error"] <-
    "closed_age_largest_absolute_error"

  if (is.null(migration_simulated)) {
    result <- closed[, c(
      "time",
      "closed_total_absolute_age_error",
      "closed_mean_absolute_age_error",
      "closed_maximum_absolute_age_error",
      "closed_age_largest_absolute_error"
    )]
    result$migration_total_absolute_age_error <- NA_real_
    result$migration_mean_absolute_age_error <- NA_real_
    result$migration_maximum_absolute_age_error <- NA_real_
    result$migration_age_largest_absolute_error <- NA_character_
    return(result)
  }

  migration <- age_distribution_diagnostics(
    migration_simulated,
    wpp,
    years,
    tiny_denominator
  )
  names(migration)[names(migration) == "total_absolute_age_distribution_error"] <-
    "migration_total_absolute_age_error"
  names(migration)[names(migration) == "mean_absolute_age_error"] <-
    "migration_mean_absolute_age_error"
  names(migration)[names(migration) == "maximum_absolute_age_error"] <-
    "migration_maximum_absolute_age_error"
  names(migration)[names(migration) == "age_largest_absolute_error"] <-
    "migration_age_largest_absolute_error"

  merge(
    closed[, c(
      "time",
      "closed_total_absolute_age_error",
      "closed_mean_absolute_age_error",
      "closed_maximum_absolute_age_error",
      "closed_age_largest_absolute_error"
    )],
    migration[, c(
      "time",
      "migration_total_absolute_age_error",
      "migration_mean_absolute_age_error",
      "migration_maximum_absolute_age_error",
      "migration_age_largest_absolute_error"
    )],
    by = "time",
    sort = FALSE
  )
}

top_age_errors <- function(simulated, wpp, year, tiny_denominator, n = 5) {
  comparison <- age_error_table(
    simulated = simulated,
    wpp = wpp,
    year = year,
    tiny_denominator = tiny_denominator
  )
  columns <- c(
    "age_group",
    "population_simulated",
    "population_wpp",
    "signed_error",
    "absolute_error",
    "relative_error"
  )

  top_absolute <- head(
    comparison[order(-comparison$absolute_error), columns, drop = FALSE],
    n
  )
  relative_comparison <- comparison[!is.na(comparison$relative_error), , drop = FALSE]
  top_relative <- head(
    relative_comparison[order(-relative_comparison$relative_error), columns, drop = FALSE],
    n
  )

  list(absolute = top_absolute, relative = top_relative)
}

simulate_annual_cohort <- function(process, initial_population, simulation_years) {
  simulate_demography(
    process = process,
    initial_state = initial_population,
    times = simulation_years,
    time_policy = "linear",
    method = "euler",
    ageing_policy = "annual_cohort"
  )
}

build_fertility_schedule_from_single_year_wpp <- function(fertility_weights_data,
                                                          tfr_data,
                                                          country,
                                                          schedule_years,
                                                          age_structure) {
  fertility_weights_input <- fertility_weights_data[
    fertility_weights_data$name == country &
      fertility_weights_data$year %in% schedule_years,
    c("year", "age", "pasfr"),
    drop = FALSE
  ]
  fertility_weights_input$year <- as.numeric(fertility_weights_input$year)
  fertility_weights_input$pasfr <- as.numeric(fertility_weights_input$pasfr)

  inferred_maternal_age_groups <- unique(as.character(fertility_weights_input$age))
  inferred_maternal_age_groups <- age_structure$age_groups[
    match(
      sort(as.numeric(inferred_maternal_age_groups)),
      age_structure$lower_bounds
    )
  ]

  tfr_input <- tfr_data[
    tfr_data$name == country & tfr_data$year %in% schedule_years,
    c("year", "tfr"),
    drop = FALSE
  ]
  tfr_input$year <- as.numeric(tfr_input$year)
  tfr_input$tfr <- as.numeric(tfr_input$tfr)

  fertility_schedule <- fertility_from_wpp_percent_asfr(
    data = fertility_weights_input,
    age_structure = age_structure,
    time_col = "year",
    age_col = "age",
    weight_col = "pasfr",
    tfr_data = tfr_input,
    tfr_time_col = "year",
    tfr_col = "tfr",
    weight_type = "percent",
    maternal_age_groups = inferred_maternal_age_groups,
    tolerance = 1e-4
  )

  list(
    schedule = fertility_schedule,
    maternal_age_groups = inferred_maternal_age_groups
  )
}

population_matrix_by_year <- function(simulated, simulation_years, age_groups) {
  matrix <- t(vapply(
    simulation_years,
    function(year) {
      year_population <- simulated[simulated$time == year, ]
      stats::setNames(
        year_population$population[match(age_groups, year_population$age_group)],
        age_groups
      )
    },
    numeric(length(age_groups))
  ))
  rownames(matrix) <- as.character(simulation_years)
  colnames(matrix) <- age_groups
  matrix
}

age_band_lookup <- function(age_structure, bands) {
  band_names <- names(bands)
  band <- rep(NA_character_, age_structure$n_age_groups)

  for (band_name in band_names) {
    limits <- bands[[band_name]]
    in_band <- age_structure$lower_bounds >= limits[1] &
      age_structure$lower_bounds <= limits[2]
    band[in_band] <- band_name
  }

  data.frame(
    age_group = age_structure$age_groups,
    age_band = factor(band, levels = band_names),
    stringsAsFactors = FALSE
  )
}

broad_age_band_residuals <- function(age_residuals, lookup) {
  merged <- merge(age_residuals, lookup, by = "age_group", all.x = TRUE, sort = FALSE)
  aggregate(
    cbind(signed_error, absolute_error) ~ time + age_band,
    merged,
    sum,
    na.rm = TRUE
  )
}

broad_age_band_population_comparison <- function(closed_simulated,
                                                 migration_simulated,
                                                 wpp,
                                                 migration_schedule,
                                                 years,
                                                 interval_years,
                                                 age_groups,
                                                 lookup) {
  aggregate_band_population <- function(data, value_name) {
    merged <- merge(data[, c("time", "age_group", "population")], lookup, by = "age_group", sort = FALSE)
    result <- aggregate(population ~ time + age_band, merged[merged$time %in% years, ], sum)
    names(result)[names(result) == "population"] <- value_name
    result
  }

  wpp_bands <- aggregate_band_population(wpp, "wpp_population")
  closed_bands <- aggregate_band_population(closed_simulated, "closed_simulated_population")
  if (is.null(migration_simulated)) {
    migration_bands <- closed_bands[, c("time", "age_band"), drop = FALSE]
    migration_bands$migration_simulated_population <- NA_real_
  } else {
    migration_bands <- aggregate_band_population(migration_simulated, "migration_simulated_population")
  }

  cumulative_migration <- do.call(
    rbind,
    lapply(
      years,
      function(year) {
        migration_by_age <- data.frame(
          time = year,
          age_group = age_groups,
          cumulative_net_migration = as.numeric(cumulative_migration_by_age(
            migration_schedule,
            year,
            interval_years,
            age_groups
          )),
          stringsAsFactors = FALSE
        )
        migration_by_age <- merge(migration_by_age, lookup, by = "age_group", sort = FALSE)
        aggregate(cumulative_net_migration ~ time + age_band, migration_by_age, sum)
      }
    )
  )

  result <- merge(wpp_bands, closed_bands, by = c("time", "age_band"), sort = FALSE)
  result <- merge(result, migration_bands, by = c("time", "age_band"), sort = FALSE)
  result <- merge(result, cumulative_migration, by = c("time", "age_band"), sort = FALSE)
  result$closed_signed_error <- result$closed_simulated_population - result$wpp_population
  result$migration_signed_error <- result$migration_simulated_population - result$wpp_population
  result[, c(
    "time",
    "age_band",
    "wpp_population",
    "closed_simulated_population",
    "closed_signed_error",
    "migration_simulated_population",
    "migration_signed_error",
    "cumulative_net_migration"
  )]
}

broad_age_band_2050_summary <- function(broad_age_band_comparison, year) {
  result <- broad_age_band_comparison[
    broad_age_band_comparison$time == year,
    c(
      "age_band",
      "wpp_population",
      "closed_signed_error",
      "migration_signed_error",
      "cumulative_net_migration"
    ),
    drop = FALSE
  ]
  result$absolute_error_reduction <-
    abs(result$closed_signed_error) - abs(result$migration_signed_error)
  result
}

cumulative_migration_by_age <- function(migration_schedule, year, interval_years, age_groups) {
  values <- stats::setNames(numeric(length(age_groups)), age_groups)
  if (is.null(migration_schedule)) {
    return(values)
  }

  elapsed_years <- interval_years[interval_years < year]
  for (interval_year in elapsed_years) {
    values <- values + schedule_rates_at(
      migration_schedule,
      interval_year,
      "migration_count",
      age_groups,
      fill_value = 0
    )
  }
  values
}

cumulative_migration_band_diagnostics <- function(migration_schedule,
                                                  closed_simulated,
                                                  migration_simulated,
                                                  wpp,
                                                  years,
                                                  interval_years,
                                                  age_groups,
                                                  lookup,
                                                  tiny_denominator) {
  diagnostics <- lapply(
    years,
    function(year) {
      migration_by_age <- data.frame(
        time = year,
        age_group = age_groups,
        cumulative_net_migration = as.numeric(cumulative_migration_by_age(
          migration_schedule,
          year,
          interval_years,
          age_groups
        )),
        stringsAsFactors = FALSE
      )
      migration_by_age <- merge(migration_by_age, lookup, by = "age_group", sort = FALSE)
      migration_bands <- aggregate(
        cumulative_net_migration ~ time + age_band,
        migration_by_age,
        sum,
        na.rm = TRUE
      )

      closed_errors <- age_error_table(closed_simulated, wpp, year, tiny_denominator)
      closed_errors$time <- year
      closed_bands <- broad_age_band_residuals(
        closed_errors[, c("time", "age_group", "signed_error", "absolute_error")],
        lookup
      )
      names(closed_bands)[names(closed_bands) == "signed_error"] <- "closed_signed_error"
      names(closed_bands)[names(closed_bands) == "absolute_error"] <- "closed_absolute_error"

      if (!is.null(migration_simulated)) {
        adjusted_errors <- age_error_table(migration_simulated, wpp, year, tiny_denominator)
        adjusted_errors$time <- year
        adjusted_bands <- broad_age_band_residuals(
          adjusted_errors[, c("time", "age_group", "signed_error", "absolute_error")],
          lookup
        )
        names(adjusted_bands)[names(adjusted_bands) == "signed_error"] <-
          "migration_adjusted_signed_error"
        names(adjusted_bands)[names(adjusted_bands) == "absolute_error"] <-
          "migration_adjusted_absolute_error"
      } else {
        adjusted_bands <- closed_bands[, c("time", "age_band"), drop = FALSE]
        adjusted_bands$migration_adjusted_signed_error <-
          closed_bands$closed_signed_error + migration_bands$cumulative_net_migration
        adjusted_bands$migration_adjusted_absolute_error <-
          abs(adjusted_bands$migration_adjusted_signed_error)
      }

      result <- merge(migration_bands, closed_bands, by = c("time", "age_band"), sort = FALSE)
      result <- merge(result, adjusted_bands, by = c("time", "age_band"), sort = FALSE)
      result$error_reduction <- result$closed_absolute_error -
        result$migration_adjusted_absolute_error
      result$migration_explains_positive_residual <-
        result$closed_signed_error > 0 & result$cumulative_net_migration < 0
      result
    }
  )

  do.call(rbind, diagnostics)
}

annual_balance_table <- function(simulated,
                                 process,
                                 fertility_schedule,
                                 mortality_schedule,
                                 simulation_years,
                                 age_groups,
                                 wpp = NULL) {
  population_matrix <- population_matrix_by_year(simulated, simulation_years, age_groups)
  annual_interval_years <- simulation_years[-length(simulation_years)]

  annual_births <- vapply(
    annual_interval_years,
    function(year) {
      fertility_rates <- schedule_rates_at(
        fertility_schedule,
        year,
        "fertility_rate",
        age_groups,
        fill_value = 0
      )
      sum(
        fertility_rates *
          process$fertility_exposure_fraction *
          population_matrix[as.character(year), ]
      )
    },
    numeric(1)
  )

  annual_deaths <- vapply(
    annual_interval_years,
    function(year) {
      mortality_rates <- schedule_rates_at(
        mortality_schedule,
        year,
        "mortality_rate",
        age_groups,
        fill_value = 0
      )
      sum((1 - exp(-mortality_rates)) * population_matrix[as.character(year), ])
    },
    numeric(1)
  )

  simulated_totals <- aggregate(population ~ time, simulated, sum)
  annual_totals <- simulated_totals$population[
    match(annual_interval_years, simulated_totals$time)
  ]
  annual_next_totals <- simulated_totals$population[
    match(annual_interval_years + 1, simulated_totals$time)
  ]

  balance <- data.frame(
    year = annual_interval_years,
    starting_simulated_total = annual_totals,
    simulated_births = annual_births,
    simulated_deaths = annual_deaths,
    simulated_net_change = annual_next_totals - annual_totals,
    crude_birth_rate_per_1000 = annual_births / annual_totals * 1000,
    crude_death_rate_per_1000 = annual_deaths / annual_totals * 1000
  )

  if (!is.null(wpp)) {
    wpp_totals <- aggregate(population ~ time, wpp, sum)
    wpp_start <- wpp_totals$population[match(annual_interval_years, wpp_totals$time)]
    wpp_end <- wpp_totals$population[match(annual_interval_years + 1, wpp_totals$time)]
    balance$wpp_total_net_change <- wpp_end - wpp_start
  } else {
    balance$wpp_total_net_change <- NA_real_
  }
  balance$simulated_minus_wpp_net_change <-
    balance$simulated_net_change - balance$wpp_total_net_change
  balance
}

find_wpp_interval_columns <- function(data, central_col = "pop") {
  if (is.null(data)) {
    return(NULL)
  }

  lower_candidates <- c(
    paste0(central_col, "_low"),
    paste0(central_col, "_lower"),
    paste0(central_col, "Low"),
    paste0(central_col, "Lower")
  )
  upper_candidates <- c(
    paste0(central_col, "_high"),
    paste0(central_col, "_upper"),
    paste0(central_col, "High"),
    paste0(central_col, "Upper")
  )
  lower_col <- lower_candidates[lower_candidates %in% names(data)][1]
  upper_col <- upper_candidates[upper_candidates %in% names(data)][1]

  if (is.na(lower_col) || is.na(upper_col)) {
    return(NULL)
  }
  if (!is.numeric(data[[lower_col]]) || !is.numeric(data[[upper_col]])) {
    return(NULL)
  }

  list(lower = lower_col, upper = upper_col)
}

build_wpp_interval_table <- function(data,
                                     age_structure,
                                     country,
                                     years,
                                     lower_col,
                                     upper_col) {
  interval_input <- data[
    data$name == country & data$year %in% years,
    c("name", "year", "age", lower_col, upper_col),
    drop = FALSE
  ]
  names(interval_input)[names(interval_input) == lower_col] <- "lower"
  names(interval_input)[names(interval_input) == upper_col] <- "upper"
  interval_input$year <- as.numeric(interval_input$year)
  interval_input$lower <- as.numeric(interval_input$lower)
  interval_input$upper <- as.numeric(interval_input$upper)

  lower_demography <- population_from_wpp(
    data = interval_input[, c("name", "year", "age", "lower")],
    age_structure = age_structure,
    time_col = "year",
    age_group_col = "age",
    population_col = "lower",
    location = country,
    location_col = "name"
  )
  upper_demography <- population_from_wpp(
    data = interval_input[, c("name", "year", "age", "upper")],
    age_structure = age_structure,
    time_col = "year",
    age_group_col = "age",
    population_col = "upper",
    location = country,
    location_col = "name"
  )

  lower_table <- demography_population_table(lower_demography)
  upper_table <- demography_population_table(upper_demography)
  names(lower_table)[names(lower_table) == "population"] <- "lower_population"
  names(upper_table)[names(upper_table) == "population"] <- "upper_population"
  merge(lower_table, upper_table, by = c("time", "age_group"), sort = FALSE)
}

interval_distance <- function(value, lower, upper) {
  pmax(lower - value, value - upper, 0)
}

total_interval_coverage <- function(simulated, interval_table) {
  simulated_totals <- aggregate(population ~ time, simulated, sum)
  lower_totals <- aggregate(lower_population ~ time, interval_table, sum)
  upper_totals <- aggregate(upper_population ~ time, interval_table, sum)
  coverage <- merge(simulated_totals, lower_totals, by = "time", sort = FALSE)
  coverage <- merge(coverage, upper_totals, by = "time", sort = FALSE)
  names(coverage)[names(coverage) == "population"] <- "simulated_population"
  coverage$inside_interval <- coverage$simulated_population >= coverage$lower_population &
    coverage$simulated_population <= coverage$upper_population
  coverage$distance_outside_interval <- interval_distance(
    coverage$simulated_population,
    coverage$lower_population,
    coverage$upper_population
  )
  coverage
}

side_by_side_total_interval_coverage <- function(closed_simulated,
                                                 migration_simulated,
                                                 interval_table) {
  closed <- total_interval_coverage(closed_simulated, interval_table)
  closed <- closed[, c(
    "time",
    "simulated_population",
    "lower_population",
    "upper_population",
    "inside_interval",
    "distance_outside_interval"
  )]
  names(closed) <- c(
    "time",
    "closed_simulated_population",
    "lower_population",
    "upper_population",
    "closed_inside_interval",
    "closed_distance_outside_interval"
  )

  if (is.null(migration_simulated)) {
    closed$migration_simulated_population <- NA_real_
    closed$migration_inside_interval <- NA
    closed$migration_distance_outside_interval <- NA_real_
    return(closed[, c(
      "time",
      "lower_population",
      "upper_population",
      "closed_simulated_population",
      "closed_inside_interval",
      "closed_distance_outside_interval",
      "migration_simulated_population",
      "migration_inside_interval",
      "migration_distance_outside_interval"
    )])
  }

  migration <- total_interval_coverage(migration_simulated, interval_table)
  migration <- migration[, c(
    "time",
    "simulated_population",
    "inside_interval",
    "distance_outside_interval"
  )]
  names(migration) <- c(
    "time",
    "migration_simulated_population",
    "migration_inside_interval",
    "migration_distance_outside_interval"
  )

  merge(closed, migration, by = "time", sort = FALSE)[, c(
    "time",
    "lower_population",
    "upper_population",
    "closed_simulated_population",
    "closed_inside_interval",
    "closed_distance_outside_interval",
    "migration_simulated_population",
    "migration_inside_interval",
    "migration_distance_outside_interval"
  )]
}

age_interval_coverage <- function(simulated, interval_table, years) {
  comparison <- merge(
    simulated[, c("time", "age_group", "population")],
    interval_table,
    by = c("time", "age_group"),
    sort = FALSE
  )
  comparison <- comparison[comparison$time %in% years, , drop = FALSE]
  comparison$inside_interval <- comparison$population >= comparison$lower_population &
    comparison$population <= comparison$upper_population
  comparison$below_lower <- comparison$population < comparison$lower_population
  comparison$above_upper <- comparison$population > comparison$upper_population
  comparison$distance_below <- pmax(comparison$lower_population - comparison$population, 0)
  comparison$distance_above <- pmax(comparison$population - comparison$upper_population, 0)
  comparison$distance_outside_interval <- pmax(
    comparison$distance_below,
    comparison$distance_above
  )

  summaries <- lapply(
    years,
    function(year) {
      year_rows <- comparison[comparison$time == year, , drop = FALSE]
      furthest_below <- year_rows[which.max(year_rows$distance_below), , drop = FALSE]
      furthest_above <- year_rows[which.max(year_rows$distance_above), , drop = FALSE]
      data.frame(
        time = year,
        ages_inside_interval = sum(year_rows$inside_interval),
        proportion_ages_inside_interval = mean(year_rows$inside_interval),
        ages_below_lower = sum(year_rows$below_lower),
        ages_above_upper = sum(year_rows$above_upper),
        largest_distance_outside_interval = max(year_rows$distance_outside_interval),
        age_furthest_below_interval = if (max(year_rows$distance_below) > 0) {
          furthest_below$age_group
        } else {
          NA_character_
        },
        age_furthest_above_interval = if (max(year_rows$distance_above) > 0) {
          furthest_above$age_group
        } else {
          NA_character_
        },
        stringsAsFactors = FALSE
      )
    }
  )

  list(details = comparison, summary = do.call(rbind, summaries))
}

side_by_side_age_interval_coverage <- function(closed_simulated,
                                               migration_simulated,
                                               interval_table,
                                               years) {
  closed <- age_interval_coverage(closed_simulated, interval_table, years)$summary
  names(closed)[names(closed) == "ages_inside_interval"] <- "closed_ages_inside_interval"
  names(closed)[names(closed) == "proportion_ages_inside_interval"] <-
    "closed_proportion_ages_inside_interval"
  names(closed)[names(closed) == "ages_below_lower"] <- "closed_ages_below_lower"
  names(closed)[names(closed) == "ages_above_upper"] <- "closed_ages_above_upper"
  names(closed)[names(closed) == "largest_distance_outside_interval"] <-
    "closed_largest_distance_outside_interval"

  closed <- closed[, c(
    "time",
    "closed_ages_inside_interval",
    "closed_proportion_ages_inside_interval",
    "closed_ages_below_lower",
    "closed_ages_above_upper",
    "closed_largest_distance_outside_interval"
  )]

  if (is.null(migration_simulated)) {
    closed$migration_ages_inside_interval <- NA_integer_
    closed$migration_proportion_ages_inside_interval <- NA_real_
    closed$migration_ages_below_lower <- NA_integer_
    closed$migration_ages_above_upper <- NA_integer_
    closed$migration_largest_distance_outside_interval <- NA_real_
    return(closed)
  }

  migration <- age_interval_coverage(migration_simulated, interval_table, years)$summary
  names(migration)[names(migration) == "ages_inside_interval"] <- "migration_ages_inside_interval"
  names(migration)[names(migration) == "proportion_ages_inside_interval"] <-
    "migration_proportion_ages_inside_interval"
  names(migration)[names(migration) == "ages_below_lower"] <- "migration_ages_below_lower"
  names(migration)[names(migration) == "ages_above_upper"] <- "migration_ages_above_upper"
  names(migration)[names(migration) == "largest_distance_outside_interval"] <-
    "migration_largest_distance_outside_interval"

  merge(
    closed,
    migration[, c(
      "time",
      "migration_ages_inside_interval",
      "migration_proportion_ages_inside_interval",
      "migration_ages_below_lower",
      "migration_ages_above_upper",
      "migration_largest_distance_outside_interval"
    )],
    by = "time",
    sort = FALSE
  )
}

broad_age_band_interval_coverage <- function(simulated, interval_table, lookup, years) {
  comparison <- merge(
    simulated[, c("time", "age_group", "population")],
    interval_table,
    by = c("time", "age_group"),
    sort = FALSE
  )
  comparison <- merge(comparison, lookup, by = "age_group", all.x = TRUE, sort = FALSE)
  comparison <- comparison[comparison$time %in% years, , drop = FALSE]
  band_totals <- aggregate(
    cbind(population, lower_population, upper_population) ~ time + age_band,
    comparison,
    sum,
    na.rm = TRUE
  )
  band_totals$inside_interval <- band_totals$population >= band_totals$lower_population &
    band_totals$population <= band_totals$upper_population
  band_totals$below_lower <- band_totals$population < band_totals$lower_population
  band_totals$above_upper <- band_totals$population > band_totals$upper_population
  band_totals$distance_outside_interval <- interval_distance(
    band_totals$population,
    band_totals$lower_population,
    band_totals$upper_population
  )
  band_totals
}

side_by_side_broad_age_band_interval_coverage <- function(closed_simulated,
                                                          migration_simulated,
                                                          interval_table,
                                                          lookup,
                                                          years) {
  closed <- broad_age_band_interval_coverage(closed_simulated, interval_table, lookup, years)
  names(closed)[names(closed) == "population"] <- "closed_population"
  names(closed)[names(closed) == "inside_interval"] <- "closed_inside_interval"
  names(closed)[names(closed) == "distance_outside_interval"] <-
    "closed_distance_outside_interval"
  closed <- closed[, c(
    "time",
    "age_band",
    "lower_population",
    "upper_population",
    "closed_population",
    "closed_inside_interval",
    "closed_distance_outside_interval"
  )]

  if (is.null(migration_simulated)) {
    closed$migration_population <- NA_real_
    closed$migration_inside_interval <- NA
    closed$migration_distance_outside_interval <- NA_real_
    return(closed)
  }

  migration <- broad_age_band_interval_coverage(
    migration_simulated,
    interval_table,
    lookup,
    years
  )
  names(migration)[names(migration) == "population"] <- "migration_population"
  names(migration)[names(migration) == "inside_interval"] <- "migration_inside_interval"
  names(migration)[names(migration) == "distance_outside_interval"] <-
    "migration_distance_outside_interval"

  merge(
    closed,
    migration[, c(
      "time",
      "age_band",
      "migration_population",
      "migration_inside_interval",
      "migration_distance_outside_interval"
    )],
    by = c("time", "age_band"),
    sort = FALSE
  )
}

print_limitations <- function() {
  cat("This benchmark validates agepi machinery, not exact WPP reproduction.\n")
  cat("No sex structure: fertility exposure is approximated with a scalar fraction.\n")
  cat("WPP fertility and mortality schedules are interpolated between available schedule years.\n")
  cat("WPP age-specific migration is treated as additive annual net migration counts when available.\n")
  cat("Mortality uses WPP mx as annual hazards within agepi's schedule convention.\n")
  cat("The script does not implement WPP's full projection machinery or assumptions.\n")
  cat("WPP uncertainty intervals are sanity checks, not pass/fail thresholds.\n")
}

# Phase 3: configuration ----------------------------------------------------

country <- "Kiribati"
simulation_years <- 2025:2050
simulation_interval_years <- simulation_years[-length(simulation_years)]
schedule_years <- seq(2023, 2053, by = 5)
diagnostic_years <- c(2030, 2040, 2050)
headline_years <- c(2025, 2030, 2040, 2050)
tiny_denominator <- 1e-8
female_exposure_fraction <- 0.5
verbose_output <- identical(tolower(Sys.getenv("AGEPI_WPP_VALIDATION_VERBOSE")), "true")
diagnostic_age_bands <- list(
  "0-4" = c(0, 4),
  "5-14" = c(5, 14),
  "15-24" = c(15, 24),
  "25-44" = c(25, 44),
  "45-64" = c(45, 64),
  "65+" = c(65, Inf)
)

internal_age_structure <- wpp_age_structure_1year(max_age = 100)
internal_age_band_lookup <- age_band_lookup(internal_age_structure, diagnostic_age_bands)

# Phase 4: WPP data loading -------------------------------------------------

population_data_1year <- load_wpp_dataset("popprojAge1dt")
fertility_weights_data <- load_wpp_dataset("percentASFR1dt")
tfr_data <- load_wpp_dataset("tfrproj5dt")
mortality_data_1year <- load_wpp_dataset("mx1dt")
migration_data_1year_result <- try_load_wpp_dataset("migprojAge1dt")

wpp_population_1year <- try_build_wpp_population(
  data = population_data_1year,
  age_structure = internal_age_structure,
  country = country,
  years = simulation_years,
  dataset_name = "popprojAge1dt"
)

wpp_migration_1year <- try_build_wpp_migration_schedule(
  data = migration_data_1year_result$data,
  age_structure = internal_age_structure,
  country = country,
  interval_years = simulation_interval_years
)

# Phase 5: single-year fertility and mortality schedule construction --------

fertility_result <- build_fertility_schedule_from_single_year_wpp(
  fertility_weights_data = fertility_weights_data,
  tfr_data = tfr_data,
  country = country,
  schedule_years = schedule_years,
  age_structure = internal_age_structure
)
fertility_schedule <- fertility_result$schedule
maternal_age_groups <- fertility_result$maternal_age_groups

mortality_input <- mortality_data_1year[
  mortality_data_1year$name == country &
    mortality_data_1year$year %in% schedule_years,
  c("year", "age", "mxB"),
  drop = FALSE
]
mortality_input$year <- as.numeric(mortality_input$year)
mortality_input$mxB <- as.numeric(mortality_input$mxB)

mortality_schedule <- mortality_from_wpp(
  data = mortality_input,
  age_structure = internal_age_structure,
  time_col = "year",
  age_col = "age",
  mortality_col = "mxB",
  quantity = "mx"
)

process <- build_demographic_process(
  age_structure = internal_age_structure,
  fertility_schedule = fertility_schedule,
  fertility_exposure_fraction = female_exposure_fraction,
  mortality_schedule = mortality_schedule,
  mode = "closed"
)

migration_process <- NULL
if (isTRUE(wpp_migration_1year$found)) {
  migration_process <- build_demographic_process(
    age_structure = internal_age_structure,
    fertility_schedule = fertility_schedule,
    fertility_exposure_fraction = female_exposure_fraction,
    mortality_schedule = mortality_schedule,
    migration_schedule = wpp_migration_1year$schedule,
    mode = "migration"
  )
}

# Phase 6: single-year benchmark simulation ---------------------------------

if (!isTRUE(wpp_population_1year$found)) {
  stop(
    "The single-year benchmark requires WPP dataset popprojAge1dt. ",
    wpp_population_1year$reason,
    call. = FALSE
  )
}

primary_initial_population <- demography_population_vector(
  wpp_population_1year$demography,
  time = min(simulation_years)
)

primary_simulated_internal <- simulate_annual_cohort(
  process = process,
  initial_population = primary_initial_population,
  simulation_years = simulation_years
)

migration_simulated_internal <- NULL
migration_result_label <- "migration-adjusted diagnostic"
migration_result_description <- "post-hoc total adjustment using cumulative WPP net migration"
if (!is.null(migration_process)) {
  migration_simulation_result <- tryCatch(
    simulate_annual_cohort(
      process = migration_process,
      initial_population = primary_initial_population,
      simulation_years = simulation_years
    ),
    error = function(error) error
  )
  if (inherits(migration_simulation_result, "error")) {
    migration_result_description <- paste(
      "post-hoc total adjustment using cumulative WPP net migration; WPP-migration simulation failed:",
      conditionMessage(migration_simulation_result)
    )
  } else {
    migration_simulated_internal <- migration_simulation_result
    migration_result_label <- "WPP-migration simulation"
    migration_result_description <- "simulate_demography() with WPP age-specific migration_count"
  }
}

primary_wpp_table <- wpp_population_1year$table
primary_label <- "WPP popprojAge1dt single-year population"

# Phase 7: single-year diagnostics and optional broad-band summaries ---------

single_year_total_comparison <- compare_total_population(
  simulated = primary_simulated_internal,
  wpp = primary_wpp_table
)

migration_total_comparison <- migration_adjusted_total_comparison(
  closed_simulated = primary_simulated_internal,
  wpp = primary_wpp_table,
  migration_schedule = wpp_migration_1year$schedule,
  migration_simulated = migration_simulated_internal,
  simulation_years = simulation_years
)
headline_total_table <- headline_total_comparison(
  migration_total_comparison,
  headline_years
)
headline_metrics <- headline_improvement_metrics(
  migration_total_comparison,
  final_year = max(simulation_years)
)

single_year_signed_residuals <- signed_age_residual_diagnostics(
  simulated = primary_simulated_internal,
  wpp = primary_wpp_table,
  years = diagnostic_years,
  tiny_denominator = tiny_denominator
)
single_year_signed_residual_summaries <- signed_age_residual_summaries(
  simulated = primary_simulated_internal,
  wpp = primary_wpp_table,
  years = diagnostic_years,
  tiny_denominator = tiny_denominator
)
single_year_broad_age_band_residuals <- broad_age_band_residuals(
  single_year_signed_residuals,
  internal_age_band_lookup
)
side_by_side_age_diagnostics <- side_by_side_age_distribution_diagnostics(
  closed_simulated = primary_simulated_internal,
  migration_simulated = migration_simulated_internal,
  wpp = primary_wpp_table,
  years = diagnostic_years,
  tiny_denominator = tiny_denominator
)
broad_age_band_comparison <- broad_age_band_population_comparison(
  closed_simulated = primary_simulated_internal,
  migration_simulated = migration_simulated_internal,
  wpp = primary_wpp_table,
  migration_schedule = wpp_migration_1year$schedule,
  years = diagnostic_years,
  interval_years = simulation_interval_years,
  age_groups = internal_age_structure$age_groups,
  lookup = internal_age_band_lookup
)
broad_age_band_2050 <- broad_age_band_2050_summary(
  broad_age_band_comparison,
  year = max(simulation_years)
)
migration_band_diagnostics <- cumulative_migration_band_diagnostics(
  migration_schedule = wpp_migration_1year$schedule,
  closed_simulated = primary_simulated_internal,
  migration_simulated = migration_simulated_internal,
  wpp = primary_wpp_table,
  years = diagnostic_years,
  interval_years = simulation_interval_years,
  age_groups = internal_age_structure$age_groups,
  lookup = internal_age_band_lookup,
  tiny_denominator = tiny_denominator
)
most_negative_migration_bands <- do.call(
  rbind,
  lapply(
    split(migration_band_diagnostics, migration_band_diagnostics$time),
    function(year_rows) {
      head(
        year_rows[order(year_rows$cumulative_net_migration), c(
          "time",
          "age_band",
          "cumulative_net_migration",
          "closed_signed_error",
          "migration_adjusted_signed_error",
          "error_reduction",
          "migration_explains_positive_residual"
        )],
        3
      )
    }
  )
)
names(most_negative_migration_bands)[
  names(most_negative_migration_bands) == "migration_adjusted_signed_error"
] <- "migration_signed_error"
single_year_age_diagnostics <- age_distribution_diagnostics(
  simulated = primary_simulated_internal,
  wpp = primary_wpp_table,
  years = diagnostic_years,
  tiny_denominator = tiny_denominator
)
single_year_top_errors_2050 <- top_age_errors(
  simulated = primary_simulated_internal,
  wpp = primary_wpp_table,
  year = 2050,
  tiny_denominator = tiny_denominator,
  n = 5
)
maximum_total_relative_difference <- max(
  single_year_total_comparison$relative_difference,
  na.rm = TRUE
)

initial_simulated_total <- sum(
  primary_simulated_internal$population[
    primary_simulated_internal$time == min(simulation_years)
  ]
)
initial_source_total <- sum(primary_initial_population)
initial_check <- isTRUE(all.equal(initial_simulated_total, initial_source_total, tolerance = 1e-10))
age_group_check <- identical(
  unique(primary_simulated_internal$age_group[
    primary_simulated_internal$time == min(simulation_years)
  ]),
  internal_age_structure$age_groups
)
finite_nonnegative_check <- all(is.finite(primary_simulated_internal$population)) &&
  all(primary_simulated_internal$population >= 0)
min_population <- min(primary_simulated_internal$population)

annual_balance <- annual_balance_table(
  simulated = primary_simulated_internal,
  process = process,
  fertility_schedule = fertility_schedule,
  mortality_schedule = mortality_schedule,
  simulation_years = simulation_years,
  age_groups = internal_age_structure$age_groups,
  wpp = primary_wpp_table
)

wpp_interval_columns <- find_wpp_interval_columns(population_data_1year, central_col = "pop")
wpp_interval_diagnostics_available <- !is.null(wpp_interval_columns)

if (wpp_interval_diagnostics_available) {
  wpp_interval_table <- build_wpp_interval_table(
    data = population_data_1year,
    age_structure = internal_age_structure,
    country = country,
    years = simulation_years,
    lower_col = wpp_interval_columns$lower,
    upper_col = wpp_interval_columns$upper
  )
  total_interval_diagnostics <- total_interval_coverage(
    simulated = primary_simulated_internal,
    interval_table = wpp_interval_table
  )
  side_by_side_total_interval_diagnostics <- side_by_side_total_interval_coverage(
    closed_simulated = primary_simulated_internal,
    migration_simulated = migration_simulated_internal,
    interval_table = wpp_interval_table
  )
  age_interval_diagnostics <- age_interval_coverage(
    simulated = primary_simulated_internal,
    interval_table = wpp_interval_table,
    years = diagnostic_years
  )
  side_by_side_age_interval_diagnostics <- side_by_side_age_interval_coverage(
    closed_simulated = primary_simulated_internal,
    migration_simulated = migration_simulated_internal,
    interval_table = wpp_interval_table,
    years = diagnostic_years
  )
  broad_age_band_interval_diagnostics <- broad_age_band_interval_coverage(
    simulated = primary_simulated_internal,
    interval_table = wpp_interval_table,
    lookup = internal_age_band_lookup,
    years = diagnostic_years
  )
  side_by_side_broad_age_band_interval_diagnostics <-
    side_by_side_broad_age_band_interval_coverage(
      closed_simulated = primary_simulated_internal,
      migration_simulated = migration_simulated_internal,
      interval_table = wpp_interval_table,
      lookup = internal_age_band_lookup,
      years = diagnostic_years
    )
} else {
  wpp_interval_table <- NULL
  total_interval_diagnostics <- NULL
  side_by_side_total_interval_diagnostics <- NULL
  age_interval_diagnostics <- NULL
  side_by_side_age_interval_diagnostics <- NULL
  broad_age_band_interval_diagnostics <- NULL
  side_by_side_broad_age_band_interval_diagnostics <- NULL
}

# Phase 8: reporting --------------------------------------------------------

section("Single-year annual-cohort WPP benchmark")
cat("Country:", country, "\n")
cat("Simulation years:", paste(range(simulation_years), collapse = "-"), "\n")
cat("Internal age grid: WPP single-year cohorts, 0 through 100+\n")
cat("Primary initial population source:", primary_label, "\n")
cat("Benchmark comparison source: WPP popprojAge1dt single-year populations\n")
cat("Fertility source: WPP percentASFR1dt single-year percent-ASFR with tfrproj5dt TFR\n")
cat("Fertile age groups inferred from percentASFR1dt:", paste(maternal_age_groups, collapse = ", "), "\n")
cat("Mortality source: WPP mx1dt single-year mxB\n")
cat("Schedule interpolation: linear between WPP schedule years\n")
cat("Fertility exposure fraction:", process$fertility_exposure_fraction, "\n")
cat("Closed baseline: closed population: births and deaths only\n")
if (isTRUE(wpp_migration_1year$found)) {
  cat("WPP migration source: migprojAge1dt age-specific mig counts, interval years ",
      paste(range(simulation_interval_years), collapse = "-"), "\n", sep = "")
  cat("Migration result type:", migration_result_label, "\n")
  cat("Migration result details:", migration_result_description, "\n")
} else {
  cat("WPP migration source: unavailable; migration diagnostics use the closed benchmark only\n")
  cat("WPP migration unavailable reason:", wpp_migration_1year$reason, "\n")
}

section("Core checks")
cat("Initial population total preserved:", initial_check, "\n")
cat("Internal age group order check:", age_group_check, "\n")
cat("Finite nonnegative population check:", finite_nonnegative_check, "\n")
cat("Minimum simulated internal-age population:", min_population, "\n")
cat("Verbose output enabled:", verbose_output, "\n")
cat("\nAnnual component-balance diagnostics:\n")
print(
  annual_balance[, c(
    "year",
    "starting_simulated_total",
    "simulated_births",
    "simulated_deaths",
    "simulated_net_change",
    "wpp_total_net_change",
    "simulated_minus_wpp_net_change"
  )],
  row.names = FALSE
)

section("Closed vs WPP-migration simulation: total population")
cat(
  "Closed baseline is births and deaths only. Migration columns report the ",
  migration_result_label,
  ". Relative errors are signed error / WPP total.\n",
  sep = ""
)
print(headline_total_table, row.names = FALSE)

section("Headline migration improvement")
print(headline_metrics, row.names = FALSE)
cat(
  "\nHeadline interpretation: the closed births/deaths-only model overpredicts ",
  "WPP totals. WPP age-specific net migration explains most of the discrepancy: ",
  "the ",
  migration_result_label,
  " reduces the 2050 total error by about ",
  round(
    headline_metrics$value[
      headline_metrics$metric == "percent_reduction_absolute_2050_total_error"
    ],
    1
  ),
  "%. Remaining differences are much smaller, but not zero.\n",
  sep = ""
)

section("Age-structure diagnostics")
cat(
  "Single-year age errors compare each model variant against WPP popprojAge1dt ",
  "for diagnostic years.\n",
  sep = ""
)
if (is.null(migration_simulated_internal)) {
  cat("Age-structure diagnostics for migration are unavailable because the migration result is post-hoc total adjustment only.\n")
}
print(side_by_side_age_diagnostics, row.names = FALSE)

section("Broad age-band residuals")
cat(
  "Signed error is simulated minus WPP. Cumulative migration is summed by band ",
  "over elapsed WPP interval years.\n",
  sep = ""
)
print(broad_age_band_comparison, row.names = FALSE)

cat("\n2050 broad age-band summary:\n")
print(broad_age_band_2050, row.names = FALSE)

section("Age-specific migration summaries")
if (!isTRUE(wpp_migration_1year$found)) {
  cat("Skipped: migprojAge1dt was unavailable or could not be standardised.\n")
} else {
  cat("Bands with the most negative cumulative WPP net migration:\n")
  print(most_negative_migration_bands, row.names = FALSE)
}

if (verbose_output) {
  section("Verbose total population comparison")
  print(migration_total_comparison, row.names = FALSE)

  section("Verbose signed single-year age residuals")
  cat(
    "Signed error is simulated minus WPP. Relative errors use signed_error / WPP ",
    "where WPP population > ",
    tiny_denominator,
    ".\n",
    sep = ""
  )
  print(single_year_signed_residuals, row.names = FALSE)
  cat("\nClosed signed residual summaries:\n")
  print(single_year_signed_residual_summaries, row.names = FALSE)

  section("Verbose largest closed age-specific errors")
  cat("Top 5 ages by absolute error in 2050:\n")
  print(single_year_top_errors_2050$absolute, row.names = FALSE)
  cat("\nTop 5 ages by relative error in 2050, excluding tiny denominators:\n")
  print(single_year_top_errors_2050$relative, row.names = FALSE)
}

section("WPP uncertainty-bound diagnostics")
if (!wpp_interval_diagnostics_available) {
  cat("WPP uncertainty-bound diagnostics skipped: lower/upper projection columns were not found.\n")
} else {
  cat(
    "Detected WPP lower/upper population columns: ",
    wpp_interval_columns$lower,
    " / ",
    wpp_interval_columns$upper,
    "\n",
    sep = ""
  )
  cat("\nTotal-population interval coverage, closed vs migration:\n")
  print(side_by_side_total_interval_diagnostics, row.names = FALSE)
  cat("\nSingle-year age interval coverage, closed vs migration:\n")
  print(side_by_side_age_interval_diagnostics, row.names = FALSE)
  cat(
    "Note: identical counts of ages inside intervals can be correct when the ",
    "migration simulation moves outside-interval ages closer to WPP without ",
    "moving them across the interval boundary; compare the distance columns.\n",
    sep = ""
  )
  cat("\nBroad age-band interval coverage, closed vs migration:\n")
  print(side_by_side_broad_age_band_interval_diagnostics, row.names = FALSE)
}

section("Interpretation")
cat(
  "The closed births/deaths-only model overpredicts WPP total population ",
  "through the projection horizon. WPP age-specific net migration explains ",
  "most of that total-population discrepancy in this benchmark.\n",
  sep = ""
)
cat(
  "Remaining discrepancies are likely due to sex aggregation, simplified ",
  "female fertility exposure, mortality conventions, schedule interpolation, ",
  "and differences from WPP's full cohort-component machinery.\n",
  sep = ""
)
cat(
  "The migration comparison is still a benchmark/decomposition, not proof of ",
  "exact WPP reproduction.\n",
  sep = ""
)

section("Limitations")
print_limitations()

# Phase 9: plotting ---------------------------------------------------------

enable_plotting <- interactive() ||
  identical(tolower(Sys.getenv("AGEPI_WPP_VALIDATION_PLOT")), "true")

if (!enable_plotting) {
  section("Plotting")
  cat(
    "Plotting skipped for this non-interactive run. Set ",
    "AGEPI_WPP_VALIDATION_PLOT=true to refresh Rplots.pdf intentionally.\n",
    sep = ""
  )
} else {
  initial_distribution <- demography_population_table(
    wpp_population_1year$demography,
    time = min(simulation_years)
  )

  final_simulated <- primary_simulated_internal[
    primary_simulated_internal$time == max(simulation_years),
  ]
  final_wpp <- primary_wpp_table[primary_wpp_table$time == max(simulation_years), ]
  plot_total_comparison <- single_year_total_comparison

  plot_device_opened <- FALSE
  if (grDevices::dev.cur() == 1) {
    grDevices::pdf("Rplots.pdf")
    plot_device_opened <- TRUE
  }

  old_par <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(8, 4, 3, 1))

  graphics::barplot(
    initial_distribution$population,
    names.arg = initial_distribution$age_group,
    las = 2,
    cex.names = 0.35,
    col = "grey70",
    border = NA,
    main = paste(country, min(simulation_years), "initial age distribution"),
    ylab = "Population"
  )

  plot(
    plot_total_comparison$time,
    plot_total_comparison$simulated_population,
    type = "o",
    pch = 16,
    col = "steelblue4",
    xlab = "Year",
    ylab = "Total population",
    main = "Annual cohort vs WPP",
    ylim = range(c(
      plot_total_comparison$simulated_population,
      plot_total_comparison$wpp_population
    ), finite = TRUE)
  )
  lines(
    plot_total_comparison$time,
    plot_total_comparison$wpp_population,
    type = "o",
    pch = 1,
    col = "firebrick3"
  )
  legend(
    "topleft",
    legend = c("agepi annual cohort", "WPP projection"),
    col = c("steelblue4", "firebrick3"),
    pch = c(16, 1),
    lty = 1,
    bty = "n",
    cex = 0.75
  )

  final_comparison <- merge(
    final_simulated[, c("age_group", "population")],
    final_wpp[, c("age_group", "population")],
    by = "age_group",
    suffixes = c("_simulated", "_wpp"),
    sort = FALSE
  )

  graphics::barplot(
    t(as.matrix(final_comparison[, c("population_simulated", "population_wpp")])),
    beside = TRUE,
    names.arg = final_comparison$age_group,
    las = 2,
    cex.names = 0.35,
    col = c("steelblue4", "firebrick3"),
    border = NA,
    main = paste(max(simulation_years), "age distribution"),
    ylab = "Population"
  )
  legend(
    "topright",
    legend = c("agepi", "WPP"),
    fill = c("steelblue4", "firebrick3"),
    bty = "n",
    cex = 0.75
  )
  graphics::par(old_par)
  if (plot_device_opened) {
    invisible(grDevices::dev.off())
  }
}
}
