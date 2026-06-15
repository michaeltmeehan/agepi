if ("package:agepi" %in% search()) {
  # Already loaded by devtools::load_all() or library(agepi).
} else if (dir.exists("R") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else if (dir.exists("R")) {
  invisible(lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), source))
} else if (requireNamespace("agepi", quietly = TRUE)) {
  library(agepi)
} else {
  stop(
    "Package agepi is not installed. Run this script from the package root ",
    "or install agepi first.",
    call. = FALSE
  )
}

# Purpose: demonstrate a lightweight WPP-connected demographic-only workflow.
# The example converts WPP-derived population, fertility, and mortality inputs
# for Kiribati into agepi demography objects, simulates a closed population
# from 2025 to 2050 using annual cohort ageing on a 1-year internal grid, and
# reports WPP 5-year-grid benchmark diagnostics.
#
# The comparison to WPP projections is intentionally qualitative. agepi's
# standalone process is not a complete WPP projection system: it does not
# reproduce WPP sex-specific exposures, migration, or full projection machinery.
# Initial 5-year WPP population counts are split uniformly to 1-year ages as an
# approximation because only 5-year population counts are used here.

required_functions <- c(
  "wpp_age_structure_5year",
  "wpp_age_structure_1year",
  "population_from_wpp",
  "fertility_from_wpp_percent_asfr",
  "mortality_from_wpp",
  "build_demographic_process",
  "simulate_demography",
  "demography_population_table",
  "demography_population_vector",
  "AgeGridMapping",
  "expand_age_counts",
  "expand_fertility_schedule_age_grid",
  "expand_mortality_schedule_age_grid",
  "aggregate_demography_trajectory_age_grid"
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

load_optional_wpp_dataset <- function(name) {
  dataset_environment <- new.env(parent = baseenv())
  tryCatch(
    {
      utils::data(list = name, package = "wpp2024", envir = dataset_environment)
      if (!exists(name, envir = dataset_environment, inherits = FALSE)) {
        return(NULL)
      }
      as.data.frame(get(name, envir = dataset_environment, inherits = FALSE))
    },
    error = function(error) NULL
  )
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

age_structure_diagnostics <- function(simulated, wpp, years, tiny_denominator = 1e-8) {
  diagnostics <- lapply(
    years,
    function(year) {
      simulated_year <- simulated[simulated$time == year, c("age_group", "population")]
      wpp_year <- wpp[wpp$time == year, c("age_group", "population")]
      comparison <- merge(
        simulated_year,
        wpp_year,
        by = "age_group",
        suffixes = c("_simulated", "_wpp"),
        sort = FALSE
      )
      comparison$absolute_error <- abs(
        comparison$population_simulated - comparison$population_wpp
      )
      relative_rows <- abs(comparison$population_wpp) > tiny_denominator
      comparison$relative_error <- NA_real_
      comparison$relative_error[relative_rows] <-
        comparison$absolute_error[relative_rows] /
          abs(comparison$population_wpp[relative_rows])

      data.frame(
        time = year,
        total_absolute_age_distribution_error = sum(comparison$absolute_error),
        maximum_relative_age_group_error = max(comparison$relative_error, na.rm = TRUE),
        age_groups_excluded_from_relative_error = sum(!relative_rows),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, diagnostics)
}

age_error_table <- function(simulated, wpp, year, tiny_denominator = 1e-8) {
  simulated_year <- simulated[simulated$time == year, c("age_group", "population")]
  wpp_year <- wpp[wpp$time == year, c("age_group", "population")]
  comparison <- merge(
    simulated_year,
    wpp_year,
    by = "age_group",
    suffixes = c("_simulated", "_wpp"),
    sort = FALSE
  )
  comparison$absolute_error <- abs(
    comparison$population_simulated - comparison$population_wpp
  )
  comparison$relative_error <- NA_real_
  relative_rows <- abs(comparison$population_wpp) > tiny_denominator
  comparison$relative_error[relative_rows] <-
    comparison$absolute_error[relative_rows] /
      abs(comparison$population_wpp[relative_rows])

  comparison
}

single_year_wpp_diagnostic <- function(
  population_data_1year,
  country,
  simulation_years,
  diagnostic_years,
  age_structure,
  process,
  tiny_denominator = 1e-8
) {
  unavailable <- function(reason) {
    list(
      found = FALSE,
      reason = reason,
      results = NULL
    )
  }

  if (is.null(population_data_1year)) {
    return(unavailable("WPP dataset popprojAge1dt is not available locally."))
  }

  required_columns <- c("name", "year", "age", "pop")
  missing_columns <- setdiff(required_columns, names(population_data_1year))
  if (length(missing_columns) > 0) {
    return(unavailable(paste(
      "WPP dataset popprojAge1dt is missing required column(s):",
      paste(missing_columns, collapse = ", ")
    )))
  }

  population_input_1year <- population_data_1year[
    population_data_1year$name == country &
      population_data_1year$year %in% simulation_years,
    required_columns,
    drop = FALSE
  ]
  population_input_1year$year <- as.numeric(population_input_1year$year)
  population_input_1year$pop <- as.numeric(population_input_1year$pop)

  if (nrow(population_input_1year) == 0) {
    return(unavailable(paste(
      "No WPP single-year population projection rows were found for",
      country,
      "over",
      paste(range(simulation_years), collapse = "-"),
      "."
    )))
  }

  wpp_population_1year <- tryCatch(
    population_from_wpp(
      data = population_input_1year,
      age_structure = age_structure,
      time_col = "year",
      age_group_col = "age",
      population_col = "pop",
      location = country,
      location_col = "name"
    ),
    error = function(error) error
  )
  if (inherits(wpp_population_1year, "error")) {
    return(unavailable(paste(
      "WPP single-year population rows could not be mapped to the internal",
      "1-year age grid:",
      conditionMessage(wpp_population_1year)
    )))
  }

  wpp_population_table_1year <- demography_population_table(wpp_population_1year)
  available_years <- sort(intersect(simulation_years, unique(wpp_population_table_1year$time)))
  if (!min(simulation_years) %in% available_years) {
    return(unavailable(paste(
      "WPP single-year projections were found, but the initial year",
      min(simulation_years),
      "is unavailable."
    )))
  }

  initial_population <- demography_population_vector(
    wpp_population_1year,
    time = min(simulation_years)
  )

  simulated_1year <- simulate_demography(
    process = process,
    initial_state = initial_population,
    times = simulation_years,
    time_policy = "linear",
    method = "euler",
    ageing_policy = "annual_cohort"
  )

  simulated_totals <- aggregate(population ~ time, simulated_1year, sum)
  names(simulated_totals)[names(simulated_totals) == "population"] <- "simulated_population"
  wpp_totals <- aggregate(population ~ time, wpp_population_table_1year, sum)
  names(wpp_totals)[names(wpp_totals) == "population"] <- "wpp_population"
  total_comparison <- merge(simulated_totals, wpp_totals, by = "time", all = FALSE)
  total_comparison <- total_comparison[
    total_comparison$time %in% available_years,
    ,
    drop = FALSE
  ]
  total_comparison$absolute_difference <- abs(
    total_comparison$simulated_population - total_comparison$wpp_population
  )
  total_comparison$relative_difference <- total_comparison$absolute_difference /
    pmax(abs(total_comparison$wpp_population), .Machine$double.eps)

  diagnostic_available_years <- intersect(diagnostic_years, available_years)
  age_diagnostics <- lapply(
    diagnostic_available_years,
    function(year) {
      comparison <- age_error_table(
        simulated_1year,
        wpp_population_table_1year,
        year,
        tiny_denominator = tiny_denominator
      )
      data.frame(
        time = year,
        total_absolute_single_year_age_distribution_error = sum(comparison$absolute_error),
        maximum_relative_single_year_age_error = max(comparison$relative_error, na.rm = TRUE),
        single_year_ages_excluded_from_relative_error = sum(is.na(comparison$relative_error)),
        stringsAsFactors = FALSE
      )
    }
  )
  age_diagnostics <- if (length(age_diagnostics) > 0) {
    do.call(rbind, age_diagnostics)
  } else {
    data.frame(
      time = numeric(0),
      total_absolute_single_year_age_distribution_error = numeric(0),
      maximum_relative_single_year_age_error = numeric(0),
      single_year_ages_excluded_from_relative_error = integer(0)
    )
  }

  top_columns <- c(
    "age_group",
    "population_simulated",
    "population_wpp",
    "absolute_error",
    "relative_error"
  )
  if (2050 %in% available_years) {
    final_comparison <- age_error_table(
      simulated_1year,
      wpp_population_table_1year,
      2050,
      tiny_denominator = tiny_denominator
    )
    top_absolute <- head(
      final_comparison[order(-final_comparison$absolute_error), top_columns],
      5
    )
    top_relative <- final_comparison[!is.na(final_comparison$relative_error), , drop = FALSE]
    top_relative <- head(
      top_relative[order(-top_relative$relative_error), top_columns],
      5
    )
  } else {
    top_absolute <- data.frame(matrix(ncol = length(top_columns), nrow = 0))
    names(top_absolute) <- top_columns
    top_relative <- top_absolute
  }

  list(
    found = TRUE,
    reason = NULL,
    results = list(
      available_years = available_years,
      total_comparison = total_comparison,
      maximum_total_relative_difference = max(total_comparison$relative_difference, na.rm = TRUE),
      age_diagnostics = age_diagnostics,
      top_absolute_2050 = top_absolute,
      top_relative_2050 = top_relative
    )
  )
}

country <- "Kiribati"
simulation_years <- 2025:2050
projection_years <- seq(2025, 2050, by = 5)
schedule_years <- seq(2023, 2053, by = 5)
diagnostic_years <- c(2030, 2040, 2050)
reporting_age_structure <- wpp_age_structure_5year(max_age = 100)
internal_age_structure <- wpp_age_structure_1year(max_age = 100)

reporting_to_internal_mapping <- AgeGridMapping(
  from_age_groups = reporting_age_structure,
  to_age_groups = internal_age_structure
)
internal_to_reporting_mapping <- AgeGridMapping(
  from_age_groups = internal_age_structure,
  to_age_groups = reporting_age_structure
)

population_data <- load_wpp_dataset("popprojAge5dt")
population_data_1year <- load_optional_wpp_dataset("popprojAge1dt")
fertility_weights_data <- load_wpp_dataset("percentASFR5dt")
tfr_data <- load_wpp_dataset("tfrproj5dt")
mortality_data <- load_wpp_dataset("mx5dt")

population_input <- population_data[
  population_data$name == country &
    population_data$year %in% projection_years,
  c("name", "year", "age", "pop"),
  drop = FALSE
]
population_input$year <- as.numeric(population_input$year)
population_input$pop <- as.numeric(population_input$pop)

if (nrow(population_input) == 0) {
  stop("No WPP population projection rows were found for ", country, ".", call. = FALSE)
}

wpp_population <- population_from_wpp(
  data = population_input,
  age_structure = reporting_age_structure,
  time_col = "year",
  age_group_col = "age",
  population_col = "pop",
  location = country,
  location_col = "name"
)
wpp_population_table <- demography_population_table(wpp_population)

maternal_age_groups <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")

fertility_weights_input <- fertility_weights_data[
  fertility_weights_data$name == country &
    fertility_weights_data$year %in% schedule_years &
    fertility_weights_data$age %in% maternal_age_groups,
  c("year", "age", "pasfr"),
  drop = FALSE
]
fertility_weights_input$year <- as.numeric(fertility_weights_input$year)
fertility_weights_input$pasfr <- as.numeric(fertility_weights_input$pasfr)

tfr_input <- tfr_data[
  tfr_data$name == country &
    tfr_data$year %in% schedule_years,
  c("year", "tfr"),
  drop = FALSE
]
tfr_input$year <- as.numeric(tfr_input$year)
tfr_input$tfr <- as.numeric(tfr_input$tfr)

fertility_schedule_5year <- fertility_from_wpp_percent_asfr(
  data = fertility_weights_input,
  age_structure = reporting_age_structure,
  time_col = "year",
  age_col = "age",
  weight_col = "pasfr",
  tfr_data = tfr_input,
  tfr_time_col = "year",
  tfr_col = "tfr",
  weight_type = "percent",
  maternal_age_groups = maternal_age_groups,
  tolerance = 1e-4
)

# WPP ASFR is defined as births per female person-year. This example simulates
# total population by age rather than sex-stratified female population by age,
# so fertility_exposure_fraction = 0.5 approximates female exposure under an
# assumed 1:1 sex ratio. This is a benchmark approximation, not a substitute
# for a full sex-structured demographic model.
female_exposure_fraction <- 0.5

# WPP mx5dt stores "0" and "1" separately, then five-year lower-bound labels.
# For this WPP 5-year reporting grid, use "0" as the first 0-4 approximation and
# omit "1"; other lower-bound labels map to their five-year age groups. The
# resulting 5-year hazards are copied to nested 1-year internal ages below.
mortality_age_labels <- c(as.character(seq(0, 95, by = 5)), "100+")
mortality_input <- mortality_data[
  mortality_data$name == country &
    mortality_data$year %in% schedule_years &
    mortality_data$age %in% mortality_age_labels,
  c("year", "age", "mxB"),
  drop = FALSE
]
mortality_input$year <- as.numeric(mortality_input$year)
mortality_input$mxB <- as.numeric(mortality_input$mxB)

mortality_schedule_5year <- mortality_from_wpp(
  data = mortality_input,
  age_structure = reporting_age_structure,
  time_col = "year",
  age_col = "age",
  mortality_col = "mxB",
  quantity = "mx"
)

# Expand WPP-derived schedules to the annual-cohort internal grid. Fertility
# rates and mortality hazards are copied to each nested 1-year age group.
fertility_schedule <- expand_fertility_schedule_age_grid(
  fertility_schedule_5year,
  reporting_to_internal_mapping
)
mortality_schedule <- expand_mortality_schedule_age_grid(
  mortality_schedule_5year,
  reporting_to_internal_mapping
)

process <- build_demographic_process(
  age_structure = internal_age_structure,
  fertility_schedule = fertility_schedule,
  fertility_exposure_fraction = female_exposure_fraction,
  mortality_schedule = mortality_schedule,
  mode = "closed"
)

initial_population_5year <- data.frame(
  age_group = reporting_age_structure$age_groups,
  population = demography_population_vector(wpp_population, time = min(simulation_years)),
  stringsAsFactors = FALSE
)

# Uniformly splitting 5-year WPP counts to 1-year internal ages is approximate:
# it preserves the source age-group total but cannot recover within-bin age
# shape from the coarser WPP population counts used by this example.
initial_population_1year_table <- expand_age_counts(
  initial_population_5year,
  mapping = reporting_to_internal_mapping,
  value_col = "population"
)
initial_population <- stats::setNames(
  initial_population_1year_table$population[
    match(internal_age_structure$age_groups, initial_population_1year_table$age_group)
  ],
  internal_age_structure$age_groups
)

simulated_internal <- simulate_demography(
  process = process,
  initial_state = initial_population,
  times = simulation_years,
  time_policy = "linear",
  method = "euler",
  ageing_policy = "annual_cohort"
)

simulated <- aggregate_demography_trajectory_age_grid(
  simulated_internal,
  internal_to_reporting_mapping
)

simulated_totals <- aggregate(population ~ time, simulated, sum)
names(simulated_totals)[names(simulated_totals) == "population"] <- "simulated_population"

wpp_totals <- aggregate(population ~ time, wpp_population_table, sum)
names(wpp_totals)[names(wpp_totals) == "population"] <- "wpp_population"

total_comparison <- merge(simulated_totals, wpp_totals, by = "time", all = FALSE)
total_comparison$absolute_difference <- abs(
  total_comparison$simulated_population - total_comparison$wpp_population
)
total_comparison$relative_difference <- total_comparison$absolute_difference /
  pmax(abs(total_comparison$wpp_population), .Machine$double.eps)

population_matrix <- t(vapply(
  simulation_years,
  function(year) {
    year_population <- simulated_internal[simulated_internal$time == year, ]
    stats::setNames(
      year_population$population[
        match(internal_age_structure$age_groups, year_population$age_group)
      ],
      internal_age_structure$age_groups
    )
  },
  numeric(length(internal_age_structure$age_groups))
))
rownames(population_matrix) <- as.character(simulation_years)
colnames(population_matrix) <- internal_age_structure$age_groups

annual_interval_years <- simulation_years[-length(simulation_years)]
annual_births <- vapply(
  annual_interval_years,
  function(year) {
    fertility_rates <- schedule_rates_at(
      fertility_schedule,
      year,
      "fertility_rate",
      internal_age_structure$age_groups,
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
      internal_age_structure$age_groups,
      fill_value = 0
    )
    sum((1 - exp(-mortality_rates)) * population_matrix[as.character(year), ])
  },
  numeric(1)
)

annual_totals <- simulated_totals$simulated_population[
  match(annual_interval_years, simulated_totals$time)
]
annual_next_totals <- simulated_totals$simulated_population[
  match(annual_interval_years + 1, simulated_totals$time)
]
annual_balance <- data.frame(
  year = annual_interval_years,
  total_population = annual_totals,
  births_approx_female_exposure = annual_births,
  deaths_from_annual_survival = annual_deaths,
  net_change = annual_next_totals - annual_totals,
  crude_birth_rate_per_1000 = annual_births / annual_totals * 1000,
  crude_death_rate_per_1000 = annual_deaths / annual_totals * 1000
)

growth_interval_start <- projection_years[-length(projection_years)]
growth_interval_end <- projection_years[-1]
implied_growth <- data.frame(
  start_year = growth_interval_start,
  end_year = growth_interval_end,
  wpp_annual_growth_rate = (
    total_comparison$wpp_population[match(growth_interval_end, total_comparison$time)] /
      total_comparison$wpp_population[match(growth_interval_start, total_comparison$time)]
  )^(1 / (growth_interval_end - growth_interval_start)) - 1,
  simulated_annual_growth_rate = (
    total_comparison$simulated_population[match(growth_interval_end, total_comparison$time)] /
      total_comparison$simulated_population[match(growth_interval_start, total_comparison$time)]
  )^(1 / (growth_interval_end - growth_interval_start)) - 1
)
implied_growth$growth_rate_difference <- implied_growth$simulated_annual_growth_rate -
  implied_growth$wpp_annual_growth_rate

raw_fertility_rate_summary <- aggregate(
  fertility_rate ~ time,
  fertility_schedule_5year$data,
  function(x) c(
    min = min(x),
    max = max(x),
    sum_5_year_age_rates = sum(x),
    implied_tfr = sum(x * 5)
  )
)
raw_fertility_rate_summary <- do.call(
  data.frame,
  raw_fertility_rate_summary
)
names(raw_fertility_rate_summary) <- c(
  "time",
  "min_rate",
  "max_rate",
  "sum_5_year_age_rates",
  "implied_tfr"
)
raw_fertility_rate_summary$adjusted_implied_tfr_for_total_population <-
  raw_fertility_rate_summary$implied_tfr * female_exposure_fraction

fertility_age_mapping <- data.frame(
  fertility_age_group = sort(unique(fertility_schedule_5year$data$age_group)),
  in_reporting_age_grid = sort(unique(fertility_schedule_5year$data$age_group)) %in%
    reporting_age_structure$age_groups,
  stringsAsFactors = FALSE
)

age_diagnostics <- age_structure_diagnostics(
  simulated = simulated,
  wpp = wpp_population_table,
  years = diagnostic_years
)

previous_exponential_diagnostics <- data.frame(
  time = c(2030, 2040, 2050),
  previous_total_absolute_age_distribution_error = c(5.988668, 9.828091, 10.632588),
  previous_maximum_relative_age_group_error = c(0.1378459, 0.2183497, 0.2374219)
)
age_diagnostics_with_baseline <- merge(
  age_diagnostics,
  previous_exponential_diagnostics,
  by = "time",
  all.x = TRUE,
  sort = FALSE
)

initial_simulated_total <- sum(simulated$population[simulated$time == min(simulation_years)])
initial_wpp_total <- sum(initial_population_5year$population)
initial_check <- isTRUE(all.equal(initial_simulated_total, initial_wpp_total, tolerance = 1e-10))

age_group_check <- identical(
  unique(simulated_internal$age_group[simulated_internal$time == min(simulation_years)]),
  internal_age_structure$age_groups
)

single_year_diagnostic <- single_year_wpp_diagnostic(
  population_data_1year = population_data_1year,
  country = country,
  simulation_years = simulation_years,
  diagnostic_years = diagnostic_years,
  age_structure = internal_age_structure,
  process = process
)

finite_nonnegative_check <- all(is.finite(simulated_internal$population)) &&
  all(simulated_internal$population >= 0)
min_population <- min(simulated_internal$population)
max_rel_diff <- max(total_comparison$relative_difference, na.rm = TRUE)

section <- function(title) {
  cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
}

section("Annual-cohort WPP benchmark")
cat("Country:", country, "\n")
cat("Simulation years:", paste(range(simulation_years), collapse = "-"), "\n")
cat("Internal age grid: 1-year\n")
cat("Reporting grid: WPP 5-year\n")
cat("Population expansion: uniform split from WPP 5-year counts to 1-year ages\n")
cat("Migration: omitted; process is closed-population except for fertility/mortality\n")
cat("Fertility exposure fraction:", process$fertility_exposure_fraction, "\n")
cat("Initial population check:", initial_check, "\n")
cat("Internal age group check:", age_group_check, "\n")
cat("Finite nonnegative population check:", finite_nonnegative_check, "\n")
cat("Minimum simulated internal-age population:", min_population, "\n")
cat("Maximum total-population relative difference:", max_rel_diff, "\n\n")
cat("Total population comparison, in WPP population units:\n")
print(total_comparison)

section("Single-year WPP diagnostic")
cat("WPP single-year population projections found:", single_year_diagnostic$found, "\n")
cat(
  "Purpose: isolate deterministic cohort ageing from 5-year population ",
  "splitting and 1-year-to-5-year aggregation.\n",
  sep = ""
)
cat(
  "Fertility and mortality still use existing WPP schedules expanded to the ",
  "1-year grid where needed; this is a diagnostic, not exact WPP validation.\n",
  sep = ""
)
if (!single_year_diagnostic$found) {
  cat(
    "The 1-year vs 1-year diagnostic cannot be run with the available local WPP data: ",
    single_year_diagnostic$reason,
    "\n",
    sep = ""
  )
} else {
  cat(
    "Available WPP single-year comparison years:",
    paste(single_year_diagnostic$results$available_years, collapse = ", "),
    "\n"
  )
  cat(
    "Maximum total-population relative difference:",
    single_year_diagnostic$results$maximum_total_relative_difference,
    "\n"
  )
  cat(
    "Single-year age-distribution errors at diagnostic years; relative errors exclude WPP denominators <=",
    1e-8,
    "\n"
  )
  # Strong total-population agreement does not guarantee exact age-structure
  # agreement because age-specific errors can offset in the total.
  print(single_year_diagnostic$results$age_diagnostics, row.names = FALSE)
  cat("\nTop 5 single-year ages by absolute error at 2050:\n")
  print(single_year_diagnostic$results$top_absolute_2050, row.names = FALSE)
  cat("\nTop 5 single-year ages by relative error at 2050, excluding tiny denominators:\n")
  # Relative errors at very old ages can be dominated by tiny WPP denominators;
  # denominators below the threshold are excluded from this ranking.
  print(single_year_diagnostic$results$top_relative_2050, row.names = FALSE)
}

section("Five-year reporting-grid comparison")
cat(
  "User-facing benchmark: annual-cohort simulation aggregated back to the WPP ",
  "5-year reporting grid.\n",
  sep = ""
)
cat("Relative age-group errors exclude WPP denominators <=", 1e-8, "\n")
print(age_diagnostics, row.names = FALSE)

section("Limitations")
cat("These are benchmark diagnostics, not exact validation against WPP.\n")
cat(
  "Strong total-population agreement can coexist with age-structure ",
  "differences because age-specific errors can offset in aggregate totals.\n",
  sep = ""
)
cat(
  "Relative errors at very old ages can be large when WPP denominators are ",
  "tiny, even when absolute differences are small.\n",
  sep = ""
)
cat(
  "Remaining age-specific discrepancies are expected: this example omits sex ",
  "structure, migration, and the full WPP projection assumptions, and expands ",
  "5-year fertility/mortality schedules onto a 1-year internal grid.\n",
  sep = ""
)

initial_distribution <- demography_population_table(wpp_population, time = min(simulation_years))
final_simulated <- simulated[simulated$time == max(simulation_years), ]
final_wpp <- demography_population_table(wpp_population, time = max(simulation_years))

old_par <- graphics::par(no.readonly = TRUE)
on.exit(graphics::par(old_par), add = TRUE)
graphics::par(mfrow = c(1, 3), mar = c(8, 4, 3, 1))

graphics::barplot(
  initial_distribution$population,
  names.arg = initial_distribution$age_group,
  las = 2,
  cex.names = 0.6,
  col = "grey70",
  border = NA,
  main = paste(country, min(simulation_years), "WPP age distribution"),
  ylab = "Population"
)

plot(
  total_comparison$time,
  total_comparison$simulated_population,
  type = "o",
  pch = 16,
  col = "steelblue4",
  xlab = "Year",
  ylab = "Total population",
  main = "Annual cohort vs WPP projection",
  ylim = range(c(
    total_comparison$simulated_population,
    total_comparison$wpp_population
  ), finite = TRUE)
)
lines(
  total_comparison$time,
  total_comparison$wpp_population,
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
  cex.names = 0.6,
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
}
