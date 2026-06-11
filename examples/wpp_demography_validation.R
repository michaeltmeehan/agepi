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
# for Kiribati into agepi's existing demography objects, simulates a closed
# population from 2025 to 2050, and reports benchmark sanity-check diagnostics.
#
# The comparison to WPP projections is intentionally qualitative. agepi's
# current demographic process is not a complete WPP projection system: it does
# not reproduce WPP sex-specific exposures, projection assumptions, or migration.

required_functions <- c(
  "wpp_age_structure_5year",
  "population_from_wpp",
  "fertility_from_wpp_percent_asfr",
  "mortality_from_wpp",
  "build_demographic_process",
  "simulate_demography",
  "demography_population_table",
  "demography_population_vector"
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

country <- "Kiribati"
simulation_years <- 2025:2050
projection_years <- seq(2025, 2050, by = 5)
schedule_years <- seq(2023, 2053, by = 5)
age_structure <- wpp_age_structure_5year(max_age = 100)

population_data <- load_wpp_dataset("popprojAge5dt")
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
  age_structure = age_structure,
  time_col = "year",
  age_group_col = "age",
  population_col = "pop",
  location = country,
  location_col = "name"
)

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
  maternal_age_groups = maternal_age_groups,
  tolerance = 1e-4
)

raw_wpp_female_asfr_schedule <- fertility_schedule

# WPP ASFR is defined as births per female person-year. This example simulates
# total population by age rather than sex-stratified female population by age,
# so applying WPP ASFR directly to the total age-group population overstates
# births. A process-level fertility_exposure_fraction of 0.5 approximates
# female exposure under an assumed 1:1 sex ratio. This is an approximate
# correction for this total-population benchmark, not a replacement for a full
# sex-structured demographic model.
female_exposure_fraction <- 0.5

# WPP mx5dt stores "0" and "1" separately, then five-year lower-bound labels.
# For this five-year agepi grid, use "0" as the first 0-4 approximation and
# omit "1"; other lower-bound labels map to their five-year age groups.
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

mortality_schedule <- mortality_from_wpp(
  data = mortality_input,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  mortality_col = "mxB",
  quantity = "mx"
)

raw_exposure_process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = raw_wpp_female_asfr_schedule,
  mortality_schedule = mortality_schedule,
  mode = "closed"
)

process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = fertility_schedule,
  fertility_exposure_fraction = female_exposure_fraction,
  mortality_schedule = mortality_schedule,
  mode = "closed"
)

initial_population <- demography_population_vector(wpp_population, time = min(simulation_years))

raw_exposure_simulated <- simulate_demography(
  process = raw_exposure_process,
  initial_state = initial_population,
  times = simulation_years,
  time_policy = "linear",
  method = "euler"
)

simulated <- simulate_demography(
  process = process,
  initial_state = initial_population,
  times = simulation_years,
  time_policy = "linear",
  method = "euler"
)

simulated_totals <- aggregate(population ~ time, simulated, sum)
names(simulated_totals)[names(simulated_totals) == "population"] <- "simulated_population"

raw_exposure_simulated_totals <- aggregate(population ~ time, raw_exposure_simulated, sum)
names(raw_exposure_simulated_totals)[
  names(raw_exposure_simulated_totals) == "population"
] <- "raw_exposure_simulated_population"

wpp_totals <- aggregate(population ~ time, demography_population_table(wpp_population), sum)
names(wpp_totals)[names(wpp_totals) == "population"] <- "wpp_population"

raw_exposure_total_comparison <- merge(
  raw_exposure_simulated_totals,
  wpp_totals,
  by = "time",
  all = FALSE
)
raw_exposure_total_comparison$relative_difference <- abs(
  raw_exposure_total_comparison$raw_exposure_simulated_population -
    raw_exposure_total_comparison$wpp_population
) / pmax(abs(raw_exposure_total_comparison$wpp_population), .Machine$double.eps)

total_comparison <- merge(simulated_totals, wpp_totals, by = "time", all = FALSE)
total_comparison$absolute_difference <- abs(
  total_comparison$simulated_population - total_comparison$wpp_population
)
total_comparison$relative_difference <- total_comparison$absolute_difference /
  pmax(abs(total_comparison$wpp_population), .Machine$double.eps)

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

population_matrix <- t(vapply(
  simulation_years,
  function(year) {
    year_population <- simulated[simulated$time == year, ]
    stats::setNames(
      year_population$population[
        match(age_structure$age_groups, year_population$age_group)
      ],
      age_structure$age_groups
    )
  },
  numeric(length(age_structure$age_groups))
))
rownames(population_matrix) <- as.character(simulation_years)
colnames(population_matrix) <- age_structure$age_groups

annual_interval_years <- simulation_years[-length(simulation_years)]
annual_births_raw_total_exposure <- vapply(
  annual_interval_years,
  function(year) {
    fertility_rates <- schedule_rates_at(
      raw_wpp_female_asfr_schedule,
      year,
      "fertility_rate",
      age_structure$age_groups,
      fill_value = 0
    )
    sum(fertility_rates * population_matrix[as.character(year), ])
  },
  numeric(1)
)
annual_births_adjusted_exposure <- vapply(
  annual_interval_years,
  function(year) {
    fertility_rates <- schedule_rates_at(
      fertility_schedule,
      year,
      "fertility_rate",
      age_structure$age_groups,
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
      age_structure$age_groups,
      fill_value = 0
    )
    sum(mortality_rates * population_matrix[as.character(year), ])
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
  raw_wpp_female_asfr_births_total_exposure = annual_births_raw_total_exposure,
  adjusted_births_approx_female_exposure = annual_births_adjusted_exposure,
  deaths = annual_deaths,
  net_change = annual_next_totals - annual_totals,
  raw_crude_birth_rate_per_1000 = annual_births_raw_total_exposure / annual_totals * 1000,
  adjusted_crude_birth_rate_per_1000 = annual_births_adjusted_exposure / annual_totals * 1000,
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
  raw_wpp_female_asfr_schedule$data,
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
  fertility_age_group = sort(unique(raw_wpp_female_asfr_schedule$data$age_group)),
  in_age_structure = sort(unique(raw_wpp_female_asfr_schedule$data$age_group)) %in% age_structure$age_groups,
  stringsAsFactors = FALSE
)

initial_simulated_total <- sum(simulated$population[simulated$time == min(simulation_years)])
initial_wpp_total <- sum(initial_population)
initial_check <- isTRUE(all.equal(initial_simulated_total, initial_wpp_total, tolerance = 1e-10))

age_group_check <- identical(
  unique(simulated$age_group[simulated$time == min(simulation_years)]),
  age_structure$age_groups
)

finite_nonnegative_check <- all(is.finite(simulated$population)) && all(simulated$population >= 0)
min_population <- min(simulated$population)
max_rel_diff <- max(total_comparison$relative_difference, na.rm = TRUE)
raw_exposure_max_rel_diff <- max(raw_exposure_total_comparison$relative_difference, na.rm = TRUE)

cat("WPP-connected demography benchmark sanity-check\n")
cat("Country:", country, "\n")
cat("Simulation years:", paste(simulation_years, collapse = ", "), "\n")
cat("Initial population check:", initial_check, "\n")
cat("Age group check:", age_group_check, "\n")
cat("Finite nonnegative population check:", finite_nonnegative_check, "\n")
cat("Minimum simulated population:", min_population, "\n")
cat("WPP projection comparison is diagnostic only; exact agreement is not expected.\n")
cat(
  "Unadjusted maximum relative difference from WPP projection:",
  raw_exposure_max_rel_diff,
  "\n"
)
cat("Adjusted maximum relative difference from WPP projection:", max_rel_diff, "\n\n")

# Differences from WPP projections are expected because this example uses a
# closed agepi population process. WPP projections include assumptions not
# reproduced here, including migration and sex-specific demographic components.
# The mortality mapping is simplified too: WPP age 0 mortality is used as an
# approximation for age group 0-4, and WPP age 1 mortality is omitted.
cat("Total population comparison, in WPP population units:\n")
print(total_comparison)

cat("\nFertility schedule unit diagnostic:\n")
cat("fertility_from_wpp_percent_asfr() returns annual births per female person-year.\n")
cat("Formula used by adapter: TFR * ASFR_fraction / 5-year age-bin width.\n")
cat("Process fertility_exposure_fraction for total-population simulation:", process$fertility_exposure_fraction, "\n")
cat("Fertility schedule rate convention:", fertility_schedule$rate_convention, "\n")
cat("Mortality schedule rate convention:", mortality_schedule$rate_convention, "\n")
print(raw_fertility_rate_summary)

cat("\nFertility age-group mapping diagnostic:\n")
print(fertility_age_mapping, row.names = FALSE)
cat(
  "Non-reproductive age groups with non-zero fertility:",
  paste(
    setdiff(unique(fertility_schedule$data$age_group), maternal_age_groups),
    collapse = ", "
  ),
  "\n"
)

cat("\nAnnual birth/death balance using adjusted simulation trajectory:\n")
cat("Raw births apply WPP female ASFR to total exposure; adjusted births use process fertility_exposure_fraction.\n")
print(annual_balance, row.names = FALSE)

cat("\nImplied annual growth rates between WPP five-year points:\n")
print(implied_growth, row.names = FALSE)

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
  main = "Closed simulation vs WPP projection",
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
  legend = c("agepi closed simulation", "WPP projection"),
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
