if ("package:agepi" %in% search()) {
  # Already loaded by library(agepi) or a development loader.
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

cat("Demographic ageing diagnostics\n")
cat("===============================\n\n")

cohort_age_structure <- function(width) {
  if (width == 1) {
    return(AgeStructure(
      age_groups = c("0", "1", "2+"),
      lower_bounds = c(0, 1, 2),
      upper_bounds = c(0, 1, Inf)
    ))
  }

  AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
}

run_single_cohort_diagnostic <- function(width, initial_population = 1) {
  ages <- cohort_age_structure(width)
  process <- DemographicProcess(age_structure = ages)
  initial_state <- rep(0, ages$n_age_groups)
  initial_state[1] <- initial_population
  times <- 0:width

  euler <- simulate_demography(
    process = process,
    initial_state = initial_state,
    times = times,
    method = "euler"
  )
  euler_final <- euler[euler$time == width, ]
  euler_remaining <- euler_final$population[euler_final$age_group == ages$age_groups[1]]

  desolve_remaining <- NA_real_
  if (requireNamespace("deSolve", quietly = TRUE)) {
    desolve <- simulate_demography(
      process = process,
      initial_state = initial_state,
      times = c(0, width),
      method = "deSolve"
    )
    desolve_final <- desolve[desolve$time == width, ]
    desolve_remaining <- desolve_final$population[
      desolve_final$age_group == ages$age_groups[1]
    ]
  }

  data.frame(
    age_bin_width = width,
    departure_rate = process$ageing_operator$departure_rate[1],
    euler_remaining_fraction = euler_remaining / initial_population,
    euler_expected_fraction = (1 - 1 / width)^width,
    desolve_remaining_fraction = desolve_remaining / initial_population,
    exponential_expected_fraction = exp(-1),
    deterministic_expected_fraction = 0,
    stringsAsFactors = FALSE
  )
}

cohort_results <- do.call(
  rbind,
  lapply(c(1, 5), run_single_cohort_diagnostic)
)

cat("Single-cohort ageing, no births/mortality/migration\n")
print(cohort_results, row.names = FALSE)
cat("\n")

if (!requireNamespace("wpp2024", quietly = TRUE)) {
  cat(
    "Skipping WPP age-distribution comparison: optional package wpp2024 ",
    "is not installed. Install it with pak::pkg_install(\"PPgp/wpp2024\") ",
    "or remotes::install_github(\"PPgp/wpp2024\").\n",
    sep = ""
  )
} else {
  load_wpp_dataset <- function(name) {
    dataset_environment <- new.env(parent = baseenv())
    utils::data(list = name, package = "wpp2024", envir = dataset_environment)
    as.data.frame(get(name, envir = dataset_environment, inherits = FALSE))
  }

  country <- "Kiribati"
  simulation_years <- 2025:2050
  projection_years <- seq(2025, 2050, by = 5)
  comparison_years <- c(2030, 2040, 2050)
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

  process <- build_demographic_process(
    age_structure = age_structure,
    fertility_schedule = fertility_schedule,
    fertility_exposure_fraction = 0.5,
    mortality_schedule = mortality_schedule,
    mode = "closed"
  )

  initial_population <- demography_population_vector(
    wpp_population,
    time = min(simulation_years)
  )

  simulated <- simulate_demography(
    process = process,
    initial_state = initial_population,
    times = simulation_years,
    time_policy = "linear",
    method = "euler"
  )

  wpp_table <- demography_population_table(wpp_population)
  compare_one_year <- function(year) {
    simulated_year <- simulated[simulated$time == year, c("age_group", "population")]
    wpp_year <- wpp_table[wpp_table$time == year, c("age_group", "population")]
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
    comparison$relative_error <- comparison$absolute_error /
      pmax(abs(comparison$population_wpp), .Machine$double.eps)

    denominator_threshold <- max(1e-6 * sum(comparison$population_wpp), 1)
    included <- comparison$population_wpp >= denominator_threshold
    data.frame(
      year = year,
      total_absolute_age_distribution_error = sum(comparison$absolute_error),
      max_relative_age_group_error = max(comparison$relative_error[included]),
      denominator_threshold = denominator_threshold,
      excluded_age_groups = sum(!included),
      stringsAsFactors = FALSE
    )
  }

  wpp_diagnostics <- do.call(rbind, lapply(comparison_years, compare_one_year))

  cat("WPP age-distribution comparison: Kiribati\n")
  cat("Simulation uses the existing WPP example setup: closed population, annual Euler steps, linear schedules, fertility exposure fraction 0.5.\n")
  print(wpp_diagnostics, row.names = FALSE)

  final_2050 <- merge(
    simulated[simulated$time == 2050, c("age_group", "population")],
    wpp_table[wpp_table$time == 2050, c("age_group", "population")],
    by = "age_group",
    suffixes = c("_simulated", "_wpp"),
    sort = FALSE
  )
  final_2050$absolute_error <- abs(
    final_2050$population_simulated - final_2050$population_wpp
  )
  final_2050$relative_error <- final_2050$absolute_error /
    pmax(abs(final_2050$population_wpp), .Machine$double.eps)

  cat("\nAge-specific simulated vs WPP population at 2050\n")
  print(final_2050, row.names = FALSE)
  cat("\n")
}

cat("Interpretation\n")
cat("--------------\n")
cat(
  "Ageing is implemented as compartmental outflow at rate 1 / age-bin width. ",
  "With deSolve this gives exponential residence times; with annual Euler ",
  "steps it gives the Euler approximation to that outflow. There is no ",
  "deterministic cohort-shift pathway in these diagnostics.\n",
  sep = ""
)
