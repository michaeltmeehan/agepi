source(file.path("development", "benchmarks", "fixtures.R"))

run_agepi_benchmarks <- function(iterations = 1) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("The bench package is required to run the benchmark suite.", call. = FALSE)
  }

  ages <- benchmark_age_structure()
  contacts <- benchmark_contact_matrix(ages)
  times <- benchmark_times()
  det_method <- benchmark_deterministic_method()

  sir_state <- benchmark_sir_state(ages)
  seir_state <- benchmark_seir_state(ages)
  generic_state <- benchmark_generic_state(ages)
  generic_model <- benchmark_generic_model(ages)
  demography <- benchmark_demographic_process(ages)
  demographic_times <- benchmark_demography_times()

  bench::mark(
    `deterministic SIR` = agepi::simulate_deterministic(
      initial_state = sir_state,
      times = times,
      model = agepi::SIRModel(gamma = 0.12, beta = 0.22),
      age_structure = ages,
      contact_matrix = contacts,
      method = det_method
    ),
    `deterministic SEIR` = agepi::simulate_deterministic(
      initial_state = seir_state,
      times = times,
      model = agepi::SEIRModel(sigma = 0.18, gamma = 0.12, beta = 0.22),
      age_structure = ages,
      contact_matrix = contacts,
      method = det_method
    ),
    `generic CompartmentModel` = agepi::simulate_deterministic(
      initial_state = generic_state,
      times = times,
      model = generic_model,
      age_structure = ages,
      contact_matrix = contacts,
      method = det_method
    ),
    `stochastic Gillespie` = agepi::simulate_stochastic(
      initial_state = sir_state,
      times = times,
      model = agepi::SIRModel(gamma = 0.12, beta = 0.22),
      age_structure = ages,
      contact_matrix = contacts,
      seed = 20260818,
      return_events = FALSE
    ),
    `deterministic SIR with demography` = agepi::simulate_deterministic(
      initial_state = sir_state,
      times = demographic_times,
      model = agepi::SIRModel(gamma = 0.12, beta = 0.22),
      age_structure = ages,
      contact_matrix = contacts,
      demographic_process = demography,
      time_policy = "linear",
      migration_policy = "proportional",
      method = det_method
    ),
    iterations = iterations,
    check = FALSE
  )
}

format_benchmark_summary <- function(results) {
  data.frame(
    benchmark = c(
      "deterministic SIR",
      "deterministic SEIR",
      "generic CompartmentModel",
      "stochastic Gillespie",
      "deterministic SIR with demography"
    ),
    median = format(results$median),
    `itr/sec` = round(as.numeric(results$`itr/sec`), 2),
    `mem alloc` = format(results$mem_alloc),
    `gc/sec` = round(as.numeric(results$`gc/sec`), 2),
    `iterations` = results$n_itr,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

if (sys.nframe() == 0L) {
  results <- run_agepi_benchmarks()
  print(format_benchmark_summary(results), row.names = FALSE)
}
