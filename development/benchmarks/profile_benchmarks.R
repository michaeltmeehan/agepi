source(file.path("development", "benchmarks", "fixtures.R"))

run_agepi_profiles <- function(
  output_dir = tempdir(),
  interval = 0.001,
  deterministic_repetitions = 20,
  stochastic_repetitions = 1
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ages <- benchmark_age_structure()
  contacts <- benchmark_contact_matrix(ages)
  times <- seq(0, 10, by = 0.5)
  det_method <- benchmark_deterministic_method()

  sir_state <- benchmark_sir_state(ages)
  generic_state <- benchmark_generic_state(ages)
  generic_model <- benchmark_generic_model(ages)

  deterministic_file <- file.path(output_dir, "agepi-deterministic.rprof")
  stochastic_file <- file.path(output_dir, "agepi-stochastic.rprof")

  utils::Rprof(deterministic_file, interval = interval, line.profiling = FALSE)
  on.exit(utils::Rprof(NULL), add = TRUE)
  for (i in seq_len(deterministic_repetitions)) {
    agepi::simulate_deterministic(
      initial_state = generic_state,
      times = times,
      model = generic_model,
      age_structure = ages,
      contact_matrix = contacts,
      method = det_method
    )
  }
  utils::Rprof(NULL)
  deterministic_summary <- summaryRprof(deterministic_file, lines = "hide")

  utils::Rprof(stochastic_file, interval = interval, line.profiling = FALSE)
  for (i in seq_len(stochastic_repetitions)) {
    agepi::simulate_stochastic(
      initial_state = sir_state,
      times = times,
      model = agepi::SIRModel(gamma = 0.12, beta = 0.22),
      age_structure = ages,
      contact_matrix = contacts,
      seed = 20260818 + i,
      return_events = FALSE
    )
  }
  utils::Rprof(NULL)
  stochastic_summary <- summaryRprof(stochastic_file, lines = "hide")

  list(
    deterministic = deterministic_summary,
    stochastic = stochastic_summary,
    files = c(deterministic = deterministic_file, stochastic = stochastic_file)
  )
}

print_profile_summary <- function(summary, title, top_n = 12) {
  cat("\n", title, "\n", sep = "")

  if (!is.null(summary$by.total) && nrow(summary$by.total) > 0) {
    cat("\nBy total time:\n")
    print(utils::head(summary$by.total[order(-summary$by.total$total.time), ], top_n))
  }

  if (!is.null(summary$by.self) && nrow(summary$by.self) > 0) {
    cat("\nBy self time:\n")
    print(utils::head(summary$by.self[order(-summary$by.self$self.time), ], top_n))
  }

  if (!is.null(summary$by.line) && nrow(summary$by.line) > 0) {
    cat("\nBy line:\n")
    print(utils::head(summary$by.line[order(-summary$by.line$self.time), ], top_n))
  }
}

if (sys.nframe() == 0L) {
  profiles <- run_agepi_profiles()
  print_profile_summary(profiles$deterministic, "Deterministic profile")
  print_profile_summary(profiles$stochastic, "Stochastic profile")
  cat("\nProfile files written to:\n")
  cat("  ", profiles$files[["deterministic"]], "\n", sep = "")
  cat("  ", profiles$files[["stochastic"]], "\n", sep = "")
}
