pkgload::load_all(".", quiet = TRUE)
source(file.path("development", "benchmarks", "fixtures.R"))

run_rate_microbenchmarks <- function(iterations = 1) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("The bench package is required to run the benchmark suite.", call. = FALSE)
  }

  ages <- benchmark_age_structure()
  contacts <- benchmark_contact_matrix(ages)

  sir_state_df <- benchmark_sir_state(ages)
  seir_state_df <- benchmark_seir_state(ages)
  generic_state_df <- benchmark_generic_state(ages)

  sir_model <- agepi::SIRModel(gamma = 0.12, beta = 0.22)
  seir_model <- agepi::SEIRModel(sigma = 0.18, gamma = 0.12, beta = 0.22)
  generic_model <- benchmark_generic_model(ages)

  sir_context <- agepi:::prepare_transition_rate_context_validated(
    model = sir_model,
    age_structure = ages,
    contact_matrix = contacts,
    include_public_template = FALSE
  )
  seir_context <- agepi:::prepare_transition_rate_context_validated(
    model = seir_model,
    age_structure = ages,
    contact_matrix = contacts,
    include_public_template = FALSE
  )
  generic_context <- agepi:::prepare_transition_rate_context_validated(
    model = generic_model,
    age_structure = ages,
    contact_matrix = contacts,
    include_public_template = FALSE
  )

  sir_state <- agepi:::transition_state_matrix_from_long(sir_state_df, sir_context)
  seir_state <- agepi:::transition_state_matrix_from_long(seir_state_df, seir_context)
  generic_state <- agepi:::transition_state_matrix_from_long(generic_state_df, generic_context)

  bench::mark(
    `specialized SIR transition rates` = agepi:::transition_rates_from_state_structure(
      sir_state,
      sir_context
    ),
    `specialized SEIR transition rates` = agepi:::transition_rates_from_state_structure(
      seir_state,
      seir_context
    ),
    `generic CompartmentModel transition rates` = agepi:::transition_rates_from_state_structure(
      generic_state,
      generic_context
    ),
    `generic infection-rate matrix` = agepi:::generic_infection_rate_matrix_from_state_matrix(
      generic_state,
      generic_context
    ),
    `force_of_infection` = agepi::force_of_infection(
      infectious = c(10, 20, 15, 12, 7, 3, 2, 1, 0),
      population = c(100, 200, 300, 400, 500, 600, 700, 800, 900),
      contact_matrix = contacts,
      beta = 0.22
    ),
    iterations = iterations,
    check = FALSE
  )
}

format_rate_benchmark_summary <- function(results) {
  data.frame(
    benchmark = c(
      "specialized SIR transition rates",
      "specialized SEIR transition rates",
      "generic CompartmentModel transition rates",
      "generic infection-rate matrix",
      "force_of_infection"
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
  results <- run_rate_microbenchmarks()
  print(format_rate_benchmark_summary(results), row.names = FALSE)
}
