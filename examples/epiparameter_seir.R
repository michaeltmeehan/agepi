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

if (!requireNamespace("epiparameter", quietly = TRUE)) {
  message("Skipping epiparameter SEIR example: epiparameter is not installed.")
} else {

  # Purpose: use an epiparameter incubation-period object to parameterise the
  # E -> I transition in an age-structured SEIR model.
  #
  # rate_from_epiparameter() converts the mean delay to 1 / mean_delay. This is
  # an exponential waiting-time approximation and does not preserve the full
  # incubation-period distribution.

  age_structure <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  incubation <- epiparameter::epiparameter_db(
    disease = "COVID-19",
    epi_name = "incubation period",
    single_epiparameter = TRUE
  )

  model <- SEIRModel(
    sigma = rate_from_epiparameter(incubation),
    gamma = 0.25
  )

  population <- c(1000, 1200)
  initial_exposed <- c(2, 1)
  initial_infected <- c(3, 2)
  initial_state <- data.frame(
    compartment = rep(c("S", "E", "I", "R"), each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = 4),
    value = c(
      population - initial_exposed - initial_infected,
      initial_exposed,
      initial_infected,
      c(0, 0)
    ),
    stringsAsFactors = FALSE
  )

  contact_matrix <- matrix(c(
    4, 2,
    2, 5
  ), nrow = age_structure$n_age_groups, byrow = TRUE)

  simulation <- simulate_deterministic(
    initial_state = initial_state,
    times = seq(0, 5, by = 1),
    model = model,
    age_structure = age_structure,
    contact_matrix = contact_matrix,
    beta = 0.04,
    method = "euler"
  )

  print(model$sigma)
  print(tail(simulation, 8))
}
