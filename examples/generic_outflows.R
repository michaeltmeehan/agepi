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

# Purpose: show a generic compartment model with explicit external outflows.
# The same outflow selector works for deterministic cumulative outputs and
# stochastic event-log-derived cumulative counts.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

initial_state <- data.frame(
  compartment = rep(c("S", "I"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 2),
  value = c(0, 0, 12, 8),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(0, nrow = age_structure$n_age_groups, ncol = age_structure$n_age_groups)

model <- CompartmentModel(
  compartments = c("S", "I"),
  outflows = data.frame(from = "I", rate = 0.4, stringsAsFactors = FALSE),
  infectious_compartments = character()
)

times <- seq(0, 4, by = 1)

deterministic_output <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0,
  method = "euler",
  cumulative_flows = list(removals = list(from = "I", to = NA_character_))
)

stochastic_output <- simulate_stochastic(
  initial_state = initial_state,
  times = times,
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0,
  seed = 123,
  return_events = TRUE,
  cumulative_flows = list(removals = list(from = "I", to = NA_character_))
)

cat("Generic outflow model transitions:\n")
print(model$transitions)

cat("\nDeterministic cumulative removals:\n")
print(deterministic_output$cumulative)

cat("\nStochastic event log:\n")
print(stochastic_output$events)

cat("\nStochastic cumulative removals:\n")
print(stochastic_output$cumulative)
