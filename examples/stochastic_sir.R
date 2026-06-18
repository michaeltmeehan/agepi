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

# Purpose: run a small age-structured stochastic SIR simulation with
# Gillespie's direct method.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10-14"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, 14)
)

population <- c(200, 250, 180)
initial_infections <- c(3, 2, 1)
initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 3),
  value = c(population - initial_infections, initial_infections, rep(0, 3)),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(c(
  4, 2, 1,
  2, 5, 2,
  1, 2, 4
), nrow = age_structure$n_age_groups, byrow = TRUE)

model <- SIRModel(gamma = 0.25)
times <- seq(0, 10, by = 1)

stochastic <- simulate_stochastic(
  initial_state = initial_state,
  times = times,
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08,
  seed = 123,
  return_events = TRUE
)

print(head(stochastic$trajectory, 12))
print(tail(stochastic$trajectory, 12))
print(head(stochastic$events, 10))

deterministic <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08
)

stochastic_totals <- aggregate(value ~ time + compartment, stochastic$trajectory, sum)
deterministic_totals <- aggregate(value ~ time + compartment, deterministic, sum)

print(stochastic_totals)
print(deterministic_totals)
