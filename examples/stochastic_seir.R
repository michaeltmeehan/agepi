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

# Purpose: run a small stochastic SEIR example and show both the realised
# trajectory and the event log returned by `simulate_stochastic()`.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

population <- c(200, 250)
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

# `return_events = TRUE` adds the realised Gillespie event log alongside the
# time-indexed trajectory. The fixed seed keeps the example reproducible.
stochastic <- simulate_stochastic(
  initial_state = initial_state,
  times = seq(0, 10, by = 1),
  model = SEIRModel(sigma = 0.35, gamma = 0.25),
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.06,
  seed = 456,
  return_events = TRUE
)

print(head(stochastic$trajectory, 12))
print(head(stochastic$events, 10))
