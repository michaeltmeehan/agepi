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

# Purpose: run the smallest age-structured deterministic SIR example.
# Main ingredients: three age groups, an age-specific contact matrix, a
# simple SIR natural-history model, and a short one-time-unit simulation.
# Expected output: the first rows of a long-format trajectory by time,
# compartment, and age group.

# AgeStructure defines the demographic strata used by every model object
# below. Here the population is split into three five-year child age bands.
age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10-14"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, 14)
)

# Initial population sizes and infections are age-specific counts. Everyone
# not initially infectious starts susceptible, and no one starts recovered.
population <- c(1000, 1200, 900)
initial_infections <- c(5, 3, 2)

# Initial states are stored in long form: each row is one
# compartment-age-group cell in the model state vector.
initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 3),
  value = c(population - initial_infections, initial_infections, rep(0, 3)),
  stringsAsFactors = FALSE
)

# Contacts are age structured. Rows are recipient age groups and columns are
# source age groups, so entry [i, j] contributes contacts from age j to age i.
contact_matrix <- matrix(c(
  4, 2, 1,
  2, 5, 2,
  1, 2, 4
), nrow = age_structure$n_age_groups, byrow = TRUE)

# gamma is the recovery rate from I to R. The force of infection is computed
# from beta, the contact matrix, and the infectious counts in each source age.
model <- SIRModel(gamma = 0.25)

# simulate_deterministic integrates the compartment model over the requested
# output times. This example has no demography, so population only changes
# through infection and recovery flows within age groups.
simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 1, by = 0.1),
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08
)

# Print a small slice of the trajectory to show the long output layout.
print(head(simulation, 12))
