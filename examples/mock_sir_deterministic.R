if (requireNamespace("agepi", quietly = TRUE)) {
  library(agepi)
} else {
  pkgload::load_all(".")
}

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10-14"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, 14)
)

population <- c(1000, 1200, 900)
initial_infections <- c(5, 3, 2)

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

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 1, by = 0.1),
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08
)

print(head(simulation, 12))
