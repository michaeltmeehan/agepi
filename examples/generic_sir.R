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

# Purpose: build the same age-structured SIR model through CompartmentModel()
# and SIRModel(), then compare their deterministic trajectories.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

contact_matrix <- matrix(c(
  4, 2, 1,
  2, 5, 2,
  1, 2, 4
), nrow = age_structure$n_age_groups, byrow = TRUE)

population <- c(1000, 1200, 900)
initial_infections <- c(5, 3, 2)
initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 3),
  value = c(population - initial_infections, initial_infections, rep(0, 3)),
  stringsAsFactors = FALSE
)

beta <- 0.08
gamma <- 0.25
susceptibility <- c(0.8, 1, 1.2)
infectiousness <- c(1.1, 1, 0.9)
times <- seq(0, 2, by = 0.1)

generic_model <- CompartmentModel(
  compartments = c("S", "I", "R"),
  infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
  transitions = data.frame(from = "I", to = "R", rate = gamma, stringsAsFactors = FALSE),
  infectious_compartments = "I"
)

specialised_model <- SIRModel(gamma = gamma)

generic_output <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = generic_model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = beta,
  susceptibility = susceptibility,
  infectiousness = infectiousness,
  method = "euler"
)

specialised_output <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = specialised_model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = beta,
  susceptibility = susceptibility,
  infectiousness = infectiousness,
  method = "euler"
)

cat("Maximum absolute trajectory difference:\n")
print(max(abs(generic_output$value - specialised_output$value)))
print(tail(generic_output, 9))
