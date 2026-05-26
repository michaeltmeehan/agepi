if ("package:agepi" %in% search()) {
  # Already loaded by devtools::load_all() or library(agepi).
} else if (dir.exists("R") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
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

# Purpose: demonstrate a custom MSIR model not available through SIRModel()
# or SEIRModel(). Maternal immunity wanes from M to S at a fixed per-capita
# rate; susceptible people can be infected; infectious people recover.

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

initial_state <- data.frame(
  compartment = rep(c("M", "S", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 4),
  value = c(
    c(100, 20, 0),
    c(890, 1170, 898),
    c(10, 10, 2),
    c(0, 0, 0)
  ),
  stringsAsFactors = FALSE
)

model <- CompartmentModel(
  compartments = c("M", "S", "I", "R"),
  infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
  transitions = data.frame(
    from = c("M", "I"),
    to = c("S", "R"),
    rate = c(0.15, 0.25),
    stringsAsFactors = FALSE
  ),
  infectious_compartments = "I",
  birth_compartment = "M",
  migration_compartment = "S"
)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 2, by = 0.1),
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08,
  method = "euler"
)

totals <- aggregate(value ~ time, simulation, sum)
cat("Maximum absolute change in total living population:\n")
print(max(abs(totals$value - totals$value[1])))
print(tail(simulation, 12))
