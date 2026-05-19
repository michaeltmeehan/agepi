if (requireNamespace("agepi", quietly = TRUE)) {
  library(agepi)
} else if (dir.exists("R")) {
  invisible(lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), source))
} else {
  stop(
    "Package agepi is not installed. Run this script from the package root ",
    "or install agepi first.",
    call. = FALSE
  )
}

# Small invented SEIR-demography example with no external data dependencies.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

contact_matrix <- matrix(
  c(4, 2, 1,
    2, 5, 2,
    1, 2, 4),
  nrow = age_structure$n_age_groups,
  byrow = TRUE
)

fertility <- FertilitySchedule(
  data.frame(
    time = c(0, 1),
    age_group = c("10+", "10+"),
    fertility_rate = c(0.025, 0.024)
  ),
  age_structure
)

mortality <- MortalitySchedule(
  data.frame(
    time = rep(c(0, 1), each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = 2),
    mortality_rate = c(0.004, 0.002, 0.015, 0.004, 0.002, 0.016)
  ),
  age_structure
)

migration <- MigrationSchedule(
  data.frame(
    time = rep(c(0, 1), each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = 2),
    migration_count = c(2, -1, 1, 2, -1, 1)
  ),
  age_structure
)

process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = fertility,
  mortality_schedule = mortality,
  migration_schedule = migration,
  mode = "migration"
)

population <- c(500, 600, 900)
initial_exposed <- c(3, 2, 1)
initial_infected <- c(2, 1, 1)

initial_state <- data.frame(
  compartment = rep(c("S", "E", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 4),
  value = c(
    population - initial_exposed - initial_infected,
    initial_exposed,
    initial_infected,
    rep(0, age_structure$n_age_groups)
  ),
  stringsAsFactors = FALSE
)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 1, by = 0.25),
  model = SEIRModel(sigma = 0.4, gamma = 0.25),
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08,
  demographic_process = process,
  time_policy = "linear",
  method = "euler"
)

print(head(simulation, 12))
print(compartment_totals(simulation))
