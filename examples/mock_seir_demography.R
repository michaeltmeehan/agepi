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

# Purpose: combine an age-structured SEIR epidemic with demographic turnover.
# Main ingredients: contact mixing, latent and infectious disease stages,
# fertility, mortality, migration, and Euler time stepping.
# Expected output: example trajectory rows plus compartment totals over time.
# All rates and counts are invented, so the example is dependency-free.

# Age groups define the shared indexing for contacts, disease compartments,
# and demographic schedules. The final open age group uses Inf as its upper
# bound to represent everyone aged 10 and older.
age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

# Contact matrices use recipient rows and source columns. These entries say
# how strongly infectious people in each source age group expose susceptible
# people in each recipient age group.
contact_matrix <- matrix(
  c(4, 2, 1,
    2, 5, 2,
    1, 2, 4),
  nrow = age_structure$n_age_groups,
  byrow = TRUE
)

# Fertility rates create new births. In this package convention, births enter
# the youngest susceptible age group rather than an exposed or infectious
# compartment.
fertility <- FertilitySchedule(
  data.frame(
    time = c(0, 1),
    age_group = c("10+", "10+"),
    fertility_rate = c(0.025, 0.024)
  ),
  age_structure
)

# Mortality rates apply across all disease compartments within each age group,
# so deaths remove susceptible, exposed, infectious, and recovered people.
mortality <- MortalitySchedule(
  data.frame(
    time = rep(c(0, 1), each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = 2),
    mortality_rate = c(0.004, 0.002, 0.015, 0.004, 0.002, 0.016)
  ),
  age_structure
)

# Migration counts are net flows by age group over time. Positive values add
# people and negative values remove people under the selected process mode.
migration <- MigrationSchedule(
  data.frame(
    time = rep(c(0, 1), each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = 2),
    migration_count = c(2, -1, 1, 2, -1, 1)
  ),
  age_structure
)

# The demographic process bundles fertility, mortality, migration, and ageing.
# With mode = "migration", the migration schedule is treated as the supplied
# net movement that remains after natural demographic change.
process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = fertility,
  mortality_schedule = mortality,
  migration_schedule = migration,
  mode = "migration"
)

# The initial epidemic state partitions each age group's population into SEIR
# compartments. Susceptibles are the remaining population after subtracting
# the seeded exposed and infectious individuals.
population <- c(500, 600, 900)
initial_exposed <- c(3, 2, 1)
initial_infected <- c(2, 1, 1)

# Long-form state layout: all susceptible rows, then exposed, infectious, and
# recovered rows, each ordered by the same age groups.
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

# SEIRModel adds sigma, the rate exposed people become infectious, and gamma,
# the recovery rate. Euler stepping is explicit and easy to inspect; deSolve
# methods can be used for smoother ODE integration when configured.
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

# Show both the detailed trajectory and an aggregated view. The totals help
# check how infection and demography jointly change compartment sizes.
print(head(simulation, 12))
print(compartment_totals(simulation))
