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

# Purpose: demonstrate a toy age-structured TB-style model with demography,
# age-assortative mixing, age-specific susceptibility/infectiousness, and
# cumulative epidemiological flows.
# These parameters are illustrative only. They are not calibrated estimates,
# do not use Kiribati data, and should not be used for policy analysis.

age_structure <- AgeStructure(
  age_groups = c("0-14", "15-29", "30-44", "45-64", "65+"),
  lower_bounds = c(0, 15, 30, 45, 65),
  upper_bounds = c(14, 29, 44, 64, Inf)
)

contact_matrix <- matrix(
  c(6.0, 2.0, 1.0, 0.6, 0.3,
    2.0, 7.0, 2.5, 1.2, 0.5,
    1.0, 2.5, 6.0, 2.0, 0.8,
    0.6, 1.2, 2.0, 4.5, 1.5,
    0.3, 0.5, 0.8, 1.5, 3.0),
  nrow = age_structure$n_age_groups,
  byrow = TRUE
)

susceptibility <- c(0.9, 1.0, 1.0, 1.1, 1.2)
infectiousness <- c(0.25, 1.0, 1.0, 1.0, 0.9)

# Annual demographic schedules. This is a closed population: births, ageing,
# and background mortality are included, but migration is omitted.
fertility <- FertilitySchedule(
  data.frame(
    time = rep(c(0, 10), each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = 2),
    fertility_rate = c(
      0.000, 0.035, 0.020, 0.001, 0.000,
      0.000, 0.033, 0.019, 0.001, 0.000
    ),
    stringsAsFactors = FALSE
  ),
  age_structure
)

mortality <- MortalitySchedule(
  data.frame(
    time = rep(c(0, 10), each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = 2),
    mortality_rate = c(
      0.004, 0.002, 0.003, 0.007, 0.040,
      0.004, 0.002, 0.003, 0.007, 0.042
    ),
    stringsAsFactors = FALSE
  ),
  age_structure
)

demographic_process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = fertility,
  mortality_schedule = mortality,
  mode = "closed"
)

# Toy TB-style transition rates in per-year units.
recent_to_remote <- 1 / 2
fast_progression <- c("0-14" = 0.06, "15-29" = 0.04, "30-44" = 0.035, "45-64" = 0.045, "65+" = 0.06)
reactivation <- c("0-14" = 0.001, "15-29" = 0.002, "30-44" = 0.003, "45-64" = 0.005, "65+" = 0.008)
treatment_initiation <- 1.5
treatment_completion <- 2.0
relapse <- c("0-14" = 0.002, "15-29" = 0.003, "30-44" = 0.004, "45-64" = 0.005, "65+" = 0.006)

transitions <- data.frame(
  from = c("Lr", "Lr", "Ld", "I", "T", "R"),
  to = c("Ld", "I", "I", "T", "R", "I"),
  stringsAsFactors = FALSE
)
transitions$rate <- I(list(
  recent_to_remote,
  fast_progression,
  reactivation,
  treatment_initiation,
  treatment_completion,
  relapse
))

tb_model <- CompartmentModel(
  compartments = c("S", "Lr", "Ld", "I", "T", "R"),
  infection_transitions = data.frame(from = "S", to = "Lr", stringsAsFactors = FALSE),
  transitions = transitions,
  infectious_compartments = "I",
  birth_compartment = "S",
  migration_compartment = "S"
)

population <- c(3500, 3000, 2600, 2300, 1200)
initial_l_recent <- c(6, 10, 8, 6, 3)
initial_l_remote <- c(40, 180, 260, 320, 210)
initial_active <- c(1, 4, 5, 6, 4)
initial_treatment <- c(0, 2, 2, 3, 2)
initial_recovered <- c(3, 20, 35, 45, 25)

initial_state <- data.frame(
  compartment = rep(c("S", "Lr", "Ld", "I", "T", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 6),
  value = c(
    population - initial_l_recent - initial_l_remote -
      initial_active - initial_treatment - initial_recovered,
    initial_l_recent,
    initial_l_remote,
    initial_active,
    initial_treatment,
    initial_recovered
  ),
  stringsAsFactors = FALSE
)

times <- seq(0, 10, by = 0.25)

tb_output <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = tb_model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.06,
  susceptibility = susceptibility,
  infectiousness = infectiousness,
  demographic_process = demographic_process,
  time_policy = "linear",
  method = "euler",
  cumulative_flows = list(
    infections = list(from = "S", to = "Lr"),
    disease_onset = list(from = c("Lr", "Ld"), to = c("I", "I")),
    treatment_initiation = list(from = "I", to = "T"),
    treatment_completion = list(from = "T", to = "R"),
    relapse = list(from = "R", to = "I")
  )
)

trajectory <- tb_output$trajectory
cumulative_flows <- tb_output$cumulative
compartment_summary <- compartment_totals(trajectory)
age_group_summary <- age_group_totals(trajectory)
population_summary <- total_population(trajectory)

disease_onset_total <- aggregate(
  value ~ time,
  cumulative_flows[cumulative_flows$cumulative_name == "disease_onset", ],
  sum
)
names(disease_onset_total)[names(disease_onset_total) == "value"] <- "cumulative_disease_onset"

print(head(trajectory, 12))
print(head(cumulative_flows, 12))
print(tail(compartment_summary, 12))
print(tail(age_group_summary, age_structure$n_age_groups))
print(tail(population_summary, 5))
print(tail(disease_onset_total, 5))
