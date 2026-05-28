source("examples/validation/validation_helpers.R")
load_local_agepi()
library(EpiDynamics)

parameters <- c(mu = 1 / (70 * 365), beta = 520 / 365, sigma = 1 / 14, gamma = 1 / 7)
initials <- c(S = 0.1, E = 1e-04, I = 1e-04, R = 1 - 0.1 - 1e-4 - 1e-4)
times <- 0:(60 * 365)

reference <- EpiDynamics::SEIR(pars = parameters, init = initials, time = times)
reference_table <- reference_results(reference)

age_structure <- AgeStructure("all", 0, Inf)
mortality <- MortalitySchedule(
  data.frame(time = range(times), age_group = "all", mortality_rate = parameters[["mu"]]),
  age_structure
)
fertility <- FertilitySchedule(
  data.frame(time = range(times), age_group = "all", fertility_rate = parameters[["mu"]]),
  age_structure
)
process <- build_demographic_process(age_structure, fertility, mortality)

initial_state <- data.frame(
  compartment = c("S", "E", "I", "R"),
  age_group = "all",
  value = unname(initials[c("S", "E", "I", "R")]),
  stringsAsFactors = FALSE
)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = SEIRModel(sigma = parameters[["sigma"]], gamma = parameters[["gamma"]]),
  age_structure = age_structure,
  contact_matrix = matrix(1, 1, 1),
  beta = parameters[["beta"]],
  demographic_process = process,
  time_policy = "linear",
  method = "deSolve"
)

agepi_table <- wide_agepi_output(simulation)
summary <- compare_shared_columns(reference_table, agepi_table, c("S", "E", "I", "R"))

cat("EpiDynamics::SEIR validation against agepi SEIRModel + demography\n")
cat("Interpretation: exact one-age validation of SEIR with equal births and deaths, up to solver tolerance.\n")
print(summary, row.names = FALSE)
