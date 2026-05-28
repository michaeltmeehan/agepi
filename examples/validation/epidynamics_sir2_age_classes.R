source("examples/validation/validation_helpers.R")
load_local_agepi()
library(EpiDynamics)

parameters <- c(
  betaCC = 100,
  betaCA = 10,
  betaAC = 10,
  betaAA = 20,
  gamma = 10,
  lC = 0.0666667,
  muC = 0.0,
  muA = 0.016667
)

initials <- c(
  SC = 0.1,
  IC = 0.0001,
  SA = 0.1,
  IA = 0.0001
)

times <- seq(0, 100, 0.01)

reference <- EpiDynamics::sir2AgeClasses(pars = parameters, init = initials, time = times)
reference_table <- reference_results(reference)
stopifnot(all(c("time", "SC", "IC", "SA", "IA") %in% names(reference_table)))

nC <- parameters[["muA"]] / (parameters[["lC"]] + parameters[["muA"]])
nA <- 1 - nC
nu <- (parameters[["lC"]] + parameters[["muA"]]) * nC

age_structure <- AgeStructure(
  age_groups = c("children", "adults"),
  lower_bounds = c(0, 15),
  upper_bounds = c(14, Inf)
)

mortality <- MortalitySchedule(
  data.frame(
    time = rep(range(times), each = 2),
    age_group = rep(age_structure$age_groups, times = 2),
    mortality_rate = rep(c(parameters[["muC"]], parameters[["muA"]]), times = 2)
  ),
  age_structure
)
fertility <- FertilitySchedule(
  data.frame(
    time = range(times),
    age_group = "adults",
    fertility_rate = nu / nA
  ),
  age_structure
)
process <- build_demographic_process(age_structure, fertility, mortality)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = 2),
  age_group = rep(age_structure$age_groups, times = 3),
  value = c(
    initials[["SC"]], initials[["SA"]],
    initials[["IC"]], initials[["IA"]],
    nC - initials[["SC"]] - initials[["IC"]],
    nA - initials[["SA"]] - initials[["IA"]]
  ),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(
  c(
    parameters[["betaCC"]] * nC, parameters[["betaCA"]] * nA,
    parameters[["betaAC"]] * nC, parameters[["betaAA"]] * nA
  ),
  nrow = 2,
  byrow = TRUE
)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = SIRModel(gamma = parameters[["gamma"]]),
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 1,
  demographic_process = process,
  time_policy = "linear",
  method = "deSolve"
)

name_map <- c(
  "S:children" = "SC",
  "I:children" = "IC",
  "S:adults" = "SA",
  "I:adults" = "IA"
)
agepi_table <- wide_agepi_output(simulation, name_map)
summary <- compare_shared_columns(reference_table, agepi_table, c("SC", "IC", "SA", "IA"))

cat("EpiDynamics::sir2AgeClasses partial reproduction with agepi SIRModel + demography\n")
cat("Reference results extracted from reference$results.\n")
cat("Interpretation: partial age-structured reproduction, not an exact validation.\n")
cat("agepi carries recovered states explicitly and applies demographic coupling to them; EpiDynamics omits recovered states from this model, so recovery exits the reported four-state system.\n")
cat("The contact matrix is column-scaled by equilibrium age totals to align the transmission convention as closely as the current public interface allows.\n")
print(summary, row.names = FALSE)
