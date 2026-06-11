source("examples/validation/validation_helpers.R")
load_local_agepi()

if (!requireNamespace("EpiDynamics", quietly = TRUE)) {
  message(
    "Skipping EpiDynamics SIR validation example: optional package ",
    "EpiDynamics is not installed."
  )
} else {

parameters <- c(beta = 1.4247, gamma = 0.14286)
initials <- c(S = 1 - 1e-06, I = 1e-06, R = 1 - (1 - 1e-06 - 1e-06))
times <- 0:70

reference <- EpiDynamics::SIR(pars = parameters, init = initials, time = times)
reference_table <- reference_results(reference)

age_structure <- AgeStructure("all", 0, Inf)
initial_state <- data.frame(
  compartment = c("S", "I", "R"),
  age_group = "all",
  value = unname(initials[c("S", "I", "R")]),
  stringsAsFactors = FALSE
)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = SIRModel(gamma = parameters[["gamma"]]),
  age_structure = age_structure,
  contact_matrix = matrix(1, 1, 1),
  beta = parameters[["beta"]],
  method = "deSolve"
)

agepi_table <- wide_agepi_output(simulation)
summary <- compare_shared_columns(reference_table, agepi_table, c("S", "I", "R"))

cat("EpiDynamics::SIR validation against agepi SIRModel\n")
cat("Interpretation: exact one-age infection-only validation, up to solver tolerance.\n")
print(summary, row.names = FALSE)
}
