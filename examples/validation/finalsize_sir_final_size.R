source("examples/validation/validation_helpers.R")
load_local_agepi()

# Optional validation against finalsize::final_size().
#
# Assumptions:
# - closed population;
# - no births, deaths, migration, waning immunity, or reinfection;
# - agepi uses recipient-row/source-column contact matrices:
#   contact_matrix[recipient_age_group, source_age_group];
# - finalsize documents contact_matrix[i, j] as contacts in group i reported by
#   group j, which is compatible with the orientation used below;
# - finalsize returns the attack rate among initially susceptible people for an
#   epidemic seeded by a vanishingly small number of infections. The agepi
#   simulation therefore uses a very small initial infectious fraction;
# - finalsize takes r0, while agepi takes beta and gamma. With unequal age
#   populations, the comparable agepi next-generation matrix is
#   (beta / gamma) * diag(demography_vector) %*% contact_matrix %*%
#   diag(1 / demography_vector). finalsize also requires
#   contact_matrix * demography_vector to have dominant real eigenvalue 1, so
#   the finalsize contact matrix is the normalized next-generation matrix
#   divided by the demography vector.

age_structure <- AgeStructure(
  age_groups = c("0-19", "20-64", "65+"),
  lower_bounds = c(0, 20, 65),
  upper_bounds = c(19, 64, Inf)
)

population <- c(900, 1700, 650)
demography_vector <- population / sum(population)
initial_infections <- population * 1e-6
initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 3),
  value = c(
    population - initial_infections,
    initial_infections,
    rep(0, age_structure$n_age_groups)
  ),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(c(
  7.0, 2.4, 0.7,
  1.2, 5.0, 1.6,
  0.5, 2.1, 3.0
), nrow = age_structure$n_age_groups, byrow = TRUE)
dimnames(contact_matrix) <- list(age_structure$age_groups, age_structure$age_groups)

gamma <- 1 / 4
beta <- 0.09
agepi_next_generation_shape <- diag(demography_vector) %*%
  contact_matrix %*%
  diag(1 / demography_vector)
shape_dominant_eigenvalue <- max(Re(eigen(
  agepi_next_generation_shape,
  only.values = TRUE
)$values))
r0_agepi <- (beta / gamma) * shape_dominant_eigenvalue
finalsize_contact_matrix <- (agepi_next_generation_shape / shape_dominant_eigenvalue) /
  demography_vector
times <- seq(0, 365, by = 1)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = times,
  model = SIRModel(gamma = gamma),
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = beta,
  method = if (requireNamespace("deSolve", quietly = TRUE)) "deSolve" else "euler",
  cumulative_flows = list(
    infections = list(from = "S", to = "I")
  )
)

final_cumulative <- simulation$cumulative[
  simulation$cumulative$time == max(times) &
    simulation$cumulative$cumulative_name == "infections",
  ,
  drop = FALSE
]
final_cumulative <- final_cumulative[
  match(age_structure$age_groups, final_cumulative$age_group),
  ,
  drop = FALSE
]

comparison <- data.frame(
  age_group = age_structure$age_groups,
  population = population,
  agepi_cumulative_infections = final_cumulative$value,
  agepi_attack_rate = final_cumulative$value / population,
  stringsAsFactors = FALSE
)

cat("agepi deterministic SIR final-size validation example\n")
cat("Closed population; no births, deaths, migration, waning, or reinfection.\n")
cat("Contact matrix orientation: recipient rows, infectious/source columns.\n")
cat("finalsize contact matrix normalization factor:", shape_dominant_eigenvalue, "\n")
cat("beta =", beta, "gamma =", gamma, "implied agepi R0 =", r0_agepi, "\n")
cat("Cumulative infection transition_id:", unique(final_cumulative$transition_id), "\n\n")

if (requireNamespace("finalsize", quietly = TRUE)) {
  finalsize_result <- finalsize::final_size(
    r0 = r0_agepi,
    contact_matrix = finalsize_contact_matrix,
    demography_vector = demography_vector,
    susceptibility = matrix(1, nrow = age_structure$n_age_groups, ncol = 1),
    p_susceptibility = matrix(1, nrow = age_structure$n_age_groups, ncol = 1)
  )

  finalsize_attack_rate <- finalsize_result$p_infected[
    match(age_structure$age_groups, as.character(finalsize_result$demo_grp))
  ]
  if (anyNA(finalsize_attack_rate)) {
    finalsize_attack_rate <- finalsize_result$p_infected
  }

  comparison$finalsize_attack_rate <- finalsize_attack_rate
  comparison$finalsize_infections <- finalsize_attack_rate * population
  comparison$attack_rate_difference <- (
    comparison$agepi_attack_rate - comparison$finalsize_attack_rate
  )

  cat("finalsize::final_size() comparison\n")
  print(comparison, row.names = FALSE)
  cat("\n")
  cat("Maximum absolute attack-rate difference:",
      max(abs(comparison$attack_rate_difference)), "\n")
} else {
  cat("Package finalsize is not installed; skipping external comparison.\n")
  cat("Install finalsize to run the optional finalsize::final_size() validation.\n")
  print(comparison, row.names = FALSE)
}
