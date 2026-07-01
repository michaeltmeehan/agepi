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

# Purpose: national public-data scaffold for a Kiribati TB model with WPP-backed
# demography. It is deterministic, age-structured, and calibration-ready, but it
# is not calibrated and must not be interpreted as a policy estimate.
#
# Public-data anchors from docs/disease_notes/kiribati_tb_modelling_inputs.md:
# - WPP 2024 national population, fertility, mortality, and migration;
# - WHO TB incidence, notifications, mortality, age/sex incidence, treatment
#   outcomes, contact/TPT, and MDR/RR-TB data as later calibration targets;
# - WUENIC BCG coverage as an optional later paediatric severe-TB modifier.
#
# TB mortality is not included as a disease transition here because agepi's
# current CompartmentModel() transitions move people between model compartments,
# while background demography already handles all-cause mortality.

# Optional Packages -----------------------------------------------------------

optional_packages <- c("wpp2024")
missing_optional_packages <- optional_packages[
  !vapply(optional_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_optional_packages) > 0) {
  message(
    "Skipping Kiribati TB realistic demography example: optional package(s) ",
    paste(missing_optional_packages, collapse = ", "),
    " are not installed. Install wpp2024 with pak::pkg_install(\"PPgp/wpp2024\") ",
    "or remotes::install_github(\"PPgp/wpp2024\") to use the public WPP inputs."
  )
} else {

# Simulation Settings ---------------------------------------------------------

country <- "Kiribati"
simulation_years <- 2025:2030
output_timestep_years <- 1 / 4
simulation_times <- seq(
  min(simulation_years),
  max(simulation_years),
  by = output_timestep_years
)
simulation_method <- "deSolve"
demographic_time_policy <- "step"

age_structure <- wpp_age_structure_1year(max_age = 95)
reporting_age_structure <- AgeStructure(
  c("0-4", "5-14", "15-24", "25-44", "45-64", "65+"),
  c(0, 5, 15, 25, 45, 65),
  c(4, 14, 24, 44, 64, Inf)
)
reporting_age_mapping <- AgeGridMapping(age_structure, reporting_age_structure, open_ended = "include")

# WPP Demography --------------------------------------------------------------

demography <- demographic_process_from_wpp(
  country = country,
  years = simulation_years,
  age_structure = age_structure,
  migration = TRUE,
  fertility_exposure_fraction = 0.5
)
wpp_population <- demography$population
demographic_process <- demography$demographic_process
population_2025 <- demography_population_vector(wpp_population, time = min(simulation_years))

# Contact Matrix --------------------------------------------------------------

# Prefer a Prem/contactdata Kiribati matrix when the installed dataset includes
# that country. Otherwise use POLYMOD UK as a documented proxy only.
if (requireNamespace("contactdata", quietly = TRUE) &&
    country %in% contactdata::list_countries()) {
  source_contacts <- load_contact_matrix_source(
    source = "prem",
    country = country,
    setting = "all"
  )
  message("Using Prem/contactdata Kiribati synthetic contact matrix.")
} else {
  if (!requireNamespace("socialmixr", quietly = TRUE)) {
    message(
      "Skipping Kiribati TB realistic demography example: neither contactdata ",
      "nor socialmixr is available for the contact-matrix fallback. Install ",
      "contactdata with Kiribati support or install socialmixr to run this ",
      "example."
    )
    quit(status = 0, save = "no")
  }
  message(
    "Prem/contactdata Kiribati contact matrix is unavailable; falling back to ",
    "POLYMOD UK as a European proxy. This is a limitation and should not be ",
    "interpreted as Kiribati-specific contact evidence."
  )
  source_contacts <- load_contact_matrix_source("polymod_uk")
  source_contacts$limitations <- c(
    source_contacts$limitations,
    "Used only because Prem/contactdata Kiribati was unavailable in this R environment."
  )
}
contact_matrix <- adapt_contact_matrix_to_age_structure(
  source_contacts,
  age_structure,
  method = "source_band"
)
contact_matrix_source <- data.frame(
  as.list(attr(contact_matrix, "contact_source")),
  stringsAsFactors = FALSE
)

# Provisional TB Assumptions --------------------------------------------------

# These TB parameters are provisional and uncalibrated scaffold assumptions.
# They are not Kiribati-specific estimates and should be replaced or calibrated
# before any policy interpretation.
age_lower <- age_structure$lower_bounds
tb_parameters <- list(
  beta = 0.13,
  recent_to_remote = 0.35,
  fast_progression = ifelse(age_lower < 5, 0.08,
    ifelse(age_lower < 15, 0.025,
      ifelse(age_lower < 65, 0.045, 0.06)
    )
  ),
  reactivation = ifelse(age_lower < 15, 0.0005,
    ifelse(age_lower < 25, 0.0015,
      ifelse(age_lower < 45, 0.0025,
        ifelse(age_lower < 65, 0.004, 0.006)
      )
    )
  ),
  treatment_initiation = 0.60,
  treatment_completion = 0.75,
  relapse = ifelse(age_lower < 15, 0.001,
    ifelse(age_lower < 45, 0.003,
      ifelse(age_lower < 65, 0.0045, 0.006)
    )
  ),
  susceptibility = ifelse(age_lower < 15, 0.75, ifelse(age_lower < 65, 1.0, 1.15)),
  infectiousness = ifelse(age_lower < 15, 0.20, ifelse(age_lower < 65, 1.0, 0.85))
)

for (parameter_name in c(
  "fast_progression",
  "reactivation",
  "relapse",
  "susceptibility",
  "infectiousness"
)) {
  names(tb_parameters[[parameter_name]]) <- age_structure$age_groups
}

tb_transitions <- data.frame(
  from = c("Lr", "Lr", "Ld", "I", "T", "R"),
  to = c("Ld", "I", "I", "T", "R", "I"),
  stringsAsFactors = FALSE
)
tb_transitions$rate <- I(list(
  tb_parameters$recent_to_remote,
  tb_parameters$fast_progression,
  tb_parameters$reactivation,
  tb_parameters$treatment_initiation,
  tb_parameters$treatment_completion,
  tb_parameters$relapse
))

tb_model <- CompartmentModel(
  compartments = c("S", "Lr", "Ld", "I", "T", "R"),
  infection_transitions = data.frame(from = "S", to = "Lr", stringsAsFactors = FALSE),
  transitions = tb_transitions,
  infectious_compartments = "I",
  birth_compartment = "S",
  migration_compartment = "S"
)

tb_cumulative_flow_specs <- list(
  infections = list(from = "S", to = "Lr"),
  progression_to_active_tb = list(from = c("Lr", "Ld"), to = c("I", "I")),
  treatment_initiation = list(from = "I", to = "T"),
  treatment_completion = list(from = "T", to = "R"),
  relapse_recurrent_tb = list(from = "R", to = "I")
)

# Initial Conditions ----------------------------------------------------------

# These initial proportions are provisional and uncalibrated. `S` is assigned
# as the residual population after latent, active, treatment, and recovered
# states are initialised.
active_prev <- ifelse(age_lower < 15, 0.00025,
  ifelse(age_lower < 25, 0.0009,
    ifelse(age_lower < 45, 0.0015,
      ifelse(age_lower < 65, 0.0022, 0.0020)
    )
  )
)
initial_proportions <- list(
  Lr = ifelse(age_lower < 15, 0.004, ifelse(age_lower < 65, 0.008, 0.006)),
  Ld = ifelse(age_lower < 5, 0.01,
    ifelse(age_lower < 15, 0.04,
      ifelse(age_lower < 25, 0.14,
        ifelse(age_lower < 45, 0.24,
          ifelse(age_lower < 65, 0.32, 0.36)
        )
      )
    )
  ),
  I = active_prev,
  T = 0.7 * active_prev,
  R = ifelse(age_lower < 15, 0.005,
    ifelse(age_lower < 45, 0.035, ifelse(age_lower < 65, 0.055, 0.065))
  )
)
for (compartment in names(initial_proportions)) {
  names(initial_proportions[[compartment]]) <- age_structure$age_groups
}

initial_state <- initialise_compartments_from_proportions(
  population = population_2025,
  proportions = initial_proportions,
  residual_compartment = "S",
  compartments = tb_model$compartments,
  age_structure = age_structure
)

# Simulation ------------------------------------------------------------------

tb_output <- simulate_deterministic(
  initial_state = initial_state,
  times = simulation_times,
  model = tb_model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = tb_parameters$beta,
  susceptibility = tb_parameters$susceptibility,
  infectiousness = tb_parameters$infectiousness,
  demographic_process = demographic_process,
  time_policy = demographic_time_policy,
  migration_policy = "proportional",
  method = simulation_method,
  cumulative_flows = tb_cumulative_flow_specs
)

trajectory <- tb_output$trajectory
cumulative_flows <- tb_output$cumulative

# Compact Output --------------------------------------------------------------

age_group_summary <- age_group_totals(trajectory)
population_summary <- total_population(trajectory)
age_group_summary_broad <- aggregate_population_summary_age_grid(age_group_summary, reporting_age_mapping)
cumulative_flow_wide <- cumulative_flow_totals_wide(cumulative_flows)

active_tb <- trajectory[trajectory$compartment == "I", ]
on_treatment <- trajectory[trajectory$compartment == "T", ]
latent_tb <- trajectory[trajectory$compartment %in% c("Lr", "Ld"), ]
tb_burden_summary <- merge(
  stats::aggregate(value ~ time, active_tb, sum),
  stats::aggregate(value ~ time, on_treatment, sum),
  by = "time",
  suffixes = c("_active_tb", "_on_treatment")
)
names(tb_burden_summary) <- c("time", "active_tb_prevalence_count", "on_treatment_count")
latent_summary <- stats::aggregate(value ~ time, latent_tb, sum)
names(latent_summary)[names(latent_summary) == "value"] <- "latent_tb_infection_count"
tb_burden_summary <- merge(tb_burden_summary, latent_summary, by = "time")
population_for_burden <- population_summary
names(population_for_burden)[names(population_for_burden) == "value"] <- "population"
tb_burden_summary <- merge(tb_burden_summary, population_for_burden, by = "time")
tb_burden_summary$active_tb_prevalence_per_100k <-
  1e5 * tb_burden_summary$active_tb_prevalence_count / tb_burden_summary$population

public_data_target_summary <- data.frame(
  target = c(
    "WHO estimated TB incidence count/rate",
    "WHO TB notifications count/rate",
    "WHO estimated TB mortality count/rate",
    "WHO treatment outcomes",
    "WHO age/sex incidence distribution",
    "WHO contact/TPT indicators",
    "WHO MDR/RR-TB estimates",
    "WUENIC BCG coverage"
  ),
  model_quantity = c(
    "annual progression_to_active_tb flow",
    "annual treatment_initiation flow as a crude notification proxy",
    "not modelled as a state transition in this scaffold",
    "annual treatment_completion among treatment starts",
    "progression_to_active_tb flow by reporting_age_group",
    "not modelled; scenario/calibration extension",
    "not modelled; scenario/reporting extension",
    "not modelled; optional paediatric severe-TB modifier"
  ),
  calibration_status = "placeholder; no formal calibration implemented",
  stringsAsFactors = FALSE
)

cat("\nKiribati TB public-data demography scaffold\n")
cat("Country:", country, "\n")
cat("Simulation years:", min(simulation_years), "-", max(simulation_years), "\n", sep = "")
cat("Output timestep: quarterly; solver: ", simulation_method, "\n", sep = "")
cat("Demography: annual WPP 2024 population, fertility, mortality, and net migration used directly with stepwise schedule lookup\n")
cat("Contact matrix:", contact_matrix_source$source_label, "\n")
cat("Status: deterministic calibration scaffold; not calibrated; not a policy estimate\n\n")

cat("Population summary:\n")
print(tail(population_summary, 6), row.names = FALSE)

cat("\nTB burden summary:\n")
print(tail(tb_burden_summary, 6), row.names = FALSE)

cat("\nBroad age-group summary at final year:\n")
print(
  age_group_summary_broad[age_group_summary_broad$time == max(simulation_years), ],
  row.names = FALSE
)

cat("\nCumulative flow summary:\n")
print(tail(cumulative_flow_wide, 6), row.names = FALSE)

cat("\nCalibration-facing public-data target map:\n")
print(public_data_target_summary, row.names = FALSE)
}
