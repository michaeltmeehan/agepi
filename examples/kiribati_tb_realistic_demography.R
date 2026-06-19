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
# Main public-data anchors from docs/disease_notes/kiribati_tb_modelling_inputs.md:
# - WPP 2024 national population, fertility, mortality, and migration;
# - WHO TB incidence, notifications, mortality, age/sex incidence, treatment
#   outcomes, contact/TPT, and MDR/RR-TB data as later calibration targets;
# - WUENIC BCG coverage as an optional later paediatric severe-TB modifier.
#
# TB mortality is not included as a disease transition here because agepi's
# current CompartmentModel() transitions move people between model compartments,
# while background demography already handles all-cause mortality. The mortality
# target is exposed below as a calibration-facing quantity to add once an
# explicit disease-death convention is available.

if (!requireNamespace("wpp2024", quietly = TRUE) ||
    !requireNamespace("socialmixr", quietly = TRUE) ||
    !requireNamespace("deSolve", quietly = TRUE)) {
  message(
    "Skipping Kiribati TB realistic demography example: optional package ",
    "wpp2024, socialmixr, or deSolve is not installed. Install these ",
    "suggested packages to use public WPP inputs, published contact data, ",
    "and the deSolve backend."
  )
} else {

simulation_years <- 2025:2030
simulation_times <- seq(2025, 2030, length.out = (2030 - 2025) * 4 + 1)
process_years <- 2025:2030
simulation_method <- "deSolve"
demographic_time_policy <- "step"
output_timestep_years <- 1 / 4
country <- "Kiribati"

age_structure <- wpp_age_structure_1year(max_age = 95)
reporting_age_structure <- AgeStructure(
  age_groups = c("0-4", "5-14", "15-24", "25-44", "45-64", "65+"),
  lower_bounds = c(0, 5, 15, 25, 45, 65),
  upper_bounds = c(4, 14, 24, 44, 64, Inf)
)
reporting_age_mapping <- AgeGridMapping(
  from_age_groups = age_structure,
  to_age_groups = reporting_age_structure,
  open_ended = "include"
)
age_groups <- age_structure$age_groups
age_lower <- age_structure$lower_bounds

population_data <- wpp_dataset("popprojAge1dt")
fertility_weight_data <- wpp_dataset("percentASFR1dt")
tfr_data <- wpp_dataset("tfrproj1dt")
mortality_data <- wpp_dataset("mx1dt")
migration_data <- wpp_dataset("migprojAge1dt")

population_input <- wpp_location_rows(
  population_data,
  location = country,
  columns = c("name", "year", "age", "pop")
)
population_input$year <- as.numeric(population_input$year)
population_input <- population_input[population_input$year %in% simulation_years, ]
population_input$pop <- 1000 * as.numeric(population_input$pop)
population_input <- collapse_wpp_open_age_counts(population_input, "pop", open_age = 95)

wpp_population <- population_from_wpp(
  data = population_input,
  age_structure = age_structure,
  time_col = "year",
  age_group_col = "age",
  population_col = "pop",
  location = country,
  location_col = "name"
)

fertility_weights <- wpp_location_rows(
  fertility_weight_data,
  location = country,
  columns = c("year", "age", "pasfr")
)
fertility_weights$year <- as.numeric(fertility_weights$year)
fertility_weights <- fertility_weights[fertility_weights$year %in% process_years, ]

tfr_source <- wpp_location_rows(tfr_data, location = country, columns = c("year", "tfr"))
tfr_source$year <- as.numeric(tfr_source$year)
tfr_source$tfr <- as.numeric(tfr_source$tfr)
tfr_years <- sort(unique(fertility_weights$year))
tfr_input <- tfr_source[tfr_source$year %in% tfr_years, , drop = FALSE]
missing_tfr_years <- setdiff(tfr_years, tfr_input$year)
if (length(missing_tfr_years) > 0) {
  stop(
    "Annual WPP TFR input is missing year(s): ",
    paste(missing_tfr_years, collapse = ", "),
    call. = FALSE
  )
}

fertility_schedule <- fertility_from_wpp_percent_asfr(
  data = fertility_weights,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  weight_col = "pasfr",
  tfr_data = tfr_input,
  tfr_time_col = "year",
  tfr_col = "tfr",
  weight_type = "percent",
  maternal_age_groups = as.character(10:54),
  tolerance = 1e-5
)

mortality_input <- wpp_location_rows(
  mortality_data,
  location = country,
  columns = c("year", "age", "mxB")
)
mortality_input$year <- as.numeric(mortality_input$year)
mortality_input <- mortality_input[mortality_input$year %in% process_years, ]
mortality_input <- collapse_wpp_open_age_rates(mortality_input, "mxB", open_age = 95)

mortality_schedule <- mortality_from_wpp(
  data = mortality_input,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  mortality_col = "mxB",
  quantity = "mx"
)

migration_input <- wpp_location_rows(
  migration_data,
  location = country,
  columns = c("year", "age", "mig")
)
migration_input$year <- as.numeric(migration_input$year)
migration_input <- migration_input[migration_input$year %in% process_years, ]
migration_input$mig <- 1000 * as.numeric(migration_input$mig)
migration_input <- collapse_wpp_open_age_counts(migration_input, "mig", open_age = 95)

migration_schedule <- standardise_wpp_migration(
  data = migration_input,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  migration_col = "mig",
  migration_type = "count"
)

demographic_process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = fertility_schedule,
  fertility_exposure_fraction = 0.5,
  mortality_schedule = mortality_schedule,
  migration_schedule = migration_schedule,
  mode = "migration"
)

# Published proxy contact matrix. No Kiribati empirical or synthetic matrix is
# bundled here, so this scaffold uses POLYMOD United Kingdom contacts through
# socialmixr as a transparent proxy. This should be replaced by Prem/conmat or
# Kiribati/Pacific-specific contacts when available.
contact_matrix <- contact_matrix_for_age_structure(age_structure, source = "polymod_uk")
contact_matrix_metadata <- attr(contact_matrix, "contact_source")
contact_matrix_source <- data.frame(
  source_label = contact_matrix_metadata$source_label,
  source_reference = contact_matrix_metadata$source_reference,
  source_age_grid = contact_matrix_metadata$source_age_grid,
  model_age_grid = contact_matrix_metadata$model_age_grid,
  expansion_note = contact_matrix_metadata$expansion_note,
  limitation = contact_matrix_metadata$limitation,
  stringsAsFactors = FALSE
)

# Epidemiological parameters. These are plausible starting values from TB
# natural-history reasoning and the modelling note's proposed structure. Treat
# them as priors or calibration placeholders, not definitive Kiribati estimates.
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
  susceptibility = ifelse(age_lower < 15, 0.75,
    ifelse(age_lower < 65, 1.0, 1.15)
  ),
  infectiousness = ifelse(age_lower < 15, 0.20,
    ifelse(age_lower < 65, 1.0, 0.85)
  )
)
names(tb_parameters$fast_progression) <- age_groups
names(tb_parameters$reactivation) <- age_groups
names(tb_parameters$relapse) <- age_groups
names(tb_parameters$susceptibility) <- age_groups
names(tb_parameters$infectiousness) <- age_groups

transitions <- data.frame(
  from = c("Lr", "Lr", "Ld", "I", "T", "R"),
  to = c("Ld", "I", "I", "T", "R", "I"),
  stringsAsFactors = FALSE
)
transitions$rate <- I(list(
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
  transitions = transitions,
  infectious_compartments = "I",
  birth_compartment = "S",
  migration_compartment = "S"
)

population_2025 <- demography_population_vector(wpp_population, time = 2025)

remote_latent_prev <- ifelse(age_lower < 5, 0.01,
  ifelse(age_lower < 15, 0.04,
    ifelse(age_lower < 25, 0.14,
      ifelse(age_lower < 45, 0.24,
        ifelse(age_lower < 65, 0.32, 0.36)
      )
    )
  )
)
recent_latent_prev <- ifelse(age_lower < 15, 0.004,
  ifelse(age_lower < 65, 0.008, 0.006)
)
active_prev <- ifelse(age_lower < 15, 0.00025,
  ifelse(age_lower < 25, 0.0009,
    ifelse(age_lower < 45, 0.0015,
      ifelse(age_lower < 65, 0.0022, 0.0020)
    )
  )
)
treatment_prev <- 0.7 * active_prev
recovered_prev <- ifelse(age_lower < 15, 0.005,
  ifelse(age_lower < 45, 0.035,
    ifelse(age_lower < 65, 0.055, 0.065)
  )
)

initial_state <- initialise_compartments_from_proportions(
  population = population_2025,
  proportions = list(
    Lr = recent_latent_prev,
    Ld = remote_latent_prev,
    I = active_prev,
    T = treatment_prev,
    R = recovered_prev
  ),
  residual_compartment = "S",
  compartments = tb_model$compartments,
  age_structure = age_structure
)

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
  cumulative_flows = list(
    infections = list(from = "S", to = "Lr"),
    progression_to_active_tb = list(from = c("Lr", "Ld"), to = c("I", "I")),
    treatment_initiation = list(from = "I", to = "T"),
    treatment_completion = list(from = "T", to = "R"),
    relapse_recurrent_tb = list(from = "R", to = "I")
  )
)

trajectory <- tb_output$trajectory
cumulative_flows <- tb_output$cumulative
compartment_summary <- compartment_totals(trajectory)
age_group_summary <- age_group_totals(trajectory)
population_summary <- total_population(trajectory)

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
population_for_rates <- population_summary
names(population_for_rates)[names(population_for_rates) == "value"] <- "population"
tb_burden_summary <- merge(tb_burden_summary, population_for_rates, by = "time")
tb_burden_summary$active_tb_prevalence_per_100k <-
  1e5 * tb_burden_summary$active_tb_prevalence_count /
    tb_burden_summary$population

cumulative_flow_summary <- cumulative_flow_totals(cumulative_flows)
cumulative_flow_wide <- cumulative_flow_totals_wide(cumulative_flows)
annual_flow_summary <- cumulative_flow_increments(
  cumulative_flows,
  times = simulation_years,
  value_col = "annual_value"
)

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

age_group_summary_broad <- aggregate_population_summary_age_grid(
  age_group_summary,
  reporting_age_mapping
)
trajectory_broad <- aggregate_epidemic_trajectory_age_grid(
  trajectory,
  reporting_age_mapping
)
cumulative_flows_broad <- aggregate_cumulative_flows_age_grid(
  cumulative_flows,
  reporting_age_mapping
)

cat("\nKiribati TB public-data demography scaffold\n")
cat("Country:", country, "\n")
cat("Simulation years:", min(simulation_years), "-", max(simulation_years), "\n", sep = "")
cat("Output timestep: quarterly; solver: deSolve\n")
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
