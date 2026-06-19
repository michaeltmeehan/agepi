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

if (!exists("kiribati_tb_parameters", mode = "function")) {
  for (helper in c(
    "kiribati_tb_parameters",
    "kiribati_tb_model",
    "kiribati_tb_initial_proportions",
    "tb_cumulative_flows",
    "summarise_tb_burden",
    "kiribati_tb_target_map",
    "print_kiribati_tb_summary"
  )) {
    assign(helper, getFromNamespace(helper, "agepi"))
  }
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

country <- "Kiribati"
simulation_years <- 2025:2030
simulation_times <- seq(2025, 2030, length.out = (2030 - 2025) * 4 + 1)
simulation_method <- "deSolve"
demographic_time_policy <- "step"
output_timestep_years <- 1 / 4

age_structure <- wpp_age_structure_1year(max_age = 95)
reporting_age_structure <- AgeStructure(
  c("0-4", "5-14", "15-24", "25-44", "45-64", "65+"),
  c(0, 5, 15, 25, 45, 65),
  c(4, 14, 24, 44, 64, Inf)
)
reporting_age_mapping <- AgeGridMapping(age_structure, reporting_age_structure, open_ended = "include")

demography <- demographic_process_from_wpp(
  country = country,
  years = simulation_years,
  age_structure = age_structure,
  migration = TRUE,
  fertility_exposure_fraction = 0.5
)
wpp_population <- demography$population
demographic_process <- demography$demographic_process
fertility_schedule <- demography$fertility_schedule
mortality_schedule <- demography$mortality_schedule
migration_schedule <- demography$migration_schedule
tfr_input <- demography$inputs$tfr
tfr_years <- sort(unique(demography$inputs$fertility_weights$year))
process_years <- simulation_years

contact_matrix <- contact_matrix_for_age_structure(age_structure, source = "polymod_uk")
contact_matrix_metadata <- attr(contact_matrix, "contact_source")
contact_matrix_source <- data.frame(as.list(contact_matrix_metadata), stringsAsFactors = FALSE)

tb_parameters <- kiribati_tb_parameters(age_structure)
tb_model <- kiribati_tb_model(tb_parameters)
population_2025 <- demography_population_vector(wpp_population, time = min(simulation_years))
initial_state <- initialise_compartments_from_proportions(
  population = population_2025,
  proportions = kiribati_tb_initial_proportions(age_structure),
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
  cumulative_flows = tb_cumulative_flows()
)

trajectory <- tb_output$trajectory
cumulative_flows <- tb_output$cumulative
compartment_summary <- compartment_totals(trajectory)
age_group_summary <- age_group_totals(trajectory)
population_summary <- total_population(trajectory)
tb_burden_summary <- summarise_tb_burden(tb_output)
cumulative_flow_summary <- cumulative_flow_totals(cumulative_flows)
cumulative_flow_wide <- cumulative_flow_totals_wide(cumulative_flows)
annual_flow_summary <- cumulative_flow_increments(
  cumulative_flows,
  times = simulation_years,
  value_col = "annual_value"
)
public_data_target_summary <- kiribati_tb_target_map()

age_group_summary_broad <- aggregate_population_summary_age_grid(age_group_summary, reporting_age_mapping)
trajectory_broad <- aggregate_epidemic_trajectory_age_grid(trajectory, reporting_age_mapping)
cumulative_flows_broad <- aggregate_cumulative_flows_age_grid(cumulative_flows, reporting_age_mapping)

print_kiribati_tb_summary(
  country = country,
  simulation_years = simulation_years,
  output_timestep_years = output_timestep_years,
  simulation_method = simulation_method,
  contact_matrix_source = contact_matrix_source,
  population_summary = population_summary,
  tb_burden_summary = tb_burden_summary,
  age_group_summary_broad = age_group_summary_broad,
  cumulative_flow_wide = cumulative_flow_wide,
  public_data_target_summary = public_data_target_summary
)
}
