# Kiribati TB scaffold parameters
#
# Returns age-specific transition rates and transmission modifiers used by the
# public-data Kiribati TB example. These values are provisional, uncalibrated
# scaffold assumptions, not Kiribati-specific estimates.
#
# @param age_structure Valid age structure.
#
# @return Named list of TB parameters.
kiribati_tb_parameters <- function(age_structure) {
  validate_age_structure(age_structure)
  age_lower <- age_structure$lower_bounds
  age_groups <- age_structure$age_groups

  parameters <- list(
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

  age_specific_parameters <- c(
    "fast_progression",
    "reactivation",
    "relapse",
    "susceptibility",
    "infectiousness"
  )
  for (parameter_name in age_specific_parameters) {
    names(parameters[[parameter_name]]) <- age_groups
  }

  parameters
}

# Kiribati TB scaffold compartment model
#
# The transition structure is part of the provisional, uncalibrated Kiribati TB
# example scaffold.
#
# @param parameters List returned by `kiribati_tb_parameters()`.
#
# @return A CompartmentModel object.
kiribati_tb_model <- function(parameters) {
  validate_kiribati_tb_parameters(parameters)

  transitions <- data.frame(
    from = c("Lr", "Lr", "Ld", "I", "T", "R"),
    to = c("Ld", "I", "I", "T", "R", "I"),
    stringsAsFactors = FALSE
  )
  transitions$rate <- I(list(
    parameters$recent_to_remote,
    parameters$fast_progression,
    parameters$reactivation,
    parameters$treatment_initiation,
    parameters$treatment_completion,
    parameters$relapse
  ))

  CompartmentModel(
    compartments = c("S", "Lr", "Ld", "I", "T", "R"),
    infection_transitions = data.frame(from = "S", to = "Lr", stringsAsFactors = FALSE),
    transitions = transitions,
    infectious_compartments = "I",
    birth_compartment = "S",
    migration_compartment = "S"
  )
}

# Kiribati TB scaffold initial compartment proportions
#
# These values are provisional, uncalibrated scaffold assumptions. `S` is
# intended to be assigned as the residual compartment.
#
# @param age_structure Valid age structure.
#
# @return Named list of non-susceptible initial proportions.
kiribati_tb_initial_proportions <- function(age_structure) {
  validate_age_structure(age_structure)
  age_lower <- age_structure$lower_bounds

  active_prev <- ifelse(age_lower < 15, 0.00025,
    ifelse(age_lower < 25, 0.0009,
      ifelse(age_lower < 45, 0.0015,
        ifelse(age_lower < 65, 0.0022, 0.0020)
      )
    )
  )

  proportions <- list(
    Lr = ifelse(age_lower < 15, 0.004,
      ifelse(age_lower < 65, 0.008, 0.006)
    ),
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
      ifelse(age_lower < 45, 0.035,
        ifelse(age_lower < 65, 0.055, 0.065)
      )
    )
  )

  for (compartment in names(proportions)) {
    names(proportions[[compartment]]) <- age_structure$age_groups
  }

  proportions
}

# TB cumulative-flow specification for the Kiribati TB scaffold
#
# @return Named list suitable for `simulate_deterministic(cumulative_flows =)`.
tb_cumulative_flows <- function() {
  list(
    infections = list(from = "S", to = "Lr"),
    progression_to_active_tb = list(from = c("Lr", "Ld"), to = c("I", "I")),
    treatment_initiation = list(from = "I", to = "T"),
    treatment_completion = list(from = "T", to = "R"),
    relapse_recurrent_tb = list(from = "R", to = "I")
  )
}

# Summarise TB burden from the Kiribati TB scaffold output
#
# @param output Simulation output data frame, or a list with a `trajectory`
#   component returned by `simulate_deterministic()`.
#
# @return Data frame with active TB, treatment, latent infection, population,
#   and active TB prevalence per 100,000.
summarise_tb_burden <- function(output) {
  trajectory <- if (is.list(output) && !is.data.frame(output) && !is.null(output$trajectory)) {
    output$trajectory
  } else {
    output
  }
  validate_simulation_output_summary_input(trajectory)

  active_tb <- trajectory[trajectory$compartment == "I", ]
  on_treatment <- trajectory[trajectory$compartment == "T", ]
  latent_tb <- trajectory[trajectory$compartment %in% c("Lr", "Ld"), ]

  burden <- merge(
    stats::aggregate(value ~ time, active_tb, sum),
    stats::aggregate(value ~ time, on_treatment, sum),
    by = "time",
    suffixes = c("_active_tb", "_on_treatment")
  )
  names(burden) <- c("time", "active_tb_prevalence_count", "on_treatment_count")

  latent_summary <- stats::aggregate(value ~ time, latent_tb, sum)
  names(latent_summary)[names(latent_summary) == "value"] <- "latent_tb_infection_count"
  burden <- merge(burden, latent_summary, by = "time")

  population <- total_population(trajectory)
  names(population)[names(population) == "value"] <- "population"
  burden <- merge(burden, population, by = "time")
  burden$active_tb_prevalence_per_100k <-
    1e5 * burden$active_tb_prevalence_count / burden$population

  burden
}

# Kiribati TB public-data target map
#
# This map is part of the example scaffold. It documents later calibration
# targets but does not implement calibration.
#
# @return Calibration-facing target map used by the Kiribati TB example.
kiribati_tb_target_map <- function() {
  data.frame(
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
}

# Print the Kiribati TB example summary
#
# @param country Country name.
# @param simulation_years Numeric simulation years.
# @param output_timestep_years Output timestep size in years.
# @param simulation_method Solver method label.
# @param contact_matrix_source Contact source metadata data frame.
# @param population_summary Population summary data frame.
# @param tb_burden_summary TB burden summary data frame.
# @param age_group_summary_broad Broad age-group summary data frame.
# @param cumulative_flow_wide Wide cumulative-flow summary.
# @param public_data_target_summary Public-data target map.
#
# @return Invisibly returns `NULL`.
print_kiribati_tb_summary <- function(country,
                                      simulation_years,
                                      output_timestep_years,
                                      simulation_method,
                                      contact_matrix_source,
                                      population_summary,
                                      tb_burden_summary,
                                      age_group_summary_broad,
                                      cumulative_flow_wide,
                                      public_data_target_summary) {
  cat("\nKiribati TB public-data demography scaffold\n")
  cat("Country:", country, "\n")
  cat("Simulation years:", min(simulation_years), "-", max(simulation_years), "\n", sep = "")
  timestep_label <- if (abs(output_timestep_years - 1 / 4) < 1e-12) {
    "quarterly"
  } else {
    paste(output_timestep_years, "years")
  }
  cat("Output timestep: ", timestep_label, "; solver: ", simulation_method, "\n", sep = "")
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

  invisible(NULL)
}

validate_kiribati_tb_parameters <- function(parameters) {
  required <- c(
    "beta",
    "recent_to_remote",
    "fast_progression",
    "reactivation",
    "treatment_initiation",
    "treatment_completion",
    "relapse",
    "susceptibility",
    "infectiousness"
  )
  missing_parameters <- setdiff(required, names(parameters))
  if (length(missing_parameters) > 0) {
    stop(
      "parameters is missing required field(s): ",
      paste(missing_parameters, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(parameters)
}
