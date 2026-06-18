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

# Purpose: demonstrate a demographic-only workflow using mock WPP-like tables.
# Main ingredients: age groups, standardised fertility/mortality/migration
# schedules, a demographic process, simulation, observed comparison, and
# residual diagnostics.
# Expected output: process mode, simulated populations, comparison summaries,
# selected age-specific errors, and residual migration-style flows.
# The data are invented and do not represent or reproduce WPP projections.

# The same AgeStructure is used to validate and align all demographic inputs.
# The final open age interval collects everyone aged 10 and older.
age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

# times are the population years to simulate. Schedule entries describe the
# intervals starting at 2020 and 2021, so the final time is not itself a new
# schedule start.
times <- c(2020, 2021, 2022)
schedule_times <- times[-length(times)]

# These mock input tables mimic the shape of external demographic data while
# staying small and local to the example.
mortality_like <- data.frame(
  year = rep(schedule_times, each = age_structure$n_age_groups),
  age = rep(age_structure$age_groups, times = length(schedule_times)),
  mx = c(
    0.006, 0.004, 0.020,
    0.006, 0.004, 0.021
  ),
  stringsAsFactors = FALSE
)

# Fertility is shown in the middle age group only to keep the toy example
# compact. Births produced by fertility enter the youngest age group.
fertility_like <- data.frame(
  year = schedule_times,
  age = rep("5-9", length(schedule_times)),
  asfr = c(0.030, 0.028),
  stringsAsFactors = FALSE
)

# Net migration can add or remove people by age group independently of births,
# deaths, and ageing.
migration_like <- data.frame(
  year = rep(schedule_times, each = age_structure$n_age_groups),
  age = rep(age_structure$age_groups, times = length(schedule_times)),
  net_migration = c(
    1.0, -0.5, 0.2,
    0.8, -0.4, 0.1
  ),
  stringsAsFactors = FALSE
)

# Standardisation adapters translate external column names into agepi schedule
# objects, checking that times and age groups match the model age structure.
mortality_schedule <- standardise_wpp_mortality(
  mortality_like,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  mortality_col = "mx"
)

# Fertility schedules are age-specific rates used by the demographic process
# to generate newborns during each simulated interval.
fertility_schedule <- standardise_wpp_fertility(
  fertility_like,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  fertility_col = "asfr"
)

# This migration schedule is interpreted as counts, not rates, because the
# source column represents net people per interval.
migration_schedule <- standardise_wpp_migration(
  migration_like,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  migration_col = "net_migration",
  migration_type = "count"
)

# The process combines natural change and migration. Ageing moves people
# between age groups, mortality removes people, fertility adds newborns, and
# migration supplies the remaining net movement.
process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = fertility_schedule,
  mortality_schedule = mortality_schedule,
  migration_schedule = migration_schedule,
  mode = "migration"
)

# Initial population is a vector ordered exactly as age_structure$age_groups.
initial_population <- c(500, 450, 800)

# simulate_demography returns projected age-specific population counts at each
# requested time, without any disease compartments.
simulated <- simulate_demography(
  process = process,
  initial_state = initial_population,
  times = times
)

# The observed object is an invented benchmark used to demonstrate diagnostics.
# It has the same long format expected from real demographic observations.
observed <- Demography(
  data.frame(
    time = rep(times, each = age_structure$n_age_groups),
    age_group = rep(age_structure$age_groups, times = length(times)),
    population = c(
      500.0, 450.0, 800.0,
      414.0, 488.0, 871.0,
      344.0, 504.0, 942.0
    ),
    stringsAsFactors = FALSE
  ),
  age_structure
)

# Comparison rows preserve age and time detail; the summary collapses those
# residuals into concise diagnostics. The implied residual asks what extra flow
# would be needed to reconcile the observed data with the specified process.
comparison <- compare_demography_to_observed(simulated, observed)
summary <- summarise_demography_comparison(comparison)
residual <- implied_demographic_residual(observed, process)

# Print progressively richer diagnostics: first the process setting, then the
# simulated population, aggregate error summaries, selected detailed errors,
# and finally residuals that can be converted back into migration schedules.
cat("Demographic process mode:\n")
print(process$mode)

cat("\nSimulated population:\n")
print(simulated)

cat("\nComparison summary:\n")
print(summary)

cat("\nSelected age-specific comparison rows:\n")
print(comparison[comparison$time %in% c(2021, 2022), c(
  "time",
  "age_group",
  "simulated_population",
  "observed_population",
  "absolute_error"
) ])

cat("\nImplied remaining residual by interval:\n")
print(residual[, c(
  "time_start",
  "time_end",
  "age_group",
  "residual_count",
  "residual_rate"
) ])

residual_migration <- residual_to_migration_schedule(
  residual,
  age_structure = age_structure,
  use = "count"
)

# Converting the residual to migration_count shows how diagnostics can be fed
# back into a migration-style schedule for calibration experiments.
cat("\nResidual converted to per-time migration_count flow:\n")
print(residual_migration$data)
