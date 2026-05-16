if (requireNamespace("agepi", quietly = TRUE)) {
  library(agepi)
} else if (dir.exists("R")) {
  invisible(lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), source))
} else {
  stop(
    "Package agepi is not installed. Run this script from the package root ",
    "or install agepi first.",
    call. = FALSE
  )
}

# This is a tiny dependency-free workflow using mock WPP-like tables.
# It demonstrates agepi's demographic adapters and diagnostics only.
# The data are invented and do not represent or reproduce WPP projections.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

times <- c(2020, 2021, 2022)
schedule_times <- times[-length(times)]

mortality_like <- data.frame(
  year = rep(schedule_times, each = age_structure$n_age_groups),
  age = rep(age_structure$age_groups, times = length(schedule_times)),
  mx = c(
    0.006, 0.004, 0.020,
    0.006, 0.004, 0.021
  ),
  stringsAsFactors = FALSE
)

fertility_like <- data.frame(
  year = schedule_times,
  age = rep("5-9", length(schedule_times)),
  asfr = c(0.030, 0.028),
  stringsAsFactors = FALSE
)

migration_like <- data.frame(
  year = rep(schedule_times, each = age_structure$n_age_groups),
  age = rep(age_structure$age_groups, times = length(schedule_times)),
  net_migration = c(
    1.0, -0.5, 0.2,
    0.8, -0.4, 0.1
  ),
  stringsAsFactors = FALSE
)

mortality_schedule <- standardise_wpp_mortality(
  mortality_like,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  mortality_col = "mx"
)

fertility_schedule <- standardise_wpp_fertility(
  fertility_like,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  fertility_col = "asfr"
)

migration_schedule <- standardise_wpp_migration(
  migration_like,
  age_structure = age_structure,
  time_col = "year",
  age_col = "age",
  migration_col = "net_migration",
  migration_type = "count"
)

process <- build_demographic_process(
  age_structure = age_structure,
  fertility_schedule = fertility_schedule,
  mortality_schedule = mortality_schedule,
  migration_schedule = migration_schedule,
  mode = "migration"
)

initial_population <- c(500, 450, 800)

simulated <- simulate_demography(
  process = process,
  initial_state = initial_population,
  times = times
)

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

comparison <- compare_demography_to_observed(simulated, observed)
summary <- summarise_demography_comparison(comparison)
residual <- implied_demographic_residual(observed, process)

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

cat("\nResidual converted to per-time migration_count flow:\n")
print(residual_migration$data)
