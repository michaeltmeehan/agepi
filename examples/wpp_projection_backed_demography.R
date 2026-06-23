# Projection-backed demography from WPP age-specific population projections.
#
# This pathway replays WPP's pre-computed age-specific population trajectory as
# an external background input. WPP has already combined its fertility,
# mortality, migration, sex structure, and cohort-component assumptions.
#
# This is different from mechanistic agepi demography, where schedules are used
# to simulate births, deaths, ageing, and migration through time.

# The example uses `population_trajectory_from_wpp()` to replay the projected
# age-specific population table, then `projection_population_vector()` to pull
# one exact year back out of that trajectory.

if (!"package:agepi" %in% search()) {
  if (dir.exists("R")) {
    invisible(lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), source))
  } else if (requireNamespace("agepi", quietly = TRUE)) {
    library(agepi)
  }
}

if (!requireNamespace("wpp2024", quietly = TRUE)) {
  message("Install wpp2024 to run this example with popprojAge1dt.")
} else {
  data("popprojAge1dt", package = "wpp2024")

  ages <- wpp_age_structure_1year(max_age = 100)

  # The public API here is the projection replay helper plus exact-time
  # accessors over the resulting age-specific table.
  projection <- population_trajectory_from_wpp(
    popprojAge1dt,
    age_structure = ages,
    location = "Kiribati",
    years = 2025:2030,
    location_col = "name",
    time_col = "year",
    age_group_col = "age",
    population_col = "pop"
  )

  print(head(projection))

  # Exact-time lookup returns the age-specific vector for one projection year.
  population_2027 <- projection_population_vector(
    projection,
    time = 2027,
    time_policy = "exact"
  )
  print(sum(population_2027))
}
