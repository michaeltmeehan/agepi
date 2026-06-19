kiribati_workflow_helper <- function(name) {
  if (exists(name, mode = "function")) {
    return(get(name, mode = "function"))
  }
  getFromNamespace(name, "agepi")
}

kiribati_tb_parameters <- kiribati_workflow_helper("kiribati_tb_parameters")
kiribati_tb_model <- kiribati_workflow_helper("kiribati_tb_model")
kiribati_tb_initial_proportions <- kiribati_workflow_helper("kiribati_tb_initial_proportions")
tb_cumulative_flows <- kiribati_workflow_helper("tb_cumulative_flows")
summarise_tb_burden <- kiribati_workflow_helper("summarise_tb_burden")
kiribati_tb_target_map <- kiribati_workflow_helper("kiribati_tb_target_map")
print_kiribati_tb_summary <- kiribati_workflow_helper("print_kiribati_tb_summary")

test_that("kiribati_tb_parameters preserves expected age-specific values", {
  ages <- AgeStructure(
    age_groups = c("0-4", "15-24", "65+"),
    lower_bounds = c(0, 15, 65),
    upper_bounds = c(4, 24, Inf)
  )

  parameters <- kiribati_tb_parameters(ages)

  expect_equal(parameters$beta, 0.13)
  expect_equal(parameters$recent_to_remote, 0.35)
  expect_equal(unname(parameters$fast_progression), c(0.08, 0.045, 0.06))
  expect_equal(unname(parameters$susceptibility), c(0.75, 1.0, 1.15))
  expect_identical(names(parameters$infectiousness), ages$age_groups)
})

test_that("kiribati_tb_model preserves compartments and transitions", {
  ages <- wpp_age_structure_1year(max_age = 2)
  parameters <- kiribati_tb_parameters(ages)

  model <- kiribati_tb_model(parameters)

  expect_identical(model$compartments, c("S", "Lr", "Ld", "I", "T", "R"))
  expect_identical(model$birth_compartment, "S")
  expect_identical(model$migration_compartment, "S")
  expect_identical(model$infectious_compartments, "I")
  expect_identical(model$infection_transitions$from, "S")
  expect_identical(model$infection_transitions$to, "Lr")
})

test_that("kiribati_tb_initial_proportions uses S as residual input", {
  ages <- AgeStructure(
    age_groups = c("0-4", "15-24", "65+"),
    lower_bounds = c(0, 15, 65),
    upper_bounds = c(4, 24, Inf)
  )

  proportions <- kiribati_tb_initial_proportions(ages)

  expect_identical(names(proportions), c("Lr", "Ld", "I", "T", "R"))
  expect_equal(unname(proportions$Ld), c(0.01, 0.14, 0.36))
  expect_equal(unname(proportions$T), 0.7 * unname(proportions$I))
  expect_true(all(Reduce(`+`, proportions) < 1))
})

test_that("tb_cumulative_flows returns the expected flow names", {
  flows <- tb_cumulative_flows()

  expect_identical(
    names(flows),
    c(
      "infections",
      "progression_to_active_tb",
      "treatment_initiation",
      "treatment_completion",
      "relapse_recurrent_tb"
    )
  )
  expect_identical(flows$progression_to_active_tb$from, c("Lr", "Ld"))
  expect_identical(flows$progression_to_active_tb$to, c("I", "I"))
})

test_that("summarise_tb_burden calculates TB burden quantities", {
  trajectory <- data.frame(
    time = rep(c(0, 1), each = 12),
    compartment = rep(c("S", "Lr", "Ld", "I", "T", "R"), each = 2, times = 2),
    age_group = rep(c("0-4", "5-9"), times = 12),
    value = c(
      80, 180, 5, 10, 10, 20, 1, 2, 2, 3, 2, 5,
      70, 170, 7, 12, 11, 21, 2, 4, 3, 4, 4, 8
    ),
    stringsAsFactors = FALSE
  )

  burden <- summarise_tb_burden(trajectory)

  expect_identical(
    names(burden),
    c(
      "time",
      "active_tb_prevalence_count",
      "on_treatment_count",
      "latent_tb_infection_count",
      "population",
      "active_tb_prevalence_per_100k"
    )
  )
  expect_equal(burden$active_tb_prevalence_count, c(3, 6))
  expect_equal(burden$on_treatment_count, c(5, 7))
  expect_equal(burden$latent_tb_infection_count, c(45, 51))
  expect_equal(burden$population, c(320, 316))
})

test_that("kiribati_tb_target_map and print helper keep public labels", {
  target_map <- kiribati_tb_target_map()
  output <- capture.output(print_kiribati_tb_summary(
    country = "Kiribati",
    simulation_years = 2025:2026,
    output_timestep_years = 1 / 4,
    simulation_method = "deSolve",
    contact_matrix_source = data.frame(source_label = "POLYMOD proxy"),
    population_summary = data.frame(time = 2025:2026, value = c(1, 2)),
    tb_burden_summary = data.frame(
      time = 2025:2026,
      active_tb_prevalence_count = c(1, 2),
      on_treatment_count = c(1, 2),
      latent_tb_infection_count = c(1, 2),
      population = c(100, 200),
      active_tb_prevalence_per_100k = c(1000, 1000)
    ),
    age_group_summary_broad = data.frame(time = 2026, age_group = "0-4", value = 1),
    cumulative_flow_wide = data.frame(time = 2025:2026, cumulative_infections = c(1, 2)),
    public_data_target_summary = target_map
  ))

  expect_true(any(grepl("not calibrated; not a policy estimate", output, fixed = TRUE)))
  expect_true(any(grepl("Output timestep: quarterly; solver: deSolve", output, fixed = TRUE)))
  expect_true(any(grepl("WHO estimated TB incidence", target_map$target)))
})

test_that("demographic_process_from_wpp builds WPP-backed inputs when WPP is available", {
  testthat::skip_if_not_installed("wpp2024")

  ages <- wpp_age_structure_1year(max_age = 95)
  demography <- demographic_process_from_wpp(
    country = "Kiribati",
    years = 2025:2026,
    age_structure = ages,
    migration = TRUE,
    fertility_exposure_fraction = 0.5
  )

  expect_s3_class(demography, "agepi_wpp_demographic_process_inputs")
  expect_s3_class(demography$population, "agepi_demography")
  expect_s3_class(demography$demographic_process, "agepi_demographic_process")
  expect_equal(demography$fertility_schedule$times, 2025:2026)
  expect_equal(demography$mortality_schedule$times, 2025:2026)
  expect_equal(demography$migration_schedule$times, 2025:2026)
})
