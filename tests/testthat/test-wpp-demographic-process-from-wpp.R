test_that("demographic_process_from_wpp returns a reusable structured object when WPP is available", {
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
  expect_identical(demography$age_structure, ages)
  expect_s3_class(demography$population, "agepi_demography")
  expect_s3_class(demography$demographic_process, "agepi_demographic_process")
  expect_named(demography$schedules, c("fertility", "mortality", "migration"))
  expect_s3_class(demography$schedules$fertility, "agepi_fertility_schedule")
  expect_s3_class(demography$schedules$mortality, "agepi_mortality_schedule")
  expect_s3_class(demography$schedules$migration, "agepi_migration_schedule")
  expect_named(demography$inputs, c("population", "fertility_weights", "tfr", "mortality", "migration"))

  expect_false("fertility_schedule" %in% names(demography))
  expect_false("mortality_schedule" %in% names(demography))
  expect_false("migration_schedule" %in% names(demography))
  expect_equal(demography$schedules$fertility$times, 2025:2026)
  expect_equal(demography$schedules$mortality$times, 2025:2026)
  expect_equal(demography$schedules$migration$times, 2025:2026)
})
