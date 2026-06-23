test_that("Kiribati TB realistic demography example runs in a clean user process", {
  skip_if(
    any(!vapply(c("wpp2024", "socialmixr", "deSolve"), requireNamespace, logical(1), quietly = TRUE)),
    "optional package wpp2024, socialmixr, or deSolve is not available"
  )

  example_file <- normalizePath(
    test_path("../../examples/kiribati_tb_realistic_demography.R"),
    mustWork = TRUE
  )

  check_script <- paste(
    "example_file <- normalizePath(commandArgs(TRUE)[1], mustWork = TRUE)",
    "package_root <- normalizePath(file.path(dirname(example_file), '..'), mustWork = FALSE)",
    "oldwd <- getwd()",
    "setwd(package_root)",
    "on.exit(setwd(oldwd), add = TRUE)",
    "invisible(lapply(list.files(file.path(package_root, 'R'), pattern = '[.]R$', full.names = TRUE), source))",
    "env <- new.env(parent = globalenv())",
    "invisible(capture.output(source(example_file, local = env)))",
    "example_lines <- readLines(example_file, warn = FALSE)",
    "stopifnot(!any(grepl('stats::approx|approxfun|time_policy = \"linear\"|method = \"euler\"|ageing_policy = \"annual_cohort\"', example_lines)))",
    "stopifnot(identical(names(env$tb_output), c('trajectory', 'cumulative')))",
    "stopifnot(is.data.frame(env$trajectory))",
    "stopifnot(is.data.frame(env$cumulative_flows))",
    "stopifnot(is.data.frame(env$contact_matrix_source))",
    "stopifnot(is.data.frame(env$population_summary))",
    "stopifnot(is.data.frame(env$age_group_summary_broad))",
    "stopifnot(is.data.frame(env$tb_burden_summary))",
    "stopifnot(is.data.frame(env$public_data_target_summary))",
    "stopifnot(all(is.finite(env$trajectory$value)))",
    "stopifnot(all(env$trajectory$value >= -1e-8))",
    "stopifnot(all(is.finite(env$cumulative_flows$value)))",
    "stopifnot(all(env$cumulative_flows$value >= -1e-8))",
    "stopifnot(identical(env$simulation_method, 'deSolve'))",
    "stopifnot(identical(env$demographic_time_policy, 'step'))",
    "stopifnot(abs(env$output_timestep_years - 1 / 4) < 1e-12)",
    "output_times <- sort(unique(env$trajectory$time))",
    "stopifnot(length(output_times) == 21)",
    "stopifnot(max(abs(diff(output_times) - 1 / 4)) < 1e-8)",
    "stopifnot(all(env$simulation_years %in% output_times))",
    "stopifnot(!exists('fertility_schedule', envir = env, inherits = FALSE))",
    "stopifnot(!exists('mortality_schedule', envir = env, inherits = FALSE))",
    "stopifnot(!exists('migration_schedule', envir = env, inherits = FALSE))",
    "stopifnot(!exists('tfr_input', envir = env, inherits = FALSE))",
    "stopifnot(!exists('tfr_years', envir = env, inherits = FALSE))",
    "stopifnot(!exists('process_years', envir = env, inherits = FALSE))",
    "stopifnot(all(env$demography$schedules$fertility$times == env$simulation_years))",
    "stopifnot(all(env$demography$schedules$mortality$times == env$simulation_years))",
    "stopifnot(all(env$demography$schedules$migration$times == env$simulation_years))",
    "stopifnot(setequal(unique(env$cumulative_flows$cumulative_name), c('infections', 'progression_to_active_tb', 'treatment_initiation', 'treatment_completion', 'relapse_recurrent_tb')))",
    "cumulative_totals <- stats::aggregate(value ~ time + cumulative_name, env$cumulative_flows, sum)",
    "for (flow_name in unique(cumulative_totals$cumulative_name)) { rows <- cumulative_totals[cumulative_totals$cumulative_name == flow_name, ]; stopifnot(all(diff(rows$value) >= -1e-8)) }",
    "stopifnot(min(env$population_summary$time) == 2025)",
    "stopifnot(max(env$population_summary$time) == 2030)",
    "stopifnot(all(env$population_summary$value > 0))",
    "stopifnot(env$population_summary$value[env$population_summary$time == 2030] > env$population_summary$value[env$population_summary$time == 2025])",
    "stopifnot(is.matrix(env$contact_matrix))",
    "stopifnot(identical(dim(env$contact_matrix), c(env$age_structure$n_age_groups, env$age_structure$n_age_groups)))",
    "stopifnot(identical(rownames(env$contact_matrix), env$age_structure$age_groups))",
    "stopifnot(identical(colnames(env$contact_matrix), env$age_structure$age_groups))",
    "stopifnot(all(is.finite(env$contact_matrix)))",
    "stopifnot(all(env$contact_matrix >= 0))",
    "stopifnot(env$contact_matrix_source$source %in% c('prem', 'polymod_uk'))",
    "if (identical(env$contact_matrix_source$source, 'prem')) stopifnot(identical(env$contact_matrix_source$country, 'Kiribati'))",
    "if (identical(env$contact_matrix_source$source, 'polymod_uk')) stopifnot(grepl('POLYMOD|socialmixr', env$contact_matrix_source$source_label))",
    "stopifnot(grepl('constant contacts within each source band', env$contact_matrix_source$adaptation_note))",
    "stopifnot(any(grepl('WHO', env$public_data_target_summary$target)))",
    "stopifnot(any(grepl('not modelled', env$public_data_target_summary$model_quantity)))",
    "cat('KIRIBATI_TB_EXAMPLE_CHECKS_OK\\n')",
    "quit(status = 0)",
    sep = "; "
  )

  check_file <- tempfile("kiribati-tb-example-check-", fileext = ".R")
  writeLines(strsplit(check_script, "; ", fixed = TRUE)[[1]], check_file)
  on.exit(unlink(check_file), add = TRUE)

  result <- suppressWarnings(
    system2("Rscript", c(check_file, example_file), stdout = TRUE, stderr = TRUE)
  )
  expect_true(any(grepl("KIRIBATI_TB_EXAMPLE_CHECKS_OK", result, fixed = TRUE)))
})
