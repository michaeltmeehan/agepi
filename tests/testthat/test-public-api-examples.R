load_agepi_for_public_api_tests <- function() {
  if ("package:agepi" %in% search()) {
    return(invisible(TRUE))
  }

  package_root <- normalizePath(test_path("../.."), mustWork = FALSE)
  r_dir <- file.path(package_root, "R")
  if (dir.exists(r_dir) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(package_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (dir.exists(r_dir)) {
    invisible(lapply(list.files(r_dir, pattern = "[.]R$", full.names = TRUE), source))
    return(invisible(TRUE))
  }
  if (requireNamespace("agepi", quietly = TRUE)) {
    library(agepi)
    return(invisible(TRUE))
  }

  stop(
    "Package agepi is not installed. Run this test from the package root ",
    "or install agepi first.",
    call. = FALSE
  )
}

source_example_command <- paste(
  "example_file <- normalizePath(commandArgs(TRUE)[1], mustWork = TRUE)",
  "package_root <- normalizePath(file.path(dirname(example_file), '..'), mustWork = FALSE)",
  "if (dir.exists(file.path(package_root, 'R')) && requireNamespace('pkgload', quietly = TRUE)) {",
  "  pkgload::load_all(package_root, quiet = TRUE)",
  "} else if (dir.exists(file.path(package_root, 'R'))) {",
  "  invisible(lapply(list.files(file.path(package_root, 'R'), pattern = '[.]R$', full.names = TRUE), source))",
  "} else if (requireNamespace('agepi', quietly = TRUE)) {",
  "  library(agepi)",
  "}",
  "source(example_file)",
  sep = "; "
)

load_agepi_for_public_api_tests()

test_that("documented public functions are exported and implemented", {
  namespace_file <- test_path("../../NAMESPACE")
  if (file.exists(namespace_file)) {
    namespace_lines <- readLines(namespace_file, warn = FALSE)
    exported <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace_lines, value = TRUE))
  } else {
    exported <- getNamespaceExports("agepi")
  }

  rd_files <- list.files(test_path("../../man"), pattern = "[.]Rd$", full.names = TRUE)
  if (length(rd_files) > 0) {
    rd_names <- vapply(rd_files, function(path) {
      name_line <- grep("^\\\\name\\{", readLines(path, warn = FALSE), value = TRUE)[1]
      sub("^\\\\name\\{(.*)\\}$", "\\1", name_line)
    }, character(1))
  } else {
    rd_names <- names(tools::Rd_db("agepi"))
    rd_names <- sub("[.]Rd$", "", rd_names)
  }

  missing_exports <- setdiff(rd_names, exported)
  missing_implementations <- exported[!vapply(exported, exists, logical(1), mode = "function")]

  expect_equal(missing_exports, character())
  expect_equal(missing_implementations, character())
})

test_that("reported public WPP and schedule helpers are available", {
  public_helpers <- c(
    "standardise_wpp_mortality",
    "FertilitySchedule",
    "wpp_age_structure_5year",
    "population_from_wpp"
  )

  expect_true(all(vapply(public_helpers, exists, logical(1), mode = "function")))
})

test_that("repository examples run after loading the package like a user", {
  source_examples_dir <- test_path("../../examples")
  installed_examples_dir <- system.file("examples", package = "agepi")
  examples_dir <- if (dir.exists(source_examples_dir)) {
    source_examples_dir
  } else if (nzchar(installed_examples_dir) && dir.exists(installed_examples_dir)) {
    installed_examples_dir
  } else {
    NA_character_
  }

  if (is.na(examples_dir)) {
    expect_false(dir.exists(source_examples_dir))
    expect_identical(installed_examples_dir, "")
    return(invisible(NULL))
  }

  core_example_names <- c(
    "demography_plots.R",
    "deterministic_cumulative_flows.R",
    "generic_msir.R",
    "generic_seir.R",
    "generic_sir.R",
    "mock_demographic_workflow.R",
    "mock_seir_demography.R",
    "mock_sir_deterministic.R",
    "observed_case_age_groups.R",
    "stochastic_cumulative_flows.R",
    "stochastic_seir.R",
    "stochastic_sir.R",
    "tb_age_structured_demography.R"
  )
  example_files <- file.path(examples_dir, core_example_names)
  example_files <- example_files[file.exists(example_files)]
  expect_gt(length(example_files), 0)

  for (example_file in example_files) {
    result <- system2(
      "Rscript",
      c(
        "-e",
        shQuote(source_example_command),
        shQuote(example_file)
      ),
      stdout = TRUE,
      stderr = TRUE
    )

    status <- attr(result, "status")
    if (is.null(status)) {
      status <- 0
    }

    expect_equal(
      status,
      0,
      info = paste(c("Example failed:", example_file, result), collapse = "\n")
    )
  }
})

test_that("optional integration examples are guarded and non-blocking", {
  source_examples_dir <- test_path("../../examples")
  installed_examples_dir <- system.file("examples", package = "agepi")
  examples_dir <- if (dir.exists(source_examples_dir)) {
    source_examples_dir
  } else if (nzchar(installed_examples_dir) && dir.exists(installed_examples_dir)) {
    installed_examples_dir
  } else {
    NA_character_
  }

  skip_if(is.na(examples_dir))

  optional_example_names <- c(
    "epiparameter_seir.R",
    "wpp_demography_validation.R"
  )
  optional_example_files <- file.path(examples_dir, optional_example_names)
  optional_example_files <- optional_example_files[file.exists(optional_example_files)]

  for (example_file in optional_example_files) {
    result <- system2(
      "Rscript",
      c(
        "-e",
        shQuote(source_example_command),
        shQuote(example_file)
      ),
      stdout = TRUE,
      stderr = TRUE
    )

    status <- attr(result, "status")
    if (is.null(status)) {
      status <- 0
    }

    expect_equal(
      status,
      0,
      info = paste(c("Optional example failed:", example_file, result), collapse = "\n")
    )
  }
})

test_that("toy age-structured TB demography example exposes expected outputs", {
  example_file <- test_path("../../examples/tb_age_structured_demography.R")
  skip_if_not(file.exists(example_file))

  env <- new.env(parent = globalenv())
  invisible(capture.output(source(example_file, local = env)))

  expect_named(env$tb_output, c("trajectory", "cumulative"))
  expect_s3_class(env$trajectory, "data.frame")
  expect_s3_class(env$cumulative_flows, "data.frame")
  expect_s3_class(env$compartment_summary, "data.frame")
  expect_s3_class(env$age_group_summary, "data.frame")
  expect_s3_class(env$population_summary, "data.frame")
  expect_s3_class(env$disease_onset_total, "data.frame")

  expect_setequal(
    unique(env$cumulative_flows$cumulative_name),
    c(
      "infections",
      "disease_onset",
      "treatment_initiation",
      "treatment_completion",
      "relapse"
    )
  )
  disease_onset_rows <- env$cumulative_flows[
    env$cumulative_flows$cumulative_name == "disease_onset",
    ,
    drop = FALSE
  ]
  expect_true(any(grepl("Lr", disease_onset_rows$from, fixed = TRUE)))
  expect_true(any(grepl("Ld", disease_onset_rows$from, fixed = TRUE)))
  expect_true(any(grepl("I", disease_onset_rows$to, fixed = TRUE)))
  expect_true(all(env$trajectory$value >= 0))
  expect_true(all(diff(env$disease_onset_total$cumulative_disease_onset) >= -1e-8))
})
