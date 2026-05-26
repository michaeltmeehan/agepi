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

  example_files <- list.files(examples_dir, pattern = "[.]R$", full.names = TRUE)
  expect_gt(length(example_files), 0)

  for (example_file in example_files) {
    result <- system2(
      "Rscript",
      c(
        "-e",
        shQuote("devtools::load_all(quiet = TRUE); source(commandArgs(TRUE)[1])"),
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
