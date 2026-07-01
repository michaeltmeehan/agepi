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

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message(
    "Skipping plotting example: optional package ggplot2 is not installed. ",
    "Install ggplot2 separately to run this example."
  )
  run_example <- FALSE
} else {
  run_example <- TRUE
}

ggplot2_built <- suppressWarnings(packageDescription("ggplot2", fields = "Built"))
ggplot2_built_r <- sub("^R ([^;]+).*$", "\\1", ggplot2_built)
if (
  run_example &&
    !interactive() &&
    !is.na(ggplot2_built_r) &&
    numeric_version(ggplot2_built_r) > getRversion()
) {
  message(
    "Plotting example skipped for this non-interactive run because ggplot2 ",
    "was built under R ", ggplot2_built_r, " and this session is R ",
    getRversion(), "."
  )
  run_example <- FALSE
}

if (run_example) {
  # Purpose: demonstrate lightweight exploratory plots for synthetic demographic
  # input tables before using them in age-structured epidemic models.
  # Expected output: three ggplot objects in a named list.

  population <- data.frame(
    time = rep(c(2020, 2025, 2030), each = 4),
    age_group = rep(c("0-4", "5-14", "15-64", "65+"), times = 3),
    population = c(
      1000, 1800, 5200, 600,
      980, 1750, 5400, 720,
      960, 1700, 5550, 860
    ),
    stringsAsFactors = FALSE
  )

  plots <- plot_demography(
    population,
    year = 2020,
    compare_year = 2030,
    age_groups = c("0-4", "15-64", "65+")
  )

  render_plots <- interactive() ||
    identical(tolower(Sys.getenv("AGEPI_DEMOGRAPHY_PLOTS_RENDER")), "true")

  if (render_plots) {
    print(plots$population_pyramid)
    print(plots$population_projection)
    print(plots$age_structure)
  } else {
    message(
      "Plot rendering skipped for this non-interactive run. ",
      "Set AGEPI_DEMOGRAPHY_PLOTS_RENDER=true to render plots intentionally."
    )
  }
}

invisible(NULL)
