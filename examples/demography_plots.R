if ("package:agepi" %in% search()) {
  # Already loaded by devtools::load_all() or library(agepi).
} else if (dir.exists("R") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
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
  stop("Install ggplot2 to run this plotting example.", call. = FALSE)
}

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

print(plots$population_pyramid)
print(plots$population_projection)
print(plots$age_structure)
