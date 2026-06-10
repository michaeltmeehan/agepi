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

# Purpose: map observed case records into the same age groups used by an agepi
# aggregate model before deriving initial conditions or calibration targets.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

observed_cases <- data.frame(
  case_id = paste0("case_", 1:6),
  onset_day = c(0, 0, 1, 1, 2, 2),
  age = c(2, 8, 11, 36, 4, 7),
  stringsAsFactors = FALSE
)

agepi_cases <- as_agepi_cases(
  observed_cases,
  age_structure = age_structure,
  age_col = "age"
)

incident_targets <- aggregate(
  case_id ~ onset_day + age_group,
  agepi_cases,
  length
)
names(incident_targets)[names(incident_targets) == "case_id"] <- "cases"

initial_infections <- table(
  factor(
    agepi_cases$age_group[agepi_cases$onset_day == 0],
    levels = age_structure$age_groups
  )
)

print(agepi_cases)
print(incident_targets)
print(as.numeric(initial_infections))
