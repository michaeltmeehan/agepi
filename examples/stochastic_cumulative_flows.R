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

# Purpose: track stochastic cumulative flows from realised Gillespie events.
# The cumulative table is derived from the event log; it is not part of the
# stochastic state vector and does not create extra propensities.

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

population <- c(200, 250)
initial_exposed <- c(4, 3)
initial_clinical <- c(2, 2)
initial_subclinical <- c(1, 1)
initial_state <- data.frame(
  compartment = rep(c("S", "E", "IP", "IC", "IS", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 6),
  value = c(
    population - initial_exposed - initial_clinical - initial_subclinical,
    initial_exposed,
    initial_clinical,
    c(0, 0),
    initial_subclinical,
    c(0, 0)
  ),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(c(
  4, 2,
  2, 5
), nrow = age_structure$n_age_groups, byrow = TRUE)

transitions <- data.frame(
  from = c("E", "E", "IP", "IC", "IS"),
  to = c("IP", "IS", "IC", "R", "R"),
  stringsAsFactors = FALSE
)
transitions$rate <- I(list(
  c("0-4" = 0.25, "5-9" = 0.35),
  c("0-4" = 0.20, "5-9" = 0.15),
  0.20,
  0.10,
  0.25
))

model <- CompartmentModel(
  compartments = c("S", "E", "IP", "IC", "IS", "R"),
  infection_transitions = data.frame(from = "S", to = "E", stringsAsFactors = FALSE),
  transitions = transitions,
  infectious_compartments = c("IP", "IC", "IS"),
  infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
)

times <- seq(0, 10, by = 1)
output <- simulate_stochastic(
  initial_state = initial_state,
  times = times,
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.06,
  seed = 789,
  return_events = TRUE,
  cumulative_flows = data.frame(
    name = c("exposures", "clinical_onsets", "subclinical_onsets"),
    from = c("S", "E", "E"),
    to = c("E", "IP", "IS"),
    stringsAsFactors = FALSE
  )
)

print(head(output$trajectory, 12))
print(head(output$events, 10))
print(head(output$cumulative, 12))

final_cumulative <- output$cumulative[output$cumulative$time == max(times), ]
print(final_cumulative)

event_counts <- aggregate(
  time ~ transition_id + age_group,
  output$events,
  length
)
names(event_counts)[names(event_counts) == "time"] <- "realised_events"
print(event_counts)
