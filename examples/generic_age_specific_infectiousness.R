library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

contact_matrix <- matrix(
  c(1.6, 0.8,
    0.8, 1.4),
  nrow = ages$n_age_groups,
  byrow = TRUE
)

initial_state <- data.frame(
  compartment = rep(c("S", "IP", "IC", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 4),
  value = c(
    c(900, 1800),
    c(20, 30),
    c(5, 10),
    c(0, 0)
  ),
  stringsAsFactors = FALSE
)

scalar_model <- CompartmentModel(
  compartments = c("S", "IP", "IC", "R"),
  infection_transitions = data.frame(from = "S", to = "IP", stringsAsFactors = FALSE),
  transitions = data.frame(
    from = c("IP", "IC"),
    to = c("IC", "R"),
    rate = c(0.15, 0.05),
    stringsAsFactors = FALSE
  ),
  infectious_compartments = c("IP", "IC"),
  infectiousness_weights = c(IP = 0.25, IC = 1)
)

age_specific_model <- CompartmentModel(
  compartments = c("S", "IP", "IC", "R"),
  infection_transitions = data.frame(from = "S", to = "IP", stringsAsFactors = FALSE),
  transitions = data.frame(
    from = c("IP", "IC"),
    to = c("IC", "R"),
    rate = c(0.15, 0.05),
    stringsAsFactors = FALSE
  ),
  infectious_compartments = c("IP", "IC"),
  infectiousness_weights = list(
    IP = c("0-4" = 0.25, "5-9" = 0.25),
    IC = c("0-4" = 1, "5-9" = 1)
  )
)

scalar_rates <- transition_rates(
  state = initial_state,
  model = scalar_model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.02
)

age_specific_rates <- transition_rates(
  state = initial_state,
  model = age_specific_model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.02
)

stopifnot(all.equal(scalar_rates$rate, age_specific_rates$rate))

child_zero_model <- CompartmentModel(
  compartments = c("S", "IP", "IC", "R"),
  infection_transitions = data.frame(from = "S", to = "IP", stringsAsFactors = FALSE),
  transitions = data.frame(
    from = c("IP", "IC"),
    to = c("IC", "R"),
    rate = c(0.15, 0.05),
    stringsAsFactors = FALSE
  ),
  infectious_compartments = c("IP", "IC"),
  infectiousness_weights = list(
    IP = c("0-4" = 0, "5-9" = 0),
    IC = c("0-4" = 1, "5-9" = 1)
  )
)

child_zero_rates <- transition_rates(
  state = initial_state,
  model = child_zero_model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.02
)

effective_infectious_by_age <- c(
  "0-4" = 0 * 20 + 1 * 5,
  "5-9" = 0 * 30 + 1 * 10
)

print(effective_infectious_by_age)
print(child_zero_rates[child_zero_rates$from == "S", c("age_group", "rate")])

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 2, by = 1),
  model = child_zero_model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.02,
  method = "euler"
)

print(compartment_totals(simulation))
