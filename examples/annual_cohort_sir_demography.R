library(agepi)

ages <- wpp_age_structure_1year(max_age = 5)

population <- c(1000, 950, 900, 850, 800, 1200)
infected <- c(2, 3, 4, 4, 3, 5)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 3),
  value = c(population - infected, infected, rep(0, ages$n_age_groups)),
  stringsAsFactors = FALSE
)

fertility <- FertilitySchedule(
  data.frame(
    time = 2020,
    age_group = "3",
    fertility_rate = 0.08,
    stringsAsFactors = FALSE
  ),
  ages
)

mortality <- MortalitySchedule(
  data.frame(
    time = 2020,
    age_group = ages$age_groups,
    mortality_rate = c(0.01, 0.002, 0.001, 0.001, 0.002, 0.03),
    stringsAsFactors = FALSE
  ),
  ages
)

process <- build_demographic_process(
  age_structure = ages,
  fertility_schedule = fertility,
  mortality_schedule = mortality
)

contacts <- matrix(1, nrow = ages$n_age_groups, ncol = ages$n_age_groups)
diag(contacts) <- 3

method <- if (requireNamespace("deSolve", quietly = TRUE)) "deSolve" else "euler"

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = c(2020, 2021),
  model = SIRModel(gamma = 0.25),
  age_structure = ages,
  contact_matrix = contacts,
  beta = 0.03,
  method = method,
  demographic_process = process,
  ageing_policy = "annual_cohort",
  time_policy = "step"
)

head(simulation)
compartment_totals(simulation)
total_population(simulation)
