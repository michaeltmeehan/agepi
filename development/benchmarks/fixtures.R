benchmark_age_structure <- function() {
  agepi::wpp_age_structure_5year(max_age = 40)
}

benchmark_contact_matrix <- function(age_structure) {
  n <- age_structure$n_age_groups
  index <- seq_len(n)
  matrix(
    outer(index, index, function(recipient, source) {
      0.15 + exp(-abs(recipient - source) / 2) + (recipient / n) * 0.03 + (source / n) * 0.02
    }),
    nrow = n,
    ncol = n,
    dimnames = list(age_structure$age_groups, age_structure$age_groups)
  )
}

benchmark_population <- function(age_structure) {
  round(seq(4000, 2500, length.out = age_structure$n_age_groups))
}

benchmark_sir_state <- function(age_structure) {
  total <- benchmark_population(age_structure)
  infectious <- pmax(2, round(total * seq(0.012, 0.006, length.out = length(total))))
  recovered <- round(total * 0.025)
  susceptible <- total - infectious - recovered

  data.frame(
    compartment = rep(c("S", "I", "R"), each = length(total)),
    age_group = rep(age_structure$age_groups, times = 3),
    value = c(susceptible, infectious, recovered),
    stringsAsFactors = FALSE
  )
}

benchmark_seir_state <- function(age_structure) {
  total <- benchmark_population(age_structure)
  exposed <- round(total * seq(0.012, 0.008, length.out = length(total)))
  infectious <- round(total * seq(0.010, 0.005, length.out = length(total)))
  recovered <- round(total * 0.03)
  susceptible <- total - exposed - infectious - recovered

  data.frame(
    compartment = rep(c("S", "E", "I", "R"), each = length(total)),
    age_group = rep(age_structure$age_groups, times = 4),
    value = c(susceptible, exposed, infectious, recovered),
    stringsAsFactors = FALSE
  )
}

benchmark_generic_state <- function(age_structure) {
  total <- benchmark_population(age_structure)
  compartments <- c("S", "E", "IP", "IC", "IS", "R")
  shares <- c(0.86, 0.04, 0.03, 0.02, 0.03, 0.02)
  values <- lapply(shares, function(share) round(total * share))
  values[[1]] <- total - Reduce(`+`, values[-1])

  data.frame(
    compartment = rep(compartments, each = length(total)),
    age_group = rep(age_structure$age_groups, times = length(compartments)),
    value = unlist(values, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

benchmark_generic_model <- function(age_structure) {
  ages <- age_structure$age_groups
  rate_growth <- setNames(seq(0.08, 0.18, length.out = length(ages)), ages)
  rate_progression <- setNames(seq(0.05, 0.14, length.out = length(ages)), ages)
  rate_clinical <- setNames(seq(0.18, 0.28, length.out = length(ages)), ages)
  rate_subclinical <- setNames(seq(0.10, 0.20, length.out = length(ages)), ages)

  transitions <- data.frame(
    from = c("E", "E", "IP", "IC", "IS"),
    to = c("IP", "IS", "IC", "R", "R"),
    stringsAsFactors = FALSE
  )
  transitions$rate <- I(list(
    rate_progression,
    rate_growth,
    rate_clinical,
    0.08,
    rate_subclinical
  ))

  agepi::CompartmentModel(
    compartments = c("S", "E", "IP", "IC", "IS", "R"),
    infection_transitions = data.frame(
      from = "S",
      to = "E",
      susceptibility = I(list(setNames(seq(0.9, 1.1, length.out = length(ages)), ages))),
      stringsAsFactors = FALSE
    ),
    transitions = transitions,
    infectious_compartments = c("IP", "IC", "IS"),
    infectiousness_weights = list(
      IS = setNames(seq(0.4, 0.6, length.out = length(ages)), ages),
      IP = setNames(seq(1.0, 1.1, length.out = length(ages)), ages),
      IC = 0.8
    ),
    beta = 0.18
  )
}

benchmark_demographic_process <- function(age_structure) {
  ages <- age_structure$age_groups
  times <- c(2000, 2005, 2010)
  fertility_ages <- ages[age_structure$lower_bounds >= 15 & age_structure$lower_bounds <= 45]
  fertility_rates <- seq(0.02, 0.08, length.out = length(fertility_ages))
  mortality_rates <- seq(0.002, 0.02, length.out = length(ages))
  migration_rates <- seq(-4, 4, length.out = length(ages))

  fertility <- agepi::FertilitySchedule(
    data.frame(
      time = rep(times, each = length(fertility_ages)),
      age_group = rep(fertility_ages, times = length(times)),
      fertility_rate = rep(fertility_rates, times = length(times)),
      stringsAsFactors = FALSE
    ),
    age_structure
  )
  mortality <- agepi::MortalitySchedule(
    data.frame(
      time = rep(times, each = length(ages)),
      age_group = rep(ages, times = length(times)),
      mortality_rate = rep(mortality_rates, times = length(times)),
      stringsAsFactors = FALSE
    ),
    age_structure
  )
  migration <- agepi::MigrationSchedule(
    data.frame(
      time = rep(times, each = length(ages)),
      age_group = rep(ages, times = length(times)),
      migration_count = rep(migration_rates, times = length(times)),
      stringsAsFactors = FALSE
    ),
    age_structure
  )

  agepi::DemographicProcess(
    age_structure = age_structure,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = "migration"
  )
}

benchmark_times <- function(start = 0, end = 3, step = 1) {
  seq(start, end, by = step)
}

benchmark_demography_times <- function() {
  seq(2000, 2003, by = 1)
}

benchmark_deterministic_method <- function() {
  "euler"
}
