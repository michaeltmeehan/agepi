validation_ages <- function(n = 3) {
  upper_bounds <- seq_len(n) - 1
  upper_bounds[n] <- Inf
  AgeStructure(
    age_groups = paste0(seq_len(n) - 1, "-", seq_len(n) - 1),
    lower_bounds = seq_len(n) - 1,
    upper_bounds = upper_bounds
  )
}

validation_sir_state <- function(S, I, R, ages) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = ages$n_age_groups),
    age_group = rep(ages$age_groups, times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

validation_generic_sir_model <- function(gamma) {
  CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = gamma, stringsAsFactors = FALSE),
    infectious_compartments = "I"
  )
}

validation_msir_state <- function(M, S, I, R, ages) {
  data.frame(
    compartment = rep(c("M", "S", "I", "R"), each = ages$n_age_groups),
    age_group = rep(ages$age_groups, times = 4),
    value = c(M, S, I, R),
    stringsAsFactors = FALSE
  )
}

validation_msir_model <- function(omega = 0.15, gamma = 0.25) {
  CompartmentModel(
    compartments = c("M", "S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(
      from = c("M", "I"),
      to = c("S", "R"),
      rate = c(omega, gamma),
      stringsAsFactors = FALSE
    ),
    infectious_compartments = "I",
    birth_compartment = "M",
    migration_compartment = "S"
  )
}

aggregate_sir_output <- function(output) {
  aggregated <- aggregate(value ~ time + compartment, output, sum)
  aggregated[order(
    aggregated$time,
    match(aggregated$compartment, c("S", "I", "R"))
  ), ]
}

scalar_sir_euler <- function(S0, I0, R0, beta, gamma, times) {
  state <- matrix(NA_real_, nrow = length(times), ncol = 3)
  colnames(state) <- c("S", "I", "R")
  state[1, ] <- c(S0, I0, R0)
  total <- S0 + I0 + R0

  for (i in seq_len(length(times) - 1)) {
    dt <- times[i + 1] - times[i]
    infection <- beta * state[i, "S"] * state[i, "I"] / total
    recovery <- gamma * state[i, "I"]
    state[i + 1, ] <- state[i, ] + dt * c(-infection, infection - recovery, recovery)
  }

  data.frame(
    time = rep(times, each = 3),
    compartment = rep(c("S", "I", "R"), times = length(times)),
    value = as.numeric(t(state)),
    stringsAsFactors = FALSE
  )
}

scalar_final_size <- function(R0, initial_infected_fraction) {
  susceptible_initial <- 1 - initial_infected_fraction
  root <- stats::uniroot(
    function(s_final) {
      log(s_final / susceptible_initial) +
        R0 * (1 - s_final)
    },
    interval = c(.Machine$double.eps, susceptible_initial)
  )
  root$root
}

vector_final_size <- function(S0, I0, R0, beta, gamma, contact_matrix, susceptibility, infectiousness) {
  population <- S0 + I0 + R0
  objective <- function(logit_s_fraction) {
    s_fraction <- stats::plogis(logit_s_fraction)
    S_final <- S0 * s_fraction
    cumulative_infectious <- (population - S_final - R0) / (gamma * population)
    expected <- S0 * exp(
      -beta * susceptibility * as.numeric(contact_matrix %*% (infectiousness * cumulative_infectious))
    )
    sum((S_final - expected)^2)
  }

  fit <- stats::optim(
    par = rep(0, length(S0)),
    fn = objective,
    method = "BFGS",
    control = list(reltol = 1e-14, maxit = 10000)
  )
  S0 * stats::plogis(fit$par)
}

test_that("age-structured SIR reduces to homogeneous scalar SIR under symmetric mixing", {
  ages <- validation_ages(4)
  population <- c(1000, 800, 1200, 600)
  infected_fraction <- 0.01
  S0 <- population * (1 - infected_fraction)
  I0 <- population * infected_fraction
  R0 <- rep(0, ages$n_age_groups)
  beta <- 0.8
  gamma <- 0.25
  times <- seq(0, 30, by = 0.05)

  output <- simulate_deterministic(
    initial_state = validation_sir_state(S0, I0, R0, ages),
    times = times,
    model = SIRModel(gamma = gamma),
    age_structure = ages,
    contact_matrix = matrix(1 / ages$n_age_groups, ages$n_age_groups, ages$n_age_groups),
    beta = beta,
    susceptibility = rep(1, ages$n_age_groups),
    infectiousness = rep(1, ages$n_age_groups),
    method = "euler"
  )

  aggregated <- aggregate_sir_output(output)
  reference <- scalar_sir_euler(sum(S0), sum(I0), sum(R0), beta, gamma, times)

  expect_equal(aggregated$value, reference$value, tolerance = 1e-8)
})

test_that("homogeneous SIR final size matches the scalar final-size equation", {
  skip_if_not_installed("deSolve")

  for (basic_reproduction_number in c(1.5, 2.5, 4)) {
    ages <- validation_ages(1)
    gamma <- 0.2
    beta <- basic_reproduction_number * gamma
    initial_infected_fraction <- 1e-5
    output <- simulate_deterministic(
      initial_state = validation_sir_state(
        S = 1 - initial_infected_fraction,
        I = initial_infected_fraction,
        R = 0,
        ages = ages
      ),
      times = seq(0, 500, by = 1),
      model = SIRModel(gamma = gamma),
      age_structure = ages,
      contact_matrix = matrix(1, 1, 1),
      beta = beta,
      method = "deSolve"
    )
    final <- output[output$time == max(output$time), ]
    simulated_s_final <- final$value[final$compartment == "S"]
    expected_s_final <- scalar_final_size(
      R0 = basic_reproduction_number,
      initial_infected_fraction = initial_infected_fraction
    )

    expect_lt(final$value[final$compartment == "I"], 1e-7)
    expect_lt(abs(simulated_s_final - expected_s_final), 5e-5)
  }
})

test_that("generic homogeneous SIR final size matches the scalar final-size equation", {
  skip_if_not_installed("deSolve")

  ages <- validation_ages(1)
  gamma <- 0.2
  basic_reproduction_number <- 2.5
  beta <- basic_reproduction_number * gamma
  initial_infected_fraction <- 1e-5

  output <- simulate_deterministic(
    initial_state = validation_sir_state(
      S = 1 - initial_infected_fraction,
      I = initial_infected_fraction,
      R = 0,
      ages = ages
    ),
    times = seq(0, 500, by = 1),
    model = validation_generic_sir_model(gamma = gamma),
    age_structure = ages,
    contact_matrix = matrix(1, 1, 1),
    beta = beta,
    method = "deSolve"
  )
  final <- output[output$time == max(output$time), ]
  simulated_s_final <- final$value[final$compartment == "S"]
  expected_s_final <- scalar_final_size(
    R0 = basic_reproduction_number,
    initial_infected_fraction = initial_infected_fraction
  )

  expect_lt(final$value[final$compartment == "I"], 1e-7)
  expect_lt(abs(simulated_s_final - expected_s_final), 5e-5)
})

test_that("generic one-age-group SIR agrees with scalar SIR reference", {
  ages <- validation_ages(1)
  times <- seq(0, 3, by = 0.1)
  beta <- 0.4
  gamma <- 0.1
  initial_state <- validation_sir_state(S = 99, I = 1, R = 0, ages = ages)

  output <- simulate_deterministic(
    initial_state = initial_state,
    times = times,
    model = validation_generic_sir_model(gamma = gamma),
    age_structure = ages,
    contact_matrix = matrix(1, 1, 1),
    beta = beta,
    method = "euler"
  )
  reference <- scalar_sir_euler(99, 1, 0, beta, gamma, times)

  expect_equal(output[, c("time", "compartment", "value")], reference, tolerance = 1e-10)
})

test_that("age-structured SIR final size matches vector final-size equations", {
  skip_if_not_installed("deSolve")

  ages <- validation_ages(2)
  S0 <- c(999, 1498)
  I0 <- c(1, 2)
  R0 <- c(0, 0)
  gamma <- 0.25
  beta <- 0.6
  susceptibility <- c(0.8, 1.3)
  infectiousness <- c(1.1, 0.7)
  contact_matrix <- matrix(c(
    5, 1,
    2, 4
  ), nrow = 2, byrow = TRUE)

  output <- simulate_deterministic(
    initial_state = validation_sir_state(S0, I0, R0, ages),
    times = seq(0, 500, by = 1),
    model = SIRModel(gamma = gamma),
    age_structure = ages,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    method = "deSolve"
  )
  final <- output[output$time == max(output$time), ]
  simulated_s_final <- final$value[final$compartment == "S"]
  expected_s_final <- vector_final_size(
    S0 = S0,
    I0 = I0,
    R0 = R0,
    beta = beta,
    gamma = gamma,
    contact_matrix = contact_matrix,
    susceptibility = susceptibility,
    infectiousness = infectiousness
  )

  expect_lt(max(final$value[final$compartment == "I"]), 1e-5)
  expect_equal(simulated_s_final, expected_s_final, tolerance = 1e-2)
})

test_that("generic age-structured SIR final size matches vector final-size equations", {
  skip_if_not_installed("deSolve")

  ages <- validation_ages(2)
  S0 <- c(999, 1498)
  I0 <- c(1, 2)
  R0 <- c(0, 0)
  gamma <- 0.25
  beta <- 0.6
  susceptibility <- c(0.8, 1.3)
  infectiousness <- c(1.1, 0.7)
  contact_matrix <- matrix(c(
    5, 1,
    2, 4
  ), nrow = 2, byrow = TRUE)

  output <- simulate_deterministic(
    initial_state = validation_sir_state(S0, I0, R0, ages),
    times = seq(0, 500, by = 1),
    model = validation_generic_sir_model(gamma = gamma),
    age_structure = ages,
    contact_matrix = contact_matrix,
    beta = beta,
    susceptibility = susceptibility,
    infectiousness = infectiousness,
    method = "deSolve"
  )
  final <- output[output$time == max(output$time), ]
  simulated_s_final <- final$value[final$compartment == "S"]
  expected_s_final <- vector_final_size(
    S0 = S0,
    I0 = I0,
    R0 = R0,
    beta = beta,
    gamma = gamma,
    contact_matrix = contact_matrix,
    susceptibility = susceptibility,
    infectiousness = infectiousness
  )

  expect_lt(max(final$value[final$compartment == "I"]), 1e-5)
  expect_equal(simulated_s_final, expected_s_final, tolerance = 1e-2)
})

test_that("generic MSIR preserves total living population without demography", {
  ages <- validation_ages(3)
  initial_state <- validation_msir_state(
    M = c(100, 20, 0),
    S = c(890, 1170, 898),
    I = c(10, 10, 2),
    R = c(0, 0, 0),
    ages = ages
  )

  output <- simulate_deterministic(
    initial_state = initial_state,
    times = seq(0, 4, by = 0.1),
    model = validation_msir_model(),
    age_structure = ages,
    contact_matrix = matrix(1, ages$n_age_groups, ages$n_age_groups),
    beta = 0.05,
    method = "euler"
  )
  totals <- aggregate(value ~ time, output, sum)

  expect_equal(totals$value, rep(totals$value[1], nrow(totals)), tolerance = 1e-8)
})

test_that("generic MSIR has no infections without initial infectious people", {
  ages <- validation_ages(3)
  initial_state <- validation_msir_state(
    M = c(100, 20, 0),
    S = c(900, 1180, 900),
    I = c(0, 0, 0),
    R = c(0, 0, 0),
    ages = ages
  )

  output <- simulate_deterministic(
    initial_state = initial_state,
    times = seq(0, 4, by = 0.1),
    model = validation_msir_model(),
    age_structure = ages,
    contact_matrix = matrix(1, ages$n_age_groups, ages$n_age_groups),
    beta = 0.5,
    method = "euler"
  )
  infectious <- output[output$compartment == "I", ]

  expect_equal(infectious$value, rep(0, nrow(infectious)), tolerance = 1e-12)
})

test_that("force of infection edge cases obey recipient and source conventions", {
  contact_matrix <- matrix(c(
    2, 4,
    3, 5
  ), nrow = 2, byrow = TRUE)

  expect_equal(
    force_of_infection(
      infectious = c(10, 20),
      population = c(100, 100),
      contact_matrix = contact_matrix,
      susceptibility = c(0, 1)
    ),
    c(0, 1.3),
    ignore_attr = TRUE
  )

  expect_equal(
    force_of_infection(
      infectious = c(10, 20),
      population = c(100, 100),
      contact_matrix = contact_matrix,
      infectiousness = c(1, 0)
    ),
    c(0.2, 0.3),
    ignore_attr = TRUE
  )

  expect_equal(
    force_of_infection(
      infectious = c(10, 20),
      population = c(100, 100),
      contact_matrix = diag(c(2, 5))
    ),
    c(0.2, 1),
    ignore_attr = TRUE
  )

  expect_equal(
    force_of_infection(
      infectious = 10,
      population = 200,
      contact_matrix = matrix(3, 1, 1),
      beta = 0.4
    ),
    0.4 * 3 * 10 / 200,
    ignore_attr = TRUE
  )
})

test_that("ageing-only demographic dynamics conserve total population", {
  ages <- validation_ages(3)
  process <- DemographicProcess(age_structure = ages)
  state <- c(100, 80, 60)

  derivative <- demographic_derivative(state, time = 0, process = process)
  output <- simulate_demography(process, initial_state = state, times = seq(0, 3, by = 0.25), method = "euler")
  totals <- aggregate(population ~ time, output, sum)

  expect_equal(sum(derivative), 0)
  expect_equal(totals$population, rep(sum(state), nrow(totals)), tolerance = 1e-10)
})

test_that("demographic derivative total change equals births minus deaths plus migration", {
  ages <- validation_ages(3)
  ageing <- AgeingOperator(ages)
  fertility_rates <- c(0, 0.05, 0.02)
  mortality_rates <- c(0.01, 0.02, 0.03)
  migration_counts <- c(4, -3, 2)
  fertility <- FertilitySchedule(
    data.frame(time = 0, age_group = ages$age_groups, fertility_rate = fertility_rates),
    ages
  )
  mortality <- MortalitySchedule(
    data.frame(time = 0, age_group = ages$age_groups, mortality_rate = mortality_rates),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(time = 0, age_group = ages$age_groups, migration_count = migration_counts),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    ageing_operator = ageing,
    fertility_schedule = fertility,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = "migration"
  )
  state <- c(100, 80, 60)

  derivative <- demographic_derivative(state, time = 0, process = process)

  expect_equal(
    sum(derivative),
    sum(fertility_rates * state) - sum(mortality_rates * state) + sum(migration_counts),
    tolerance = 1e-12
  )
})

test_that("simulate_demography agrees with deterministic SIR demography when transmission is zero", {
  ages <- validation_ages(3)
  times <- seq(0, 2, by = 0.25)
  schedule_times <- times[-length(times)]
  mortality <- MortalitySchedule(
    data.frame(
      time = rep(schedule_times, each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = length(schedule_times)),
      mortality_rate = rep(c(0.01, 0.02, 0.03), times = length(schedule_times))
    ),
    ages
  )
  migration <- MigrationSchedule(
    data.frame(
      time = rep(schedule_times, each = ages$n_age_groups),
      age_group = rep(ages$age_groups, times = length(schedule_times)),
      migration_count = rep(c(1, -2, 3), times = length(schedule_times))
    ),
    ages
  )
  process <- DemographicProcess(
    age_structure = ages,
    mortality_schedule = mortality,
    migration_schedule = migration,
    mode = "migration"
  )
  initial_population <- c(100, 80, 60)

  demographic_output <- simulate_demography(
    process = process,
    initial_state = initial_population,
    times = times,
    method = "euler",
    time_policy = "step"
  )
  deterministic_output <- simulate_deterministic(
    initial_state = validation_sir_state(initial_population, rep(0, ages$n_age_groups), rep(0, ages$n_age_groups), ages),
    times = times,
    model = SIRModel(gamma = 0),
    age_structure = ages,
    contact_matrix = matrix(0, ages$n_age_groups, ages$n_age_groups),
    beta = 0,
    method = "euler",
    demographic_process = process,
    time_policy = "step"
  )
  deterministic_totals <- aggregate(value ~ time + age_group, deterministic_output, sum)
  deterministic_totals <- deterministic_totals[order(
    deterministic_totals$time,
    match(deterministic_totals$age_group, ages$age_groups)
  ), ]

  expect_equal(deterministic_totals$value, demographic_output$population, tolerance = 1e-10)
})

test_that("implied residual migration reproduces one-interval WPP-style population table under Euler", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5+"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, Inf)
  )
  baseline_process <- DemographicProcess(age_structure = ages)
  observed <- Demography(
    data.frame(
      time = c(2020, 2020, 2021, 2021),
      age_group = rep(ages$age_groups, times = 2),
      population = c(100, 200, 90, 230),
      stringsAsFactors = FALSE
    ),
    ages
  )
  residual <- implied_demographic_residual(observed, baseline_process)
  migration <- residual_to_migration_schedule(residual, ages, use = "count")
  residual_process <- DemographicProcess(
    age_structure = ages,
    migration_schedule = migration,
    mode = "migration"
  )

  output <- simulate_demography(
    process = residual_process,
    initial_state = demography_population_vector(observed, 2020),
    times = c(2020, 2021),
    method = "euler"
  )

  expect_equal(
    output$population[output$time == 2021],
    as.numeric(demography_population_vector(observed, 2021)),
    tolerance = 1e-12
  )
})
