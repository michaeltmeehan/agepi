scenario_test_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

scenario_test_state <- function(S = c(90, 180), I = c(10, 20), R = c(0, 0)) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

scenario_test_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

scenario_test_trajectory <- function() {
  simulate_deterministic(
    initial_state = scenario_test_state(),
    times = c(0, 0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = scenario_test_ages(),
    contact_matrix = scenario_test_contacts(),
    method = "euler"
  )
}

test_that("initial_state_from_simulation works on a plain trajectory data frame", {
  trajectory <- scenario_test_trajectory()

  state <- initial_state_from_simulation(trajectory, time = 0.1)
  expected <- trajectory[trajectory$time == 0.1, c("compartment", "age_group", "value")]
  row.names(expected) <- NULL

  expect_identical(names(state), c("compartment", "age_group", "value"))
  expect_equal(state, expected)
})

test_that("initial_state_from_simulation works on a simulation result list", {
  trajectory <- scenario_test_trajectory()
  result <- list(
    trajectory = trajectory,
    cumulative = data.frame(
      time = c(0, 0.1),
      cumulative_name = "ignored",
      value = c(0, 1)
    )
  )

  state <- initial_state_from_simulation(result, time = 0.2)
  expected <- trajectory[trajectory$time == 0.2, c("compartment", "age_group", "value")]
  row.names(expected) <- NULL

  expect_equal(state, expected)
})

test_that("extracted state can be passed back into simulate_deterministic", {
  trajectory <- scenario_test_trajectory()
  state <- initial_state_from_simulation(trajectory, time = 0.1)

  output <- simulate_deterministic(
    initial_state = state,
    times = c(0.1, 0.2),
    model = SIRModel(gamma = 0.2),
    age_structure = scenario_test_ages(),
    contact_matrix = scenario_test_contacts(),
    method = "euler"
  )

  expect_equal(
    output[output$time == 0.1, c("compartment", "age_group", "value")],
    state,
    ignore_attr = TRUE
  )
})

test_that("initial_state_from_simulation errors clearly when time is absent", {
  expect_error(
    initial_state_from_simulation(scenario_test_trajectory(), time = 0.15),
    "does not contain requested time"
  )
})

test_that("initial_state_from_simulation errors clearly when required columns are missing", {
  trajectory <- scenario_test_trajectory()
  trajectory$value <- NULL

  expect_error(
    initial_state_from_simulation(trajectory, time = 0.1),
    "missing required column"
  )
})

test_that("initial_state_from_simulation rejects ambiguous rows", {
  trajectory <- scenario_test_trajectory()
  duplicate <- trajectory[trajectory$time == 0.1, ][1, ]
  trajectory <- rbind(trajectory, duplicate)

  expect_error(
    initial_state_from_simulation(trajectory, time = 0.1),
    "multiple rows"
  )
})

test_that("Scenario validates name, overrides, modifier, and metadata", {
  modifier <- function(args, scenario) args

  scenario <- Scenario(
    "baseline",
    description = "Current assumptions",
    overrides = list(beta = 0.1),
    modifier = modifier,
    metadata = list(source = "test")
  )

  expect_s3_class(scenario, "agepi_scenario")
  expect_identical(scenario$name, "baseline")
  expect_identical(scenario$overrides$beta, 0.1)
  expect_identical(scenario$modifier, modifier)
  expect_identical(scenario$metadata$source, "test")

  expect_error(Scenario(""), "name")
  expect_error(Scenario("bad", overrides = list(0.1)), "overrides")
  expect_error(Scenario("bad", modifier = "not a function"), "modifier")
  expect_error(Scenario("bad", metadata = list("test")), "metadata")
})

test_that("ScenarioSet accepts multiple scenarios", {
  baseline <- Scenario("baseline")
  intervention <- Scenario("intervention", overrides = list(beta = 0.05))

  scenario_set <- ScenarioSet(
    baseline,
    intervention,
    baseline_name = "baseline",
    branch_time = 10,
    cumulative_policy = "continue"
  )

  expect_s3_class(scenario_set, "agepi_scenario_set")
  expect_identical(names(scenario_set$scenarios), c("baseline", "intervention"))
  expect_identical(scenario_set$baseline_name, "baseline")
  expect_identical(scenario_set$branch_time, 10)
  expect_identical(scenario_set$cumulative_policy, "continue")
})

test_that("ScenarioSet rejects duplicate names", {
  expect_error(
    ScenarioSet(Scenario("baseline"), Scenario("baseline")),
    "unique"
  )
})

test_that("ScenarioSet validates baseline_name", {
  expect_error(
    ScenarioSet(Scenario("baseline"), baseline_name = "missing"),
    "baseline_name"
  )
})

test_that("ScenarioSet validates cumulative_policy", {
  expect_error(
    ScenarioSet(Scenario("baseline"), cumulative_policy = "restart"),
    "cumulative_policy"
  )
})
