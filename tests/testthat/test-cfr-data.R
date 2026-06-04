test_that("as_cfr_data converts deterministic cumulative flows to increments", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 3),
    value = c(90, 180, 10, 20, 0, 0),
    stringsAsFactors = FALSE
  )

  output <- simulate_deterministic(
    initial_state = initial_state,
    times = c(0, 1, 2),
    model = SIRModel(gamma = 0.1),
    age_structure = ages,
    contact_matrix = matrix(c(2, 1, 1, 2), nrow = 2, byrow = TRUE),
    beta = 0.01,
    method = "euler",
    cumulative_flows = list(
      infections = list(from = "S", to = "I"),
      removals = list(from = "I", to = "R")
    )
  )

  cfr_data <- as_cfr_data(output, cases = "infections", deaths = "removals")
  cumulative_cases <- aggregate(
    value ~ time,
    output$cumulative[output$cumulative$cumulative_name == "infections", ],
    sum
  )
  cumulative_deaths <- aggregate(
    value ~ time,
    output$cumulative[output$cumulative$cumulative_name == "removals", ],
    sum
  )

  expect_identical(names(cfr_data), c("time", "cases", "deaths"))
  expect_equal(cfr_data$time, c(0, 1, 2))
  expect_equal(cfr_data$cases, c(cumulative_cases$value[1], diff(cumulative_cases$value)))
  expect_equal(cfr_data$deaths, c(cumulative_deaths$value[1], diff(cumulative_deaths$value)))
})

test_that("as_cfr_data can preserve cumulative values and group by age", {
  cumulative <- data.frame(
    time = rep(c(0, 1, 2), each = 4),
    cumulative_name = rep(c("cases", "deaths"), each = 2, times = 3),
    transition_id = rep(c("infection:S->I", "transition:I->R"), each = 2, times = 3),
    from = rep(c("S", "I"), each = 2, times = 3),
    to = rep(c("I", "R"), each = 2, times = 3),
    age_group = rep(c("0-4", "5-9"), times = 6),
    value = c(0, 0, 0, 0, 5, 6, 1, 2, 8, 9, 3, 4),
    stringsAsFactors = FALSE
  )

  cfr_data <- as_cfr_data(
    list(trajectory = data.frame(), cumulative = cumulative),
    cases = list(transition_id = "infection:S->I"),
    deaths = list(cumulative_name = "deaths"),
    cumulative = FALSE,
    by_age_group = TRUE
  )

  expect_identical(names(cfr_data), c("time", "age_group", "cases", "deaths"))
  expect_equal(cfr_data$cases, c(0, 0, 5, 6, 8, 9))
  expect_equal(cfr_data$deaths, c(0, 0, 1, 2, 3, 4))
})

test_that("as_cfr_data converts stochastic cumulative outputs", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  initial_state <- data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(ages$age_groups, times = 3),
    value = c(3, 3, 5, 5, 0, 0),
    stringsAsFactors = FALSE
  )

  output <- simulate_stochastic(
    initial_state = initial_state,
    times = c(0, 1, 5),
    model = SIRModel(gamma = 0.5),
    age_structure = ages,
    contact_matrix = matrix(c(2, 1, 1, 2), nrow = 2, byrow = TRUE),
    beta = 1,
    seed = 101,
    cumulative_flows = list(
      infections = list(from = "S", to = "I"),
      removals = list(from = "I", to = "R")
    )
  )

  cfr_data <- as_cfr_data(
    output,
    cases = list(cumulative_name = "infections"),
    deaths = "transition:I->R"
  )

  expect_identical(names(cfr_data), c("time", "cases", "deaths"))
  expect_equal(cfr_data$time, c(0, 1, 5))
  expect_true(all(cfr_data$cases >= 0))
  expect_true(all(cfr_data$deaths >= 0))
})

test_that("as_cfr_data validates selectors and cumulative output shape", {
  cumulative <- data.frame(
    time = rep(c(0, 1), each = 4),
    cumulative_name = rep(c("shared", "shared"), each = 2, times = 2),
    transition_id = rep(c("infection:S->I", "transition:I->R"), each = 2, times = 2),
    age_group = rep(c("0-4", "5-9"), times = 4),
    value = c(0, 0, 0, 0, 1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  expect_error(
    as_cfr_data(data.frame(time = 0, value = 0), cases = "x", deaths = "y"),
    "missing required column"
  )
  expect_error(
    as_cfr_data(cumulative, cases = "missing", deaths = "shared"),
    "did not match"
  )
  expect_error(
    as_cfr_data(cumulative, cases = "shared", deaths = "transition:I->R"),
    "ambiguous"
  )
  expect_error(
    as_cfr_data(cumulative, cases = list(foo = "bar"), deaths = "transition:I->R"),
    "unknown field"
  )
})
