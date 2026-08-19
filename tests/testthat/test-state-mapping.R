test_that("long-format state converts to a compartment-major vector", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")
  state_long <- data.frame(
    age_group = c("0-4", "0-4", "0-4", "5-9", "5-9", "5-9"),
    compartment = c("S", "I", "R", "S", "I", "R"),
    value = c(100, 2, 0, 120, 3, 1)
  )

  state_vector <- state_long_to_vector(state_long, ages, compartments)

  expect_equal(state_vector, c(100, 120, 2, 3, 0, 1), ignore_attr = TRUE)
  expect_identical(names(state_vector), c("S_0-4", "S_5-9", "I_0-4", "I_5-9", "R_0-4", "R_5-9"))
})

test_that("state vectors convert back to long format in the same ordering", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")

  state_long <- state_vector_to_long(c(100, 120, 2, 3, 0, 1), ages, compartments)

  expect_identical(state_long$compartment, c("S", "S", "I", "I", "R", "R"))
  expect_identical(state_long$age_group, c("0-4", "5-9", "0-4", "5-9", "0-4", "5-9"))
  expect_equal(state_long$value, c(100, 120, 2, 3, 0, 1))
})

test_that("round trip preserves totals by age and compartment", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")
  state_long <- data.frame(
    age_group = c("5-9", "0-4", "5-9", "0-4", "5-9", "0-4"),
    compartment = c("R", "S", "I", "R", "S", "I"),
    value = c(1, 100, 3, 0, 120, 2)
  )

  round_trip <- state_vector_to_long(
    state_long_to_vector(state_long, ages, compartments),
    ages,
    compartments
  )

  expect_equal(
    tapply(round_trip$value, round_trip$age_group, sum),
    tapply(state_long$value, state_long$age_group, sum)
  )
  expect_equal(
    tapply(round_trip$value, round_trip$compartment, sum),
    tapply(state_long$value, state_long$compartment, sum)
  )
})

test_that("round trip works with more than two compartments and age groups", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-14"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, 14)
  )
  compartments <- c("S", "E", "I", "R")
  state_long <- expand.grid(
    age_group = c("10-14", "0-4", "5-9"),
    compartment = c("I", "S", "R", "E"),
    stringsAsFactors = FALSE
  )
  state_long$value <- seq_len(nrow(state_long))

  state_vector <- state_long_to_vector(state_long, ages, compartments)
  round_trip <- state_vector_to_long(state_vector, ages, compartments)

  expect_identical(
    names(state_vector),
    c(
      "S_0-4", "S_5-9", "S_10-14",
      "E_0-4", "E_5-9", "E_10-14",
      "I_0-4", "I_5-9", "I_10-14",
      "R_0-4", "R_5-9", "R_10-14"
    )
  )
  expect_equal(sum(round_trip$value), sum(state_long$value))
  expect_identical(round_trip$compartment, rep(compartments, each = 3))
  expect_identical(round_trip$age_group, rep(ages$age_groups, times = 4))
})

test_that("state mapping requires exactly one row per compartment-age pair", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I")

  missing_row <- data.frame(
    age_group = c("0-4", "5-9", "0-4"),
    compartment = c("S", "S", "I"),
    value = c(1, 2, 3)
  )
  duplicate_row <- data.frame(
    age_group = c("0-4", "5-9", "0-4", "0-4"),
    compartment = c("S", "S", "I", "I"),
    value = c(1, 2, 3, 4)
  )

  expect_error(state_long_to_vector(missing_row, ages, compartments), "missing")
  expect_error(state_long_to_vector(duplicate_row, ages, compartments), "duplicate")
})

test_that("state mapping rejects non-numeric and missing state values", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I")

  non_numeric <- data.frame(
    age_group = c("0-4", "5-9", "0-4", "5-9"),
    compartment = c("S", "S", "I", "I"),
    value = c("1", "2", "3", "4")
  )
  missing_value <- data.frame(
    age_group = c("0-4", "5-9", "0-4", "5-9"),
    compartment = c("S", "S", "I", "I"),
    value = c(1, 2, NA, 4)
  )

  expect_error(state_long_to_vector(non_numeric, ages, compartments), "numeric")
  expect_error(state_long_to_vector(missing_value, ages, compartments), "missing")
})

test_that("state vectors must match the expected length", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    state_vector_to_long(c(1, 2, 3), ages, c("S", "I")),
    "length"
  )
})

test_that("state vector names are ignored when converting to long format", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  state_vector <- c(foo = 100, bar = 120, baz = 2, qux = 3)

  state_long <- state_vector_to_long(state_vector, ages, c("S", "I"))

  expect_identical(state_long$age_group, c("0-4", "5-9", "0-4", "5-9"))
  expect_identical(state_long$compartment, c("S", "S", "I", "I"))
  expect_equal(state_long$value, c(100, 120, 2, 3))
})

test_that("cached output templates and numeric state matrices match public state mappings", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")
  state_vector <- c(100, 120, 2, 3, 0, 1)

  context <- agepi:::prepare_transition_rate_context_validated(
    model = SIRModel(gamma = 0.2),
    age_structure = ages,
    contact_matrix = matrix(1, nrow = 2, ncol = 2),
    beta = 0.1,
    include_public_template = FALSE
  )
  expected_template <- agepi:::state_order(ages, compartments)[, c("compartment", "age_group"), drop = FALSE]
  expected_matrix <- matrix(
    state_vector,
    nrow = ages$n_age_groups,
    ncol = length(compartments),
    dimnames = list(ages$age_groups, compartments)
  )

  expect_identical(context$state_output_template, expected_template)
  expect_identical(
    agepi:::transition_state_matrix_from_vector(state_vector, context),
    expected_matrix
  )
  expect_equal(
    agepi:::simulation_state_output(
      state_vector,
      time = 7,
      age_structure = ages,
      compartments = compartments
    ),
    agepi:::simulation_state_output(
      state_vector,
      time = 7,
      age_structure = ages,
      compartments = compartments,
      state_template = expected_template
    )
  )
})

test_that("state trajectory materialisation matches repeated single-time outputs", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  compartments <- c("S", "I", "R")
  state_template <- agepi:::state_order(ages, compartments)[, c("compartment", "age_group"), drop = FALSE]
  trajectory_matrix <- matrix(
    c(100, 120, 2, 3, 0, 1, 90, 110, 5, 7, 1, 2),
    nrow = 2,
    byrow = TRUE
  )
  expected <- do.call(
    rbind,
    list(
      agepi:::simulation_state_output(
        trajectory_matrix[1, ],
        time = 0,
        age_structure = ages,
        compartments = compartments,
        state_template = state_template
      ),
      agepi:::simulation_state_output(
        trajectory_matrix[2, ],
        time = 1,
        age_structure = ages,
        compartments = compartments,
        state_template = state_template
      )
    )
  )
  actual <- agepi:::state_trajectory_to_data_frame(
    trajectory_matrix,
    times = c(0, 1),
    state_template = state_template
  )

  expect_identical(actual, expected)
})

test_that("initialise_compartments_from_proportions allocates residual population", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  population <- c(100, 200)
  proportions <- list(
    I = c(0.1, 0.2),
    R = c(0.3, 0.1)
  )

  state <- initialise_compartments_from_proportions(
    population = population,
    proportions = proportions,
    residual_compartment = "S",
    compartments = c("S", "I", "R"),
    age_structure = ages
  )

  expect_identical(state$compartment, c("S", "S", "I", "I", "R", "R"))
  expect_identical(state$age_group, rep(ages$age_groups, times = 3))
  expect_equal(state$value, c(60, 140, 10, 40, 30, 20))
})

test_that("initialise_compartments_from_proportions rejects over-allocation", {
  population <- c(a = 100, b = 200)

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(I = c(0.8, 0.2), R = c(0.3, 0.1)),
      residual_compartment = "S"
    ),
    "exceed 1"
  )
})

initialisation_test_age_structure <- function() {
  AgeStructure(
    age_groups = c("child", "adult"),
    lower_bounds = c(0, 15),
    upper_bounds = c(14, Inf)
  )
}

initialisation_test_compartments <- c(
  "M.tb",
  "Incipient",
  "Contained",
  "Cleared",
  "Sub.clin.lowinf",
  "Sub.clin.inf",
  "Clin.lowinf",
  "Clin.inf",
  "Treatment",
  "Recovered"
)

test_that("initialise_compartments_from_proportions supports concise scalar inputs and omitted compartments", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)
  concise <- initialise_compartments_from_proportions(
    population = population,
    proportions = list(
      Clin.inf = 0.01
    ),
    residual_compartment = "M.tb",
    compartments = initialisation_test_compartments,
    age_structure = ages
  )

  expect_identical(concise$compartment, rep(initialisation_test_compartments, each = ages$n_age_groups))
  expect_identical(concise$age_group, rep(ages$age_groups, times = length(initialisation_test_compartments)))
  expect_equal(
    concise$value,
    c(
      990, 1980,
      rep(0, 2 * 6),
      10, 20,
      rep(0, 2),
      0, 0
    )
  )
})

test_that("initialise_compartments_from_proportions treats an empty named list as all residual", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)

  state <- initialise_compartments_from_proportions(
    population = population,
    proportions = setNames(list(), character()),
    residual_compartment = "M.tb",
    compartments = initialisation_test_compartments,
    age_structure = ages
  )

  expect_equal(unname(state$value[state$compartment == "M.tb"]), unname(population))
  expect_equal(
    unname(state$value[state$compartment != "M.tb"]),
    rep(0, sum(state$compartment != "M.tb"))
  )
})

test_that("initialise_compartments_from_proportions realigns named age-specific vectors", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)

  state <- initialise_compartments_from_proportions(
    population = population,
    proportions = list(
      Clin.inf = c(adult = 0.02, child = 0.01),
      Recovered = c(0, 0)
    ),
    residual_compartment = "M.tb",
    compartments = c("M.tb", "Clin.inf", "Recovered"),
    age_structure = ages
  )

  expect_equal(
    state$value[state$compartment == "Clin.inf"],
    c(10, 40)
  )
  expect_equal(
    state$value[state$compartment == "Recovered"],
    c(0, 0)
  )
  expect_equal(
    state$value[state$compartment == "M.tb"],
    c(990, 1960)
  )
})

test_that("initialise_compartments_from_proportions accepts unnamed full-length age-specific vectors", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)

  state <- initialise_compartments_from_proportions(
    population = population,
    proportions = list(
      Clin.inf = c(0.01, 0.02)
    ),
    residual_compartment = "M.tb",
    compartments = c("M.tb", "Clin.inf"),
    age_structure = ages
  )

  expect_equal(state$value[state$compartment == "Clin.inf"], c(10, 40))
  expect_equal(state$value[state$compartment == "M.tb"], c(990, 1960))
})

test_that("initialise_compartments_from_proportions rejects malformed proportion lists", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(c(0.01, 0.02)),
      residual_compartment = "M.tb",
      compartments = initialisation_test_compartments,
      age_structure = ages
    ),
    "named list"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(),
      residual_compartment = "M.tb",
      compartments = initialisation_test_compartments,
      age_structure = ages
    ),
    "named list"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = structure(list(c(0.01, 0.02), c(0.03, 0.04)), names = c("Clin.inf", "")),
      residual_compartment = "M.tb",
      compartments = initialisation_test_compartments,
      age_structure = ages
    ),
    "missing or empty compartment names"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = structure(list(c(0.01, 0.02), c(0.03, 0.04)), names = c("Clin.inf", "Clin.inf")),
      residual_compartment = "M.tb",
      compartments = initialisation_test_compartments,
      age_structure = ages
    ),
    "Duplicate compartment"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clincal.inf = c(0.01, 0.02)),
      residual_compartment = "M.tb",
      compartments = initialisation_test_compartments,
      age_structure = ages
    ),
    "Unknown compartment"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(M.tb = c(0.01, 0.02)),
      residual_compartment = "M.tb",
      compartments = initialisation_test_compartments,
      age_structure = ages
    ),
    "must not be supplied"
  )
})

test_that("initialise_compartments_from_proportions rejects malformed age-specific vectors", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = c(0.01, 0.02, 0.03)),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "length 1 or 2"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = c(child = 0.01, 0.02)),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "each model age group exactly once"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = c(child = 0.01, child = 0.02)),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "each model age group exactly once"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = c(child = 0.01, unknown = 0.02)),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "each model age group exactly once"
  )
})

test_that("initialise_compartments_from_proportions rejects invalid proportion values", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = -0.01),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "between 0 and 1"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = 1.01),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "between 0 and 1"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = NA_real_),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "finite and non-missing"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = NaN),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "finite and non-missing"
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(Clin.inf = Inf),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf"),
      age_structure = ages
    ),
    "finite and non-missing"
  )
})

test_that("initialise_compartments_from_proportions enforces age-specific totals and tolerance", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)

  exact_total <- initialise_compartments_from_proportions(
    population = population,
    proportions = list(
      Clin.inf = c(0.3, 0.2),
      Recovered = c(0.7, 0.8)
    ),
    residual_compartment = "M.tb",
    compartments = c("M.tb", "Clin.inf", "Recovered"),
    age_structure = ages
  )

  expect_equal(exact_total$value[exact_total$compartment == "M.tb"], c(0, 0))
  expect_equal(exact_total$value[exact_total$compartment == "Clin.inf"], c(300, 400))
  expect_equal(exact_total$value[exact_total$compartment == "Recovered"], c(700, 1600))

  tolerance <- sqrt(.Machine$double.eps) / 2
  near_total <- initialise_compartments_from_proportions(
    population = population,
    proportions = list(
      Clin.inf = c(0.6, 0.4),
      Recovered = c(0.4 + tolerance, 0.6)
    ),
    residual_compartment = "M.tb",
    compartments = c("M.tb", "Clin.inf", "Recovered"),
    age_structure = ages
  )

  expect_equal(near_total$value[near_total$compartment == "M.tb"], c(0, 0))
  expect_equal(
    near_total$value[near_total$compartment == "Clin.inf"],
    c(600, 800)
  )
  expect_equal(
    near_total$value[near_total$compartment == "Recovered"],
    c(400 + tolerance * 1000, 1200)
  )

  expect_error(
    initialise_compartments_from_proportions(
      population = population,
      proportions = list(
        Clin.inf = c(0.8, 0.2),
        Recovered = c(0.3, 0.1)
      ),
      residual_compartment = "M.tb",
      compartments = c("M.tb", "Clin.inf", "Recovered"),
      age_structure = ages
    ),
    "exceed 1"
  )
})

test_that("concise and expanded state initialisation are identical", {
  ages <- initialisation_test_age_structure()
  population <- c(child = 1000, adult = 2000)
  tb_compartments <- initialisation_test_compartments

  concise <- initialise_compartments_from_proportions(
    population = population,
    proportions = list(
      Clin.inf = 0.01
    ),
    residual_compartment = "M.tb",
    compartments = tb_compartments,
    age_structure = ages
  )

  expanded <- initialise_compartments_from_proportions(
    population = population,
    proportions = list(
      Incipient = rep(0, ages$n_age_groups),
      Contained = rep(0, ages$n_age_groups),
      Cleared = rep(0, ages$n_age_groups),
      Sub.clin.lowinf = rep(0, ages$n_age_groups),
      Sub.clin.inf = rep(0, ages$n_age_groups),
      Clin.lowinf = rep(0, ages$n_age_groups),
      Clin.inf = rep(0.01, ages$n_age_groups),
      Treatment = rep(0, ages$n_age_groups),
      Recovered = rep(0, ages$n_age_groups)
    ),
    residual_compartment = "M.tb",
    compartments = tb_compartments,
    age_structure = ages
  )

  expect_identical(concise, expanded)
})
