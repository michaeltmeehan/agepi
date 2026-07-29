test_inspection_ages <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
}

test_inspection_contacts <- function() {
  matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)
}

test_inspection_state <- function(S = c(95, 90), I = c(5, 10), R = c(0, 0)) {
  data.frame(
    compartment = rep(c("S", "I", "R"), each = 2),
    age_group = rep(c("0-4", "5-9"), times = 3),
    value = c(S, I, R),
    stringsAsFactors = FALSE
  )
}

test_that("inspect_transitions returns logical transition metadata for SIR and generic models", {
  sir <- inspect_transitions(SIRModel(gamma = 0.2))

  expect_equal(
    sir[, c(
      "transition_id",
      "transition_type",
      "source",
      "destination",
      "rate_definition_type",
      "rate_source",
      "is_infection_driven"
    )],
    data.frame(
      transition_id = c("infection:S->I", "transition:I->R"),
      transition_type = c("infection", "transition"),
      source = c("S", "I"),
      destination = c("I", "R"),
      rate_definition_type = c("derived", "scalar"),
      rate_source = c("force_of_infection()", "model$gamma"),
      is_infection_driven = c(TRUE, FALSE),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(unclass(sir$rate_definition), list(NA_real_, 0.2))

  generic <- inspect_transitions(
    CompartmentModel(
      compartments = c("S", "I", "R"),
      infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
      transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
      outflows = data.frame(from = "I", rate = 0.05, stringsAsFactors = FALSE),
      infectious_compartments = "I"
    )
  )

  expect_equal(generic$transition_id, c("infection:S->I", "transition:I->R", "outflow:I"))
  expect_equal(generic$transition_type, c("infection", "transition", "outflow"))
  expect_equal(generic$rate_definition_type, c("scalar", "scalar", "scalar"))
  expect_equal(unclass(generic$rate_definition), list(1, 0.2, 0.05))
})

test_that("inspect_compartment_flows summarizes compartment connectivity", {
  model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    outflows = data.frame(from = "I", rate = 0.05, stringsAsFactors = FALSE),
    infectious_compartments = "I"
  )

  summary <- inspect_compartment_flows(model)
  expect_equal(summary$compartment, c("S", "I", "R"))
  expect_equal(summary$n_inflows, c(0, 1, 1))
  expect_equal(summary$n_outflows, c(1, 2, 0))
  expect_equal(
    lapply(summary$inflow_transition_ids, unclass),
    list(character(), "infection:S->I", "transition:I->R")
  )
  expect_equal(
    lapply(summary$outflow_transition_ids, unclass),
    list("infection:S->I", c("transition:I->R", "outflow:I"), character())
  )

  selected <- inspect_compartment_flows(model, "I")
  expect_equal(selected$direction, c("inflow", "outflow", "outflow"))
  expect_equal(selected$transition_id, c("infection:S->I", "transition:I->R", "outflow:I"))
  expect_equal(selected$transition_type, c("infection", "transition", "outflow"))
})

test_that("diagnose_model_structure reports reachability and duplicate transition issues", {
  model <- CompartmentModel(
    compartments = c("S", "I", "R"),
    infection_transitions = data.frame(from = "S", to = "I", stringsAsFactors = FALSE),
    transitions = data.frame(from = "I", to = "R", rate = 0.2, stringsAsFactors = FALSE),
    outflows = data.frame(from = "I", rate = 0.05, stringsAsFactors = FALSE),
    infectious_compartments = "I"
  )

  report <- diagnose_model_structure(
    model,
    initial_state = data.frame(
      compartment = rep(c("S", "I", "R"), each = 2),
      age_group = rep(c("0-4", "5-9"), times = 3),
      value = c(95, 90, 5, 10, 0, 0),
      stringsAsFactors = FALSE
    )
  )

  expect_true(any(report$check == "reachability" & report$severity == "pass"))
  expect_true(any(report$check == "reachable compartments" & report$severity == "pass"))

  duplicated_model <- model
  duplicated_model$transitions <- rbind(duplicated_model$transitions, duplicated_model$transitions[1, ])
  duplicated_model$transitions$transition_id[2] <- duplicated_model$transitions$transition_id[1]

  duplicate_report <- diagnose_model_structure(duplicated_model)
  expect_true(any(duplicate_report$check == "duplicate logical transition id" & duplicate_report$severity == "error"))
  expect_true(any(duplicate_report$check == "duplicate logical transition pair" & duplicate_report$severity == "error"))
})

test_that("inspect_transition_rates enriches evaluated rates with source-population diagnostics", {
  ages <- test_inspection_ages()
  state <- test_inspection_state()
  model <- SIRModel(gamma = 0.2)

  result <- inspect_transition_rates(
    state = state,
    model = model,
    age_structure = ages,
    contact_matrix = test_inspection_contacts(),
    beta = 0.1
  )

  expect_identical(
    names(result),
    c(
      "from",
      "to",
      "age_group",
      "transition_id",
      "transition_type",
      "rate_definition",
      "rate_definition_type",
      "rate_source",
      "is_infection_driven",
      "rate",
      "source_population",
      "per_capita_source_rate",
      "is_finite",
      "is_negative",
      "source_empty",
      "nonzero_from_empty_source"
    )
  )
  expect_equal(result$rate, c(1.9, 1, 4.95, 2))
  expect_equal(result$source_population, c(95, 5, 90, 10))
  expect_equal(result$per_capita_source_rate, c(0.02, 0.2, 0.055, 0.2))
  expect_true(all(result$is_finite))
  expect_false(any(result$is_negative))
  expect_false(any(result$source_empty))
  expect_false(any(result$nonzero_from_empty_source))
})

test_that("check_population_balance validates closed systems and external outflows", {
  ages <- test_inspection_ages()
  state <- test_inspection_state()
  sir_model <- SIRModel(gamma = 0.2)

  closed_balance <- check_population_balance(
    state = state,
    model = sir_model,
    age_structure = ages,
    contact_matrix = test_inspection_contacts(),
    beta = 0.1
  )

  expect_true(closed_balance$summary$passes_balance_check)
  expect_equal(closed_balance$summary$residual_balance_error, 0)
  expect_equal(closed_balance$summary$total_external_outflow, 0)

  outflow_model <- CompartmentModel(
    compartments = c("I"),
    transitions = data.frame(from = character(), to = character(), rate = numeric(), stringsAsFactors = FALSE),
    outflows = data.frame(from = "I", rate = 0.1, stringsAsFactors = FALSE),
    infectious_compartments = character()
  )
  outflow_state <- data.frame(
    compartment = "I",
    age_group = "0-4",
    value = 10,
    stringsAsFactors = FALSE
  )
  outflow_ages <- AgeStructure("0-4", 0, 4)

  outflow_balance <- check_population_balance(
    state = outflow_state,
    model = outflow_model,
    age_structure = outflow_ages,
    contact_matrix = matrix(0, nrow = 1, ncol = 1)
  )

  expect_true(outflow_balance$summary$passes_balance_check)
  expect_equal(outflow_balance$summary$total_external_outflow, 1)
  expect_equal(outflow_balance$summary$total_derivative, -1)
  expect_equal(outflow_balance$by_age$total_external_outflow, 1)
  expect_equal(outflow_balance$by_age$residual_balance_error, 0)
})
