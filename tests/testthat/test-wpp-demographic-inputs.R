test_that("population_from_wpp creates sorted Demography from tidy WPP-style input", {
  ages <- wpp_age_structure_5year(max_age = 10)
  population_like <- data.frame(
    location = "Exampleland",
    year = rep(c(2025, 2020), each = ages$n_age_groups),
    age = rep(rev(ages$age_groups), times = 2),
    pop = c(130, 120, 110, 100, 90, 80),
    stringsAsFactors = FALSE
  )

  demography <- population_from_wpp(
    data = population_like,
    age_structure = ages,
    time_col = "year",
    age_group_col = "age",
    population_col = "pop",
    location = "Exampleland",
    location_col = "location"
  )

  expect_s3_class(demography, "agepi_demography")
  expect_identical(demography$demography$time, rep(c(2020, 2025), each = ages$n_age_groups))
  expect_identical(demography$demography$age_group, rep(ages$age_groups, times = 2))
  expect_identical(
    demography_population_vector(demography, time = 2020),
    stats::setNames(c(80, 90, 100), ages$age_groups)
  )
})

test_that("population_from_wpp maps WPP lower-bound and open-ended age labels", {
  ages <- wpp_age_structure_5year(max_age = 10)
  population_like <- data.frame(
    year = 2020,
    age = c("0", "5", "10+"),
    pop = c(100, 90, 80),
    stringsAsFactors = FALSE
  )

  demography <- population_from_wpp(
    data = population_like,
    age_structure = ages,
    time_col = "year",
    age_group_col = "age",
    population_col = "pop"
  )

  expect_identical(demography$demography$age_group, ages$age_groups)
})

test_that("population_from_wpp rejects missing duplicate and invalid population inputs", {
  ages <- wpp_age_structure_5year(max_age = 10)
  population_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    pop = c(100, 90, 80),
    stringsAsFactors = FALSE
  )

  expect_error(
    population_from_wpp(
      data = population_like[-1, ],
      age_structure = ages,
      time_col = "year",
      age_group_col = "age",
      population_col = "pop"
    ),
    "missing age_group"
  )
  expect_error(
    population_from_wpp(
      data = rbind(population_like, population_like[1, ]),
      age_structure = ages,
      time_col = "year",
      age_group_col = "age",
      population_col = "pop"
    ),
    "duplicate time-age_group"
  )

  negative_population <- population_like
  negative_population$pop[1] <- -1
  expect_error(
    population_from_wpp(
      data = negative_population,
      age_structure = ages,
      time_col = "year",
      age_group_col = "age",
      population_col = "pop"
    ),
    "population cannot be negative"
  )

  infinite_population <- population_like
  infinite_population$pop[1] <- Inf
  expect_error(
    population_from_wpp(
      data = infinite_population,
      age_structure = ages,
      time_col = "year",
      age_group_col = "age",
      population_col = "pop"
    ),
    "population must be finite"
  )
})

test_that("population_from_wpp does not silently mix multiple locations", {
  ages <- wpp_age_structure_5year(max_age = 10)
  one_location <- data.frame(
    location = "A",
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    pop = c(100, 90, 80),
    stringsAsFactors = FALSE
  )
  population_like <- rbind(one_location, transform(one_location, location = "B"))

  expect_error(
    population_from_wpp(
      data = population_like,
      age_structure = ages,
      time_col = "year",
      age_group_col = "age",
      population_col = "pop",
      location_col = "location"
    ),
    "multiple location values"
  )
})

test_that("standardise_wpp_mortality creates sorted MortalitySchedule", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality_like <- data.frame(
    year = rep(c(2025, 2020), each = ages$n_age_groups),
    age = rep(rev(ages$age_groups), times = 2),
    mx = rep(0.01, 2 * ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  mortality <- standardise_wpp_mortality(
    mortality_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    mortality_col = "mx"
  )

  expect_s3_class(mortality, "agepi_mortality_schedule")
  expect_identical(mortality$rate_convention, "annual_hazard")
  expect_identical(mortality$data$time[1:ages$n_age_groups], rep(2020, ages$n_age_groups))
  expect_identical(mortality$data$age_group[1:ages$n_age_groups], ages$age_groups)
})

test_that("mortality_from_wpp_mx converts annual central death rates", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality_like <- data.frame(
    year = rep(c(2025, 2020), each = ages$n_age_groups),
    age = rep(rev(ages$age_groups), times = 2),
    mx = c(0.04, 0.02, 0.01, 0.03, 0.015, 0.005),
    stringsAsFactors = FALSE
  )

  mortality <- mortality_from_wpp_mx(
    mortality_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    mx_col = "mx"
  )

  expect_s3_class(mortality, "agepi_mortality_schedule")
  expect_identical(mortality$rate_convention, "annual_hazard")
  expect_identical(mortality$data$age_group[1:ages$n_age_groups], ages$age_groups)
  expect_equal(mortality$data$mortality_rate[1:ages$n_age_groups], c(0.005, 0.015, 0.03))
})

test_that("mortality_from_wpp_mx is compatible with demographic derivatives", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality <- mortality_from_wpp_mx(
    data.frame(
      year = rep(2020, ages$n_age_groups),
      age = ages$age_groups,
      mx = c(0.01, 0.02, 0.03),
      stringsAsFactors = FALSE
    ),
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    mx_col = "mx"
  )
  process <- DemographicProcess(ages, mortality_schedule = mortality)

  derivative <- demographic_derivative(c(100, 100, 100), time = 2020, process = process)

  expect_equal(unname(derivative), c(-21, -2, 17))
})

test_that("mortality_from_wpp_mx rejects missing duplicate and invalid mortality inputs", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    mx = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )

  expect_error(
    mortality_from_wpp_mx(mortality_like[-1, ], ages, "year", "age", "mx"),
    "missing age_group"
  )
  expect_error(
    mortality_from_wpp_mx(rbind(mortality_like, mortality_like[1, ]), ages, "year", "age", "mx"),
    "duplicate time-age_group"
  )

  negative_mortality <- mortality_like
  negative_mortality$mx[1] <- -0.01
  expect_error(
    mortality_from_wpp_mx(negative_mortality, ages, "year", "age", "mx"),
    "negative"
  )

  non_finite_mortality <- mortality_like
  non_finite_mortality$mx[1] <- Inf
  expect_error(
    mortality_from_wpp_mx(non_finite_mortality, ages, "year", "age", "mx"),
    "finite non-missing"
  )
})

test_that("mortality_from_wpp rejects unsupported ambiguous mortality quantities", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    qx = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )

  expect_error(
    mortality_from_wpp(mortality_like, ages, "year", "age", "qx", quantity = "qx"),
    "'arg' should be"
  )
})

test_that("standardise_wpp_mortality validates input and delegates rate checks", {
  ages <- wpp_age_structure_5year(max_age = 10)
  mortality_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    mx = rep(0.01, ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  missing_mx <- mortality_like
  missing_mx$mx <- NULL
  expect_error(
    standardise_wpp_mortality(missing_mx, ages, "year", "age", "mx"),
    "missing required column"
  )

  invalid_age <- mortality_like
  invalid_age$age[1] <- "not-an-age"
  expect_error(
    standardise_wpp_mortality(invalid_age, ages, "year", "age", "mx"),
    "Unsupported WPP age label"
  )

  negative_mortality <- mortality_like
  negative_mortality$mx[1] <- -0.01
  expect_error(
    standardise_wpp_mortality(negative_mortality, ages, "year", "age", "mx"),
    "negative"
  )
})

test_that("standardise_wpp_fertility creates FertilitySchedule with partial age coverage", {
  ages <- wpp_age_structure_5year()
  fertility_like <- data.frame(
    year = c(2025, 2020, 2020),
    age = c("25-29", "20-24", "15-19"),
    asfr = c(0.1, 0.08, 0.04),
    stringsAsFactors = FALSE
  )

  fertility <- standardise_wpp_fertility(
    fertility_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    fertility_col = "asfr"
  )

  expect_s3_class(fertility, "agepi_fertility_schedule")
  expect_identical(fertility$rate_convention, "births_per_female_person_year")
  expect_identical(fertility$data$time, c(2020, 2020, 2025))
  expect_identical(fertility$data$age_group, c("15-19", "20-24", "25-29"))
})

test_that("standardise_wpp_fertility validates inputs and non-negative rates", {
  ages <- wpp_age_structure_5year()
  fertility_like <- data.frame(
    year = 2020,
    age = "20-24",
    asfr = 0.08,
    stringsAsFactors = FALSE
  )

  missing_asfr <- fertility_like
  missing_asfr$asfr <- NULL
  expect_error(
    standardise_wpp_fertility(missing_asfr, ages, "year", "age", "asfr"),
    "missing required column"
  )

  negative_fertility <- fertility_like
  negative_fertility$asfr <- -0.01
  expect_error(
    standardise_wpp_fertility(negative_fertility, ages, "year", "age", "asfr"),
    "negative"
  )
})

wpp_percent_asfr_weights <- function(value_type = c("percent", "fraction")) {
  value_type <- match.arg(value_type)
  values <- c(5, 20, 30, 25, 15, 4, 1)
  if (identical(value_type, "fraction")) {
    values <- values / 100
  }

  data.frame(
    year = rep(c(2025, 2020), each = 7),
    age = rep(c("45-49", "40-44", "35-39", "30-34", "25-29", "20-24", "15-19"), times = 2),
    weight = rep(rev(values), times = 2),
    stringsAsFactors = FALSE
  )
}

test_that("fertility_from_wpp_percent_asfr converts fraction schedules", {
  ages <- wpp_age_structure_5year()
  weights <- wpp_percent_asfr_weights("fraction")
  tfr <- data.frame(year = c(2020, 2025), tfr = c(2.0, 2.5))

  fertility <- fertility_from_wpp_percent_asfr(
    weights,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    weight_col = "weight",
    tfr_data = tfr,
    tfr_time_col = "year",
    tfr_col = "tfr",
    weight_type = "fraction"
  )

  expect_s3_class(fertility, "agepi_fertility_schedule")
  expect_identical(fertility$rate_convention, "births_per_female_person_year")
  expect_identical(fertility$data$time, rep(c(2020, 2025), each = 7))
  expect_identical(
    fertility$data$age_group,
    rep(c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"), times = 2)
  )
  expect_equal(
    fertility$data$fertility_rate[1:7],
    2.0 * c(0.05, 0.20, 0.30, 0.25, 0.15, 0.04, 0.01) / 5
  )
})

test_that("fertility_from_wpp_percent_asfr converts percent schedules", {
  ages <- wpp_age_structure_5year()
  weights <- wpp_percent_asfr_weights("percent")
  tfr <- data.frame(year = c(2020, 2025), tfr = c(2.0, 2.5))

  fertility <- fertility_from_wpp_percent_asfr(
    weights,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    weight_col = "weight",
    tfr_data = tfr,
    tfr_time_col = "year",
    tfr_col = "tfr",
    weight_type = "percent"
  )

  expect_equal(
    fertility$data$fertility_rate[1:7],
    2.0 * c(5, 20, 30, 25, 15, 4, 1) / 100 / 5
  )
})

test_that("fertility_from_wpp_percent_asfr rejects missing age groups duplicates and bad sums", {
  ages <- wpp_age_structure_5year()
  weights <- subset(wpp_percent_asfr_weights("percent"), year == 2020)
  tfr <- data.frame(year = 2020, tfr = 2.0)

  expect_error(
    fertility_from_wpp_percent_asfr(
      weights[-1, ],
      ages,
      "year",
      "age",
      "weight",
      tfr,
      "year",
      "tfr"
    ),
    "missing age_group"
  )

  expect_error(
    fertility_from_wpp_percent_asfr(
      rbind(weights, weights[1, ]),
      ages,
      "year",
      "age",
      "weight",
      tfr,
      "year",
      "tfr"
    ),
    "duplicate time-age_group"
  )

  bad_sum <- weights
  bad_sum$weight[1] <- bad_sum$weight[1] + 1
  expect_error(
    fertility_from_wpp_percent_asfr(
      bad_sum,
      ages,
      "year",
      "age",
      "weight",
      tfr,
      "year",
      "tfr"
    ),
    "must sum to 100"
  )
})

test_that("fertility_from_wpp_percent_asfr validates TFR values and age-bin widths", {
  finite_ages <- AgeStructure(
    age_groups = c("15-19", "20-24"),
    lower_bounds = c(15, 20),
    upper_bounds = c(19, 24)
  )
  finite_weights <- data.frame(
    year = c(2020, 2020),
    age = c("15-19", "20-24"),
    weight = c(50, 50)
  )

  expect_error(
    fertility_from_wpp_percent_asfr(
      finite_weights,
      finite_ages,
      "year",
      "age",
      "weight",
      data.frame(year = 2020, tfr = -1),
      "year",
      "tfr",
      maternal_age_groups = c("15-19", "20-24")
    ),
    "TFR values cannot be negative"
  )

  open_ended_ages <- AgeStructure(
    age_groups = c("15-19", "20+"),
    lower_bounds = c(15, 20),
    upper_bounds = c(19, Inf)
  )
  weights <- data.frame(
    year = c(2020, 2020),
    age = c("15-19", "20+"),
    weight = c(50, 50)
  )

  expect_error(
    fertility_from_wpp_percent_asfr(
      weights,
      open_ended_ages,
      "year",
      "age",
      "weight",
      data.frame(year = 2020, tfr = 2),
      "year",
      "tfr",
      maternal_age_groups = c("15-19", "20+")
    ),
    "age-bin widths"
  )
})

test_that("fertility_from_wpp_percent_asfr rejects missing TFR times", {
  ages <- wpp_age_structure_5year()
  weights <- wpp_percent_asfr_weights("percent")
  tfr <- data.frame(year = 2020, tfr = 2.0)

  expect_error(
    fertility_from_wpp_percent_asfr(
      weights,
      ages,
      "year",
      "age",
      "weight",
      tfr,
      "year",
      "tfr"
    ),
    "missing TFR"
  )
})

test_that("fertility_from_wpp_percent_asfr does not silently mix full multi-location inputs", {
  ages <- wpp_age_structure_5year()
  one_location <- subset(wpp_percent_asfr_weights("percent"), year == 2020)
  weights <- rbind(
    transform(one_location, location = "A"),
    transform(one_location, location = "B")
  )

  expect_error(
    fertility_from_wpp_percent_asfr(
      weights,
      ages,
      "year",
      "age",
      "weight",
      data.frame(year = 2020, tfr = 2.0),
      "year",
      "tfr"
    ),
    "duplicate time-age_group"
  )

  expect_error(
    fertility_from_wpp_percent_asfr(
      one_location,
      ages,
      "year",
      "age",
      "weight",
      data.frame(year = c(2020, 2020), tfr = c(2.0, 3.0), location = c("A", "B")),
      "year",
      "tfr"
    ),
    "duplicate time"
  )
})

test_that("WPP percent ASFR fertility schedules work with demographic derivatives", {
  ages <- wpp_age_structure_5year(max_age = 50)
  weights <- subset(wpp_percent_asfr_weights("percent"), year == 2020)
  tfr <- data.frame(year = 2020, tfr = 2.0)
  fertility <- fertility_from_wpp_percent_asfr(
    weights,
    ages,
    "year",
    "age",
    "weight",
    tfr,
    "year",
    "tfr"
  )
  process <- DemographicProcess(ages, fertility_schedule = fertility)
  state <- rep(100, ages$n_age_groups)

  derivative <- demographic_derivative(state, time = 2020, process = process)

  expect_equal(derivative[["0-4"]], sum(fertility$data$fertility_rate * 100) - 20)
})

test_that("WPP population fertility and mortality adapters compose into demographic process", {
  ages <- wpp_age_structure_5year(max_age = 50)
  population_like <- data.frame(
    location = "Exampleland",
    year = rep(c(2020, 2021), each = ages$n_age_groups),
    age = rep(ages$age_groups, times = 2),
    pop = rep(seq(1000, 500, length.out = ages$n_age_groups), times = 2),
    stringsAsFactors = FALSE
  )
  fertility_weights <- subset(wpp_percent_asfr_weights("percent"), year == 2020)
  mortality_like <- data.frame(
    year = rep(c(2020, 2021), each = ages$n_age_groups),
    age = rep(ages$age_groups, times = 2),
    mx = rep(0.01, 2 * ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  demography <- population_from_wpp(
    population_like,
    age_structure = ages,
    time_col = "year",
    age_group_col = "age",
    population_col = "pop",
    location = "Exampleland",
    location_col = "location"
  )
  fertility <- fertility_from_wpp_percent_asfr(
    fertility_weights,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    weight_col = "weight",
    tfr_data = data.frame(year = 2020, tfr = 2.0),
    tfr_time_col = "year",
    tfr_col = "tfr"
  )
  mortality <- mortality_from_wpp_mx(
    mortality_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    mx_col = "mx"
  )

  process <- build_demographic_process(
    age_structure = ages,
    fertility_schedule = fertility,
    mortality_schedule = mortality
  )
  initial_state <- demography_population_at(demography, time = 2020)

  expect_s3_class(process, "agepi_demographic_process")
  expect_null(process$times)
  expect_silent(validate_demographic_process(process))
  expect_silent(demographic_derivative(initial_state, time = 2020, process = process))
  expect_silent(
    simulate_demography(
      process = process,
      initial_state = initial_state,
      times = c(2020, 2021),
      time_policy = "step"
    )
  )
})

test_that("WPP-derived demographic process supports missing optional schedules", {
  ages <- wpp_age_structure_5year(max_age = 10)
  population_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    pop = c(100, 90, 80),
    stringsAsFactors = FALSE
  )
  demography <- population_from_wpp(
    population_like,
    age_structure = ages,
    time_col = "year",
    age_group_col = "age",
    population_col = "pop"
  )

  process <- build_demographic_process(age_structure = ages)
  derivative <- demographic_derivative(
    demography_population_at(demography, time = 2020),
    time = 2020,
    process = process
  )

  expect_s3_class(process, "agepi_demographic_process")
  expect_null(process$fertility_schedule)
  expect_null(process$mortality_schedule)
  expect_equal(unname(derivative), c(-20, 2, 18))
})

test_that("WPP-derived demographic process rejects incompatible age structures", {
  ages <- wpp_age_structure_5year(max_age = 10)
  other_ages <- wpp_age_structure_5year(max_age = 15)
  mortality <- mortality_from_wpp_mx(
    data.frame(
      year = rep(2020, other_ages$n_age_groups),
      age = other_ages$age_groups,
      mx = rep(0.01, other_ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    age_structure = other_ages,
    time_col = "year",
    age_col = "age",
    mx_col = "mx"
  )

  expect_error(
    build_demographic_process(
      age_structure = ages,
      mortality_schedule = mortality
    ),
    "mortality_schedule must use the same age_structure"
  )
})

test_that("WPP-derived demographic process preserves existing time-policy behavior", {
  ages <- wpp_age_structure_5year(max_age = 50)
  fertility <- fertility_from_wpp_percent_asfr(
    subset(wpp_percent_asfr_weights("percent"), year == 2020),
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    weight_col = "weight",
    tfr_data = data.frame(year = 2020, tfr = 2.0),
    tfr_time_col = "year",
    tfr_col = "tfr"
  )
  mortality <- mortality_from_wpp_mx(
    data.frame(
      year = rep(c(2020, 2025), each = ages$n_age_groups),
      age = rep(ages$age_groups, times = 2),
      mx = rep(0.01, 2 * ages$n_age_groups),
      stringsAsFactors = FALSE
    ),
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    mx_col = "mx"
  )
  process <- build_demographic_process(
    age_structure = ages,
    fertility_schedule = fertility,
    mortality_schedule = mortality
  )

  expect_null(process$times)
  expect_error(
    demographic_derivative(rep(100, ages$n_age_groups), time = 2021, process = process),
    "Exact time 2021 is not available.*no interpolation"
  )
  expect_error(
    demographic_derivative(
      rep(100, ages$n_age_groups),
      time = 2021,
      process = process,
      time_policy = "step"
    ),
    "after the final available schedule time 2020"
  )
})

test_that("standardise_wpp_migration creates count and rate schedules", {
  ages <- wpp_age_structure_5year(max_age = 10)
  migration_like <- data.frame(
    year = rep(c(2025, 2020), each = ages$n_age_groups),
    age = rep(rev(ages$age_groups), times = 2),
    net = c(-10, 0, 3, 5, -2, 1),
    stringsAsFactors = FALSE
  )

  migration_count <- standardise_wpp_migration(
    migration_like,
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    migration_col = "net",
    migration_type = "count"
  )
  expect_s3_class(migration_count, "agepi_migration_schedule")
  expect_identical(migration_count$migration_type, "count")
  expect_true(any(migration_count$data$migration_count < 0))
  expect_identical(migration_count$data$age_group[1:ages$n_age_groups], ages$age_groups)

  migration_rate <- standardise_wpp_migration(
    transform(migration_like, net = net / 1000),
    age_structure = ages,
    time_col = "year",
    age_col = "age",
    migration_col = "net",
    migration_type = "rate"
  )
  expect_s3_class(migration_rate, "agepi_migration_schedule")
  expect_identical(migration_rate$migration_type, "rate")
  expect_true(any(migration_rate$data$migration_rate < 0))
})

test_that("standardise_wpp_migration validates migration_type", {
  ages <- wpp_age_structure_5year(max_age = 10)
  migration_like <- data.frame(
    year = rep(2020, ages$n_age_groups),
    age = ages$age_groups,
    net = rep(0, ages$n_age_groups),
    stringsAsFactors = FALSE
  )

  expect_error(
    standardise_wpp_migration(
      migration_like,
      age_structure = ages,
      time_col = "year",
      age_col = "age",
      migration_col = "net",
      migration_type = "invalid"
    ),
    "'arg' should be one of"
  )
})

test_that("WPP age labels map to one-year age structures", {
  ages <- wpp_age_structure_1year(max_age = 3)
  mortality_like <- data.frame(
    year = 2020,
    age = c(0, 1, 2, 3),
    mx = rep(0.01, ages$n_age_groups)
  )

  mortality <- standardise_wpp_mortality(mortality_like, ages, "year", "age", "mx")

  expect_identical(mortality$data$age_group, c("0", "1", "2", "3+"))
})

test_that("WPP age labels map to five-year and lower-bound age structures", {
  ages <- wpp_age_structure_5year(max_age = 10)
  five_year_like <- data.frame(
    year = 2020,
    age = c("0-4", "5-9", "10+"),
    mx = rep(0.01, ages$n_age_groups)
  )
  lower_bound_like <- data.frame(
    year = 2020,
    age = c("0", "5", "10+"),
    mx = rep(0.01, ages$n_age_groups)
  )

  five_year <- standardise_wpp_mortality(five_year_like, ages, "year", "age", "mx")
  lower_bound <- standardise_wpp_mortality(lower_bound_like, ages, "year", "age", "mx")

  expect_identical(five_year$data$age_group, ages$age_groups)
  expect_identical(lower_bound$data$age_group, ages$age_groups)
})

test_that("wpp2024 remains optional for demographic input standardisation", {
  testthat::skip_if_not_installed("wpp2024")
  expect_true(requireNamespace("wpp2024", quietly = TRUE))
})

test_that("wpp_location_rows filters one location and keeps requested columns", {
  data <- data.frame(
    name = c("A", "B", "A"),
    year = c(2020, 2020, 2021),
    age = c("0", "0", "1"),
    value = 1:3,
    stringsAsFactors = FALSE
  )

  rows <- wpp_location_rows(data, "A", columns = c("name", "year", "value"))

  expect_identical(rows$name, c("A", "A"))
  expect_identical(names(rows), c("name", "year", "value"))
})

test_that("collapse_wpp_open_age_counts sums terminal open ages by grouping columns", {
  data <- data.frame(
    name = "Exampleland",
    year = c(2020, 2020, 2020, 2021, 2021, 2021),
    age = c("94", "95", "96", "94", "95", "96"),
    pop = c(10, 2, 3, 20, 4, 5),
    stringsAsFactors = FALSE
  )

  collapsed <- collapse_wpp_open_age_counts(data, "pop", open_age = 95)

  expect_equal(collapsed$pop[collapsed$year == 2020 & collapsed$age == "95+"], 5)
  expect_equal(collapsed$pop[collapsed$year == 2021 & collapsed$age == "95+"], 9)
  expect_true(any(collapsed$age == "94"))
})

test_that("collapse_wpp_open_age_rates keeps the terminal open-age rate row", {
  data <- data.frame(
    year = 2020,
    age = c("94", "95", "96"),
    mx = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )

  collapsed <- collapse_wpp_open_age_rates(data, "mx", open_age = 95)

  expect_identical(collapsed$age, c("94", "95+"))
  expect_equal(collapsed$mx, c(0.1, 0.2))
})
