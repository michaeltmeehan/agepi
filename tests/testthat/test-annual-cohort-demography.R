annual_cohort_test_ages <- function() {
  wpp_age_structure_1year(max_age = 2)
}

annual_cohort_zero_mortality <- function(ages = annual_cohort_test_ages()) {
  stats::setNames(rep(0, ages$n_age_groups), ages$age_groups)
}

annual_cohort_zero_fertility <- function(ages = annual_cohort_test_ages()) {
  stats::setNames(rep(0, ages$n_age_groups), ages$age_groups)
}

test_that("annual cohort step shifts a single cohort from age 0 to age 1", {
  ages <- annual_cohort_test_ages()

  out <- annual_cohort_demographic_step(
    population = c("0" = 10, "1" = 0, "2+" = 0),
    age_structure = ages,
    fertility = annual_cohort_zero_fertility(ages),
    mortality = annual_cohort_zero_mortality(ages)
  )

  expect_equal(out$population, c(0, 10, 0))
})

test_that("annual cohort step retains and receives survivors in final open-ended group", {
  ages <- annual_cohort_test_ages()

  out <- annual_cohort_demographic_step(
    population = c("0" = 0, "1" = 7, "2+" = 11),
    age_structure = ages,
    fertility = annual_cohort_zero_fertility(ages),
    mortality = annual_cohort_zero_mortality(ages)
  )

  expect_equal(out$population, c(0, 0, 18))
})

test_that("annual cohort ageing alone conserves population", {
  ages <- annual_cohort_test_ages()

  out <- annual_cohort_demographic_step(
    population = c("0" = 4, "1" = 6, "2+" = 8),
    age_structure = ages,
    fertility = annual_cohort_zero_fertility(ages),
    mortality = annual_cohort_zero_mortality(ages)
  )

  expect_equal(sum(out$population), 18)
})

test_that("annual cohort mortality uses annual hazards converted to survival", {
  ages <- AgeStructure("0+", lower_bounds = 0, upper_bounds = Inf)

  out <- annual_cohort_demographic_step(
    population = c("0+" = 100),
    age_structure = ages,
    fertility = c("0+" = 0),
    mortality = c("0+" = log(2))
  )

  expect_equal(out$population, 50)
  expect_equal(attr(out, "diagnostics")$deaths, 50)
})

test_that("annual cohort births enter age 0 after survival and ageing", {
  ages <- annual_cohort_test_ages()

  out <- annual_cohort_demographic_step(
    population = c("0" = 0, "1" = 20, "2+" = 0),
    age_structure = ages,
    fertility = c("0" = 0, "1" = 0.25, "2+" = 0),
    mortality = annual_cohort_zero_mortality(ages)
  )

  expect_equal(out$population, c(5, 0, 20))
  expect_equal(attr(out, "diagnostics")$births, 5)
})

test_that("annual cohort fertility exposure fraction scales births", {
  ages <- annual_cohort_test_ages()

  out <- annual_cohort_demographic_step(
    population = c("0" = 0, "1" = 20, "2+" = 0),
    age_structure = ages,
    fertility = c("0" = 0, "1" = 0.25, "2+" = 0),
    mortality = annual_cohort_zero_mortality(ages),
    fertility_exposure_fraction = 0.5
  )

  expect_equal(out$population, c(2.5, 0, 20))
  expect_equal(attr(out, "diagnostics")$births, 2.5)
})

test_that("annual cohort count-based migration adds and subtracts destination-age counts", {
  ages <- annual_cohort_test_ages()

  out <- annual_cohort_demographic_step(
    population = c("0" = 10, "1" = 20, "2+" = 30),
    age_structure = ages,
    fertility = annual_cohort_zero_fertility(ages),
    mortality = annual_cohort_zero_mortality(ages),
    migration = data.frame(
      age_group = ages$age_groups,
      migration_count = c(3, -4, 5),
      stringsAsFactors = FALSE
    )
  )

  expect_equal(out$population, c(3, 6, 55))
  expect_equal(sum(out$population), 60 + 4)
  expect_equal(attr(out, "diagnostics")$net_migration, 4)
})

test_that("annual cohort rate-based migration multiplies destination population after ageing and births", {
  ages <- annual_cohort_test_ages()

  out <- annual_cohort_demographic_step(
    population = c("0" = 10, "1" = 20, "2+" = 30),
    age_structure = ages,
    fertility = annual_cohort_zero_fertility(ages),
    mortality = annual_cohort_zero_mortality(ages),
    migration = data.frame(
      age_group = ages$age_groups,
      migration_rate = c(0.5, -0.1, 0),
      stringsAsFactors = FALSE
    )
  )

  expect_equal(out$population, c(0, 9, 50))
  expect_equal(attr(out, "diagnostics")$migration_type, "rate")
})

test_that("annual cohort step errors when migration makes final population negative", {
  ages <- annual_cohort_test_ages()

  expect_error(
    annual_cohort_demographic_step(
      population = c("0" = 0, "1" = 1, "2+" = 0),
      age_structure = ages,
      fertility = annual_cohort_zero_fertility(ages),
      mortality = annual_cohort_zero_mortality(ages),
      migration = data.frame(
        age_group = ages$age_groups,
        migration_count = c(0, -2, 0),
        stringsAsFactors = FALSE
      )
    ),
    "negative population"
  )
})

test_that("annual cohort step errors clearly for missing age groups", {
  ages <- annual_cohort_test_ages()

  expect_error(
    annual_cohort_demographic_step(
      population = data.frame(
        age_group = c("0", "2+"),
        population = c(10, 5),
        stringsAsFactors = FALSE
      ),
      age_structure = ages,
      fertility = annual_cohort_zero_fertility(ages),
      mortality = annual_cohort_zero_mortality(ages)
    ),
    "population is missing required age_group value"
  )
})
