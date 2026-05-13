test_age_structure <- function() {
  AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
}

test_demography_table <- function() {
  data.frame(
    time = c(1, 1, 1, 0, 0, 0),
    age_group = c("0-4", "5-9", "10+", "0-4", "5-9", "10+"),
    population = c(100, 90, 80, 110, 95, 85),
    stringsAsFactors = FALSE
  )
}

test_that("valid demography tables pass validation", {
  ages <- test_age_structure()
  demography <- test_demography_table()

  expect_silent(validate_demography_table(demography, ages))
  expect_invisible(validate_demography_table(demography, ages))
})

test_that("Demography constructs a valid inspectable object", {
  ages <- test_age_structure()
  demography <- Demography(test_demography_table(), ages)

  expect_s3_class(demography, "agepi_demography")
  expect_identical(demography$age_structure, ages)
  expect_identical(demography$times, c(0, 1))
  expect_identical(demography$n_times, 2L)
  expect_identical(demography$age_groups, ages$age_groups)
  expect_identical(demography$n_age_groups, ages$n_age_groups)
  expect_identical(
    names(demography),
    c(
      "demography",
      "age_structure",
      "times",
      "n_times",
      "age_groups",
      "n_age_groups"
    )
  )
  expect_identical(names(demography$demography), c("time", "age_group", "population"))
})

test_that("demography validation rejects non-data-frame input", {
  expect_error(
    validate_demography_table(list(), test_age_structure()),
    "demography must be a data frame"
  )
})

test_that("demography validation requires time, age_group, and population columns", {
  demography <- test_demography_table()
  demography$population <- NULL

  expect_error(
    validate_demography_table(demography, test_age_structure()),
    "missing required column"
  )

  demography <- test_demography_table()
  demography$time <- NULL
  expect_error(
    validate_demography_table(demography, test_age_structure()),
    "missing required column"
  )

  demography <- test_demography_table()
  demography$age_group <- NULL
  expect_error(
    validate_demography_table(demography, test_age_structure()),
    "missing required column"
  )
})

test_that("demography validation rejects invalid time values", {
  ages <- test_age_structure()
  demography <- test_demography_table()
  demography$time <- as.character(demography$time)
  expect_error(validate_demography_table(demography, ages), "time must be numeric")

  demography <- test_demography_table()
  demography$time[1] <- NA_real_
  expect_error(validate_demography_table(demography, ages), "time must be finite and non-missing")

  demography <- test_demography_table()
  demography$time[1] <- Inf
  expect_error(validate_demography_table(demography, ages), "time must be finite and non-missing")
})

test_that("demography validation rejects invalid population values", {
  ages <- test_age_structure()
  demography <- test_demography_table()
  demography$population <- as.character(demography$population)
  expect_error(validate_demography_table(demography, ages), "population must be numeric")

  demography <- test_demography_table()
  demography$population[1] <- -1
  expect_error(validate_demography_table(demography, ages), "population cannot be negative")

  demography <- test_demography_table()
  demography$population[1] <- NA_real_
  expect_error(validate_demography_table(demography, ages), "population must be finite and non-missing")

  demography <- test_demography_table()
  demography$population[1] <- Inf
  expect_error(validate_demography_table(demography, ages), "population must be finite and non-missing")
})

test_that("demography validation rejects duplicate time-age_group rows", {
  demography <- rbind(test_demography_table(), test_demography_table()[1, ])

  expect_error(
    validate_demography_table(demography, test_age_structure()),
    "duplicate time-age_group"
  )
})

test_that("demography validation rejects missing age groups at a time point", {
  demography <- test_demography_table()
  demography <- demography[-1, ]

  expect_error(
    validate_demography_table(demography, test_age_structure()),
    "missing age_group"
  )
})

test_that("demography validation rejects extra age groups", {
  demography <- test_demography_table()
  demography$age_group[1] <- "15+"

  expect_error(
    validate_demography_table(demography, test_age_structure()),
    "not in age_structure"
  )
})

test_that("demography validation rejects invalid age structures", {
  ages <- test_age_structure()
  ages$n_age_groups <- 4

  expect_error(
    validate_demography_table(test_demography_table(), ages),
    "n_age_groups"
  )
})

test_that("Demography sorts rows by increasing time and age-structure order", {
  demography <- Demography(test_demography_table(), test_age_structure())

  expect_identical(demography$times, c(0, 1))
  expect_identical(demography$demography$time, c(0, 0, 0, 1, 1, 1))
  expect_identical(
    demography$demography$age_group,
    c("0-4", "5-9", "10+", "0-4", "5-9", "10+")
  )
})
