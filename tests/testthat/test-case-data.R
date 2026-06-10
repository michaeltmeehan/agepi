test_that("as_agepi_cases maps exact numeric ages to age groups", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
  cases <- data.frame(
    case_id = c("a", "b", "c"),
    age = c(0, 9, 25),
    onset_day = c(1, 1, 2),
    stringsAsFactors = FALSE
  )

  mapped <- as_agepi_cases(cases, ages, age_col = "age")

  expect_s3_class(mapped, "data.frame")
  expect_identical(names(mapped), c("case_id", "age", "onset_day", "age_group"))
  expect_identical(as.character(mapped$age_group), c("0-4", "5-9", "10+"))
  expect_identical(levels(mapped$age_group), ages$age_groups)
  expect_true(is.ordered(mapped$age_group))
})

test_that("as_agepi_cases validates existing age-group labels", {
  ages <- AgeStructure(c("0-4", "5+"), c(0, 5), c(4, Inf))
  cases <- data.frame(
    case_id = 1:3,
    reported_group = c("5+", "0-4", "5+"),
    stringsAsFactors = FALSE
  )

  mapped <- as_agepi_cases(cases, ages, age_group_col = "reported_group")

  expect_identical(names(mapped), c("case_id", "reported_group", "age_group"))
  expect_identical(as.character(mapped$age_group), c("5+", "0-4", "5+"))
  expect_identical(levels(mapped$age_group), ages$age_groups)
})

test_that("as_agepi_cases requires one age source", {
  ages <- AgeStructure(c("0+"), 0, Inf)
  cases <- data.frame(age = 1, age_group = "0+")

  expect_error(as_agepi_cases(cases, ages), "exactly one")
  expect_error(
    as_agepi_cases(cases, ages, age_col = "age", age_group_col = "age_group"),
    "exactly one"
  )
})

test_that("as_agepi_cases errors clearly for invalid ages", {
  ages <- AgeStructure(c("0-4", "5-9"), c(0, 5), c(4, 9))

  expect_error(
    as_agepi_cases(data.frame(age = c(1, NA)), ages, age_col = "age"),
    "missing age values"
  )
  expect_error(
    as_agepi_cases(data.frame(age = c(1, Inf)), ages, age_col = "age"),
    "non-finite age values"
  )
  expect_error(
    as_agepi_cases(data.frame(age = c(1, 10)), ages, age_col = "age"),
    "outside age_structure bounds: 10"
  )
  expect_error(
    as_agepi_cases(data.frame(age = "1"), ages, age_col = "age"),
    "must be numeric"
  )
})

test_that("as_agepi_cases errors clearly for invalid age-group labels", {
  ages <- AgeStructure(c("0-4", "5+"), c(0, 5), c(4, Inf))

  expect_error(
    as_agepi_cases(data.frame(group = c("0-4", NA)), ages, age_group_col = "group"),
    "missing age-group values"
  )
  expect_error(
    as_agepi_cases(data.frame(group = c("0-4", "10+")), ages, age_group_col = "group"),
    "invalid age-group label"
  )
})

test_that("as_agepi_cases preserves linelist tags when linelist is available", {
  skip_if_not_installed("linelist")

  ages <- AgeStructure(c("0-4", "5+"), c(0, 5), c(4, Inf))
  cases <- data.frame(
    id = c("a", "b"),
    age = c(2, 8),
    onset = as.Date(c("2020-01-01", "2020-01-03"))
  )
  linelist_cases <- linelist::make_linelist(cases, id = "id", age = "age")

  mapped <- as_agepi_cases(linelist_cases, ages, age_col = "age")

  expect_s3_class(mapped, "linelist")
  expect_identical(linelist::tags(mapped)$id, "id")
  expect_identical(linelist::tags(mapped)$age, "age")
  expect_identical(as.character(mapped$age_group), c("0-4", "5+"))
})
