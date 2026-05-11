test_that("valid age structures pass validation", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )

  expect_s3_class(ages, "AgeStructure")
  expect_identical(ages$labels, c("0-4", "5-9", "10+"))
  expect_silent(validate_age_structure(ages))
})

test_that("age structures require matching field lengths", {
  expect_error(
    AgeStructure(
      age_groups = c("0-4", "5-9"),
      lower_bounds = c(0, 5),
      upper_bounds = c(4, 9, Inf)
    ),
    "same length"
  )
})

test_that("age structures require at least one age group", {
  expect_error(
    AgeStructure(
      age_groups = character(),
      lower_bounds = numeric(),
      upper_bounds = numeric()
    ),
    "at least one"
  )
})

test_that("age group labels must be character values", {
  expect_error(
    AgeStructure(
      age_groups = c(1, 2),
      lower_bounds = c(0, 5),
      upper_bounds = c(4, 9)
    ),
    "character"
  )
})

test_that("overlapping age bins fail validation", {
  expect_error(
    AgeStructure(
      age_groups = c("0-5", "5-9"),
      lower_bounds = c(0, 5),
      upper_bounds = c(5, 9)
    ),
    "non-overlapping"
  )
})

test_that("age bins must have upper bounds greater than lower bounds", {
  expect_error(
    AgeStructure(
      age_groups = c("0"),
      lower_bounds = c(0),
      upper_bounds = c(0)
    ),
    "greater than"
  )
})

test_that("unsorted age bins fail validation", {
  expect_error(
    AgeStructure(
      age_groups = c("5-9", "0-4"),
      lower_bounds = c(5, 0),
      upper_bounds = c(9, 4)
    ),
    "strictly increasing"
  )
})

test_that("duplicate age labels fail validation", {
  expect_error(
    AgeStructure(
      age_groups = c("0-4", "0-4"),
      lower_bounds = c(0, 5),
      upper_bounds = c(4, 9)
    ),
    "unique"
  )
})

test_that("only the final age bin may be open-ended", {
  expect_error(
    AgeStructure(
      age_groups = c("0+", "5-9"),
      lower_bounds = c(0, 5),
      upper_bounds = c(Inf, 9)
    ),
    "final"
  )
})
