test_that("aggregate_age_vector performs simple exact aggregation", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-14", "15-19"),
    lower_bounds = c(0, 5, 10, 15),
    upper_bounds = c(4, 9, 14, 19)
  )
  to_ages <- AgeStructure(
    age_groups = c("0-9", "10-19"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, 19)
  )

  expect_identical(
    aggregate_age_vector(c(100, 200, 300, 400), from_ages, to_ages),
    c("0-9" = 300, "10-19" = 700)
  )
})

test_that("aggregate_age_vector supports identity transformations", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )

  expect_identical(
    aggregate_age_vector(c(10, 20, 30), ages, ages),
    c("0-4" = 10, "5-9" = 20, "10+" = 30)
  )
})

test_that("aggregate_age_vector validates value length", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    aggregate_age_vector(1, from_ages, from_ages),
    "values length"
  )
})

test_that("aggregate_age_vector requires numeric finite non-missing values", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    aggregate_age_vector(c("1", "2"), ages, ages),
    "numeric"
  )
  expect_error(
    aggregate_age_vector(c(1, NA_real_), ages, ages),
    "finite"
  )
  expect_error(
    aggregate_age_vector(c(1, Inf), ages, ages),
    "finite"
  )
})

test_that("aggregate_age_vector rejects non-exact target bins", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = "0-8",
    lower_bounds = 0,
    upper_bounds = 8
  )

  expect_error(
    aggregate_age_vector(c(100, 200), from_ages, to_ages),
    "exact union"
  )
})

test_that("aggregate_age_vector rejects partial overlaps", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = "1-9",
    lower_bounds = 1,
    upper_bounds = 9
  )

  expect_error(
    aggregate_age_vector(c(100, 200), from_ages, to_ages),
    "exact union"
  )
})

test_that("aggregate_age_vector rejects target bins outside source range", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = "0-14",
    lower_bounds = 0,
    upper_bounds = 14
  )

  expect_error(
    aggregate_age_vector(c(100, 200), from_ages, to_ages),
    "exact union"
  )

  to_ages <- AgeStructure(
    age_groups = "under-9",
    lower_bounds = -1,
    upper_bounds = 9
  )

  expect_error(
    aggregate_age_vector(c(100, 200), from_ages, to_ages),
    "exact union"
  )
})

test_that("aggregate_age_vector rejects target bins spanning source gaps", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "10-14"),
    lower_bounds = c(0, 10),
    upper_bounds = c(4, 14)
  )
  to_ages <- AgeStructure(
    age_groups = "0-14",
    lower_bounds = 0,
    upper_bounds = 14
  )

  expect_error(
    aggregate_age_vector(c(100, 200), from_ages, to_ages),
    "exact union"
  )
})

test_that("aggregate_age_vector validates unordered or invalid age structures", {
  valid_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  unordered_ages <- list(
    age_groups = c("5-9", "0-4"),
    n_age_groups = 2,
    lower_bounds = c(5, 0),
    upper_bounds = c(9, 4)
  )
  class(unordered_ages) <- "AgeStructure"

  expect_error(
    aggregate_age_vector(c(100, 200), unordered_ages, valid_ages),
    "strictly increasing"
  )
  expect_error(
    aggregate_age_vector(c(100, 200), valid_ages, unordered_ages),
    "strictly increasing"
  )
})

test_that("segregate_age_vector performs simple exact segregation", {
  from_ages <- AgeStructure(
    age_groups = c("0-9", "10-19"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, 19)
  )
  to_ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-14", "15-19"),
    lower_bounds = c(0, 5, 10, 15),
    upper_bounds = c(4, 9, 14, 19)
  )

  expect_identical(
    segregate_age_vector(c(100, 300), from_ages, to_ages, weights = c(1, 3, 2, 4)),
    c("0-4" = 25, "5-9" = 75, "10-14" = 100, "15-19" = 200)
  )
})

test_that("segregate_age_vector supports identity transformations", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )

  expect_identical(
    segregate_age_vector(c(10, 20, 30), ages, ages, weights = c(1, 2, 3)),
    c("0-4" = 10, "5-9" = 20, "10+" = 30)
  )
})

test_that("segregate_age_vector validates value length", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    segregate_age_vector(1, ages, ages, weights = c(1, 1)),
    "values length"
  )
})

test_that("segregate_age_vector requires numeric finite non-missing values", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    segregate_age_vector(c("1", "2"), ages, ages, weights = c(1, 1)),
    "values must be numeric"
  )
  expect_error(
    segregate_age_vector(c(1, NA_real_), ages, ages, weights = c(1, 1)),
    "values must be finite"
  )
  expect_error(
    segregate_age_vector(c(1, Inf), ages, ages, weights = c(1, 1)),
    "values must be finite"
  )
})

test_that("segregate_age_vector validates weight length", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    segregate_age_vector(c(1, 2), ages, ages, weights = 1),
    "weights length"
  )
})

test_that("segregate_age_vector requires numeric finite non-missing non-negative weights", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    segregate_age_vector(c(1, 2), ages, ages, weights = c("1", "1")),
    "weights must be numeric"
  )
  expect_error(
    segregate_age_vector(c(1, 2), ages, ages, weights = c(1, NA_real_)),
    "weights must be finite"
  )
  expect_error(
    segregate_age_vector(c(1, 2), ages, ages, weights = c(1, Inf)),
    "weights must be finite"
  )
  expect_error(
    segregate_age_vector(c(1, 2), ages, ages, weights = c(1, -1)),
    "negative"
  )
})

test_that("segregate_age_vector rejects zero total weights within a source bin", {
  from_ages <- AgeStructure(
    age_groups = "0-9",
    lower_bounds = 0,
    upper_bounds = 9
  )
  to_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    segregate_age_vector(100, from_ages, to_ages, weights = c(0, 0)),
    "positive"
  )
})

test_that("segregate_age_vector rejects target bins partially overlapping two source bins", {
  from_ages <- AgeStructure(
    age_groups = c("0-9", "10-19"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, 19)
  )
  to_ages <- AgeStructure(
    age_groups = "5-14",
    lower_bounds = 5,
    upper_bounds = 14
  )

  expect_error(
    segregate_age_vector(c(100, 300), from_ages, to_ages, weights = 1),
    "fully contained"
  )
})

test_that("segregate_age_vector rejects target bins outside upper source range", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = c("0-4", "5-14"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 14)
  )

  expect_error(
    segregate_age_vector(c(100, 200), from_ages, to_ages, weights = c(1, 1)),
    "fully contained"
  )
})

test_that("segregate_age_vector rejects target bins outside lower source range", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = c("under-4", "5-9"),
    lower_bounds = c(-1, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    segregate_age_vector(c(100, 200), from_ages, to_ages, weights = c(1, 1)),
    "fully contained"
  )
})

test_that("segregate_age_vector rejects target gaps within a source bin", {
  from_ages <- AgeStructure(
    age_groups = "0-14",
    lower_bounds = 0,
    upper_bounds = 14
  )
  to_ages <- AgeStructure(
    age_groups = c("0-4", "10-14"),
    lower_bounds = c(0, 10),
    upper_bounds = c(4, 14)
  )

  expect_error(
    segregate_age_vector(100, from_ages, to_ages, weights = c(1, 1)),
    "exactly covered"
  )
})

test_that("segregate_age_vector rejects source bins not fully covered by target bins", {
  from_ages <- AgeStructure(
    age_groups = c("0-9", "10-19"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, 19)
  )
  to_ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-14"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, 14)
  )

  expect_error(
    segregate_age_vector(c(100, 300), from_ages, to_ages, weights = c(1, 1, 1)),
    "exactly covered"
  )
})

test_that("segregate_age_vector preserves target age-group names", {
  from_ages <- AgeStructure(
    age_groups = "children",
    lower_bounds = 0,
    upper_bounds = 9
  )
  to_ages <- AgeStructure(
    age_groups = c("early", "late"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_identical(
    names(segregate_age_vector(100, from_ages, to_ages, weights = c(2, 3))),
    c("early", "late")
  )
})
