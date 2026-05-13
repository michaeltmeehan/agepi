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

test_that("transform_age_vector returns a named copy for identical age structures", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  result <- transform_age_vector(c("ignored" = 10, "also_ignored" = 20), ages, ages)

  expect_identical(result, c("0-4" = 10, "5-9" = 20))
})

test_that("transform_age_vector auto-detects exact aggregation", {
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
    transform_age_vector(c(100, 200, 300, 400), from_ages, to_ages),
    c("0-9" = 300, "10-19" = 700)
  )
})

test_that("transform_age_vector auto-detects exact segregation", {
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
    transform_age_vector(c(100, 300), from_ages, to_ages),
    c("0-4" = 50, "5-9" = 50, "10-14" = 150, "15-19" = 150)
  )
  expect_identical(
    transform_age_vector(
      c(100, 300),
      from_ages,
      to_ages,
      weights = c(1, 3, 2, 4),
      split_method = "weights"
    ),
    c("0-4" = 25, "5-9" = 75, "10-14" = 100, "15-19" = 200)
  )
})

test_that("transform_age_vector can reject auto-detected splitting without weights", {
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
    transform_age_vector(100, from_ages, to_ages, split_method = "error"),
    "rejects source age bins"
  )
  expect_error(
    transform_age_vector(100, from_ages, to_ages, split_method = "weights"),
    "weights must be supplied"
  )
})

test_that("transform_age_vector rejects incompatible structures in auto mode", {
  from_ages <- AgeStructure(
    age_groups = "0-9",
    lower_bounds = 0,
    upper_bounds = 9
  )
  to_ages <- AgeStructure(
    age_groups = c("0-3", "6-9"),
    lower_bounds = c(0, 6),
    upper_bounds = c(3, 9)
  )

  expect_error(
    transform_age_vector(100, from_ages, to_ages),
    "not compatible|compatible with an exact age transformation"
  )
})

test_that("transform_age_vector handles explicit directions", {
  fine_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  coarse_ages <- AgeStructure(
    age_groups = "0-9",
    lower_bounds = 0,
    upper_bounds = 9
  )

  expect_identical(
    transform_age_vector(c(100, 200), fine_ages, coarse_ages, direction = "aggregate"),
    c("0-9" = 300)
  )
  expect_identical(
    transform_age_vector(300, coarse_ages, fine_ages, weights = c(1, 2), direction = "segregate"),
    c("0-4" = 100, "5-9" = 200)
  )
  expect_error(
    transform_age_vector(c(100, 200), fine_ages, coarse_ages, weights = 1, direction = "aggregate"),
    "weights must not be supplied"
  )
  expect_error(
    transform_age_vector(300, coarse_ages, fine_ages, direction = "segregate"),
    "weights must be supplied"
  )
  expect_error(
    transform_age_vector(c(100, 200), fine_ages, coarse_ages, direction = "sideways"),
    "'arg' should be one of"
  )
})

test_that("transform_age_vector preserves target age-group names", {
  from_ages <- AgeStructure(
    age_groups = c("a", "b"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = "combined",
    lower_bounds = 0,
    upper_bounds = 9
  )

  expect_identical(
    names(transform_age_vector(c(5, 10), from_ages, to_ages)),
    "combined"
  )
})

test_that("transform_age_vector supports mixed exact aggregation and splitting", {
  from_ages <- AgeStructure(
    age_groups = c("0-9", "10-14", "15-19", "20+"),
    lower_bounds = c(0, 10, 15, 20),
    upper_bounds = c(9, 14, 19, Inf)
  )
  to_ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-19", "20+"),
    lower_bounds = c(0, 5, 10, 20),
    upper_bounds = c(4, 9, 19, Inf)
  )

  expect_identical(
    transform_age_vector(c(100, 200, 300, 400), from_ages, to_ages),
    c("0-4" = 50, "5-9" = 50, "10-19" = 500, "20+" = 400)
  )
})

test_that("transform_age_vector supports split methods for mixed exact transformations", {
  from_ages <- AgeStructure(
    age_groups = c("0-9", "10-19"),
    lower_bounds = c(0, 10),
    upper_bounds = c(9, 19)
  )
  to_ages <- AgeStructure(
    age_groups = c("0-3", "4-9", "10-19"),
    lower_bounds = c(0, 4, 10),
    upper_bounds = c(3, 9, 19)
  )

  expect_identical(
    transform_age_vector(c(100, 200), from_ages, to_ages, split_method = "width"),
    c("0-3" = 40, "4-9" = 60, "10-19" = 200)
  )
  expect_identical(
    transform_age_vector(c(100, 200), from_ages, to_ages, split_method = "equal"),
    c("0-3" = 50, "4-9" = 50, "10-19" = 200)
  )
  expect_identical(
    transform_age_vector(
      c(100, 200),
      from_ages,
      to_ages,
      weights = c(1, 3, 1),
      split_method = "weights"
    ),
    c("0-3" = 25, "4-9" = 75, "10-19" = 200)
  )
})

test_that("transform_age_vector rejects width splitting of open-ended bins", {
  from_ages <- AgeStructure(
    age_groups = "20+",
    lower_bounds = 20,
    upper_bounds = Inf
  )
  to_ages <- AgeStructure(
    age_groups = c("20-29", "30+"),
    lower_bounds = c(20, 30),
    upper_bounds = c(29, Inf)
  )

  expect_error(
    transform_age_vector(100, from_ages, to_ages, split_method = "width"),
    "open-ended"
  )
  expect_identical(
    transform_age_vector(100, from_ages, to_ages, split_method = "equal"),
    c("20-29" = 50, "30+" = 50)
  )
})
