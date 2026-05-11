test_that("force of infection computes a basic two-age-group calculation", {
  infectious <- c(10, 20)
  population <- c(100, 200)
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(infectious, population, contact_matrix)

  expect_equal(lambda, c(0.3, 0.7), ignore_attr = TRUE)
})

test_that("force of infection returns an unnamed vector without age_structure", {
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(c(10, 20), c(100, 200), contact_matrix)

  expect_type(lambda, "double")
  expect_length(lambda, 2)
  expect_null(names(lambda))
})

test_that("contact matrix rows are recipients and columns are sources", {
  infectious <- c(10, 40)
  population <- c(100, 200)
  contact_matrix <- matrix(c(
    0, 10,
    5, 0
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(infectious, population, contact_matrix)

  expect_equal(lambda, c(2, 0.5), ignore_attr = TRUE)
})

test_that("force of infection applies beta scaling", {
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(c(10, 20), c(100, 200), contact_matrix, beta = 2)

  expect_equal(lambda, c(0.6, 1.4), ignore_attr = TRUE)
})

test_that("force of infection applies susceptibility by recipient age group", {
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(
    infectious = c(10, 20),
    population = c(100, 200),
    contact_matrix = contact_matrix,
    susceptibility = c(0.5, 2)
  )

  expect_equal(lambda, c(0.15, 1.4), ignore_attr = TRUE)
})

test_that("force of infection applies infectiousness by source age group", {
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(
    infectious = c(10, 20),
    population = c(100, 200),
    contact_matrix = contact_matrix,
    infectiousness = c(2, 0.5)
  )

  expect_equal(lambda, c(0.45, 0.8), ignore_attr = TRUE)
})

test_that("force of infection combines beta, susceptibility, and infectiousness", {
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(
    infectious = c(10, 20),
    population = c(100, 200),
    contact_matrix = contact_matrix,
    beta = 3,
    susceptibility = c(0.5, 2),
    infectiousness = c(2, 0.5)
  )

  expect_equal(lambda, c(0.675, 4.8), ignore_attr = TRUE)
})

test_that("zero infectious counts return zero force of infection", {
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(c(0, 0), c(100, 200), contact_matrix)

  expect_equal(lambda, c(0, 0), ignore_attr = TRUE)
})

test_that("force of infection names output from age_structure age_groups", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  lambda <- force_of_infection(c(10, 20), c(100, 200), contact_matrix, age_structure = ages)

  expect_equal(lambda, c("0-4" = 0.3, "5-9" = 0.7))
  expect_length(lambda, ages$n_age_groups)
  expect_identical(names(lambda), ages$age_groups)
})

test_that("force of infection rejects mismatched infectious length", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  contact_matrix <- diag(2)

  expect_error(
    force_of_infection(c(1, 2, 3), c(100, 200), contact_matrix, age_structure = ages),
    "infectious length"
  )
})

test_that("force of infection rejects mismatched population length", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  contact_matrix <- diag(2)

  expect_error(
    force_of_infection(c(1, 2), c(100, 200, 300), contact_matrix, age_structure = ages),
    "population length"
  )
})

test_that("force of infection rejects non-square contact matrices", {
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), matrix(1, nrow = 2, ncol = 3)),
    "square"
  )
})

test_that("force of infection rejects contact matrix dimension mismatch", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(3), age_structure = ages),
    "dimensions"
  )
})

test_that("force of infection rejects age_structure length mismatch", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10-14"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, 14)
  )

  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), age_structure = ages),
    "infectious length"
  )
})

test_that("force of infection rejects negative infectious counts", {
  expect_error(
    force_of_infection(c(1, -2), c(100, 200), diag(2)),
    "infectious cannot contain negative"
  )
})

test_that("force of infection rejects zero population values", {
  expect_error(
    force_of_infection(c(1, 2), c(100, 0), diag(2)),
    "positive"
  )
})

test_that("force of infection rejects negative population values", {
  expect_error(
    force_of_infection(c(1, 2), c(100, -200), diag(2)),
    "population cannot contain negative"
  )
})

test_that("force of infection rejects negative contact matrix entries", {
  contact_matrix <- matrix(c(
    1, -1,
    0, 1
  ), nrow = 2, byrow = TRUE)

  expect_error(
    force_of_infection(c(1, 2), c(100, 200), contact_matrix),
    "contact_matrix cannot contain negative"
  )
})

test_that("force of infection rejects negative susceptibility and infectiousness", {
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), susceptibility = c(1, -1)),
    "susceptibility cannot contain negative"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), infectiousness = c(1, -1)),
    "infectiousness cannot contain negative"
  )
})

test_that("force of infection rejects missing values", {
  expect_error(
    force_of_infection(c(1, NA), c(100, 200), diag(2)),
    "missing"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, NA), diag(2)),
    "missing"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), matrix(c(1, NA, 0, 1), nrow = 2)),
    "missing"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), susceptibility = c(1, NA)),
    "missing"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), infectiousness = c(1, NA)),
    "missing"
  )
})

test_that("force of infection rejects non-numeric inputs", {
  expect_error(
    force_of_infection(c("1", "2"), c(100, 200), diag(2)),
    "infectious must be a numeric vector"
  )
  expect_error(
    force_of_infection(c(1, 2), c("100", "200"), diag(2)),
    "population must be a numeric vector"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), matrix(c("1", "0", "0", "1"), nrow = 2)),
    "contact_matrix must be a numeric matrix"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), susceptibility = c("1", "1")),
    "susceptibility must be a numeric vector"
  )
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), infectiousness = c("1", "1")),
    "infectiousness must be a numeric vector"
  )
})

test_that("force of infection rejects negative beta", {
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), beta = -1),
    "beta cannot be negative"
  )
})

test_that("force of infection rejects non-finite beta", {
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), beta = Inf),
    "finite"
  )
})

test_that("force of infection rejects invalid susceptibility length", {
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), susceptibility = c(1, 1, 1)),
    "susceptibility length"
  )
})

test_that("force of infection rejects invalid infectiousness length", {
  expect_error(
    force_of_infection(c(1, 2), c(100, 200), diag(2), infectiousness = c(1, 1, 1)),
    "infectiousness length"
  )
})
