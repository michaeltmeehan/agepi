test_that("transform_contact_matrix identity returns the original matrix with age-group names", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  contact_matrix <- matrix(c(
    2, 1,
    3, 4
  ), nrow = 2, byrow = TRUE)

  result <- transform_contact_matrix(contact_matrix, ages, ages, population = c(100, 200))

  expected <- contact_matrix
  dimnames(expected) <- list(ages$age_groups, ages$age_groups)
  expect_equal(result, expected)
})

test_that("transform_contact_matrix performs simple exact aggregation", {
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
  contact_matrix <- matrix(
    c(
      1, 2, 3, 4,
      5, 6, 7, 8,
      9, 10, 11, 12,
      13, 14, 15, 16
    ),
    nrow = 4,
    byrow = TRUE
  )

  result <- transform_contact_matrix(
    contact_matrix,
    from_ages,
    to_ages,
    population = c(10, 10, 10, 10)
  )

  expected <- matrix(
    c(
      7, 11,
      23, 27
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(to_ages$age_groups, to_ages$age_groups)
  )
  expect_equal(result, expected)
})

test_that("transform_contact_matrix aggregation is recipient-population weighted", {
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
  contact_matrix <- matrix(
    c(
      1, 2, 10, 20,
      5, 6, 50, 60,
      7, 8, 70, 80,
      9, 10, 90, 100
    ),
    nrow = 4,
    byrow = TRUE
  )

  result <- transform_contact_matrix(
    contact_matrix,
    from_ages,
    to_ages,
    population = c(1, 3, 2, 8)
  )

  expected <- matrix(
    c(
      (1 * 3 + 3 * 11) / 4, (1 * 30 + 3 * 110) / 4,
      (2 * 15 + 8 * 19) / 10, (2 * 150 + 8 * 190) / 10
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(to_ages$age_groups, to_ages$age_groups)
  )
  expect_equal(result, expected)
})

test_that("transform_contact_matrix rejects invalid population length", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    transform_contact_matrix(diag(2), ages, ages, population = 100),
    "population length"
  )
})

test_that("transform_contact_matrix rejects missing infinite and negative population", {
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )

  expect_error(
    transform_contact_matrix(diag(2), ages, ages, population = c(100, NA_real_)),
    "finite"
  )
  expect_error(
    transform_contact_matrix(diag(2), ages, ages, population = c(100, Inf)),
    "finite"
  )
  expect_error(
    transform_contact_matrix(diag(2), ages, ages, population = c(100, -1)),
    "negative"
  )
})

test_that("transform_contact_matrix rejects zero aggregate recipient population denominators", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = "0-9",
    lower_bounds = 0,
    upper_bounds = 9
  )

  expect_error(
    transform_contact_matrix(diag(2), from_ages, to_ages, population = c(0, 0)),
    "positive"
  )
})

test_that("transform_contact_matrix rejects incompatible age structures", {
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
    transform_contact_matrix(diag(2), from_ages, to_ages, population = c(100, 200)),
    "exact union"
  )
})

test_that("transform_contact_matrix rejects transformations requiring source-bin splitting", {
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
    transform_contact_matrix(matrix(1, nrow = 1, ncol = 1), from_ages, to_ages, population = 100),
    "exact union"
  )
})

test_that("transform_contact_matrix returns target row and column names", {
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

  result <- transform_contact_matrix(diag(2), from_ages, to_ages, population = c(100, 200))

  expect_identical(rownames(result), to_ages$age_groups)
  expect_identical(colnames(result), to_ages$age_groups)
})

test_that("transform_contact_matrix result passes contact matrix validation", {
  from_ages <- AgeStructure(
    age_groups = c("0-4", "5-9"),
    lower_bounds = c(0, 5),
    upper_bounds = c(4, 9)
  )
  to_ages <- AgeStructure(
    age_groups = "0-9",
    lower_bounds = 0,
    upper_bounds = 9
  )

  result <- transform_contact_matrix(diag(2), from_ages, to_ages, population = c(100, 200))

  expect_silent(validate_contact_matrix(result, to_ages))
})

test_that("force_of_infection accepts transformed contact matrices", {
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
  contact_matrix <- matrix(
    c(
      1, 2, 3, 4,
      5, 6, 7, 8,
      9, 10, 11, 12,
      13, 14, 15, 16
    ),
    nrow = 4,
    byrow = TRUE
  )
  population <- c(100, 200, 300, 400)
  transformed_contact_matrix <- transform_contact_matrix(
    contact_matrix,
    from_ages,
    to_ages,
    population = population
  )
  transformed_population <- transform_age_vector(population, from_ages, to_ages)

  lambda <- force_of_infection(
    infectious = c(10, 20),
    population = transformed_population,
    contact_matrix = transformed_contact_matrix,
    age_structure = to_ages
  )

  expect_type(lambda, "double")
  expect_length(lambda, to_ages$n_age_groups)
  expect_identical(names(lambda), to_ages$age_groups)
})
