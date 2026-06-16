synthetic_wpp_projection <- function() {
  ages <- wpp_age_structure_1year(max_age = 3)
  data <- data.frame(
    location = rep(c("Kiribati", "Exampleland"), each = 2 * ages$n_age_groups),
    year = rep(rep(c(2020, 2021), each = ages$n_age_groups), times = 2),
    age = rep(c(0, 1, 2, "3+"), times = 4),
    pop = c(
      100, 90, 80, 70,
      110, 95, 82, 72,
      200, 180, 160, 140,
      210, 185, 162, 142
    ),
    stringsAsFactors = FALSE
  )

  list(ages = ages, data = data)
}

test_that("population_trajectory_from_wpp creates a complete projection table", {
  fixture <- synthetic_wpp_projection()

  trajectory <- population_trajectory_from_wpp(
    fixture$data,
    age_structure = fixture$ages,
    location = "Kiribati",
    years = 2020:2021,
    location_col = "location",
    population_col = "pop"
  )

  expect_s3_class(trajectory, "data.frame")
  expect_identical(names(trajectory), c("time", "age_group", "population"))
  expect_equal(trajectory$time, rep(2020:2021, each = fixture$ages$n_age_groups))
  expect_identical(trajectory$age_group, rep(fixture$ages$age_groups, times = 2))
  expect_equal(trajectory$population, c(100, 90, 80, 70, 110, 95, 82, 72))
})

test_that("population_trajectory_from_wpp rejects missing required columns", {
  fixture <- synthetic_wpp_projection()
  missing_population <- fixture$data
  missing_population$pop <- NULL

  expect_error(
    population_trajectory_from_wpp(
      missing_population,
      age_structure = fixture$ages,
      location = "Kiribati",
      years = 2020,
      location_col = "location",
      population_col = "pop"
    ),
    "missing required column"
  )
})

test_that("population_trajectory_from_wpp rejects missing locations and years", {
  fixture <- synthetic_wpp_projection()

  expect_error(
    population_trajectory_from_wpp(
      fixture$data,
      age_structure = fixture$ages,
      location = "Missingland",
      years = 2020,
      location_col = "location",
      population_col = "pop"
    ),
    "Requested location is not present"
  )

  expect_error(
    population_trajectory_from_wpp(
      fixture$data,
      age_structure = fixture$ages,
      location = "Kiribati",
      years = 2022,
      location_col = "location",
      population_col = "pop"
    ),
    "missing requested year"
  )
})

test_that("population_trajectory_from_wpp rejects incomplete age grids", {
  fixture <- synthetic_wpp_projection()
  incomplete <- fixture$data[-which(
    fixture$data$location == "Kiribati" &
      fixture$data$year == 2021 &
      fixture$data$age == 1
  )[1], ]

  expect_error(
    population_trajectory_from_wpp(
      incomplete,
      age_structure = fixture$ages,
      location = "Kiribati",
      years = 2020:2021,
      location_col = "location",
      population_col = "pop"
    ),
    "missing age_group"
  )
})

test_that("population_trajectory_from_wpp rejects invalid age labels", {
  fixture <- synthetic_wpp_projection()
  bad_age <- fixture$data
  bad_age$age[bad_age$location == "Kiribati" & bad_age$year == 2020][1] <- "not-an-age"

  expect_error(
    population_trajectory_from_wpp(
      bad_age,
      age_structure = fixture$ages,
      location = "Kiribati",
      years = 2020:2021,
      location_col = "location",
      population_col = "pop"
    ),
    "Unsupported WPP age label"
  )
})

test_that("population_trajectory_from_wpp rejects negative or nonfinite populations", {
  fixture <- synthetic_wpp_projection()
  negative <- fixture$data
  negative$pop[1] <- -1

  expect_error(
    population_trajectory_from_wpp(
      negative,
      age_structure = fixture$ages,
      location = "Kiribati",
      years = 2020,
      location_col = "location",
      population_col = "pop"
    ),
    "population cannot be negative"
  )

  nonfinite <- fixture$data
  nonfinite$pop[1] <- Inf

  expect_error(
    population_trajectory_from_wpp(
      nonfinite,
      age_structure = fixture$ages,
      location = "Kiribati",
      years = 2020,
      location_col = "location",
      population_col = "pop"
    ),
    "population must be finite"
  )
})

test_that("projection_population_vector supports exact step and linear policies", {
  fixture <- synthetic_wpp_projection()
  trajectory <- population_trajectory_from_wpp(
    fixture$data,
    age_structure = fixture$ages,
    location = "Kiribati",
    years = 2020:2021,
    location_col = "location",
    population_col = "pop"
  )

  expect_equal(
    projection_population_vector(trajectory, 2020, time_policy = "exact"),
    stats::setNames(c(100, 90, 80, 70), fixture$ages$age_groups)
  )
  expect_equal(
    projection_population_vector(trajectory, 2020.5, time_policy = "step"),
    stats::setNames(c(100, 90, 80, 70), fixture$ages$age_groups)
  )
  expect_equal(
    projection_population_vector(trajectory, 2020.5, time_policy = "linear"),
    stats::setNames(c(105, 92.5, 81, 71), fixture$ages$age_groups)
  )
  expect_error(
    projection_population_vector(trajectory, 2020.5, time_policy = "exact"),
    "time is not available"
  )
})
