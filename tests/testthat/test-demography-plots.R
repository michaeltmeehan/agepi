test_plot_population <- function() {
  data.frame(
    time = rep(c(2020, 2025, 2030), each = 3),
    age_group = rep(c("0-4", "5-9", "10+"), times = 3),
    population = c(100, 90, 80, 105, 88, 82, 110, 86, 85),
    stringsAsFactors = FALSE
  )
}

test_that("demography plotting helpers return ggplot objects", {
  testthat::skip_if_not_installed("ggplot2")
  population <- test_plot_population()

  expect_s3_class(plot_population_pyramid(population, year = 2020), "ggplot")
  expect_s3_class(
    plot_population_pyramid(population, year = 2020, compare_year = 2030),
    "ggplot"
  )
  expect_s3_class(plot_population_projection(population), "ggplot")
  expect_s3_class(plot_age_structure(population), "ggplot")
})

test_that("demography plotting helpers accept agepi_demography objects", {
  testthat::skip_if_not_installed("ggplot2")
  ages <- AgeStructure(
    age_groups = c("0-4", "5-9", "10+"),
    lower_bounds = c(0, 5, 10),
    upper_bounds = c(4, 9, Inf)
  )
  demography <- Demography(test_plot_population(), ages)

  expect_s3_class(plot_population_projection(demography), "ggplot")
})

test_that("demography plotting helpers report missing required columns", {
  population <- test_plot_population()
  population$population <- NULL

  expect_error(
    plot_population_projection(population),
    "missing required column"
  )
})

test_that("population pyramid reports missing requested years", {
  testthat::skip_if_not_installed("ggplot2")

  expect_error(
    plot_population_pyramid(test_plot_population(), year = 2040),
    "requested year\\(s\\) not available"
  )
  expect_error(
    plot_population_pyramid(test_plot_population(), year = 2020, compare_year = 2040),
    "requested year\\(s\\) not available"
  )
})

test_that("age structure validates requested age groups", {
  testthat::skip_if_not_installed("ggplot2")

  expect_error(
    plot_age_structure(test_plot_population(), age_groups = c("0-4", "80+")),
    "requested age_groups not available"
  )
})

test_that("plot_demography returns the standard named plot list", {
  testthat::skip_if_not_installed("ggplot2")
  plots <- plot_demography(test_plot_population(), year = 2020)

  expect_type(plots, "list")
  expect_identical(
    names(plots),
    c("population_pyramid", "population_projection", "age_structure")
  )
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
})

test_that("age-group proportions sum to one within each time point", {
  proportions <- calculate_age_group_proportions(test_plot_population())
  sums <- stats::aggregate(proportion ~ time, data = proportions, FUN = sum)

  expect_equal(sums$proportion, rep(1, nrow(sums)))
})

test_that("age-group proportions sum to one after age-group filtering", {
  proportions <- calculate_age_group_proportions(
    test_plot_population(),
    age_groups = c("0-4", "10+")
  )
  sums <- stats::aggregate(proportion ~ time, data = proportions, FUN = sum)

  expect_equal(sums$proportion, rep(1, nrow(sums)))
})
