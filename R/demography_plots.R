#' Plot an age-specific population distribution
#'
#' Creates a lightweight exploratory age-distribution plot for one population
#' year, with optional comparison to a second year. The input may be a
#' `Demography()`/`population_from_wpp()` object or a data frame with age, time,
#' and population columns.
#'
#' This is a plotting convenience only. It does not download data, interpolate
#' population values, or modify demographic or epidemic simulation semantics.
#'
#' @param population Population input as an `agepi_demography` object or a data
#'   frame.
#' @param year Year or time point to plot.
#' @param compare_year Optional second year or time point to compare.
#' @param age_col,time_col,value_col Column names containing age groups, times,
#'   and population values.
#'
#' @return A `ggplot` object.
#' @export
plot_population_pyramid <- function(
  population,
  year,
  compare_year = NULL,
  age_col = "age_group",
  time_col = "time",
  value_col = "population"
) {
  population <- demography_plot_population_table(
    population,
    age_col = age_col,
    time_col = time_col,
    value_col = value_col
  )
  selected <- filter_demography_plot_years(
    population,
    years = c(year, compare_year),
    time_col = time_col
  )
  selected <- aggregate_population_by_time_age(
    selected,
    age_col = age_col,
    time_col = time_col,
    value_col = value_col
  )
  require_ggplot2()
  selected$.year <- factor(selected[[time_col]], levels = unique(selected[[time_col]]))
  selected[[age_col]] <- factor(selected[[age_col]], levels = unique(population[[age_col]]))

  ggplot2::ggplot(
    selected,
    ggplot2::aes(x = get(age_col), y = get(value_col), fill = get(".year"))
  ) +
    ggplot2::geom_col(position = if (is.null(compare_year)) "stack" else "dodge") +
    ggplot2::labs(
      x = "Age group",
      y = "Population",
      fill = "Year"
    ) +
    ggplot2::theme_minimal()
}

#' Plot total population over time
#'
#' Aggregates an age-specific population table over age groups and plots total
#' population size by time.
#'
#' @inheritParams plot_population_pyramid
#'
#' @return A `ggplot` object.
#' @export
plot_population_projection <- function(
  population,
  age_col = "age_group",
  time_col = "time",
  value_col = "population"
) {
  population <- demography_plot_population_table(
    population,
    age_col = age_col,
    time_col = time_col,
    value_col = value_col
  )
  totals <- aggregate_total_population(
    population,
    time_col = time_col,
    value_col = value_col
  )
  require_ggplot2()

  ggplot2::ggplot(
    totals,
    ggplot2::aes(x = get(time_col), y = get(value_col))
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Year",
      y = "Total population"
    ) +
    ggplot2::theme_minimal()
}

#' Plot age-group population proportions over time
#'
#' Aggregates population by time and age group, converts counts to proportions
#' within each time point, and plots the age structure through time.
#'
#' @inheritParams plot_population_pyramid
#' @param age_groups Optional character vector of age groups to include.
#'
#' @return A `ggplot` object.
#' @export
plot_age_structure <- function(
  population,
  age_groups = NULL,
  age_col = "age_group",
  time_col = "time",
  value_col = "population"
) {
  population <- demography_plot_population_table(
    population,
    age_col = age_col,
    time_col = time_col,
    value_col = value_col
  )
  proportions <- calculate_age_group_proportions(
    population,
    age_groups = age_groups,
    age_col = age_col,
    time_col = time_col,
    value_col = value_col
  )
  require_ggplot2()
  proportions[[age_col]] <- factor(proportions[[age_col]], levels = unique(population[[age_col]]))

  ggplot2::ggplot(
    proportions,
    ggplot2::aes(x = get(time_col), y = get("proportion"), colour = get(age_col))
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Year",
      y = "Population proportion",
      colour = "Age group"
    ) +
    ggplot2::theme_minimal()
}

#' Create a standard set of exploratory demography plots
#'
#' Returns a named list of lightweight demography plots for a supplied
#' population table or `Demography()` object. The plots are not arranged into a
#' dashboard, avoiding any additional layout dependency.
#'
#' @inheritParams plot_age_structure
#' @param year Year or time point for the population pyramid. If `NULL`, the
#'   earliest available year greater than or equal to the current calendar year
#'   is used when present; otherwise the first available year is used.
#' @param compare_year Optional second year or time point for the population
#'   pyramid.
#'
#' @return A named list of `ggplot` objects: `population_pyramid`,
#'   `population_projection`, and `age_structure`.
#' @export
plot_demography <- function(
  population,
  year = NULL,
  compare_year = NULL,
  age_groups = NULL,
  age_col = "age_group",
  time_col = "time",
  value_col = "population"
) {
  population_table <- demography_plot_population_table(
    population,
    age_col = age_col,
    time_col = time_col,
    value_col = value_col
  )

  if (is.null(year)) {
    available_years <- sort(unique(population_table[[time_col]]))
    current_year <- as.integer(format(Sys.Date(), "%Y"))
    future_years <- available_years[available_years >= current_year]
    year <- if (length(future_years) > 0) future_years[1] else available_years[1]
  }

  list(
    population_pyramid = plot_population_pyramid(
      population_table,
      year = year,
      compare_year = compare_year,
      age_col = age_col,
      time_col = time_col,
      value_col = value_col
    ),
    population_projection = plot_population_projection(
      population_table,
      age_col = age_col,
      time_col = time_col,
      value_col = value_col
    ),
    age_structure = plot_age_structure(
      population_table,
      age_groups = age_groups,
      age_col = age_col,
      time_col = time_col,
      value_col = value_col
    )
  )
}

require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package ggplot2 is required for demography plotting utilities. ",
      "Install ggplot2 or use the non-plotting demographic helpers.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

demography_plot_population_table <- function(
  population,
  age_col = "age_group",
  time_col = "time",
  value_col = "population"
) {
  if (inherits(population, "agepi_demography")) {
    population <- population$demography
  }

  check_demography_plot_columns(
    population,
    required_columns = c(age_col, time_col, value_col)
  )

  if (!is.numeric(population[[time_col]])) {
    stop("population ", time_col, " must be numeric.", call. = FALSE)
  }
  if (anyNA(population[[time_col]]) || any(!is.finite(population[[time_col]]))) {
    stop("population ", time_col, " must be finite and non-missing.", call. = FALSE)
  }
  if (!is.numeric(population[[value_col]])) {
    stop("population ", value_col, " must be numeric.", call. = FALSE)
  }
  if (anyNA(population[[value_col]]) || any(!is.finite(population[[value_col]]))) {
    stop("population ", value_col, " must be finite and non-missing.", call. = FALSE)
  }
  if (any(population[[value_col]] < 0)) {
    stop("population ", value_col, " cannot be negative.", call. = FALSE)
  }
  if (anyNA(population[[age_col]])) {
    stop("population ", age_col, " cannot contain missing values.", call. = FALSE)
  }

  population[[age_col]] <- as.character(population[[age_col]])
  population
}

check_demography_plot_columns <- function(population, required_columns) {
  if (!is.data.frame(population)) {
    stop("population must be a data frame or agepi_demography object.", call. = FALSE)
  }

  missing_columns <- setdiff(required_columns, names(population))
  if (length(missing_columns) > 0) {
    stop(
      "population is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(population)
}

filter_demography_plot_years <- function(population, years, time_col = "time") {
  if (!is.numeric(years) || anyNA(years) || any(!is.finite(years))) {
    stop("requested year(s) must be finite numeric values.", call. = FALSE)
  }

  years <- unique(years)
  available_years <- sort(unique(population[[time_col]]))
  missing_years <- setdiff(years, available_years)
  if (length(missing_years) > 0) {
    stop(
      "requested year(s) not available in population: ",
      paste(missing_years, collapse = ", "),
      ". Available year(s): ",
      paste(available_years, collapse = ", "),
      call. = FALSE
    )
  }

  selected <- population[population[[time_col]] %in% years, , drop = FALSE]
  row.names(selected) <- NULL
  selected
}

aggregate_total_population <- function(
  population,
  time_col = "time",
  value_col = "population"
) {
  totals <- stats::aggregate(
    population[value_col],
    by = population[time_col],
    FUN = sum
  )
  totals <- totals[order(totals[[time_col]]), , drop = FALSE]
  row.names(totals) <- NULL
  totals
}

aggregate_population_by_time_age <- function(
  population,
  age_col = "age_group",
  time_col = "time",
  value_col = "population"
) {
  totals <- stats::aggregate(
    population[value_col],
    by = population[c(time_col, age_col)],
    FUN = sum
  )
  totals <- totals[order(totals[[time_col]], match(totals[[age_col]], unique(population[[age_col]]))), ]
  row.names(totals) <- NULL
  totals
}

calculate_age_group_proportions <- function(
  population,
  age_groups = NULL,
  age_col = "age_group",
  time_col = "time",
  value_col = "population"
) {
  age_totals <- aggregate_population_by_time_age(
    population,
    age_col = age_col,
    time_col = time_col,
    value_col = value_col
  )

  if (!is.null(age_groups)) {
    age_groups <- as.character(age_groups)
    missing_age_groups <- setdiff(age_groups, unique(age_totals[[age_col]]))
    if (length(missing_age_groups) > 0) {
      stop(
        "requested age_groups not available in population: ",
        paste(missing_age_groups, collapse = ", "),
        call. = FALSE
      )
    }
    age_totals <- age_totals[age_totals[[age_col]] %in% age_groups, , drop = FALSE]
  }

  time_totals <- aggregate_total_population(
    age_totals,
    time_col = time_col,
    value_col = value_col
  )
  names(time_totals)[names(time_totals) == value_col] <- ".total_population"
  proportions <- merge(age_totals, time_totals, by = time_col, sort = FALSE)
  proportions$proportion <- ifelse(
    proportions$.total_population > 0,
    proportions[[value_col]] / proportions$.total_population,
    NA_real_
  )
  proportions <- proportions[
    order(proportions[[time_col]], match(proportions[[age_col]], unique(population[[age_col]]))),
    c(time_col, age_col, value_col, "proportion")
  ]
  row.names(proportions) <- NULL
  proportions
}
