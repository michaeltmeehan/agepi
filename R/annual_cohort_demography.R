annual_cohort_demographic_step <- function(population,
                                           age_structure,
                                           fertility,
                                           mortality,
                                           migration = NULL,
                                           fertility_exposure_fraction = 1,
                                           migration_type = NULL,
                                           birth_age_group = NULL,
                                           open_ended = c("terminal"),
                                           tolerance = sqrt(.Machine$double.eps)) {
  open_ended <- match.arg(open_ended)
  validate_annual_cohort_age_structure(age_structure)
  fertility_exposure_fraction <- validate_fertility_exposure_fraction(fertility_exposure_fraction)
  validate_annual_cohort_tolerance(tolerance)

  age_groups <- age_structure$age_groups
  n_age_groups <- age_structure$n_age_groups

  if (is.null(birth_age_group)) {
    birth_age_group <- age_groups[1]
  }
  birth_index <- match_annual_cohort_birth_age_group(birth_age_group, age_groups)

  population <- annual_cohort_values_by_age(
    population,
    age_structure,
    value_name = "population",
    allow_negative = FALSE,
    require_full_coverage = TRUE
  )
  fertility <- annual_cohort_values_by_age(
    fertility,
    age_structure,
    value_name = "fertility_rate",
    allow_negative = FALSE,
    require_full_coverage = FALSE,
    fill_value = 0
  )
  mortality <- annual_cohort_values_by_age(
    mortality,
    age_structure,
    value_name = "mortality_rate",
    allow_negative = FALSE,
    require_full_coverage = TRUE
  )
  migration_info <- annual_cohort_migration_by_age(
    migration,
    age_structure,
    migration_type = migration_type
  )

  survivors <- population * exp(-mortality)
  deaths <- population - survivors

  aged_population <- numeric(n_age_groups)
  if (n_age_groups == 1) {
    aged_population[1] <- survivors[1]
  } else {
    aged_population[seq.int(2, n_age_groups)] <- survivors[seq_len(n_age_groups - 1)]
    aged_population[n_age_groups] <- aged_population[n_age_groups] + survivors[n_age_groups]
  }
  names(aged_population) <- age_groups

  births <- sum(fertility * fertility_exposure_fraction * population)
  aged_population[birth_index] <- aged_population[birth_index] + births

  if (identical(migration_info$type, "rate")) {
    net_migration <- migration_info$values * aged_population
  } else {
    net_migration <- migration_info$values
  }

  next_population <- aged_population + net_migration
  validate_annual_cohort_next_population(next_population, tolerance)
  next_population[next_population < 0] <- 0

  result <- data.frame(
    age_group = age_groups,
    population = as.numeric(next_population),
    stringsAsFactors = FALSE
  )
  attr(result, "diagnostics") <- list(
    births = births,
    deaths = sum(deaths),
    net_migration = sum(net_migration),
    total_population_before = sum(population),
    total_population_after = sum(next_population),
    mortality_convention = "annual_hazard_survival_exp_minus_m",
    fertility_exposure_fraction = fertility_exposure_fraction,
    migration_type = migration_info$type,
    open_ended_policy = open_ended
  )

  result
}

validate_annual_cohort_age_structure <- function(age_structure) {
  validate_age_structure(age_structure)
  validate_contiguous_age_structure(age_structure)

  n_age_groups <- age_structure$n_age_groups
  if (n_age_groups > 1) {
    finite_widths <- age_structure$upper_bounds[-n_age_groups] -
      age_structure$lower_bounds[-n_age_groups] + 1
    if (any(finite_widths != 1)) {
      stop(
        "annual cohort demographic step requires 1-year finite age groups.",
        call. = FALSE
      )
    }
  }

  invisible(age_structure)
}

validate_annual_cohort_tolerance <- function(tolerance) {
  if (!is.numeric(tolerance) || length(tolerance) != 1 || anyNA(tolerance) || !is.finite(tolerance)) {
    stop("tolerance must be a single finite numeric value.", call. = FALSE)
  }

  if (tolerance < 0) {
    stop("tolerance must be non-negative.", call. = FALSE)
  }

  invisible(tolerance)
}

match_annual_cohort_birth_age_group <- function(birth_age_group, age_groups) {
  if (!is.character(birth_age_group) || length(birth_age_group) != 1 || is.na(birth_age_group) || !nzchar(birth_age_group)) {
    stop("birth_age_group must be a single non-empty age_group value.", call. = FALSE)
  }

  birth_index <- match(birth_age_group, age_groups)
  if (is.na(birth_index)) {
    stop("birth_age_group must be present in age_structure.", call. = FALSE)
  }

  birth_index
}

annual_cohort_values_by_age <- function(x,
                                        age_structure,
                                        value_name,
                                        allow_negative,
                                        require_full_coverage,
                                        fill_value = NA_real_) {
  age_groups <- age_structure$age_groups

  if (is.data.frame(x)) {
    return(annual_cohort_table_values_by_age(
      x = x,
      age_structure = age_structure,
      value_name = value_name,
      allow_negative = allow_negative,
      require_full_coverage = require_full_coverage,
      fill_value = fill_value
    ))
  }

  if (!is.numeric(x) || is.matrix(x)) {
    stop(value_name, " must be a numeric vector or data frame.", call. = FALSE)
  }

  if (length(x) != age_structure$n_age_groups) {
    stop(value_name, " length must equal the number of age groups.", call. = FALSE)
  }

  if (!is.null(names(x)) && length(names(x)) > 0 && !identical(names(x), age_groups)) {
    stop(value_name, " names must be complete and ordered as age_structure$age_groups.", call. = FALSE)
  }

  values <- as.numeric(x)
  validate_annual_cohort_numeric_values(values, value_name, allow_negative)
  stats::setNames(values, age_groups)
}

annual_cohort_table_values_by_age <- function(x,
                                              age_structure,
                                              value_name,
                                              allow_negative,
                                              require_full_coverage,
                                              fill_value) {
  required_columns <- c("age_group", value_name)
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0) {
    stop(
      value_name,
      " data is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  age_group <- as.character(x$age_group)
  values <- x[[value_name]]
  validate_annual_cohort_numeric_values(values, value_name, allow_negative)

  unknown_age_groups <- setdiff(unique(age_group), age_structure$age_groups)
  if (length(unknown_age_groups) > 0) {
    stop(
      value_name,
      " contains age_group value(s) not in age_structure: ",
      paste(unknown_age_groups, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(duplicated(age_group))) {
    stop(value_name, " contains duplicate age_group rows.", call. = FALSE)
  }

  if (require_full_coverage) {
    missing_age_groups <- setdiff(age_structure$age_groups, age_group)
    if (length(missing_age_groups) > 0) {
      stop(
        value_name,
        " is missing required age_group value(s): ",
        paste(missing_age_groups, collapse = ", "),
        call. = FALSE
      )
    }

    if (!identical(age_group, age_structure$age_groups)) {
      stop(value_name, " age_group rows must be complete and ordered as age_structure$age_groups.", call. = FALSE)
    }
  } else {
    age_positions <- match(age_group, age_structure$age_groups)
    if (is.unsorted(age_positions, strictly = TRUE)) {
      stop(value_name, " age_group rows must be ordered as age_structure$age_groups.", call. = FALSE)
    }
  }

  result <- rep(fill_value, age_structure$n_age_groups)
  result[match(age_group, age_structure$age_groups)] <- values
  stats::setNames(result, age_structure$age_groups)
}

validate_annual_cohort_numeric_values <- function(values, value_name, allow_negative) {
  if (!is.numeric(values)) {
    stop(value_name, " values must be numeric.", call. = FALSE)
  }

  if (anyNA(values) || any(!is.finite(values))) {
    stop(value_name, " values must be finite and non-missing.", call. = FALSE)
  }

  if (!allow_negative && any(values < 0)) {
    stop(value_name, " values must be non-negative.", call. = FALSE)
  }

  invisible(values)
}

annual_cohort_migration_by_age <- function(migration, age_structure, migration_type = NULL) {
  if (is.null(migration)) {
    return(list(
      type = "count",
      values = stats::setNames(numeric(age_structure$n_age_groups), age_structure$age_groups)
    ))
  }

  migration_type <- infer_annual_cohort_migration_type(migration, migration_type)
  value_name <- paste0("migration_", migration_type)

  values <- annual_cohort_values_by_age(
    migration,
    age_structure,
    value_name = value_name,
    allow_negative = TRUE,
    require_full_coverage = TRUE
  )

  list(type = migration_type, values = values)
}

infer_annual_cohort_migration_type <- function(migration, migration_type) {
  if (is.null(migration_type)) {
    if (is.data.frame(migration)) {
      has_rate <- "migration_rate" %in% names(migration)
      has_count <- "migration_count" %in% names(migration)
      if (has_rate && has_count) {
        stop("migration must supply exactly one of migration_rate or migration_count, not both.", call. = FALSE)
      }
      if (has_rate) {
        return("rate")
      }
      if (has_count) {
        return("count")
      }
    }

    return("count")
  }

  if (!is.character(migration_type) || length(migration_type) != 1 || is.na(migration_type)) {
    stop("migration_type must be a single character value.", call. = FALSE)
  }

  migration_type <- match.arg(migration_type, c("count", "rate"))
  migration_type
}

validate_annual_cohort_next_population <- function(next_population, tolerance) {
  if (any(!is.finite(next_population))) {
    stop("annual cohort demographic step produced non-finite population value.", call. = FALSE)
  }

  negative <- next_population < -tolerance
  if (any(negative)) {
    stop(
      "annual cohort demographic step produced negative population value(s): ",
      paste(names(next_population)[negative], collapse = ", "),
      call. = FALSE
    )
  }

  invisible(next_population)
}
