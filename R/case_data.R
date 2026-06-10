#' Map observed case data to agepi age groups
#'
#' Prepares a user-supplied case or event data frame for aggregate agepi
#' workflows by validating or deriving an `age_group` column against an
#' [AgeStructure()]. This helper only maps observed records into model age bins;
#' it does not create individual-level line lists from agepi compartment
#' trajectories.
#'
#' @param cases A `data.frame` with one row per observed case or event.
#' @param age_structure An [AgeStructure()] object.
#' @param age_col Optional name of a numeric exact-age column in `cases`.
#' @param age_group_col Optional name of an existing age-group column in `cases`.
#'
#' @return A `data.frame` containing all original columns plus a validated
#'   `age_group` column ordered as `age_structure$age_groups`. If `cases` is a
#'   `linelist` object and the `linelist` package is installed, existing
#'   linelist tags are preserved where possible.
#' @export
#'
#' @examples
#' ages <- AgeStructure(c("0-4", "5-9", "10+"), c(0, 5, 10), c(4, 9, Inf))
#' observed_cases <- data.frame(
#'   case_id = c("a", "b", "c", "d"),
#'   onset_day = c(1, 1, 2, 2),
#'   age = c(2, 7, 12, 41)
#' )
#'
#' mapped_cases <- as_agepi_cases(
#'   observed_cases,
#'   age_structure = ages,
#'   age_col = "age"
#' )
#'
#' aggregate(case_id ~ onset_day + age_group, mapped_cases, length)
as_agepi_cases <- function(
  cases,
  age_structure,
  age_col = NULL,
  age_group_col = NULL
) {
  if (!is.data.frame(cases)) {
    stop("cases must be a data.frame.", call. = FALSE)
  }
  validate_age_structure(age_structure)

  selectors <- c(!is.null(age_col), !is.null(age_group_col))
  if (sum(selectors) != 1) {
    stop("Specify exactly one of age_col or age_group_col.", call. = FALSE)
  }

  input_was_linelist <- inherits(cases, "linelist")
  linelist_tags <- NULL
  if (input_was_linelist && requireNamespace("linelist", quietly = TRUE)) {
    linelist_tags <- linelist::tags(cases)
  }

  result <- as.data.frame(cases, stringsAsFactors = FALSE)

  if (!is.null(age_col)) {
    validate_case_column(result, age_col, "age_col")
    result$age_group <- age_groups_from_numeric_age(
      result[[age_col]],
      age_structure,
      age_col
    )
  } else {
    validate_case_column(result, age_group_col, "age_group_col")
    result$age_group <- validate_case_age_groups(
      result[[age_group_col]],
      age_structure,
      age_group_col
    )
  }

  result$age_group <- factor(
    result$age_group,
    levels = age_structure$age_groups,
    ordered = TRUE
  )

  restore_linelist_cases(result, input_was_linelist, linelist_tags)
}

validate_case_column <- function(cases, column, argument) {
  if (!is.character(column) || length(column) != 1 || is.na(column) || column == "") {
    stop(argument, " must be a non-empty character scalar.", call. = FALSE)
  }
  if (!column %in% names(cases)) {
    stop(argument, " '", column, "' was not found in cases.", call. = FALSE)
  }
  invisible(column)
}

age_groups_from_numeric_age <- function(age, age_structure, age_col) {
  if (!is.numeric(age)) {
    stop("age_col '", age_col, "' must be numeric.", call. = FALSE)
  }
  if (anyNA(age)) {
    stop("age_col '", age_col, "' contains missing age values.", call. = FALSE)
  }
  if (any(!is.finite(age))) {
    stop("age_col '", age_col, "' contains non-finite age values.", call. = FALSE)
  }

  age_group <- rep(NA_character_, length(age))
  for (i in seq_along(age_structure$age_groups)) {
    in_group <- age >= age_structure$lower_bounds[[i]] &
      age <= age_structure$upper_bounds[[i]]
    age_group[in_group] <- age_structure$age_groups[[i]]
  }

  if (anyNA(age_group)) {
    invalid <- unique(age[is.na(age_group)])
    stop(
      "age_col '", age_col, "' contains age value(s) outside age_structure bounds: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }

  age_group
}

validate_case_age_groups <- function(age_group, age_structure, age_group_col) {
  if (anyNA(age_group)) {
    stop("age_group_col '", age_group_col, "' contains missing age-group values.", call. = FALSE)
  }

  age_group <- as.character(age_group)
  invalid <- unique(age_group[!age_group %in% age_structure$age_groups])
  if (length(invalid) > 0) {
    stop(
      "age_group_col '", age_group_col,
      "' contains invalid age-group label(s): ",
      paste(invalid, collapse = ", "),
      ". Expected labels are: ",
      paste(age_structure$age_groups, collapse = ", "),
      call. = FALSE
    )
  }

  age_group
}

restore_linelist_cases <- function(result, input_was_linelist, linelist_tags) {
  if (!input_was_linelist || !requireNamespace("linelist", quietly = TRUE)) {
    return(result)
  }

  if (is.null(linelist_tags)) {
    return(result)
  }

  linelist_tags <- linelist_tags[unlist(linelist_tags, use.names = FALSE) %in% names(result)]
  args <- c(list(x = result, allow_extra = TRUE), linelist_tags)
  do.call(linelist::make_linelist, args)
}
