#' Build an exact age-grid mapping
#'
#' @param from_age_groups Source age structure.
#' @param to_age_groups Target age structure.
#' @param open_ended How to handle open-ended age groups.
#'
#' @return An `AgeGridMapping` object.
#' @export
AgeGridMapping <- function(from_age_groups, to_age_groups, open_ended = c("include", "error")) {
  open_ended <- match.arg(open_ended)

  validate_age_structure(from_age_groups)
  validate_age_structure(to_age_groups)
  validate_age_grid_mapping_labels(from_age_groups, "from_age_groups")
  validate_age_grid_mapping_labels(to_age_groups, "to_age_groups")

  if (open_ended == "error" && any(is.infinite(c(from_age_groups$upper_bounds, to_age_groups$upper_bounds)))) {
    stop(
      "AgeGridMapping() requires open_ended = 'include' when either age grid has an open-ended age group.",
      call. = FALSE
    )
  }

  from_to_target_index <- age_grid_source_to_target_index(from_age_groups, to_age_groups)
  target_to_source_index <- age_grid_source_to_target_index(to_age_groups, from_age_groups)

  source_indices_by_target <- age_grid_indices_by_target(
    from_age_groups,
    to_age_groups,
    from_to_target_index,
    "target"
  )
  target_indices_by_source <- age_grid_indices_by_target(
    to_age_groups,
    from_age_groups,
    target_to_source_index,
    "source"
  )

  can_aggregate <- all(!is.na(from_to_target_index)) &&
    all(lengths(source_indices_by_target) > 0)
  can_expand <- all(!is.na(target_to_source_index)) &&
    all(lengths(target_indices_by_source) > 0)

  if (!can_aggregate && !can_expand) {
    stop(
      "AgeGridMapping() requires age grids that support exact nested aggregation or expansion.",
      call. = FALSE
    )
  }

  mapping <- list(
    from_age_groups = from_age_groups,
    to_age_groups = to_age_groups,
    open_ended = open_ended,
    from_to_target_index = from_to_target_index,
    source_indices_by_target = source_indices_by_target,
    target_to_source_index = target_to_source_index,
    target_indices_by_source = target_indices_by_source,
    can_aggregate = can_aggregate,
    can_expand = can_expand
  )

  class(mapping) <- "AgeGridMapping"
  mapping
}

#' Aggregate age-specific counts across an age-grid mapping
#'
#' @param x Data frame containing `age_group` and count values.
#' @param mapping An `AgeGridMapping` object.
#' @param value_col Count column name.
#' @param allow_negative Whether negative count values are allowed.
#'
#' @return Aggregated count table.
#' @export
aggregate_age_counts <- function(x, mapping, value_col = "value", allow_negative = FALSE) {
  validate_age_grid_mapping(mapping)
  if (!isTRUE(mapping$can_aggregate)) {
    stop("mapping does not support aggregation from source age groups to target age groups.", call. = FALSE)
  }

  x <- validate_age_count_table(x, mapping$from_age_groups, value_col, allow_negative = allow_negative)
  source_index <- match(as.character(x$age_group), mapping$from_age_groups$age_groups)
  x$age_group <- mapping$to_age_groups$age_groups[mapping$from_to_target_index[source_index]]

  sum_age_count_table(x, value_col)
}

#' Aggregate a demography trajectory to a reporting age grid
#'
#' @param x Data frame with `time`, `age_group`, and `population`.
#' @param mapping An `AgeGridMapping` object.
#'
#' @return Aggregated demography trajectory.
#' @export
aggregate_demography_trajectory_age_grid <- function(x, mapping) {
  aggregate_output_age_grid(
    x = x,
    mapping = mapping,
    required_columns = c("time", "age_group", "population"),
    value_col = "population",
    output_name = "demography trajectory"
  )
}

#' Aggregate an epidemic trajectory to a reporting age grid
#'
#' @param x Data frame with `time`, `compartment`, `age_group`, and `value`.
#' @param mapping An `AgeGridMapping` object.
#'
#' @return Aggregated epidemic trajectory.
#' @export
aggregate_epidemic_trajectory_age_grid <- function(x, mapping) {
  aggregate_output_age_grid(
    x = x,
    mapping = mapping,
    required_columns = c("time", "compartment", "age_group", "value"),
    value_col = "value",
    output_name = "epidemic trajectory"
  )
}

#' Aggregate an age-group population summary to a reporting age grid
#'
#' @param x Data frame with `time`, `age_group`, and `value`.
#' @param mapping An `AgeGridMapping` object.
#'
#' @return Aggregated population summary.
#' @export
aggregate_population_summary_age_grid <- function(x, mapping) {
  aggregate_output_age_grid(
    x = x,
    mapping = mapping,
    required_columns = c("time", "age_group", "value"),
    value_col = "value",
    output_name = "population summary"
  )
}

#' Aggregate cumulative-flow output to a reporting age grid
#'
#' @param x Cumulative-flow data frame.
#' @param mapping An `AgeGridMapping` object.
#'
#' @return Aggregated cumulative-flow data frame.
#' @export
aggregate_cumulative_flows_age_grid <- function(x, mapping) {
  aggregate_output_age_grid(
    x = x,
    mapping = mapping,
    required_columns = c(
      "time",
      "cumulative_name",
      "transition_id",
      "from",
      "to",
      "age_group",
      "value"
    ),
    value_col = "value",
    output_name = "cumulative flows"
  )
}

expand_age_counts <- function(x, mapping, value_col = "value", method = c("uniform"), allow_negative = FALSE) {
  method <- match.arg(method)
  validate_age_grid_mapping(mapping)
  if (!isTRUE(mapping$can_expand)) {
    stop("mapping does not support expansion from source age groups to target age groups.", call. = FALSE)
  }

  x <- validate_age_count_table(x, mapping$from_age_groups, value_col, allow_negative = allow_negative)
  grouping_cols <- setdiff(names(x), c("age_group", value_col))

  expanded_rows <- vector("list", nrow(x))
  for (row_index in seq_len(nrow(x))) {
    source_index <- match(as.character(x$age_group[row_index]), mapping$from_age_groups$age_groups)
    target_indices <- mapping$target_indices_by_source[[source_index]]
    row <- x[rep(row_index, length(target_indices)), , drop = FALSE]
    row$age_group <- mapping$to_age_groups$age_groups[target_indices]
    row[[value_col]] <- row[[value_col]] / length(target_indices)
    expanded_rows[[row_index]] <- row
  }

  expanded <- do.call(rbind, expanded_rows)
  row.names(expanded) <- NULL
  expanded[, c(grouping_cols, "age_group", value_col), drop = FALSE]
}

expand_age_rates <- function(x, mapping, value_col = "value") {
  expand_age_rate_values(
    x = x,
    mapping = mapping,
    value_col = value_col,
    quantity = "rate"
  )
}

expand_age_hazards <- function(x, mapping, value_col = "value") {
  expand_age_rate_values(
    x = x,
    mapping = mapping,
    value_col = value_col,
    quantity = "hazard"
  )
}

expand_age_interval_probabilities <- function(x, mapping, value_col = "value") {
  expand_age_rate_values(
    x = x,
    mapping = mapping,
    value_col = value_col,
    quantity = "interval_probability"
  )
}

expand_age_rate_values <- function(x, mapping, value_col, quantity) {
  validate_age_grid_mapping(mapping)
  if (!isTRUE(mapping$can_expand)) {
    stop("mapping does not support expansion from source age groups to target age groups.", call. = FALSE)
  }

  x <- validate_age_rate_value_table(x, mapping$from_age_groups, value_col, quantity)
  grouping_cols <- setdiff(names(x), c("age_group", value_col))

  expanded_rows <- vector("list", nrow(x))
  for (row_index in seq_len(nrow(x))) {
    source_index <- match(as.character(x$age_group[row_index]), mapping$from_age_groups$age_groups)
    target_indices <- mapping$target_indices_by_source[[source_index]]
    row <- x[rep(row_index, length(target_indices)), , drop = FALSE]
    row$age_group <- mapping$to_age_groups$age_groups[target_indices]
    row[[value_col]] <- expand_one_age_rate_value(
      value = x[[value_col]][row_index],
      mapping = mapping,
      source_index = source_index,
      target_indices = target_indices,
      quantity = quantity
    )
    expanded_rows[[row_index]] <- row
  }

  expanded <- do.call(rbind, expanded_rows)
  row.names(expanded) <- NULL
  expanded[, c(grouping_cols, "age_group", value_col), drop = FALSE]
}

expand_one_age_rate_value <- function(value, mapping, source_index, target_indices, quantity) {
  if (quantity %in% c("rate", "hazard")) {
    return(rep(value, length(target_indices)))
  }

  if (length(target_indices) == 1) {
    return(value)
  }

  source_width <- age_grid_interval_width(mapping$from_age_groups, source_index)
  target_widths <- vapply(
    target_indices,
    function(target_index) age_grid_interval_width(mapping$to_age_groups, target_index),
    numeric(1)
  )

  1 - (1 - value)^(target_widths / source_width)
}

age_grid_interval_width <- function(age_structure, age_index) {
  width <- age_structure$upper_bounds[age_index] - age_structure$lower_bounds[age_index] + 1
  if (!is.finite(width) || width <= 0) {
    stop(
      "interval_probability expansion requires finite age-bin widths when a source age group is split.",
      call. = FALSE
    )
  }

  width
}

validate_age_grid_mapping <- function(mapping) {
  if (!inherits(mapping, "AgeGridMapping")) {
    stop("mapping must be an AgeGridMapping object.", call. = FALSE)
  }

  required_fields <- c(
    "from_age_groups",
    "to_age_groups",
    "open_ended",
    "from_to_target_index",
    "source_indices_by_target",
    "target_to_source_index",
    "target_indices_by_source",
    "can_aggregate",
    "can_expand"
  )
  missing_fields <- setdiff(required_fields, names(mapping))
  if (length(missing_fields) > 0) {
    stop(
      "mapping is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  validate_age_structure(mapping$from_age_groups)
  validate_age_structure(mapping$to_age_groups)
  invisible(mapping)
}

validate_age_grid_mapping_labels <- function(age_structure, argument_name) {
  parsed <- lapply(age_structure$age_groups, parse_age_grid_label)
  parsed_lower <- vapply(parsed, `[[`, numeric(1), "lower")
  parsed_upper <- vapply(parsed, `[[`, numeric(1), "upper")

  if (
    any(parsed_lower != age_structure$lower_bounds) ||
      any(!age_grid_upper_bounds_equal(parsed_upper, age_structure$upper_bounds))
  ) {
    stop(
      argument_name,
      " age_group labels must parse to their explicit lower_bounds and upper_bounds.",
      call. = FALSE
    )
  }

  invisible(age_structure)
}

parse_age_grid_label <- function(label) {
  if (!is.character(label) || length(label) != 1 || is.na(label) || !nzchar(label)) {
    stop("age_group labels must be non-empty character values.", call. = FALSE)
  }

  if (grepl("^[0-9]+$", label)) {
    lower <- as.numeric(label)
    return(list(lower = lower, upper = lower))
  }

  if (grepl("^[0-9]+\\+$", label)) {
    return(list(
      lower = as.numeric(sub("\\+$", "", label)),
      upper = Inf
    ))
  }

  if (grepl("^[0-9]+-[0-9]+$", label)) {
    bounds <- strsplit(label, "-", fixed = TRUE)[[1]]
    return(list(
      lower = as.numeric(bounds[1]),
      upper = as.numeric(bounds[2])
    ))
  }

  stop("age_group labels must use forms like '0', '0-4', or '100+'.", call. = FALSE)
}

age_grid_upper_bounds_equal <- function(x, y) {
  (is.infinite(x) & is.infinite(y)) | x == y
}

age_grid_source_to_target_index <- function(source_age_groups, target_age_groups) {
  vapply(
    seq_len(source_age_groups$n_age_groups),
    function(source_index) {
      target_indices <- which(
        target_age_groups$lower_bounds <= source_age_groups$lower_bounds[source_index] &
          target_age_groups$upper_bounds >= source_age_groups$upper_bounds[source_index]
      )

      if (length(target_indices) == 0) {
        return(NA_integer_)
      }

      if (length(target_indices) != 1) {
        stop(
          "source age group '",
          source_age_groups$age_groups[source_index],
          "' maps to more than one target age group.",
          call. = FALSE
        )
      }

      target_indices
    },
    integer(1)
  )
}

age_grid_indices_by_target <- function(source_age_groups, target_age_groups, source_to_target_index, target_kind) {
  indices <- vector("list", target_age_groups$n_age_groups)
  for (target_index in seq_len(target_age_groups$n_age_groups)) {
    source_indices <- which(source_to_target_index == target_index)
    indices[[target_index]] <- source_indices

    if (length(source_indices) == 0) {
      next
    }

    source_lower <- source_age_groups$lower_bounds[source_indices]
    source_upper <- source_age_groups$upper_bounds[source_indices]
    if (
      source_lower[1] != target_age_groups$lower_bounds[target_index] ||
        !age_grid_upper_bounds_equal(source_upper[length(source_upper)], target_age_groups$upper_bounds[target_index])
    ) {
      stop_non_exact_age_grid_mapping(target_kind, target_age_groups$age_groups[target_index])
    }

    if (length(source_indices) > 1) {
      expected_next_lower <- source_upper[-length(source_upper)] + 1
      if (any(source_lower[-1] != expected_next_lower)) {
        stop_non_exact_age_grid_mapping(target_kind, target_age_groups$age_groups[target_index])
      }
    }
  }

  indices
}

stop_non_exact_age_grid_mapping <- function(target_kind, age_group) {
  stop(
    target_kind,
    " age group '",
    age_group,
    "' is not exactly covered by nested source age groups.",
    call. = FALSE
  )
}

aggregate_output_age_grid <- function(x, mapping, required_columns, value_col, output_name) {
  if (!is.data.frame(x)) {
    stop(output_name, " must be a data frame.", call. = FALSE)
  }

  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0) {
    stop(
      output_name,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  aggregate_age_counts(x, mapping, value_col = value_col)
}

validate_age_count_table <- function(x, age_structure, value_col, allow_negative = FALSE) {
  if (!is.data.frame(x)) {
    stop("x must be a data frame.", call. = FALSE)
  }

  if (!is.character(value_col) || length(value_col) != 1 || is.na(value_col) || !nzchar(value_col)) {
    stop("value_col must be a single non-empty string.", call. = FALSE)
  }

  required_columns <- c("age_group", value_col)
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0) {
    stop(
      "x is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(x[[value_col]])) {
    stop("x value column must be numeric.", call. = FALSE)
  }

  if (anyNA(x[[value_col]]) || any(!is.finite(x[[value_col]]))) {
    stop("x value column must be finite and non-missing.", call. = FALSE)
  }

  if (!is.logical(allow_negative) || length(allow_negative) != 1 || is.na(allow_negative)) {
    stop("allow_negative must be TRUE or FALSE.", call. = FALSE)
  }

  if (!allow_negative && any(x[[value_col]] < 0)) {
    stop("x value column cannot be negative.", call. = FALSE)
  }

  age_group <- as.character(x$age_group)
  extra_age_groups <- setdiff(unique(age_group), age_structure$age_groups)
  if (length(extra_age_groups) > 0) {
    stop(
      "x contains age_group value(s) not in mapping source age groups: ",
      paste(extra_age_groups, collapse = ", "),
      call. = FALSE
    )
  }

  grouping_cols <- setdiff(names(x), c("age_group", value_col))
  validate_age_count_grouping_columns(x, grouping_cols)
  duplicate_key <- age_count_group_key(x, c(grouping_cols, "age_group"))
  if (any(duplicated(duplicate_key))) {
    stop("x contains duplicate rows for the same non-age columns and age_group.", call. = FALSE)
  }

  x$age_group <- age_group
  x
}

validate_age_rate_value_table <- function(x, age_structure, value_col, quantity) {
  if (!is.data.frame(x)) {
    stop("x must be a data frame.", call. = FALSE)
  }

  if (!is.character(value_col) || length(value_col) != 1 || is.na(value_col) || !nzchar(value_col)) {
    stop("value_col must be a single non-empty string.", call. = FALSE)
  }

  required_columns <- c("age_group", value_col)
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0) {
    stop(
      "x is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(x[[value_col]])) {
    stop("x value column must be numeric.", call. = FALSE)
  }

  if (anyNA(x[[value_col]]) || any(!is.finite(x[[value_col]]))) {
    stop("x value column must be finite and non-missing.", call. = FALSE)
  }

  if (quantity == "hazard" && any(x[[value_col]] < 0)) {
    stop("hazard values cannot be negative.", call. = FALSE)
  }

  if (quantity == "interval_probability" && any(x[[value_col]] < 0 | x[[value_col]] > 1)) {
    stop("interval_probability values must be between 0 and 1.", call. = FALSE)
  }

  age_group <- as.character(x$age_group)
  extra_age_groups <- setdiff(unique(age_group), age_structure$age_groups)
  if (length(extra_age_groups) > 0) {
    stop(
      "x contains age_group value(s) not in mapping source age groups: ",
      paste(extra_age_groups, collapse = ", "),
      call. = FALSE
    )
  }

  grouping_cols <- setdiff(names(x), c("age_group", value_col))
  validate_age_count_grouping_columns(x, grouping_cols)
  duplicate_key <- age_count_group_key(x, c(grouping_cols, "age_group"))
  if (any(duplicated(duplicate_key))) {
    stop("x contains duplicate rows for the same non-age columns and age_group.", call. = FALSE)
  }

  x$age_group <- age_group
  x
}

validate_age_count_grouping_columns <- function(x, grouping_cols) {
  unsupported <- grouping_cols[vapply(x[grouping_cols], function(col) is.list(col) || is.data.frame(col), logical(1))]
  if (length(unsupported) > 0) {
    stop(
      "x contains unsupported non-age column(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(x)
}

sum_age_count_table <- function(x, value_col) {
  grouping_cols <- setdiff(names(x), c("age_group", value_col))
  output_cols <- c(grouping_cols, "age_group", value_col)
  group_key <- age_count_group_key(x, c(grouping_cols, "age_group"))
  group_levels <- unique(group_key)
  group_index <- match(group_key, group_levels)
  first_rows <- match(group_levels, group_key)

  out <- x[first_rows, output_cols, drop = FALSE]
  out[[value_col]] <- as.numeric(tapply(
    x[[value_col]],
    factor(group_index, levels = seq_along(group_levels)),
    sum
  ))
  row.names(out) <- NULL
  out
}

age_count_group_key <- function(x, cols) {
  if (length(cols) == 0) {
    return(rep("", nrow(x)))
  }

  do.call(paste, c(x[cols], list(sep = "\r")))
}
