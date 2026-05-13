#' Aggregate an age-indexed vector to a coarser age structure
#'
#' Transforms numeric values from one age structure to another when every
#' target age bin is an exact union of complete source age bins.
#'
#' @param values Numeric vector with one value per source age group.
#' @param from_age_structure Source age structure.
#' @param to_age_structure Target age structure.
#'
#' @return Numeric vector with one value per target age group, named with
#'   `to_age_structure$age_groups`.
#' @export
aggregate_age_vector <- function(values, from_age_structure, to_age_structure) {
  validate_age_structure(from_age_structure)
  validate_age_structure(to_age_structure)

  if (!is.numeric(values)) {
    stop("values must be numeric.", call. = FALSE)
  }

  if (length(values) != from_age_structure$n_age_groups) {
    stop(
      "values length must equal from_age_structure$n_age_groups: ",
      from_age_structure$n_age_groups,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(values) || any(!is.finite(values))) {
    stop("values must be finite and cannot contain missing values.", call. = FALSE)
  }

  aggregated <- vapply(
    seq_len(to_age_structure$n_age_groups),
    function(i) {
      source_indices <- exact_source_indices_for_target_bin(
        from_age_structure,
        to_age_structure$lower_bounds[i],
        to_age_structure$upper_bounds[i],
        to_age_structure$age_groups[i]
      )

      sum(values[source_indices])
    },
    numeric(1)
  )

  names(aggregated) <- to_age_structure$age_groups
  aggregated
}

#' Segregate an age-indexed vector to a finer age structure
#'
#' Transforms numeric values from one age structure to another when every
#' target age bin is fully nested inside exactly one source age bin. Source
#' values are allocated to target bins in proportion to explicit target
#' weights within each source age bin.
#'
#' @param values Numeric vector with one value per source age group.
#' @param from_age_structure Source age structure.
#' @param to_age_structure Target age structure.
#' @param weights Numeric non-negative vector with one weight per target age
#'   group.
#'
#' @return Numeric vector with one value per target age group, named with
#'   `to_age_structure$age_groups`.
#' @export
segregate_age_vector <- function(values, from_age_structure, to_age_structure, weights) {
  validate_age_structure(from_age_structure)
  validate_age_structure(to_age_structure)

  if (!is.numeric(values)) {
    stop("values must be numeric.", call. = FALSE)
  }

  if (length(values) != from_age_structure$n_age_groups) {
    stop(
      "values length must equal from_age_structure$n_age_groups: ",
      from_age_structure$n_age_groups,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(values) || any(!is.finite(values))) {
    stop("values must be finite and cannot contain missing values.", call. = FALSE)
  }

  if (!is.numeric(weights)) {
    stop("weights must be numeric.", call. = FALSE)
  }

  if (length(weights) != to_age_structure$n_age_groups) {
    stop(
      "weights length must equal to_age_structure$n_age_groups: ",
      to_age_structure$n_age_groups,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("weights must be finite and cannot contain missing values.", call. = FALSE)
  }

  if (any(weights < 0)) {
    stop("weights cannot contain negative values.", call. = FALSE)
  }

  target_source_indices <- source_index_for_nested_target_bins(from_age_structure, to_age_structure)

  segregated <- numeric(to_age_structure$n_age_groups)
  for (source_index in seq_len(from_age_structure$n_age_groups)) {
    target_indices <- which(target_source_indices == source_index)
    validate_target_coverage_for_source_bin(
      from_age_structure,
      to_age_structure,
      source_index,
      target_indices
    )

    source_weight_total <- sum(weights[target_indices])
    if (source_weight_total == 0) {
      stop(
        "weights for from_age_structure age group '",
        from_age_structure$age_groups[source_index],
        "' must sum to a positive value.",
        call. = FALSE
      )
    }

    segregated[target_indices] <- values[source_index] * weights[target_indices] / source_weight_total
  }

  names(segregated) <- to_age_structure$age_groups
  segregated
}

#' Transform an age-indexed vector between age structures
#'
#' Convenience wrapper for exact age-vector transformations.
#' Supports identity, aggregation, segregation, and mixed exact
#' transformations when source and target bins align to a common set of
#' boundaries.
#'
#' @param values Numeric vector with one value per source age group.
#' @param from_age_structure Source age structure.
#' @param to_age_structure Target age structure.
#' @param weights Optional numeric non-negative vector with one weight per
#'   target age group.
#' @param direction Transformation direction.
#' @param split_method How to allocate source bins that are split across
#'   multiple target bins.
#'
#' @return Numeric vector with one value per target age group, named with
#'   `to_age_structure$age_groups`.
#' @export
transform_age_vector <- function(
  values,
  from_age_structure,
  to_age_structure,
  weights = NULL,
  direction = c("auto", "aggregate", "segregate"),
  split_method = c("width", "weights", "equal", "error")
) {
  direction <- match.arg(direction)
  split_method <- match.arg(split_method)

  validate_age_structure(from_age_structure)
  validate_age_structure(to_age_structure)

  if (direction == "aggregate") {
    if (!is.null(weights)) {
      stop("weights must not be supplied for aggregate age transformations.", call. = FALSE)
    }

    return(aggregate_age_vector(values, from_age_structure, to_age_structure))
  }

  if (direction == "segregate") {
    if (is.null(weights)) {
      stop("weights must be supplied for segregate age transformations.", call. = FALSE)
    }

    return(segregate_age_vector(values, from_age_structure, to_age_structure, weights))
  }

  validate_age_transform_values(values, from_age_structure)
  if (!is.null(weights)) {
    validate_age_transform_weights(weights, to_age_structure)
  }
  if (split_method == "weights" && is.null(weights)) {
    stop("weights must be supplied when split_method = 'weights'.", call. = FALSE)
  }

  transform_age_vector_exact(
    values,
    from_age_structure,
    to_age_structure,
    weights,
    split_method
  )
}

exact_source_indices_for_target_bin <- function(from_age_structure, target_lower, target_upper, target_age_group) {
  source_indices <- which(
    from_age_structure$lower_bounds >= target_lower &
      from_age_structure$upper_bounds <= target_upper
  )

  if (length(source_indices) == 0) {
    stop_non_exact_target_bin(target_age_group)
  }

  selected_lower <- from_age_structure$lower_bounds[source_indices]
  selected_upper <- from_age_structure$upper_bounds[source_indices]

  if (selected_lower[1] != target_lower || selected_upper[length(selected_upper)] != target_upper) {
    stop_non_exact_target_bin(target_age_group)
  }

  if (length(source_indices) > 1) {
    expected_next_lower <- selected_upper[-length(selected_upper)] + 1
    if (any(selected_lower[-1] != expected_next_lower)) {
      stop_non_exact_target_bin(target_age_group)
    }
  }

  source_indices
}

validate_age_transform_values <- function(values, from_age_structure) {
  if (!is.numeric(values)) {
    stop("values must be numeric.", call. = FALSE)
  }

  if (length(values) != from_age_structure$n_age_groups) {
    stop(
      "values length must equal from_age_structure$n_age_groups: ",
      from_age_structure$n_age_groups,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(values) || any(!is.finite(values))) {
    stop("values must be finite and cannot contain missing values.", call. = FALSE)
  }

  invisible(values)
}

validate_age_transform_weights <- function(weights, to_age_structure) {
  if (!is.numeric(weights)) {
    stop("weights must be numeric.", call. = FALSE)
  }

  if (length(weights) != to_age_structure$n_age_groups) {
    stop(
      "weights length must equal to_age_structure$n_age_groups: ",
      to_age_structure$n_age_groups,
      ".",
      call. = FALSE
    )
  }

  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("weights must be finite and cannot contain missing values.", call. = FALSE)
  }

  if (any(weights < 0)) {
    stop("weights cannot contain negative values.", call. = FALSE)
  }

  invisible(weights)
}

can_aggregate_age_structures <- function(from_age_structure, to_age_structure) {
  tryCatch(
    {
      for (i in seq_len(to_age_structure$n_age_groups)) {
        exact_source_indices_for_target_bin(
          from_age_structure,
          to_age_structure$lower_bounds[i],
          to_age_structure$upper_bounds[i],
          to_age_structure$age_groups[i]
        )
      }

      TRUE
    },
    error = function(e) FALSE
  )
}

transform_age_vector_exact <- function(values, from_age_structure, to_age_structure, weights, split_method) {
  validate_exact_age_transform_coverage(from_age_structure, to_age_structure)

  transformed <- numeric(to_age_structure$n_age_groups)
  for (source_index in seq_len(from_age_structure$n_age_groups)) {
    target_indices <- overlapping_target_indices(from_age_structure, to_age_structure, source_index)
    split_proportions <- source_split_proportions(
      from_age_structure,
      to_age_structure,
      source_index,
      target_indices,
      weights,
      split_method
    )

    transformed[target_indices] <- transformed[target_indices] +
      values[source_index] * split_proportions
  }

  names(transformed) <- to_age_structure$age_groups
  transformed
}

validate_exact_age_transform_coverage <- function(from_age_structure, to_age_structure) {
  for (source_index in seq_len(from_age_structure$n_age_groups)) {
    validate_overlapping_bins_cover_interval(
      from_age_structure,
      to_age_structure,
      from_age_structure$lower_bounds[source_index],
      from_age_structure$upper_bounds[source_index],
      overlapping_target_indices(from_age_structure, to_age_structure, source_index),
      "from_age_structure",
      from_age_structure$age_groups[source_index]
    )
  }

  for (target_index in seq_len(to_age_structure$n_age_groups)) {
    validate_overlapping_bins_cover_interval(
      to_age_structure,
      from_age_structure,
      to_age_structure$lower_bounds[target_index],
      to_age_structure$upper_bounds[target_index],
      overlapping_source_indices(from_age_structure, to_age_structure, target_index),
      "to_age_structure",
      to_age_structure$age_groups[target_index]
    )
  }

  invisible(TRUE)
}

overlapping_target_indices <- function(from_age_structure, to_age_structure, source_index) {
  which(
    to_age_structure$lower_bounds <= from_age_structure$upper_bounds[source_index] &
      to_age_structure$upper_bounds >= from_age_structure$lower_bounds[source_index]
  )
}

overlapping_source_indices <- function(from_age_structure, to_age_structure, target_index) {
  which(
    from_age_structure$lower_bounds <= to_age_structure$upper_bounds[target_index] &
      from_age_structure$upper_bounds >= to_age_structure$lower_bounds[target_index]
  )
}

validate_overlapping_bins_cover_interval <- function(interval_structure, covering_structure, interval_lower, interval_upper, covering_indices, interval_kind, interval_age_group) {
  if (length(covering_indices) == 0) {
    stop_incompatible_exact_transform(interval_kind, interval_age_group)
  }

  covering_lower <- covering_structure$lower_bounds[covering_indices]
  covering_upper <- covering_structure$upper_bounds[covering_indices]

  if (covering_lower[1] > interval_lower || covering_upper[length(covering_upper)] < interval_upper) {
    stop_incompatible_exact_transform(interval_kind, interval_age_group)
  }

  if (length(covering_indices) > 1) {
    expected_next_lower <- covering_upper[-length(covering_upper)] + 1
    if (any(covering_lower[-1] > expected_next_lower)) {
      stop_incompatible_exact_transform(interval_kind, interval_age_group)
    }
  }

  invisible(interval_structure)
}

source_split_proportions <- function(from_age_structure, to_age_structure, source_index, target_indices, weights, split_method) {
  if (length(target_indices) == 1) {
    return(1)
  }

  if (split_method == "error" && is.null(weights)) {
    stop(
      "split_method = 'error' rejects source age bins that require splitting unless weights are supplied.",
      call. = FALSE
    )
  }

  if (split_method == "width") {
    if (is.infinite(from_age_structure$upper_bounds[source_index])) {
      stop("Cannot split open-ended age bins using split_method = 'width'.", call. = FALSE)
    }

    split_weights <- pmin(to_age_structure$upper_bounds[target_indices], from_age_structure$upper_bounds[source_index]) -
      pmax(to_age_structure$lower_bounds[target_indices], from_age_structure$lower_bounds[source_index]) + 1
  } else if (split_method == "equal") {
    split_weights <- rep(1, length(target_indices))
  } else {
    if (is.null(weights)) {
      stop("weights must be supplied for weighted age-bin splitting.", call. = FALSE)
    }

    split_weights <- weights[target_indices]
  }

  split_weight_total <- sum(split_weights)
  if (split_weight_total == 0) {
    stop(
      "split weights for from_age_structure age group '",
      from_age_structure$age_groups[source_index],
      "' must sum to a positive value.",
      call. = FALSE
    )
  }

  split_weights / split_weight_total
}

source_index_for_nested_target_bins <- function(from_age_structure, to_age_structure) {
  vapply(
    seq_len(to_age_structure$n_age_groups),
    function(i) {
      source_indices <- which(
        from_age_structure$lower_bounds <= to_age_structure$lower_bounds[i] &
          from_age_structure$upper_bounds >= to_age_structure$upper_bounds[i]
      )

      if (length(source_indices) != 1) {
        stop_non_nested_target_bin(to_age_structure$age_groups[i])
      }

      source_indices
    },
    integer(1)
  )
}

can_segregate_age_structures <- function(from_age_structure, to_age_structure) {
  tryCatch(
    {
      target_source_indices <- source_index_for_nested_target_bins(from_age_structure, to_age_structure)

      for (source_index in seq_len(from_age_structure$n_age_groups)) {
        validate_target_coverage_for_source_bin(
          from_age_structure,
          to_age_structure,
          source_index,
          which(target_source_indices == source_index)
        )
      }

      TRUE
    },
    error = function(e) FALSE
  )
}

validate_target_coverage_for_source_bin <- function(from_age_structure, to_age_structure, source_index, target_indices) {
  if (length(target_indices) == 0) {
    stop_non_exact_source_coverage(from_age_structure$age_groups[source_index])
  }

  target_lower <- to_age_structure$lower_bounds[target_indices]
  target_upper <- to_age_structure$upper_bounds[target_indices]

  if (
    target_lower[1] != from_age_structure$lower_bounds[source_index] ||
      target_upper[length(target_upper)] != from_age_structure$upper_bounds[source_index]
  ) {
    stop_non_exact_source_coverage(from_age_structure$age_groups[source_index])
  }

  if (length(target_indices) > 1) {
    expected_next_lower <- target_upper[-length(target_upper)] + 1
    if (any(target_lower[-1] != expected_next_lower)) {
      stop_non_exact_source_coverage(from_age_structure$age_groups[source_index])
    }
  }

  invisible(target_indices)
}

stop_non_exact_target_bin <- function(target_age_group) {
  stop(
    "to_age_structure age group '",
    target_age_group,
    "' cannot be represented as an exact union of from_age_structure age bins.",
    call. = FALSE
  )
}

stop_non_nested_target_bin <- function(target_age_group) {
  stop(
    "to_age_structure age group '",
    target_age_group,
    "' must be fully contained inside exactly one from_age_structure age bin.",
    call. = FALSE
  )
}

stop_non_exact_source_coverage <- function(source_age_group) {
  stop(
    "from_age_structure age group '",
    source_age_group,
    "' must be exactly covered by one or more to_age_structure age bins.",
    call. = FALSE
  )
}

stop_incompatible_exact_transform <- function(interval_kind, interval_age_group) {
  stop(
    interval_kind,
    " age group '",
    interval_age_group,
    "' is not compatible with an exact age transformation.",
    call. = FALSE
  )
}
