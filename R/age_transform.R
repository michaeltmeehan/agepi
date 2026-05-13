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
