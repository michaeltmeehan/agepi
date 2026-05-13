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

stop_non_exact_target_bin <- function(target_age_group) {
  stop(
    "to_age_structure age group '",
    target_age_group,
    "' cannot be represented as an exact union of from_age_structure age bins.",
    call. = FALSE
  )
}
