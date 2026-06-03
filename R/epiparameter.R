#' Convert an epiparameter delay to a Markov transition rate
#'
#' `rate_from_epiparameter()` extracts the mean delay from an
#' `<epiparameter>` object and returns the corresponding per-capita Markov
#' transition rate, `1 / mean_delay`.
#'
#' This is an exponential waiting-time approximation: it collapses the full
#' delay distribution represented by the `<epiparameter>` object to a single
#' constant hazard with the same mean. It does not preserve the shape, variance,
#' truncation, or other distributional features of the original delay. Erlang
#' or gamma dwell-time approximations may be added in future work.
#'
#' @param x An `<epiparameter>` object, for example one returned by
#'   `epiparameter::epiparameter_db()`.
#'
#' @return A finite positive numeric scalar equal to `1 / mean(x)`.
#' @examples
#' if (requireNamespace("epiparameter", quietly = TRUE)) {
#'   incubation <- epiparameter::epiparameter_db(
#'     disease = "COVID-19",
#'     epi_name = "incubation period",
#'     single_epiparameter = TRUE
#'   )
#'   SEIRModel(
#'     sigma = rate_from_epiparameter(incubation),
#'     gamma = 0.25
#'   )
#' }
#' @export
rate_from_epiparameter <- function(x) {
  if (!inherits(x, "epiparameter")) {
    stop("x must be an <epiparameter> object.", call. = FALSE)
  }

  if (!requireNamespace("epiparameter", quietly = TRUE)) {
    stop(
      "The epiparameter package is required to use rate_from_epiparameter(). ",
      "Install it to convert <epiparameter> objects.",
      call. = FALSE
    )
  }

  mean_delay <- mean(x)

  if (!is.numeric(mean_delay) || length(mean_delay) != 1 || anyNA(mean_delay) ||
      !is.finite(mean_delay)) {
    stop("The epiparameter mean delay must be a finite numeric scalar.", call. = FALSE)
  }

  if (mean_delay <= 0) {
    stop("The epiparameter mean delay must be positive.", call. = FALSE)
  }

  1 / as.numeric(mean_delay)
}
