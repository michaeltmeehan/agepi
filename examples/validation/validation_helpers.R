load_local_agepi <- function() {
  if ("package:agepi" %in% search()) {
    return(invisible(TRUE))
  }
  if (dir.exists("R") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }
  if (dir.exists("R")) {
    invisible(lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), source))
    return(invisible(TRUE))
  }
  if (requireNamespace("agepi", quietly = TRUE)) {
    library(agepi)
    return(invisible(TRUE))
  }
  stop(
    "Package agepi is not installed. Run this script from the package root ",
    "or install agepi first.",
    call. = FALSE
  )
}

reference_results <- function(reference) {
  if (is.list(reference) && is.data.frame(reference$results)) {
    return(reference$results)
  }
  if (is.data.frame(reference) || is.matrix(reference)) {
    return(as.data.frame(reference))
  }
  stop("Could not find a comparable trajectory table in the reference output.", call. = FALSE)
}

wide_agepi_output <- function(output, name_map = NULL) {
  if (is.null(name_map)) {
    value_name <- output$compartment
  } else {
    value_name <- name_map[paste(output$compartment, output$age_group, sep = ":")]
  }
  keep <- !is.na(value_name)
  output <- output[keep, , drop = FALSE]
  value_name <- value_name[keep]
  output_key <- paste(output$time, value_name, sep = ":")
  if (anyDuplicated(output_key)) {
    stop(
      "agepi output contains duplicate rows for the same time and variable; ",
      "supply a name_map for multi-age comparisons.",
      call. = FALSE
    )
  }

  wide <- reshape(
    data.frame(
      time = output$time,
      variable = value_name,
      value = output$value,
      stringsAsFactors = FALSE
    ),
    idvar = "time",
    timevar = "variable",
    direction = "wide"
  )
  names(wide) <- sub("^value[.]", "", names(wide))
  wide[order(wide$time), , drop = FALSE]
}

compare_shared_columns <- function(reference, candidate, variables) {
  missing_reference <- setdiff(c("time", variables), names(reference))
  missing_candidate <- setdiff(c("time", variables), names(candidate))
  if (length(missing_reference) > 0 || length(missing_candidate) > 0) {
    stop(
      "Missing comparison columns. Reference: ",
      paste(missing_reference, collapse = ", "),
      "; candidate: ",
      paste(missing_candidate, collapse = ", "),
      call. = FALSE
    )
  }

  reference <- reference[, c("time", variables), drop = FALSE]
  candidate <- candidate[, c("time", variables), drop = FALSE]
  reference <- reference[order(reference$time), , drop = FALSE]
  candidate <- candidate[order(candidate$time), , drop = FALSE]
  if (
    nrow(reference) != nrow(candidate) ||
      any(abs(reference$time - candidate$time) > sqrt(.Machine$double.eps))
  ) {
    stop("Reference and agepi outputs must use the same aligned time grid.", call. = FALSE)
  }

  data.frame(
    variable = variables,
    max_abs_difference = vapply(
      variables,
      function(variable) {
        max(abs(reference[[variable]] - candidate[[variable]]))
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}
