validate_cumulative_flows <- function(cumulative_flows, transition_rate_table) {
  if (is.null(cumulative_flows)) {
    return(data.frame(
      cumulative_name = character(),
      selector_index = integer(),
      transition_label = character(),
      transition_id = character(),
      transition_type = character(),
      from = character(),
      to = character(),
      stringsAsFactors = FALSE
    ))
  }

  validate_cumulative_transition_rate_table(transition_rate_table)
  normalized <- normalize_cumulative_flows(cumulative_flows)
  resolved <- resolve_cumulative_flow_selectors(normalized, transition_rate_table)
  row.names(resolved) <- NULL
  resolved[, c(
    "cumulative_name",
    "selector_index",
    "transition_label",
    "transition_id",
    "transition_type",
    "from",
    "to"
  ), drop = FALSE]
}

normalize_cumulative_flows <- function(cumulative_flows) {
  if (is.data.frame(cumulative_flows)) {
    return(normalize_cumulative_flow_data_frame(cumulative_flows))
  }

  if (is.list(cumulative_flows)) {
    return(normalize_cumulative_flow_list(cumulative_flows))
  }

  stop("cumulative_flows must be a named list or a data frame.", call. = FALSE)
}

normalize_cumulative_flow_list <- function(cumulative_flows) {
  flow_names <- names(cumulative_flows)
  if (is.null(flow_names) || length(flow_names) != length(cumulative_flows)) {
    stop("cumulative_flows list entries must be named.", call. = FALSE)
  }
  validate_cumulative_flow_names(flow_names)

  rows <- vector("list", length(cumulative_flows))
  for (i in seq_along(cumulative_flows)) {
    rows[[i]] <- normalize_cumulative_flow_entry(
      cumulative_name = flow_names[i],
      selector = cumulative_flows[[i]]
    )
  }

  do.call(rbind, rows)
}

normalize_cumulative_flow_data_frame <- function(cumulative_flows) {
  cumulative_name <- normalize_cumulative_flow_name_column(cumulative_flows)
  validate_cumulative_flow_names(cumulative_name)
  to_supplied <- "to" %in% names(cumulative_flows)

  selector_columns <- c("label", "transition_id", "transition_type", "from", "to")
  if (!any(selector_columns %in% names(cumulative_flows))) {
    stop(
      "cumulative_flows data frame must include at least one selector column: label, transition_id, transition_type, from, or to.",
      call. = FALSE
    )
  }

  normalized <- data.frame(
    cumulative_name = cumulative_name,
    selector_index = ave(seq_len(nrow(cumulative_flows)), cumulative_name, FUN = seq_along),
    transition_label = normalize_cumulative_flow_character_column(
      cumulative_flows,
      "label"
    ),
    transition_id = normalize_cumulative_flow_character_column(
      cumulative_flows,
      "transition_id"
    ),
    transition_type = normalize_cumulative_flow_character_column(
      cumulative_flows,
      "transition_type"
    ),
    from = normalize_cumulative_flow_character_column(
      cumulative_flows,
      "from"
    ),
    to = normalize_cumulative_flow_character_column(
      cumulative_flows,
      "to",
      allow_na = TRUE
    ),
    stringsAsFactors = FALSE
  )

  normalized$transition_type <- normalize_cumulative_flow_transition_type(normalized$transition_type)
  normalized$transition_type[is.na(normalized$transition_type) &
    is.na(normalized$transition_label) &
    is.na(normalized$transition_id) &
    to_supplied &
    is.na(normalized$to)] <- "outflow"

  if (any(!normalize_cumulative_flow_row_has_selector(normalized))) {
    bad_names <- unique(normalized$cumulative_name[!normalize_cumulative_flow_row_has_selector(normalized)])
    stop(
      "cumulative flow(s) must specify at least one selector column: ",
      paste(bad_names, collapse = ", "),
      call. = FALSE
    )
  }

  normalized
}

normalize_cumulative_flow_entry <- function(cumulative_name, selector) {
  if (!is.list(selector) || is.data.frame(selector)) {
    stop("each cumulative flow must be a list with optional selector fields.", call. = FALSE)
  }

  normalized <- normalize_cumulative_flow_selector_fields(
    cumulative_name = cumulative_name,
    selector = selector
  )
  normalized
}

normalize_cumulative_flow_selector_fields <- function(cumulative_name, selector) {
  transition_label <- selector$label
  transition_id <- selector$transition_id
  transition_type <- selector$transition_type
  from <- selector$from
  to <- selector$to
  to_supplied <- "to" %in% names(selector)

  supplied_fields <- list(
    label = transition_label,
    transition_id = transition_id,
    transition_type = transition_type,
    from = from,
    to = to
  )
  for (field in names(supplied_fields)) {
    value <- supplied_fields[[field]]
    if (!is.null(value) && length(value) == 0) {
      stop(
        "cumulative flow '",
        cumulative_name,
        "' selector field '",
        field,
        "' has length 0; expected length 1 or more.",
        call. = FALSE
      )
    }
  }

  supplied_lengths <- c(
    if (!is.null(transition_label)) length(transition_label),
    if (!is.null(transition_id)) length(transition_id),
    if (!is.null(transition_type)) length(transition_type),
    if (!is.null(from)) length(from),
    if (!is.null(to)) length(to)
  )
  supplied_lengths <- supplied_lengths[supplied_lengths > 0]
  if (length(supplied_lengths) == 0) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' must specify at least one selector field.",
      call. = FALSE
    )
  }

  n <- max(supplied_lengths)
  transition_label <- recycle_cumulative_selector_field(
    transition_label,
    "label",
    n,
    cumulative_name
  )
  transition_id <- recycle_cumulative_selector_field(
    transition_id,
    "transition_id",
    n,
    cumulative_name
  )
  transition_type <- recycle_cumulative_selector_field(
    transition_type,
    "transition_type",
    n,
    cumulative_name
  )
  from <- recycle_cumulative_selector_field(from, "from", n, cumulative_name)
  to <- recycle_cumulative_selector_field(to, "to", n, cumulative_name, allow_na = TRUE)

  normalized <- data.frame(
    cumulative_name = rep(cumulative_name, n),
    selector_index = seq_len(n),
    transition_label = transition_label,
    transition_id = transition_id,
    transition_type = normalize_cumulative_flow_transition_type(transition_type),
    from = from,
    to = to,
    stringsAsFactors = FALSE
  )

  normalized$transition_type[is.na(normalized$transition_type) &
    is.na(normalized$transition_label) &
    is.na(normalized$transition_id) &
    to_supplied &
    is.na(normalized$to)] <- "outflow"

  if (!all(normalize_cumulative_flow_row_has_selector(normalized))) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' must specify at least one selector field.",
      call. = FALSE
    )
  }

  normalized
}

recycle_cumulative_selector_field <- function(value, field, n, cumulative_name, allow_na = FALSE) {
  if (is.null(value)) {
    return(rep(NA_character_, n))
  }

  if (is.factor(value)) {
    value <- as.character(value)
  }
  if (!is.character(value)) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' selector field '",
      field,
      "' must be a character vector.",
      call. = FALSE
    )
  }
  if (length(value) == 0) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' selector field '",
      field,
      "' has length 0; expected length 1 or ",
      n,
      ".",
      call. = FALSE
    )
  }
  if (!(length(value) %in% c(1, n))) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' selector field '",
      field,
      "' has length ",
      length(value),
      "; expected length 1 or ",
      n,
      ".",
      call. = FALSE
    )
  }

  if (!allow_na && any(value == "", na.rm = TRUE)) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' selector field '",
      field,
      "' cannot contain empty values.",
      call. = FALSE
    )
  }

  if (allow_na && any(value == "", na.rm = TRUE)) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' selector field '",
      field,
      "' cannot contain empty values.",
      call. = FALSE
    )
  }

  if (length(value) == 1) {
    rep(value, n)
  } else {
    value
  }
}

normalize_cumulative_flow_name_column <- function(cumulative_flows) {
  if ("name" %in% names(cumulative_flows) && "cumulative_name" %in% names(cumulative_flows)) {
    name_values <- normalize_cumulative_flow_character_column(cumulative_flows, "name")
    cumulative_values <- normalize_cumulative_flow_character_column(cumulative_flows, "cumulative_name")
    if (!identical(name_values, cumulative_values)) {
      stop("cumulative_flows name and cumulative_name columns must match when both are supplied.", call. = FALSE)
    }
    return(name_values)
  }

  if ("name" %in% names(cumulative_flows)) {
    return(normalize_cumulative_flow_character_column(cumulative_flows, "name"))
  }

  if ("cumulative_name" %in% names(cumulative_flows)) {
    return(normalize_cumulative_flow_character_column(cumulative_flows, "cumulative_name"))
  }

  stop("cumulative_flows data frame is missing required column(s): name or cumulative_name", call. = FALSE)
}

normalize_cumulative_flow_character_column <- function(cumulative_flows, field, allow_na = FALSE) {
  if (!field %in% names(cumulative_flows)) {
    return(rep(NA_character_, nrow(cumulative_flows)))
  }

  value <- cumulative_flows[[field]]
  if (is.factor(value)) {
    value <- as.character(value)
  }
  if (!is.character(value)) {
    stop("cumulative_flows data frame ", field, " must contain character values.", call. = FALSE)
  }
  if (length(value) != nrow(cumulative_flows)) {
    stop("cumulative_flows data frame ", field, " must have one value per row.", call. = FALSE)
  }
  if (any(value == "", na.rm = TRUE)) {
    stop("cumulative_flows data frame ", field, " must contain non-empty character values.", call. = FALSE)
  }
  if (!allow_na && anyNA(value)) {
    stop("cumulative_flows data frame ", field, " must contain non-empty character values.", call. = FALSE)
  }
  value
}

normalize_cumulative_flow_transition_type <- function(transition_type) {
  if (all(is.na(transition_type))) {
    return(transition_type)
  }

  transition_type <- as.character(transition_type)
  transition_type[transition_type == "internal"] <- "transition"
  invalid_types <- setdiff(unique(transition_type[!is.na(transition_type)]), c("infection", "transition", "outflow"))
  if (length(invalid_types) > 0) {
    stop(
      "cumulative_flows transition_type must be infection, transition, outflow, or internal.",
      call. = FALSE
    )
  }
  transition_type
}

normalize_cumulative_flow_row_has_selector <- function(normalized) {
  !is.na(normalized$transition_label) |
    !is.na(normalized$transition_id) |
    !is.na(normalized$transition_type) |
    !is.na(normalized$from) |
    !is.na(normalized$to)
}

validate_cumulative_transition_rate_table <- function(transition_rate_table) {
  if (!is.data.frame(transition_rate_table)) {
    stop("transition_rate_table must be a data frame.", call. = FALSE)
  }

  required_columns <- c("transition_id", "from", "to")
  missing_columns <- setdiff(required_columns, names(transition_rate_table))
  if (length(missing_columns) > 0) {
    stop(
      "transition_rate_table is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(transition_rate_table$transition_id) ||
      anyNA(transition_rate_table$from) ||
      any(transition_rate_table$transition_id == "") ||
      any(transition_rate_table$from == "")) {
    stop("transition_rate_table transition_id and from cannot contain missing or empty values.", call. = FALSE)
  }

  if (!"transition_label" %in% names(transition_rate_table)) {
    transition_rate_table$transition_label <- rep(NA_character_, nrow(transition_rate_table))
  }
  if (!"transition_type" %in% names(transition_rate_table)) {
    transition_rate_table$transition_type <- ifelse(is.na(transition_rate_table$to), "outflow", transition_type_from_ids(transition_rate_table$transition_id))
  }

  transition_rate_table$transition_type <- normalize_cumulative_flow_transition_type(transition_rate_table$transition_type)

  if (any(!is.na(transition_rate_table$transition_label) & transition_rate_table$transition_label == "")) {
    stop("transition_rate_table transition_label cannot contain empty values.", call. = FALSE)
  }

  outflow_rows <- is.na(transition_rate_table$to)
  if (any(outflow_rows)) {
    outflow_ids <- transition_rate_table$transition_id[outflow_rows]
    if (any(!startsWith(outflow_ids, "outflow:"))) {
      stop("transition_rate_table rows with missing destinations must be outflows.", call. = FALSE)
    }
  }

  invisible(transition_rate_table)
}

transition_rate_logical_table <- function(transition_rate_table) {
  validate_cumulative_transition_rate_table(transition_rate_table)

  if (!"transition_label" %in% names(transition_rate_table)) {
    transition_rate_table$transition_label <- rep(NA_character_, nrow(transition_rate_table))
  }
  if (!"transition_type" %in% names(transition_rate_table)) {
    transition_rate_table$transition_type <- ifelse(is.na(transition_rate_table$to), "outflow", transition_type_from_ids(transition_rate_table$transition_id))
  }

  logical_columns <- c("transition_id", "transition_label", "transition_type", "from", "to")
  transition_ids <- unique(transition_rate_table$transition_id)
  rows <- vector("list", length(transition_ids))

  for (i in seq_along(transition_ids)) {
    transition_id <- transition_ids[i]
    slice <- unique(
      transition_rate_table[transition_rate_table$transition_id == transition_id, logical_columns, drop = FALSE]
    )
    if (nrow(slice) != 1) {
      stop(
        "transition_rate_table contains inconsistent logical transition metadata for transition_id '",
        transition_id,
        "'.",
        call. = FALSE
      )
    }
    rows[[i]] <- slice
  }

  do.call(rbind, rows)
}

resolve_cumulative_flow_selectors <- function(normalized, transition_rate_table) {
  logical_table <- transition_rate_logical_table(transition_rate_table)
  rows <- vector("list", nrow(normalized))

  for (i in seq_len(nrow(normalized))) {
    rows[[i]] <- match_cumulative_flow_selector(
      selector = normalized[i, , drop = FALSE],
      logical_table = logical_table
    )
  }

  do.call(rbind, rows)
}

match_cumulative_flow_selector <- function(selector, logical_table) {
  cumulative_name <- selector$cumulative_name[1]
  selector_index <- selector$selector_index[1]
  matches <- rep(TRUE, nrow(logical_table))

  selector_fields <- c("transition_label", "transition_id", "transition_type", "from", "to")
  for (field in selector_fields) {
    value <- selector[[field]][1]
    if (is.na(value)) {
      next
    }
    if (field == "transition_type" && value == "internal") {
      value <- "transition"
    }
    matches <- matches & logical_table[[field]] == value
  }

  selected <- logical_table[matches, , drop = FALSE]
  if (nrow(selected) == 0) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' selector ",
      selector_index,
      " matched no transitions.\n",
      format_cumulative_selector_message(selector),
      call. = FALSE
    )
  }

  if (nrow(selected) > 1) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' selector ",
      selector_index,
      " matched multiple transitions:\n\n",
      format_cumulative_transition_table(selected),
      "\nAdd label, transition_id, or transition_type to disambiguate.",
      call. = FALSE
    )
  }

  data.frame(
    cumulative_name = cumulative_name,
    selector_index = selector_index,
    transition_label = selected$transition_label[1],
    transition_id = selected$transition_id[1],
    transition_type = selected$transition_type[1],
    from = selected$from[1],
    to = selected$to[1],
    stringsAsFactors = FALSE
  )
}

format_cumulative_selector_message <- function(selector) {
  values <- c(
    transition_label = selector$transition_label[1],
    transition_id = selector$transition_id[1],
    transition_type = selector$transition_type[1],
    from = selector$from[1],
    to = selector$to[1]
  )

  lines <- names(values)[!is.na(values)]
  if (length(lines) == 0) {
    return("  <no selectors supplied>")
  }

  paste0(
    "  ",
    paste(names(values)[!is.na(values)], values[!is.na(values)], sep = " = ", collapse = "\n  ")
  )
}

format_cumulative_transition_table <- function(selected) {
  selected <- selected[, c("transition_label", "transition_id", "transition_type", "from", "to"), drop = FALSE]
  selected$transition_label[is.na(selected$transition_label)] <- "<NA>"
  selected$to[is.na(selected$to)] <- "<NA>"

  header <- sprintf(
    "  %-22s %-16s %-14s %-16s %-16s",
    "transition_label",
    "transition_id",
    "transition_type",
    "from",
    "to"
  )
  rows <- apply(selected, 1, function(row) {
    sprintf("  %-22s %-16s %-14s %-16s %-16s", row[["transition_label"]], row[["transition_id"]], row[["transition_type"]], row[["from"]], row[["to"]])
  })
  paste(c(header, rows), collapse = "\n")
}
