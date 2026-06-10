validate_cumulative_flows <- function(cumulative_flows, transition_rate_table) {
  if (is.null(cumulative_flows)) {
    return(data.frame(
      cumulative_name = character(),
      transition_id = character(),
      from = character(),
      to = character(),
      stringsAsFactors = FALSE
    ))
  }

  validate_cumulative_transition_rate_table(transition_rate_table)
  normalized <- normalize_cumulative_flows(cumulative_flows)

  rows <- vector("list", nrow(normalized))
  for (i in seq_len(nrow(normalized))) {
    rows[[i]] <- match_cumulative_flow_transition(
      cumulative_name = normalized$cumulative_name[i],
      from = normalized$from[i],
      to = normalized$to[i],
      transition_rate_table = transition_rate_table
    )
  }

  result <- do.call(rbind, rows)
  row.names(result) <- NULL
  result
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
  validate_cumulative_flow_list_names(flow_names)

  rows <- vector("list", length(cumulative_flows))
  for (i in seq_along(cumulative_flows)) {
    flow <- cumulative_flows[[i]]
    if (!is.list(flow) || is.data.frame(flow)) {
      stop("each cumulative flow must be a list with from and to fields.", call. = FALSE)
    }
    if (is.null(flow$from) || is.null(flow$to)) {
      stop("each cumulative flow must include from and to fields.", call. = FALSE)
    }

    from <- validate_cumulative_flow_compartments(flow$from, "from")
    to <- validate_cumulative_flow_compartments(flow$to, "to")
    if (length(from) != length(to)) {
      stop("cumulative flow from and to must have the same length.", call. = FALSE)
    }

    rows[[i]] <- data.frame(
      cumulative_name = rep(flow_names[i], length(from)),
      from = from,
      to = to,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

normalize_cumulative_flow_data_frame <- function(cumulative_flows) {
  required_columns <- c("name", "from", "to")
  missing_columns <- setdiff(required_columns, names(cumulative_flows))
  if (length(missing_columns) > 0) {
    stop(
      "cumulative_flows data frame is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  cumulative_name <- validate_cumulative_flow_character_column(
    cumulative_flows$name,
    "name"
  )
  validate_cumulative_flow_names(cumulative_name)

  data.frame(
    cumulative_name = cumulative_name,
    from = validate_cumulative_flow_character_column(
      cumulative_flows$from,
      "from"
    ),
    to = validate_cumulative_flow_character_column(
      cumulative_flows$to,
      "to"
    ),
    stringsAsFactors = FALSE
  )
}

validate_cumulative_flow_list_names <- function(cumulative_names) {
  validate_cumulative_flow_names(cumulative_names)
}

validate_cumulative_flow_names <- function(cumulative_names) {
  if (anyNA(cumulative_names) || any(cumulative_names == "")) {
    stop("cumulative flow names must be present and non-empty.", call. = FALSE)
  }

  duplicated_names <- unique(cumulative_names[duplicated(cumulative_names)])
  if (length(duplicated_names) > 0) {
    stop(
      "cumulative flow names must be unique; duplicate name(s): ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(cumulative_names)
}

validate_cumulative_flow_compartments <- function(x, field) {
  if (!is.character(x) || length(x) < 1 || anyNA(x) || any(x == "")) {
    stop("cumulative flow ", field, " must contain non-empty character value(s).", call. = FALSE)
  }

  x
}

validate_cumulative_flow_character_column <- function(x, field) {
  if (!is.character(x) || anyNA(x) || any(x == "")) {
    stop("cumulative_flows data frame ", field, " must contain non-empty character values.", call. = FALSE)
  }

  x
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
      anyNA(transition_rate_table$to) ||
      any(transition_rate_table$transition_id == "") ||
      any(transition_rate_table$from == "") ||
      any(transition_rate_table$to == "")) {
    stop("transition_rate_table transition_id, from, and to cannot contain missing or empty values.", call. = FALSE)
  }

  invisible(transition_rate_table)
}

match_cumulative_flow_transition <- function(cumulative_name, from, to, transition_rate_table) {
  known_from <- unique(transition_rate_table$from)
  if (!from %in% known_from) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' references unknown source compartment: ",
      from,
      call. = FALSE
    )
  }

  known_to <- unique(transition_rate_table$to)
  if (!to %in% known_to) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' references unknown destination compartment: ",
      to,
      call. = FALSE
    )
  }

  matches <- transition_rate_table[
    transition_rate_table$from == from & transition_rate_table$to == to,
    c("transition_id", "from", "to"),
    drop = FALSE
  ]
  if (nrow(matches) == 0) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' does not match a declared transition: ",
      from,
      "->",
      to,
      call. = FALSE
    )
  }

  logical_matches <- unique(matches)
  if (length(unique(logical_matches$transition_id)) != 1) {
    stop(
      "cumulative flow '",
      cumulative_name,
      "' is ambiguous for transition: ",
      from,
      "->",
      to,
      call. = FALSE
    )
  }

  data.frame(
    cumulative_name = cumulative_name,
    transition_id = logical_matches$transition_id[1],
    from = logical_matches$from[1],
    to = logical_matches$to[1],
    stringsAsFactors = FALSE
  )
}

cumulative_flow_state_spec <- function(flows, age_groups) {
  flow_names <- unique(flows$cumulative_name)
  age_order <- data.frame(
    age_group = age_groups,
    .age_order = seq_along(age_groups),
    stringsAsFactors = FALSE
  )

  rows <- vector("list", length(flow_names))
  transition_sets <- vector("list", length(flow_names) * length(age_groups))
  transition_index <- 1

  for (i in seq_along(flow_names)) {
    flow_rows <- flows[flows$cumulative_name == flow_names[i], , drop = FALSE]
    rows[[i]] <- data.frame(
      cumulative_name = flow_names[i],
      transition_id = paste(flow_rows$transition_id, collapse = ","),
      from = paste(flow_rows$from, collapse = ","),
      to = paste(flow_rows$to, collapse = ","),
      age_group = age_groups,
      .flow_order = i,
      .age_order = age_order$.age_order,
      stringsAsFactors = FALSE
    )

    for (j in seq_along(age_groups)) {
      transition_sets[[transition_index]] <- flow_rows$transition_id
      transition_index <- transition_index + 1
    }
  }

  state_order <- do.call(rbind, rows)
  state_order <- state_order[order(state_order$.flow_order, state_order$.age_order), ]
  row.names(state_order) <- NULL

  list(
    state_order = state_order[, c(
      "cumulative_name",
      "transition_id",
      "from",
      "to",
      "age_group"
    )],
    transition_sets = transition_sets
  )
}
