library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-17", "18-64", "65+"),
  lower_bounds = c(0, 5, 18, 65),
  upper_bounds = c(4, 17, 64, Inf)
)

contact_matrix <- diag(ages$n_age_groups)

model <- CompartmentModel(
  compartments = c("S", "V", "I", "R"),
  infection_transitions = data.frame(
    from = c("S", "V"),
    to = c("I", "I"),
    susceptibility = I(list(
      1,
      c("0-4" = 0.3, "5-17" = 0.4, "18-64" = 0.5, "65+" = 0.7)
    )),
    stringsAsFactors = FALSE
  ),
  transitions = data.frame(
    from = "I",
    to = "R",
    rate = 0.25,
    stringsAsFactors = FALSE
  ),
  infectious_compartments = "I"
)

initial_state <- data.frame(
  compartment = rep(c("S", "V", "I", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 4),
  value = c(
    c(900, 1000, 1100, 1200),
    c(300, 400, 500, 600),
    c(5, 5, 5, 5),
    c(0, 0, 0, 0)
  ),
  stringsAsFactors = FALSE
)

output <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 4, by = 0.5),
  model = model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  method = "euler",
  cumulative_flows = data.frame(
    name = c("infections_from_S", "breakthrough_infections"),
    from = c("S", "V"),
    to = c("I", "I"),
    stringsAsFactors = FALSE
  )
)

final_time <- max(output$trajectory$time)
final_cumulative <- output$cumulative[output$cumulative$time == final_time, ]
summary_table <- aggregate(value ~ cumulative_name + age_group, final_cumulative, sum)
print(summary_table, row.names = FALSE)

summary_by_source <- aggregate(value ~ time + cumulative_name, output$cumulative, sum)
summary_by_source <- summary_by_source[order(summary_by_source$time), ]
summary_by_source <- reshape(
  summary_by_source,
  idvar = "time",
  timevar = "cumulative_name",
  direction = "wide"
)

matplot(
  summary_by_source$time,
  summary_by_source[, -1],
  type = "l",
  lwd = 2,
  lty = 1,
  col = c("steelblue", "firebrick"),
  xlab = "Time",
  ylab = "Cumulative infections"
)
legend(
  "topleft",
  legend = c("from S", "from V"),
  col = c("steelblue", "firebrick"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

