# agepi

`agepi` provides tools for deterministic age-structured epidemic and
demographic modelling in R.

The package is under active development. Its current focus is on clear,
inspectable workflows for age-structured SIR, SEIR, custom compartment models,
contact matrices, and demographic schedules.

## Overview

`agepi` helps you build models where people are grouped by age and where
transmission depends on contacts between age groups. It currently supports:

- defining age groups with `AgeStructure()`;
- running deterministic SIR and SEIR simulations with `simulate_deterministic()`;
- tracking selected deterministic transition flows as cumulative outputs;
- running fixed-population Gillespie SIR, SEIR, and supported generic
  compartment simulations with `simulate_stochastic()`;
- summarising selected stochastic transition flows from realised event logs;
- using age-specific contact matrices, susceptibility, and infectiousness;
- building custom deterministic compartment models with `CompartmentModel()`;
- simulating demographic-only processes with ageing, fertility, mortality, and
  migration;
- using dependency-free adapters for WPP-style demographic tables and
  socialmixr/conmat-style contact matrix outputs;
- summarising simulation output with `compartment_totals()`,
  `age_group_totals()`, and `total_population()`.

The package is intended for researchers, analysts, students, and modellers who
want a small R toolkit for prototyping age-structured epidemic and demographic
models before moving to a larger modelling framework.

## Installation

Install the development version from a local checkout:

```r
install.packages("devtools")
devtools::install(".")
```

Or, if this repository is available on GitHub:

```r
remotes::install_github("michaeltmeehan/agepi")
```

`deSolve` is suggested, not required. When it is installed, infection-only SIR
and SEIR simulations can use `deSolve::ode()`; otherwise explicit Euler
stepping is available.

## Quick Start

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

population <- c(1000, 1200)
infected <- c(5, 3)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 3),
  value = c(population - infected, infected, rep(0, ages$n_age_groups)),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(
  c(4, 2,
    2, 5),
  nrow = ages$n_age_groups,
  byrow = TRUE
)

model <- SIRModel(gamma = 0.25)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 1, by = 0.1),
  model = model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  method = "euler"
)

head(simulation)
compartment_totals(simulation)
```

## Age-Structured Transmission

Age-structured transmission combines:

- an `AgeStructure()` describing the age groups;
- a contact matrix whose rows are recipient age groups and columns are
  infectious/source age groups;
- the number infectious in each source age group;
- optional age-specific susceptibility and infectiousness vectors;
- a transmission scaling parameter, `beta`.

A slightly richer SIR example with three age groups and age-specific
susceptibility/infectiousness is:

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

contact_matrix <- matrix(c(
  4, 2, 1,
  2, 5, 2,
  1, 2, 4
), nrow = ages$n_age_groups, byrow = TRUE)

population <- c(1000, 1200, 900)
infected <- c(5, 3, 2)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 3),
  value = c(population - infected, infected, rep(0, ages$n_age_groups)),
  stringsAsFactors = FALSE
)

sir <- SIRModel(gamma = 0.25)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 2, by = 0.1),
  model = sir,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  susceptibility = c(0.8, 1.0, 1.2),
  infectiousness = c(1.1, 1.0, 0.9),
  method = "euler"
)

tail(compartment_totals(simulation))
```

See [examples/mock_sir_deterministic.R](examples/mock_sir_deterministic.R) for
a compact script version.

## Custom Compartment Models

Use `CompartmentModel()` when you want an age-structured model that is not one
of the built-in `SIRModel()` or `SEIRModel()` constructors. Deterministic
simulation supports generic models through the transition-rate interface, and
stochastic simulation supports the fixed-population subset that can be
represented as within-age transitions through `transition_rates()`.
For example, this defines an MSIR-style model where maternal immunity wanes
from `M` to `S`, susceptible people are infected, and infectious people recover:

```r
msir <- CompartmentModel(
  compartments = c("M", "S", "I", "R"),
  infection_transitions = data.frame(from = "S", to = "I"),
  transitions = data.frame(
    from = c("M", "I"),
    to = c("S", "R"),
    rate = c(0.15, 0.25)
  ),
  infectious_compartments = "I",
  birth_compartment = "M",
  migration_compartment = "S"
)
```

Generic models can also declare scalar or named age-specific per-capita
transition rates, competing outgoing transitions from the same source
compartment, and more than one infectious source compartment. Use named
`infectiousness_weights` for relative infectiousness by compartment; the force
of infection still uses source-age totals across all compartments in the
denominator and recipient-row/source-column contact matrices.

```r
transitions <- data.frame(
  from = c("E", "E", "IP", "IC", "IS"),
  to = c("IP", "IS", "IC", "R", "R")
)
transitions$rate <- I(list(
  c("0-4" = 0.20, "5-9" = 0.35),
  c("0-4" = 0.30, "5-9" = 0.15),
  0.25,
  0.10,
  c("0-4" = 0.30, "5-9" = 0.40)
))

covid_like <- CompartmentModel(
  compartments = c("S", "E", "IP", "IC", "IS", "R"),
  infection_transitions = data.frame(from = "S", to = "E"),
  transitions = transitions,
  infectious_compartments = c("IP", "IC", "IS"),
  infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
)
```

This is only a small expressiveness fixture; it is not a full Davies-style
COVID model and does not add branching-probability, R0, or setting-specific
contact tools.

See [docs/generic_compartment_models.md](docs/generic_compartment_models.md),
[examples/generic_sir.R](examples/generic_sir.R),
[examples/generic_seir.R](examples/generic_seir.R), and
[examples/generic_msir.R](examples/generic_msir.R).

## Demography

`agepi` includes demographic-only workflows for age-specific population
dynamics. You can define fertility, mortality, migration, and ageing processes,
then simulate them with `simulate_demography()`.

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

mortality <- MortalitySchedule(
  data.frame(
    time = 0,
    age_group = ages$age_groups,
    mortality_rate = c(0.006, 0.004, 0.020)
  ),
  ages
)

process <- build_demographic_process(
  age_structure = ages,
  mortality_schedule = mortality
)

simulate_demography(
  process = process,
  initial_state = c(500, 450, 800),
  times = c(0, 1),
  time_policy = "step"
)
```

See [examples/mock_demographic_workflow.R](examples/mock_demographic_workflow.R)
and [docs/demographic_residuals.md](docs/demographic_residuals.md) for a fuller
demographic workflow and diagnostics.

## Contact Matrices And External Data

Contact matrices can be supplied directly as numeric matrices. The current
convention is:

```text
contact_matrix[recipient_age_group, source_age_group]
```

The package also includes dependency-free adapters for common external shapes:

- `contact_matrix_from_socialmixr()` for socialmixr-like objects with a numeric
  `matrix` element;
- `contact_matrix_from_conmat()` for conmat-style long data frames;
- `population_from_wpp()` and `demography_from_wpp()` for WPP-style population
  tables;
- WPP-style fertility, mortality, and migration standardisers.

These adapters reshape and validate supplied data; they do not attempt to
reproduce full external projection systems. See
[docs/external_data_adapters.md](docs/external_data_adapters.md) and
[docs/contact_matrix_integration_design.md](docs/contact_matrix_integration_design.md).

## Examples

- [examples/mock_sir_deterministic.R](examples/mock_sir_deterministic.R):
  deterministic age-structured SIR.
- [examples/generic_sir.R](examples/generic_sir.R): SIR through
  `CompartmentModel()`.
- [examples/generic_seir.R](examples/generic_seir.R): SEIR through
  `CompartmentModel()`.
- [examples/generic_msir.R](examples/generic_msir.R): custom MSIR model.
- [examples/mock_seir_demography.R](examples/mock_seir_demography.R): SEIR with
  demographic turnover.
- [examples/mock_demographic_workflow.R](examples/mock_demographic_workflow.R):
  demographic-only workflow with diagnostics.

## Documentation

- [docs/model_conventions.md](docs/model_conventions.md): state-vector,
  force-of-infection, contact-matrix, and demographic schedule conventions.
- [docs/development_status.md](docs/development_status.md): current scope,
  limitations, solver notes, adapter boundaries, and roadmap notes.
- [docs/generic_compartment_models.md](docs/generic_compartment_models.md):
  custom compartment models.
- [docs/age_structured_transmission_design.md](docs/age_structured_transmission_design.md):
  age-structured transmission design notes.
- [docs/external_data_adapters.md](docs/external_data_adapters.md): WPP-style
  demographic and contact-matrix adapter notes.
- [docs/demographic_residuals.md](docs/demographic_residuals.md): demographic
  residual diagnostics.

## Development Status

`agepi` is an early-stage package. The current release is best suited for
deterministic prototype models, teaching examples, and development of modelling
workflows. SIR, SEIR, generic deterministic compartment models,
fixed-population stochastic SIR, SEIR, and supported generic compartment
simulation, demographic processes, and adapter utilities are implemented, but
broader stochastic methods, stochastic demography, fitting/calibration, plotting
helpers, time-varying contact matrices, and full WPP projection matching are not
yet in scope.

For detailed technical limitations and implementation conventions, see
[docs/development_status.md](docs/development_status.md) and
[docs/model_conventions.md](docs/model_conventions.md).

## Contributing And Reporting Issues

Bug reports, example requests, and documentation improvements are welcome. When
reporting an issue, please include a small reproducible example, the package
version, and whether `deSolve` is installed.
