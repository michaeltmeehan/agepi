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
- coupling deterministic epidemic models to demography either with the default
  continuous derivative policy or with opt-in annual-cohort operator splitting;
- using dependency-free adapters for WPP-style demographic tables and
  socialmixr/conmat-style contact matrix outputs;
- mapping observed case/event records into agepi age groups with
  `as_agepi_cases()`;
- summarising simulation output with `compartment_totals()`,
  `age_group_totals()`, and `total_population()`.

The package is intended for researchers, analysts, students, and modellers who
want a small R toolkit for prototyping age-structured epidemic and demographic
models before moving to a larger modelling framework.

## Installation

`agepi` is split into clear installation tiers so you can install only the
pieces you need.

| Tier | Best for | Install |
| --- | --- | --- |
| Minimal core | Core simulation, demography with direct inputs, and other functionality that only needs `Imports` | `pak::pkg_install("michaeltmeehan/agepi", dependencies = c("Depends", "Imports", "LinkingTo"))` |
| Recommended | Plotting, examples, `epiparameter`, `linelist`, `contactdata`, `socialmixr`, and the package test suite | `pak::pkg_install("michaeltmeehan/agepi", dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))` |
| WPP-backed demography | `population_from_wpp()` / `demography_from_wpp()` workflows that rely on the `wpp2024` data package | `pak::pkg_install("PPgp/wpp2024")` or `remotes::install_github("PPgp/wpp2024")` |
| Advanced contact matrices | `load_contact_matrix_source(source = "conmat")` workflows | `install.packages("conmat", repos = c("https://idem-lab.r-universe.dev", "https://cloud.r-project.org"))` or `remotes::install_github("idem-lab/conmat")` |
| Developer/test | Local development, examples, docs, and tests from a checkout | `devtools::install(".")` or `pak::pkg_install(".", dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))` |

The default package metadata keeps `deSolve` in `Imports`, so the core install
stays light. Ordinary optional helpers remain in `Suggests`, while `wpp2024`
and `conmat` are intentionally treated as separate external backends so they do
not block minimal installs or ordinary dependency resolution.

If you are starting from a local checkout, `devtools::install(".")` is the
simplest path. If you are working from GitHub, `remotes::install_github("michaeltmeehan/agepi")`
also works for the core package, and you can add the optional backend installs
above when you need WPP-backed demography or synthetic conmat contacts.

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

## Optional Epiparameter Interoperability

`rate_from_epiparameter()` converts an `<epiparameter>` delay object, such as
one returned by `epiparameter::epiparameter_db()`, to a Markov per-capita
transition rate by taking `1 / mean(delay)`. The `epiparameter` package is
suggested, not imported, and is only required when this helper is called.
Epiverse integrations are optional; install their packages separately when you
want to run the corresponding examples.

This intentionally collapses the full delay distribution into an exponential
waiting-time approximation with the same mean. It does not preserve the
distribution shape, variance, truncation, or other non-Markov dwell-time
features.

```r
if (requireNamespace("epiparameter", quietly = TRUE)) {
  incubation <- epiparameter::epiparameter_db(
    disease = "COVID-19",
    epi_name = "incubation period",
    single_epiparameter = TRUE
  )

  seir <- SEIRModel(
    sigma = rate_from_epiparameter(incubation),
    gamma = 0.25
  )
}
```

See [examples/epiparameter_seir.R](examples/epiparameter_seir.R) for a compact
age-structured SEIR example. Erlang/gamma dwell-time helpers are future work.

## Cumulative Flow Outputs

Selected transition flows can be requested at simulation time. Deterministic
simulation integrates selected transition rates as auxiliary cumulative output:

```r
deterministic_output <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 10, by = 1),
  model = SIRModel(gamma = 0.25),
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  cumulative_flows = list(
    infections = list(from = "S", to = "I"),
    recoveries = list(from = "I", to = "R")
  )
)

head(deterministic_output$trajectory)
head(deterministic_output$cumulative)
```

One cumulative name can also sum multiple disease transitions by supplying
equal-length `from` and `to` vectors:

```r
cumulative_flows = list(
  disease_onset = list(from = c("Lr", "Ld"), to = c("I", "I"))
)
```

Fixed-population stochastic simulation uses the same request shape, but the
cumulative table counts realised Gillespie events at each requested output
time. These counts are derived from event logs and do not add state variables
or propensities:

```r
stochastic_output <- simulate_stochastic(
  initial_state = initial_state,
  times = seq(0, 10, by = 1),
  model = SIRModel(gamma = 0.25),
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  seed = 123,
  return_events = TRUE,
  cumulative_flows = list(
    infections = list(from = "S", to = "I")
  )
)

head(stochastic_output$events)
head(stochastic_output$cumulative)
```

For downstream CFR or severity workflows, `as_cfr_data()` selects one
cumulative flow as cases and one as deaths, then converts cumulative values to
interval increments by default:

```r
cfr_input <- as_cfr_data(
  deterministic_output,
  cases = "infections",
  deaths = "recoveries"
)

if (requireNamespace("cfr", quietly = TRUE)) {
  head(cfr_input)
}
```

See [examples/deterministic_cumulative_flows.R](examples/deterministic_cumulative_flows.R)
and [examples/stochastic_cumulative_flows.R](examples/stochastic_cumulative_flows.R)
for runnable scripts.

[examples/tb_age_structured_demography.R](examples/tb_age_structured_demography.R)
shows the same public API pieces in a toy TB-style model with susceptible,
recent latent, remote latent, infectious, treatment, and recovered states,
closed demographic turnover, age-assortative mixing, and cumulative event
accounting. Its parameters are illustrative only; it is not a calibrated
Kiribati model and should not be used for policy analysis.

## Observed Case Data

`as_agepi_cases()` maps a user-supplied case or event data frame into the
`age_group` labels used by an `AgeStructure()`. It can derive groups from a
numeric exact-age column or validate an existing age-group column. This is an
input-preparation helper for observed records; it does not turn deterministic
compartment trajectories into individual-level line lists.

```r
observed_cases <- data.frame(
  case_id = c("a", "b", "c"),
  onset_day = c(0, 1, 1),
  age = c(3, 8, 24)
)

agepi_cases <- as_agepi_cases(
  observed_cases,
  age_structure = ages,
  age_col = "age"
)

aggregate(case_id ~ onset_day + age_group, agepi_cases, length)
```

If the optional `linelist` package is installed and the input already inherits
from `linelist`, existing tags are preserved where possible. `linelist` is
suggested, not imported.

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
demographic workflow and diagnostics. For a conceptual guide to choosing
between the demography workflows, see
[docs/demography_workflows.md](docs/demography_workflows.md).

Demographic-only and deterministic epidemic-demography simulations support two
ageing policies:

- `ageing_policy = "exponential"` is the default derivative-based path and
  preserves existing behaviour.
- `ageing_policy = "annual_cohort"` applies a discrete cohort projection once
  per annual interval. In deterministic epidemic simulations, infection ODEs run
  over each one-year interval and the annual demographic operator is applied at
  the year boundary. This mode is opt-in, deterministic-only, and requires a
  complete 1-year internal age grid ending in an open-ended age group when a
  `demographic_process` is supplied.

Annual-cohort coupling keeps births in the configured birth compartment
(`S` for `SIRModel()` and `SEIRModel()`), applies background mortality and
ageing to every compartment, and allocates age-total net migration using
`migration_policy`. Use `output_age_structure` to aggregate 1-year internal
outputs back to broader reporting groups. See
[docs/coupled_epidemic_demography_operator_splitting_design.md](docs/coupled_epidemic_demography_operator_splitting_design.md)
and [examples/annual_cohort_sir_demography.R](examples/annual_cohort_sir_demography.R).

Lightweight exploratory plotting helpers are available when `ggplot2` is
installed. `plot_population_pyramid()`, `plot_population_projection()`,
`plot_age_structure()`, and `plot_demography()` inspect supplied population
tables or `Demography()` objects without downloading data, interpolating
population values, or changing simulation behaviour.

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
- `load_contact_matrix_source()` plus
  `adapt_contact_matrix_to_age_structure()` for explicitly separating contact
  source loading from model age-grid adaptation. Current source loaders cover
  POLYMOD via optional `socialmixr`, Prem et al. matrices via optional
  `contactdata`, and conmat-generated matrices from caller-supplied population
  data via optional `conmat`;
- `population_from_wpp()` and `demography_from_wpp()` for WPP-style population
  tables;
- WPP-style fertility, mortality, and migration standardisers.

These adapters reshape and validate supplied data; they do not attempt to
reproduce full external projection systems. See
[docs/external_data_adapters.md](docs/external_data_adapters.md) and
[docs/contact_matrix_integration_design.md](docs/contact_matrix_integration_design.md).
For a user-facing explanation of contact-matrix workflows, see
[docs/contact_matrix_workflows.md](docs/contact_matrix_workflows.md).

If `contactdata`, `socialmixr`, `wpp2024`, or `conmat` are unavailable, the
matching loader or example will tell you which package is missing and how to
proceed. For the WPP and conmat pathways, the simplest fallback is to supply a
plain data frame or matrix directly and use the dependency-free adapters.

The intended source workflow keeps provenance separate from adaptation:

```r
if (requireNamespace("contactdata", quietly = TRUE) &&
    "Fiji" %in% contactdata::list_countries()) {
  source_contacts <- load_contact_matrix_source(
    source = "prem",
    country = "Fiji",
    setting = "all"
  )

  contact_matrix <- adapt_contact_matrix_to_age_structure(
    source_contacts,
    age_structure,
    population = source_grid_population,
    method = "source_band"
  )
}
```

`population` is required when adapting a finer source matrix to coarser model
age groups because aggregation is source-grid recipient-population weighted.
When expanding a coarser source matrix to nested target ages, the adapter
preserves total contacts with each original source age band and splits them
across target source age groups using target-grid population weights when
supplied, or equal weights otherwise.

Prem/contactdata matrices are preferred when the requested country is available
in the installed `contactdata` dataset. Kiribati may be absent from that
dataset; nearby Pacific Prem matrices should only be used as explicit proxy
assumptions. The Kiribati TB example falls back to POLYMOD UK only as a
documented proxy, not as Kiribati-specific contact data.

## Examples

Core examples are designed to run with the package's regular imports. Optional
integration examples check for their extra packages with `requireNamespace()`
and print a message instead of failing when those packages are unavailable.

- [examples/mock_sir_deterministic.R](examples/mock_sir_deterministic.R):
  deterministic age-structured SIR.
- [examples/generic_sir.R](examples/generic_sir.R): SIR through
  `CompartmentModel()`.
- [examples/generic_seir.R](examples/generic_seir.R): SEIR through
  `CompartmentModel()`.
- [examples/generic_msir.R](examples/generic_msir.R): custom MSIR model.
- [examples/epiparameter_seir.R](examples/epiparameter_seir.R): optional
  `epiparameter` incubation-period parameterisation for SEIR.
- [examples/deterministic_cumulative_flows.R](examples/deterministic_cumulative_flows.R):
  deterministic cumulative infection and recovery flows.
- [examples/stochastic_cumulative_flows.R](examples/stochastic_cumulative_flows.R):
  stochastic cumulative flows derived from realised event logs.
- [examples/tb_age_structured_demography.R](examples/tb_age_structured_demography.R):
  toy TB-style `CompartmentModel()` with demography and cumulative flows.
- [examples/observed_case_age_groups.R](examples/observed_case_age_groups.R):
  observed case records mapped into agepi age groups.
- [examples/mock_seir_demography.R](examples/mock_seir_demography.R): SEIR with
  demographic turnover.
- [examples/annual_cohort_sir_demography.R](examples/annual_cohort_sir_demography.R):
  toy deterministic SIR with annual-cohort demographic operator splitting.
- [examples/mock_demographic_workflow.R](examples/mock_demographic_workflow.R):
  demographic-only workflow with diagnostics.
- [examples/wpp_demography_validation.R](examples/wpp_demography_validation.R):
  optional WPP-connected demography benchmark that skips cleanly unless
  `wpp2024` is installed.
- [examples/kiribati_tb_realistic_demography.R](examples/kiribati_tb_realistic_demography.R):
  optional Kiribati TB public-data scaffold using WPP 2024 demography and a
  public contact-matrix source. It prefers Prem/contactdata when Kiribati is
  available in the installed dataset; otherwise it falls back to a documented
  POLYMOD UK proxy. It skips cleanly unless the WPP and solver dependencies
  are installed. It is not calibrated and is not a policy model.
- [examples/demography_plots.R](examples/demography_plots.R): synthetic
  exploratory demography plots using optional `ggplot2`.
- [examples/validation/finalsize_sir_final_size.R](examples/validation/finalsize_sir_final_size.R):
  optional closed-population final-size comparison using `finalsize` when it is
  installed.
- `examples/validation/epidynamics_*.R`: optional comparisons with
  `EpiDynamics` when it is installed.

## Documentation

- [docs/model_conventions.md](docs/model_conventions.md): state-vector,
  force-of-infection, contact-matrix, and demographic schedule conventions.
- [docs/development_status.md](docs/development_status.md): current scope,
  limitations, solver notes, adapter boundaries, and roadmap notes.
- [docs/generic_compartment_models.md](docs/generic_compartment_models.md):
  custom compartment models.
- [docs/epidemic_model_workflows.md](docs/epidemic_model_workflows.md):
  conceptual guide to epidemic and compartmental model workflows.
- [docs/cumulative_flow_states_design.md](docs/cumulative_flow_states_design.md):
  deterministic auxiliary cumulative states and stochastic event-log-derived
  cumulative outputs.
- [docs/age_structured_transmission_design.md](docs/age_structured_transmission_design.md):
  age-structured transmission design notes.
- [docs/external_data_adapters.md](docs/external_data_adapters.md): WPP-style
  demographic and contact-matrix adapter notes.
- [docs/contact_matrix_workflows.md](docs/contact_matrix_workflows.md):
  conceptual guide to choosing and adapting contact matrices.
- [docs/demography_workflows.md](docs/demography_workflows.md): conceptual
  guide to choosing between demographic workflows.
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
