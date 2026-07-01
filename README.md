# agepi

`agepi` is a lightweight R toolkit for transparent, inspectable age-structured epidemic and demographic modelling workflows.

It is built for people who want to define age groups explicitly, work with compartment-by-age state vectors, use contact matrices, and keep the modelling steps easy to inspect.

## What is agepi for?

`agepi` supports workflows built from:

- age structures;
- compartment models by age group;
- contact matrices;
- SIR and SEIR starting models;
- custom compartment models via `CompartmentModel()`;
- age-specific susceptibility and infectiousness;
- cumulative flow outputs;
- demographic schedules and age turnover;
- optional adapters for external data sources.

The package is designed for clear modelling workflows rather than black-box fitting. SIR and SEIR are convenient starting constructors, but they are not the conceptual centre of the package.

The usual pattern is straightforward: define an age structure, create an initial state by compartment and age group, choose a contact matrix, select a disease model, and then simulate. Outputs stay in ordinary R data frames so they can be inspected, summarised, or passed into downstream analysis without special glue code.

## Installation

Install the core package from GitHub:

```r
remotes::install_github("michaeltmeehan/agepi")
```

Optional integrations such as WPP-backed demography, external contact matrices, plotting helpers, and validation examples are documented in the package docs and examples. They may require additional packages.

If you are starting from a local checkout, `devtools::install(".")` also works. The key point is that the core package stays small, while optional backends remain optional.

## A minimal workflow

This is the smallest age-structured deterministic SIR workflow:

```r
library(agepi)

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9", "10-14"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, 14)
)

population <- c(1000, 1200, 900)
initial_infections <- c(5, 3, 2)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 3),
  value = c(population - initial_infections, initial_infections, rep(0, 3)),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(
  c(
    4, 2, 1,
    2, 5, 2,
    1, 2, 4
  ),
  nrow = age_structure$n_age_groups,
  byrow = TRUE
)

model <- SIRModel(gamma = 0.25)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 10, by = 1),
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08
)

head(compartment_totals(simulation))
```

The key pieces are the age structure, the initial compartment state, the contact matrix, and the disease model. For most users, that is enough to start building and inspecting a simple age-structured model.

If you need more control, `simulate_deterministic()` also supports age-specific susceptibility and infectiousness, demographic coupling, output aggregation, and cumulative flow tracking. Those features are documented in the workflow pages below rather than expanded here.

## Learning path

If you want a guided path through the package, start with the new tutorials in
order:

1. [tutorials/01_getting_started_age_structured_models.qmd](tutorials/01_getting_started_age_structured_models.qmd)
2. [tutorials/02_custom_models_and_cumulative_flows.qmd](tutorials/02_custom_models_and_cumulative_flows.qmd)
3. [tutorials/03_demography_and_coupling.qmd](tutorials/03_demography_and_coupling.qmd)
4. [tutorials/04_stochastic_simulation.qmd](tutorials/04_stochastic_simulation.qmd)
5. [tutorials/05_demography_only_and_diagnostics.qmd](tutorials/05_demography_only_and_diagnostics.qmd)
6. [tutorials/06_external_data_and_case_study_scaffolds.qmd](tutorials/06_external_data_and_case_study_scaffolds.qmd)

The full index lives in [tutorials/README.md](tutorials/README.md).

## Main workflows

| Goal | Start here |
| --- | --- |
| Run a simple age-structured SIR or SEIR model | [tutorials/01_getting_started_age_structured_models.qmd](tutorials/01_getting_started_age_structured_models.qmd), [docs/epidemic_model_workflows.md](docs/epidemic_model_workflows.md) |
| Define a custom compartment model | [tutorials/02_custom_models_and_cumulative_flows.qmd](tutorials/02_custom_models_and_cumulative_flows.qmd), [docs/generic_compartment_models.md](docs/generic_compartment_models.md) |
| Track cumulative infections, recoveries, deaths, or other flows | [tutorials/02_custom_models_and_cumulative_flows.qmd](tutorials/02_custom_models_and_cumulative_flows.qmd), [docs/cumulative_flow_states_design.md](docs/cumulative_flow_states_design.md) |
| Add demographic turnover | [tutorials/03_demography_and_coupling.qmd](tutorials/03_demography_and_coupling.qmd), [docs/demography_workflows.md](docs/demography_workflows.md) |
| Try annual-cohort demographic coupling | [tutorials/03_demography_and_coupling.qmd](tutorials/03_demography_and_coupling.qmd), [docs/coupled_epidemic_demography_operator_splitting_design.md](docs/coupled_epidemic_demography_operator_splitting_design.md) |
| Work with contact matrices | [docs/contact_matrix_workflows.md](docs/contact_matrix_workflows.md), [docs/external_data_adapters.md](docs/external_data_adapters.md) |
| Use WPP-style demographic inputs | [examples/wpp_demography_validation.R](examples/wpp_demography_validation.R), [examples/wpp_projection_backed_demography.R](examples/wpp_projection_backed_demography.R) |
| Explore stochastic simulation workflows | [docs/stochastic_simulation_design.md](docs/stochastic_simulation_design.md), [examples/stochastic_sir.R](examples/stochastic_sir.R), [examples/stochastic_seir.R](examples/stochastic_seir.R) |
| Inspect the Kiribati TB scaffold | [examples/kiribati_tb_realistic_demography.R](examples/kiribati_tb_realistic_demography.R), [docs/disease_notes/kiribati_tb_modelling_inputs.md](docs/disease_notes/kiribati_tb_modelling_inputs.md) |

## Current scope

`agepi` is under active development. It is currently best suited for research prototypes, method development, and transparent age-structured simulation workflows where the model structure, age grouping, contact assumptions, demographic schedules, and outputs need to remain explicit and inspectable.

The package currently supports deterministic SIR, SEIR, and custom compartment models; fixed-population stochastic simulations for supported models; cumulative flow outputs; demographic-only simulations; deterministic epidemic-demography coupling; and adapters for selected external demographic and contact-matrix data shapes.

It is not yet intended to be a full calibration, inference, stochastic-demography, or production forecasting framework.

## Documentation

Key documentation pages:

- [docs/model_conventions.md](docs/model_conventions.md)
- [docs/epidemic_model_workflows.md](docs/epidemic_model_workflows.md)
- [docs/generic_compartment_models.md](docs/generic_compartment_models.md)
- [docs/demography_workflows.md](docs/demography_workflows.md)
- [docs/contact_matrix_workflows.md](docs/contact_matrix_workflows.md)
- [docs/external_data_adapters.md](docs/external_data_adapters.md)
- [docs/development_status.md](docs/development_status.md)

## Contributing and issues

Bug reports, example requests, and documentation improvements are welcome.

When reporting an issue, a small reproducible example and the package version are especially helpful.
