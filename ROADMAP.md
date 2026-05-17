# Roadmap

## Completed

### Milestone 1: Age structure and state mapping

Implemented:

- `AgeStructure()`;
- `validate_age_structure()`;
- `state_long_to_vector()`;
- `state_vector_to_long()`;
- `aggregate_age_vector()`;
- `segregate_age_vector()`;
- `transform_age_vector()`;
- tests for age validation and state mapping.

### Milestone 2A: Force of infection

Implemented:

- reusable `force_of_infection()`;
- recipient-row, source-column contact-matrix convention;
- age-specific susceptibility and infectiousness vectors.

### Milestone 2B: SIR model and transition rates

Implemented:

- `SIRModel()`;
- `validate_disease_model()`;
- `transition_rates()` for SIR infection and recovery rates;
- tests using small manually checkable examples.

### Milestone 3A: Deterministic derivative

Implemented:

- `rates_to_derivative()`;
- deterministic derivative construction from transition-rate tables.

### Milestone 3B: Deterministic simulation

Implemented:

- `simulate_deterministic()` using explicit Euler time steps;
- optional `method = "deSolve"` backend using the suggested `deSolve` package
  for the same static deterministic SIR model;
- tidy long-format simulation outputs;
- `compartment_totals()`, `age_group_totals()`, and `total_population()`;
- mock deterministic SIR example.

### Milestone 4A: Minimal demography layer

Implemented:

- `validate_demography_table()`;
- `Demography()`;
- optional `demography_from_wpp()` adapter for WPP-style tabular data;
- sorted storage of tidy `time`, `age_group`, and `population` tables;
- `demography_times()`;
- `demography_population_vector()`;
- `demography_population_table()`.

This milestone established population trajectory storage and exact-time access.
Later completed demographic-only helpers add ageing, fertility, mortality,
migration, simulation, comparison, and residual diagnostics without adding
interpolation, WPP projection matching, residual forcing, or infection-demography
coupling.

### Milestone 4B: Contact matrix utilities

Implemented:

- `validate_contact_matrix()`;
- `as_agepi_contact_matrix()`;
- optional `contact_matrix_from_socialmixr()` adapter;
- optional `contact_matrix_from_conmat()` adapter;
- dependency-free coercion for numeric matrices, numeric data frames,
  socialmixr-like lists with `x$matrix`, and conmat-style long data frames;
- `transform_contact_matrix()` for exact fine-to-coarse contact-matrix
  aggregation.

Current contact-matrix support does not include reciprocity correction,
population balancing, package dependencies, socialmixr/conmat-specific S3
methods, or contact-matrix splitting/rebinning.

### Milestone 4C: Exact-time external input accessors

Implemented:

- `demography_population_at()` for exact-time population-vector lookup;
- `ContactSchedule()` for externally supplied contact matrices indexed by time;
- `contact_matrix_at()` for exact-time contact-matrix lookup.

This is the first step toward time-varying infection inputs. Infection-simulator
integration, interpolation, reciprocity correction, and population balancing
remain future work.

### Milestone 4D: Demographic-only process components

Implemented:

- WPP-style one-year and five-year age-grid helpers;
- `AgeingOperator()` and `validate_ageing_operator()`;
- `FertilitySchedule()`, `MortalitySchedule()`, and `MigrationSchedule()`;
- `DemographicProcess()` and `build_demographic_process()`;
- `demographic_derivative()` and `simulate_demography()` using exact-time Euler
  semantics by default, with opt-in interval-start stepwise schedule lookup via
  `time_policy = "step"`;
- `compare_demography_to_observed()` and
  `summarise_demography_comparison()`;
- `standardise_wpp_fertility()`, `standardise_wpp_mortality()`, and
  `standardise_wpp_migration()` as dependency-free adapters for supplied
  WPP-like tables;
- `implied_demographic_residual()` and
  `residual_to_migration_schedule()` for diagnostic residual accounting.

These helpers are demographic-only. Stepwise schedule lookup does not implement
interpolation, automatic residual forcing, WPP projection matching, or
infection-demography coupling.

### Milestone 4E: First-pass SIR-demography coupling

Implemented:

- optional `demographic_process` argument for `simulate_deterministic()`;
- Euler-only deterministic coupling for SIR models;
- demographic schedule lookup via `time_policy = c("exact", "step")`;
- births entering the youngest susceptible age group;
- mortality applied proportionally to `S`, `I`, and `R`;
- ageing applied independently to `S`, `I`, and `R`;
- net migration applied to susceptible compartments only;
- force of infection continuing to use current `S + I + R` by age group.

This is not WPP projection matching. Coupled demography does not support
`deSolve`/`ode`, interpolation, time-varying contact matrices, SEIR, or
stochastic simulation.

## Future Work

### Milestone 5: Package-level documentation

Planned:

- package overview documentation;
- function examples aligned with the current deterministic SIR scope;
- clearer user-facing notes on current limitations, including optional external
  data adapters.

### Milestone 6: Deterministic solver backend

Planned:

- broader solver interface for deterministic simulation;
- tests that preserve the transition-rate to derivative spine.

### Milestone 7: Demography extensions

Planned:

- explicit policy for any future demographic interpolation or projection
  helpers;
- optional solver/backend refinements for demographic-only simulations;
- explicit residual application policies if residual forcing is ever added;
- richer infection-demography policies beyond the first-pass SIR coupling.

These are not part of the current deterministic SIR simulator.

### Milestone 8: Contact matrix integration

Planned:

- integration of contact schedules into future simulators;
- contact-matrix splitting and general rebinning;
- reciprocity correction or population balancing if explicitly needed;
- optional mixing/contact object constructor.

No `socialmixr` or `conmat` dependency is currently included. Projection
dynamics and reciprocity correction remain future work.

### Milestone 9: SEIR and additional disease structures

Planned:

- SEIR model support;
- broader transition-rate validation for additional compartment structures.

Only SIR is currently implemented.

### Milestone 10: Stochastic simulation infrastructure

Planned:

- stochastic event-intensity interface based on transition rates;
- stochastic simulator design and tests.

No stochastic simulation is currently implemented.
