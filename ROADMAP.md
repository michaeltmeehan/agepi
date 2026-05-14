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

### Milestone 3B: Deterministic Euler simulation

Implemented:

- `simulate_deterministic()` using explicit Euler time steps;
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

Current demography support is validation, storage, sorting, and exact-time
population access only. It does not include interpolation, fertility, mortality,
births, deaths, ageing, migration, or demographic projection dynamics.

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

## Future Work

### Milestone 5: Package-level documentation

Planned:

- package overview documentation;
- function examples aligned with the current deterministic SIR scope;
- clearer user-facing notes on current limitations, including optional external
  data adapters.

### Milestone 6: Deterministic solver backend

Planned:

- solver interface for deterministic simulation;
- optional non-Euler backend after the current explicit Euler path is stable;
- tests that preserve the transition-rate to derivative spine.

### Milestone 7: Demography extensions

Planned:

- time-varying demographic accessors;
- demographic interpolation or projection helpers;
- births, deaths, ageing, migration, fertility, and mortality dynamics;
- integration of demographic dynamics into future simulators.

These are not part of the current deterministic SIR simulator.

### Milestone 8: Contact matrix integration

Planned:

- time-varying contact matrix accessors;
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
