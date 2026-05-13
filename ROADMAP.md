# Roadmap

## Completed

### Milestone 1: Age structure and state mapping

Implemented:

- `AgeStructure()`;
- `validate_age_structure()`;
- `state_long_to_vector()`;
- `state_vector_to_long()`;
- `aggregate_age_vector()`;
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

## Future Work

### Milestone 4: Package-level documentation

Planned:

- package overview documentation;
- function examples aligned with the current deterministic SIR scope;
- clearer user-facing notes on current limitations.

### Milestone 5: Deterministic solver backend

Planned:

- solver interface for deterministic simulation;
- optional non-Euler backend after the current explicit Euler path is stable;
- tests that preserve the transition-rate to derivative spine.

### Milestone 6: Minimal demography layer

Planned:

- demography object constructor;
- demographic quantity lookup or interpolation;
- mock demography inputs that remain separate from simulation code.

This does not currently include WPP integration, ageing, births, deaths, or migration in the implemented simulator.

### Milestone 7: Contact matrix integration

Planned:

- mixing/contact object constructor;
- contact matrix validation at the package boundary;
- compatibility with externally prepared contact matrices.

This does not currently include `socialmixr` or `conmat` integration.

### Milestone 8: Stochastic simulation infrastructure

Planned:

- stochastic event-intensity interface based on transition rates;
- stochastic simulator design and tests.

No stochastic simulation is currently implemented.
