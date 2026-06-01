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
  for infection-only deterministic SIR and SEIR models, and for SIR/SEIR with
  demographic coupling;
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
  `time_policy = "step"` and bounded linear rate interpolation via
  `time_policy = "linear"`;
- `compare_demography_to_observed()` and
  `summarise_demography_comparison()`;
- `standardise_wpp_fertility()`, `standardise_wpp_mortality()`, and
  `standardise_wpp_migration()` as dependency-free adapters for supplied
  WPP-like tables;
- `population_from_wpp()`/`demography_from_wpp()` for WPP-style population
  tables;
- `fertility_from_wpp_percent_asfr()` for WPP percent-ASFR distributions plus
  TFR tables;
- `mortality_from_wpp_mx()` and `mortality_from_wpp(quantity = "mx")` for
  WPP-style central death rates;
- `implied_demographic_residual()` and
  `residual_to_migration_schedule()` for diagnostic residual accounting.

These helpers include schedule-level step and linear lookup for rate-like
fertility, mortality, and migration schedules. They do not implement population
interpolation, automatic residual forcing, or WPP projection matching. Infection
coupling is handled separately by `simulate_deterministic()`.

### Milestone 4E: First-pass SIR/SEIR-demography coupling

Implemented:

- optional `demographic_process` argument for `simulate_deterministic()`;
- deterministic coupling for SIR and SEIR models with Euler or optional `deSolve`;
- demographic schedule lookup via `time_policy = c("exact", "step", "linear")`
  for rate-like fertility, mortality, and migration schedules;
- SIR births entering the youngest susceptible age group;
- SIR mortality applied proportionally to `S`, `I`, and `R`;
- SIR ageing applied independently to `S`, `I`, and `R`;
- net migration allocation controlled by `migration_policy`, defaulting to
  susceptible-only allocation for backwards compatibility;
- force of infection continuing to use current `S + I + R` by age group for
  SIR, or `S + E + I + R` for SEIR.

This is not WPP projection matching. Coupled demography does not support
population interpolation, time-varying contact matrices, or stochastic
simulation.

The default `S`-only migration rule is an allocation convention for age-total
net migration and residual-derived migration schedules. It should not be read as
a mechanistic movement model for susceptible, infectious, or recovered people.
The explicit proportional policy allocates positive or negative net migration by
current age-specific compartment shares; the explicit error policy rejects
non-zero migration when allocation would be ambiguous.

### Milestone 4F: SEIR deterministic infection-only support

Implemented:

- `SEIRModel(sigma, gamma)`;
- SEIR transition rates for infection, progression, and recovery;
- deterministic infection-only SEIR simulation with Euler;
- optional infection-only SEIR simulation with `method = "deSolve"`;
- focused tests for SEIR model validation, transition rates, derivatives, and
  supported deterministic solver backends.

Implemented SEIR-demography policy:

- fertility exposure uses current `S + E + I + R` by age group;
- births enter the youngest susceptible compartment only;
- background mortality applies independently to `S`, `E`, `I`, and `R`;
- ageing applies independently to `S`, `E`, `I`, and `R`;
- net migration uses the same `migration_policy` choices as SIR, defaulting to
  the existing `S` convention for age-total net migration schedules;
- `E -> I` progression and `I -> R` recovery remain disease-model transitions;
- force of infection depends on `I`, not `E`;
- no disease-induced mortality, vaccination, waning immunity,
  compartment-specific demographic rates, or WPP projection matching.

## Future Work

### Milestone 5: Package-level documentation

Planned:

- package overview documentation;
- function examples aligned with the current deterministic SIR scope;
- clearer user-facing notes on current limitations, including optional external
  data adapters.

### Milestone 6: Deterministic solver backend

Implemented:

- optional `method = "deSolve"`/`method = "ode"` backend for infection-only
  SIR and SEIR models;
- optional `deSolve` backend for SIR and SEIR with demographic coupling;
- tests that preserve the transition-rate to derivative spine across Euler and
  deSolve where applicable.

Future solver work should focus on explicit event policies, time-varying
contact inputs, and any additional solver interface guarantees.

### Milestone 7: Demography extensions

Planned:

- explicit policy for any future population interpolation or projection helpers;
- optional solver/backend refinements for demographic-only simulations;
- explicit residual application policies if residual forcing is ever added;
- richer infection-demography policies beyond the first-pass SIR/SEIR coupling.

Implemented interpolation path: `time_policy = "linear"` is schedule-level and
generic, not WPP-specific. It applies first to rate-like fertility, mortality,
and migration schedules, rejects extrapolation by default, and leaves population
`Demography()` accessors exact until state-trajectory semantics are specified.

Population interpolation and residual forcing remain outside the current
deterministic simulator.

### Milestone 8: Contact matrix integration

Planned:

- integration of contact schedules into future simulators;
- contact-matrix splitting and general rebinning;
- reciprocity correction or population balancing if explicitly needed;
- optional mixing/contact object constructor.

No `socialmixr` or `conmat` dependency is currently included. Projection
dynamics and reciprocity correction remain future work.

### Milestone 9: SEIR and additional disease structures

Implemented:

- minimal deterministic infection-only SEIR support.

Planned:

- broader transition-rate validation for additional compartment structures;
- vaccination, waning immunity, and event handling only after explicit design
  decisions.

### Milestone 10: Stochastic simulation infrastructure

Implemented:

- `simulate_stochastic()`;
- fixed-population stochastic SIR, SEIR, and supported generic
  `CompartmentModel()` support;
- Gillespie/direct-method dynamics only;
- observation-time-aligned trajectory output;
- optional event log output;
- reproducible trajectories when `seed` is supplied;
- focused stochastic simulator tests;
- examples for stochastic SIR and SEIR.

Current stochastic support is deliberately narrow:

- fixed population only;
- no demography;
- no ageing;
- no mortality;
- no fertility;
- no migration;
- no tau-leaping;
- generic `CompartmentModel()` support is limited to fixed-population
  transition structures expressible through `transition_rates()`;
- SIR event semantics: `S -> I` infection and `I -> R` recovery;
- SEIR event semantics: `S -> E` infection, `E -> I` progression, and
  `I -> R` recovery;
- force of infection uses current stochastic infectious counts and the existing
  recipient-row, source-column contact-matrix convention.

Future stochastic work may include:

- broader stochastic simulator design beyond fixed-population supported
  transition-rate models;
- additional stochastic methods beyond Gillespie/direct-method dynamics;
- ageing remains deterministic;
- stochastic mortality, fertility, and migration may be later optional
  extensions with explicit policy, especially for gross versus net migration.
