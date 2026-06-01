# Development Status And Current Limitations

`agepi` is an early-stage R package for age-structured epidemic and demographic
model prototypes. The package is useful for deterministic examples and
workflow development, but its modelling scope is intentionally narrow.

## Current Support

The current implementation supports deterministic age-structured SIR, SEIR, and
generic `CompartmentModel()` workflows, plus fixed-population Gillespie
simulation for supported disease transition structures:

- define and validate age groups with `AgeStructure()` and
  `validate_age_structure()`;
- convert between long-form state data and numeric solver vectors;
- compute a reusable age-structured force of infection;
- validate and coerce contact matrices;
- construct minimal SIR and SEIR models with `SIRModel()` and `SEIRModel()`;
- construct custom deterministic compartment models with `CompartmentModel()`;
- compute transition rates and convert them to deterministic derivatives;
- run deterministic simulations through a shared solver path;
- run fixed-population stochastic SIR, SEIR, and supported generic
  `CompartmentModel()` simulations with `simulate_stochastic()`;
- optionally couple deterministic models to a first-pass demographic process;
- summarise deterministic simulation output with `compartment_totals()`,
  `age_group_totals()`, and `total_population()`.

The package also includes utility layers for age-bin transformations,
demographic-only ODE components, demography tables, residual diagnostics,
WPP-style demographic adapters, and contact-matrix handling. These utilities
help prepare or inspect inputs, can run demographic-only workflows, and can be
coupled to deterministic epidemic models in the supported narrow mode.

## Solver Notes

`simulate_deterministic()` defaults to `deSolve` for infection-only
deterministic SIR and SEIR models when the suggested `deSolve` package is
installed, and falls back to Euler otherwise.

SIR, SEIR, or generic `CompartmentModel()` simulations with demographic
coupling keep the Euler default because exact schedule lookup remains the
default demographic policy. Pass `method = "deSolve"` explicitly for supported
coupled runs. The Euler and deSolve backends use the same static
transition-rate pathway for supported models.

Euler updates are intentionally not truncated. If a step would produce negative
compartment values, simulation stops with an error.

## Current Scope

The current package scope is deliberately small:

- static contact matrix;
- static `beta`, susceptibility, and infectiousness inputs;
- deterministic simulation plus fixed-population stochastic SIR, SEIR, and
  supported generic compartment models;
- optional `deSolve` backend for documented deterministic combinations;
- first-pass epidemic-demography coupling;
- dependency-free external-data adapters;
- no plotting, fitting, or calibration layer.

## Current Limitations

Important limitations to keep in mind:

- no stochastic demography, tau-leaping, or stochastic models beyond
  fixed-population supported transition-rate structures;
- no time-varying contact matrices in epidemic simulation;
- no demographic residual forcing or WPP projection matching;
- no reciprocity correction or population balancing for contact matrices;
- no automatic population interpolation for `Demography()` objects;
- no vaccination, waning immunity, disease-induced mortality, or
  compartment-specific demographic rates;
- no event handling or additional solver backends beyond Euler and optional
  `deSolve`;
- contact-matrix transformation supports exact aggregation only, not general
  rebinning or source-bin splitting.

## Generic CompartmentModel Limitations

`CompartmentModel()` supports compartment models that fit the current
transition-rate interface:

- infection transitions use the shared `force_of_infection()` convention;
- fixed transitions are per-capita flows with static scalar or named
  age-specific rates;
- multiple outgoing fixed transitions from the same source compartment are
  supported when destinations differ;
- one or more infectious compartments can contribute to force of infection
  through named, non-negative relative infectiousness weights;
- demographic births and susceptible-policy migration can target configured
  compartments.

Its stochastic support is limited to fixed-population Gillespie simulation for
supported within-age transition structures expressible through
`transition_rates()`. It does not currently provide time-varying interventions,
arbitrary nonlinear transition functions, vaccination schedules, stochastic
demography, or model-specific demographic rates.

See [generic_compartment_models.md](generic_compartment_models.md) for details.

## External Adapter Boundaries

External adapters are optional and dependency-free where possible. The package
does not require `wpp2024`, `socialmixr`, or `conmat` for core functionality.

The adapters reshape and validate supplied data. They do not implement
projection dynamics, interpolation, reciprocity correction, population
balancing, or full external package semantics. WPP-style inputs should be
pre-filtered to one country, location, or entity unless an adapter argument
selects a single location.

See [external_data_adapters.md](external_data_adapters.md) and
[contact_matrix_integration_design.md](contact_matrix_integration_design.md).

## Roadmap Notes

Likely future directions include broader compartment structures, time-varying
inputs, richer demographic coupling, broader stochastic simulation, calibration
and inference workflows, plotting helpers, multi-country batching, and possible
porting of the simulation core to Julia. These are design goals rather than
implemented features.

For lower-level modelling conventions, see [model_conventions.md](model_conventions.md).
