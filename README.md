# agepi

`agepi` is an early-stage R package for age-structured epidemic model prototypes.

The current implementation supports a deterministic age-structured SIR workflow:

- define and validate age groups with `AgeStructure()` and `validate_age_structure()`;
- convert between long-form state data and numeric solver vectors;
- compute a reusable age-structured force of infection;
- validate and coerce contact matrices;
- construct a minimal SIR model with `SIRModel()`;
- compute SIR infection and recovery transition rates;
- convert transition rates to deterministic derivatives;
- run a deterministic SIR simulation with explicit Euler time steps;
- summarise deterministic simulation output with `compartment_totals()`, `age_group_totals()`, and `total_population()`.

It also includes small utility layers for age-bin transformations, demography
tables, and contact-matrix handling. These utilities are deliberately separate
from the current simulator: they help prepare or inspect inputs, but they do not
add demographic projection dynamics or additional disease models.

## Current limitations

The package currently supports deterministic SIR only. `simulate_deterministic()` currently supports only `method = "euler"`.

The current scope is deliberately small:

- static contact matrix;
- static `beta`, susceptibility, and infectiousness inputs;
- mock examples only;
- no SEIR models or stochastic simulation;
- no adaptive or external deterministic solver backend;
- no births, deaths, ageing, migration, fertility, mortality, or demographic projection dynamics;
- no WPP integration;
- no `socialmixr` or `conmat` dependency or package-specific import workflow;
- no reciprocity correction or population balancing for contact matrices;
- no plotting or fitting.

## State-vector convention

Numeric state vectors use compartment-major, age-group-minor ordering. For compartments `c("S", "I", "R")` and age groups `c("0-4", "5-9")`, the vector order is:

```text
S_0-4, S_5-9, I_0-4, I_5-9, R_0-4, R_5-9
```

`state_vector_to_long()` and `simulate_deterministic()` interpret numeric vectors by position only. Names on numeric state vectors are ignored when converting back to long form or simulating.

## Force-of-infection convention

The force of infection is:

```text
lambda = beta * susceptibility *
  contact_matrix %*% (infectiousness * infectious / population)
```

Equivalently:

```text
lambda_a(t) =
  beta(t) * s_a(t) *
  sum_b C_ab(t) * iota_b(t) * I_b(t) / N_b(t)
```

Rows of `contact_matrix` are recipient age groups `a`; columns are source age groups `b`.
That is, `contact_matrix[a, b]` gives contacts made by recipient group `a` with source group `b`.

## Minimal deterministic SIR example

```r
library(agepi)

age_structure <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

population <- c(1000, 1200)
initial_infections <- c(5, 3)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = age_structure$n_age_groups),
  age_group = rep(age_structure$age_groups, times = 3),
  value = c(population - initial_infections, initial_infections, 0, 0),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(
  c(4, 2,
    2, 5),
  nrow = age_structure$n_age_groups,
  byrow = TRUE
)

model <- SIRModel(gamma = 0.25)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 1, by = 0.1),
  model = model,
  age_structure = age_structure,
  contact_matrix = contact_matrix,
  beta = 0.08,
  method = "euler"
)

head(simulation)
```

See `examples/mock_sir_deterministic.R` for a slightly larger mock-only example.

## Small utilities

### Simulation summaries

`compartment_totals()`, `age_group_totals()`, and `total_population()` summarise
long-format simulation output returned by `simulate_deterministic()`.

### Age-vector transformations

`aggregate_age_vector()` aggregates a numeric vector from one age structure to a
coarser age structure when each target age bin is an exact union of complete
source age bins.

`segregate_age_vector()` splits a coarse age vector into a finer age structure
when each target age bin is nested inside one source age bin. Splits use explicit
target weights within each source age bin.

`transform_age_vector()` is a convenience wrapper for exact age transformations,
including identity, aggregation, segregation, and mixed exact transformations
when the source and target bins align to a common set of boundaries. Its
`split_method` controls how source bins split across multiple target bins.

### Demography helpers

`validate_demography_table()` checks a tidy table with `time`, `age_group`, and
`population` columns against an age structure.

`Demography()` stores a validated demography table sorted by time and age-group
order. `demography_times()`, `demography_population_vector()`, and
`demography_population_table()` provide exact-time accessors.

Current demography support is validation, storage, sorting, and population
access only. It does not implement demographic projection, interpolation,
fertility, mortality, births, deaths, ageing, migration, or WPP integration.

### Contact matrices

`validate_contact_matrix()` checks that a contact matrix is numeric, finite,
non-missing, non-negative, square, and optionally has dimensions matching an age
structure.

`as_agepi_contact_matrix()` coerces supported inputs to agepi's recipient-source
matrix convention. Current supported inputs are numeric matrices, numeric data
frames, socialmixr-like lists with a numeric `matrix` element, and conmat-style
long data frames with `age_group_from`, `age_group_to`, and `contacts` columns.
This is dependency-free coercion; agepi does not depend on `socialmixr` or
`conmat`.

`transform_contact_matrix()` aggregates contact matrices from a finer source age
structure to a coarser target age structure when every target age bin is an exact
union of source age bins. It currently supports exact aggregation only and
rejects transformations that would require source-bin splitting or general
rebinning.

## Design notes

See `docs/age_structured_transmission_design.md`.
