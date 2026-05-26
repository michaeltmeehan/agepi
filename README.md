# agepi

`agepi` is an early-stage R package for age-structured epidemic model prototypes.

The current implementation supports deterministic age-structured SIR and SEIR
workflows:

- define and validate age groups with `AgeStructure()` and `validate_age_structure()`;
- convert between long-form state data and numeric solver vectors;
- compute a reusable age-structured force of infection;
- validate and coerce contact matrices;
- construct minimal SIR and SEIR models with `SIRModel()` and `SEIRModel()`;
- compute SIR infection/recovery and SEIR infection/progression/recovery
  transition rates;
- convert transition rates to deterministic derivatives;
- run deterministic SIR or SEIR simulations with explicit Euler time steps, or
  with optional `deSolve::ode()` integration when `deSolve` is installed;
- optionally couple SIR or SEIR simulation to a first-pass demographic process with
  Euler or `deSolve`;
- summarise deterministic simulation output with `compartment_totals()`, `age_group_totals()`, and `total_population()`.

It also includes utility layers for age-bin transformations, demographic-only
ODE components, demography tables, residual diagnostics, and contact-matrix
handling. These utilities are deliberately separate from the current infection
simulator: they help prepare or inspect inputs, can run demographic-only
workflows, and can now be coupled to deterministic SIR and SEIR in a narrow
mode.

## Current limitations

The infection simulator currently supports deterministic SIR and SEIR.
`simulate_deterministic()` defaults to `method = "euler"` and can optionally
use `method = "deSolve"` for infection-only deterministic SIR and SEIR models,
and for SIR or SEIR with demographic coupling, when the suggested `deSolve`
package is installed.
SIR and SEIR support use the same static transition-rate pathway; the focused
package tests cover the Euler and deSolve backends.

The current scope is deliberately small:

- static contact matrix;
- static `beta`, susceptibility, and infectiousness inputs;
- mock examples only;
- no stochastic simulation;
- optional `deSolve` backend for infection-only SIR/SEIR and SIR/SEIR-demography
  coupling only;
- first-pass SIR/SEIR-demography coupling: births enter `S`, deaths and ageing
  apply independently to all disease compartments, and migration defaults to
  `S` with optional explicit proportional or error-on-migration policies;
- no demographic residual forcing or WPP projection matching;
- optional external-data adapters only; no required WPP, `socialmixr`, or `conmat` dependency;
- no reciprocity correction or population balancing for contact matrices;
- no plotting or fitting.

## Capability summary

| Area | Current support | Not in current scope |
|---|---|---|
| Disease models | Deterministic SIR and SEIR | Stochastic simulation, vaccination, waning immunity, disease-induced mortality |
| Solvers | Euler; optional `deSolve` for documented SIR/SEIR combinations | Event handling or additional solver backends |
| Demography coupling | Births to `S`; net migration defaults to `S` or can be proportional/error; mortality and ageing by compartment | Compartment-specific demographic rates |
| Time-varying inputs | Demographic schedules with exact, step, or linear rate lookup | Time-varying contact matrices |
| WPP-style data | Dependency-free population, fertility, mortality, and migration adapters | qx/survival conversion, residual migration fitting, projection matching |

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
See `examples/mock_seir_demography.R` for a small dependency-free SEIR example
with first-pass demographic coupling.

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
order. `demography_times()`, `demography_population_at()`/
`demography_population_vector()`, and `demography_population_table()` provide
exact-time accessors. Population interpolation is deliberately not implemented:
population tables may represent initial conditions, observed trajectories, or
projection targets, and those meanings need different interpolation policies.

Demographic-only helpers now cover age grids, ageing operators, fertility,
mortality, migration, demographic process assembly, Euler simulation,
comparison, and residual diagnostics. `simulate_demography()` uses exact-time
schedule lookup by default and offers opt-in interval-start stepwise lookup via
`time_policy = "step"` or bounded linear rate interpolation via
`time_policy = "linear"`. Stepwise schedules are left-continuous: a schedule row
at `t_i` applies from `t_i` up to the next schedule time. Linear interpolation
applies only to rate-like fertility, mortality, and migration schedules between
available schedule times; it does not extrapolate and does not interpolate
`Demography()` population access. `migration_count` remains a per-time additive
flow rather than an interval total. `simulate_deterministic()` can optionally use a
`DemographicProcess()` for first-pass SIR/SEIR-demography coupling with Euler or
`deSolve`. For adaptive `deSolve` runs with demographic schedules,
`time_policy = "linear"` is generally recommended; `time_policy = "step"` gives
piecewise-constant demographic rates, while `time_policy = "exact"` is
generally unsuitable unless all solver evaluation times are schedule times.
This coupling is not WPP projection matching and does not support time-varying
contact matrices or additional disease structures.

Current SIR-demography allocation rules are deliberately narrow. Fertility uses
the current total age-specific infection-state population `S + I + R` as its
exposure, births enter only the youngest susceptible compartment, mortality
applies independently to `S`, `I`, and `R`, and ageing moves `S`, `I`, and `R`
independently through the ageing operator. Net migration, including
residual-derived migration schedules, defaults to allocation entirely to `S`.
This `S`-only migration rule is an allocation convention for age-total net
migration inputs, not a mechanistic model of who moves while infected or
recovered. `simulate_deterministic(migration_policy = "proportional")` instead
allocates age-total net migration across compartments by current age-specific
compartment shares; `migration_policy = "error"` allows zero migration but
errors when non-zero migration would require an allocation choice.

The initial SEIR-demography policy follows the same convention: fertility
exposure uses `S + E + I + R`; births enter only the youngest `S`; mortality and
ageing apply independently to `S`, `E`, `I`, and `R`; migration is allocated
by the same `migration_policy`; and `E -> I` progression plus `I -> R` recovery
remain disease-model transitions. The force of infection will continue to
depend on `I`, not `E`. This policy does not add disease-induced mortality,
vaccination, waning immunity, compartment-specific demographic rates, or WPP
projection matching.

See `examples/mock_demographic_workflow.R` for a small dependency-free
demographic-only workflow using invented WPP-like fertility, mortality, and
migration tables. The example demonstrates standardisation, process assembly,
simulation, and diagnostic comparison; it is not a WPP projection reproduction.
See `docs/demographic_residuals.md` for the residual diagnostic and
residual-derived migration schedule conventions.

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

`ContactSchedule()` stores externally supplied contact matrices by time, and
`contact_matrix_at()` retrieves a matrix at an exact available time point. This
prepares agepi for later time-varying simulation work without adding
interpolation, reciprocity correction, population balancing, or simulator
integration.

### External data adapters

`population_from_wpp()`/`demography_from_wpp()` convert WPP-style tidy
population tables into a `Demography()` object. Population values are interpreted
as caller-supplied counts with no scaling. `standardise_wpp_fertility()` accepts
already-computed age-specific fertility rates, while
`fertility_from_wpp_percent_asfr()` converts WPP 2024 `percentASFR1dt`-style
fertility weights to agepi fertility rates using
`fertility_rate = TFR * fraction / age_bin_width`. Percent weights are divided
by 100 first, and TFR is interpreted as births per woman over the reproductive
lifetime. The result is a `FertilitySchedule()` whose values are annual births
per female person-year. `mortality_from_wpp_mx()` converts WPP-style central
death rates (`mx`) into a `MortalitySchedule()` using agepi's `annual_hazard`
convention. Death probabilities (`qx`) and survival probabilities are not
converted because their period and interval conventions require extra metadata.
These are adapter-layer helpers, not a complete WPP projection system. Inputs
should be pre-filtered to one country or location, and maternal age groups must
have finite age-bin widths, so open-ended maternal bins are rejected.

```r
ages <- wpp_age_structure_5year()
percent_asfr <- data.frame(
  year = rep(2020, 7),
  age = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
  percent_asfr = c(5, 20, 30, 25, 15, 4, 1)
)
tfr <- data.frame(year = 2020, tfr = 2.1)

fertility <- fertility_from_wpp_percent_asfr(
  percent_asfr,
  age_structure = ages,
  time_col = "year",
  age_col = "age",
  weight_col = "percent_asfr",
  tfr_data = tfr,
  tfr_time_col = "year",
  tfr_col = "tfr"
)
```

The resulting WPP-style schedule objects can be assembled with the existing
demographic process helper:

```r
population <- population_from_wpp(
  population_table,
  age_structure = ages,
  time_col = "year",
  age_group_col = "age",
  population_col = "population"
)

mortality <- mortality_from_wpp_mx(
  mortality_table,
  age_structure = ages,
  time_col = "year",
  age_col = "age",
  mx_col = "mx"
)

process <- build_demographic_process(
  age_structure = ages,
  fertility_schedule = fertility,
  mortality_schedule = mortality
)

initial_state <- demography_population_at(population, time = 2020)
simulate_demography(
  process = process,
  initial_state = initial_state,
  times = c(2020, 2021),
  time_policy = "step"
)
```

This composition uses existing demographic semantics: exact-time schedule lookup
by default, optional left-continuous step lookup with `time_policy = "step"`,
or bounded linear interpolation of rate-like fertility, mortality, and migration
schedules with `time_policy = "linear"`. WPP-derived annual or five-year rate
schedules can therefore be evaluated stepwise or linearly inside their schedule
range, but this is still not WPP projection matching. Population
`Demography()` accessors remain exact-time only.

`contact_matrix_from_socialmixr()` and
`contact_matrix_from_conmat()` convert socialmixr-like and conmat-style contact
outputs into agepi contact matrices. These adapters are optional: `wpp2024`,
`socialmixr`, and `conmat` are not required for core agepi functionality.

The adapters only reshape and validate supplied data. They do not implement
projection dynamics, interpolation, reciprocity correction, or population
balancing. Contact matrix rows are recipient age groups and columns are source
age groups.

## Design notes

See `docs/age_structured_transmission_design.md`.
