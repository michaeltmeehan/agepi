# Model Conventions

This note collects implementation conventions that are important when extending
or auditing `agepi` models. The README gives the user-facing overview; this file
keeps the sharper details close to the package documentation.

## State-Vector Convention

Numeric state vectors use compartment-major, age-group-minor ordering. For
compartments `c("S", "I", "R")` and age groups `c("0-4", "5-9")`, the vector
order is:

```text
S_0-4, S_5-9, I_0-4, I_5-9, R_0-4, R_5-9
```

`state_vector_to_long()` and `simulate_deterministic()` interpret numeric
vectors by position only. Names on numeric state vectors are ignored when
converting back to long form or simulating.

## Force-Of-Infection Convention

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

Rows of `contact_matrix` are recipient age groups `a`; columns are source age
groups `b`. That is, `contact_matrix[a, b]` gives contacts made by recipient
group `a` with source group `b`.

## Contact Matrices

`validate_contact_matrix()` checks that a contact matrix is numeric, finite,
non-missing, non-negative, square, and optionally has dimensions matching an age
structure.

`as_agepi_contact_matrix()` coerces supported inputs to agepi's
recipient-source matrix convention. Current supported inputs are numeric
matrices, numeric data frames, socialmixr-like lists with a numeric `matrix`
element, and conmat-style long data frames with `age_group_from`,
`age_group_to`, and `contacts` columns. This is dependency-free coercion; agepi
does not depend on `socialmixr` or `conmat`.

`transform_contact_matrix()` aggregates contact matrices from a finer source
age structure to a coarser target age structure when every target age bin is an
exact union of source age bins. It remains a narrow exact aggregation helper
and rejects transformations that would require source-bin splitting or general
rebinning.

Published or generated contact sources use the newer source-layer workflow:
`load_contact_matrix_source()` loads a matrix and provenance metadata on its
native/source age grid, and `adapt_contact_matrix_to_age_structure()` adapts it
to a target `AgeStructure()`. The preferred adaptation method is
`method = "source_band"`: fine-to-coarse aggregation uses
source-grid recipient-population weighting and coarse-to-fine expansion
preserves total contacts with each source band, splitting across nested target
source groups by target-grid population weights when supplied and equally
otherwise. The older `method = "exact"` spelling is
deprecated and retained only as an alias for this source-band assumption.

`ContactSchedule()` stores externally supplied contact matrices by time, and
`contact_matrix_at()` retrieves a matrix at an exact available time point. This
prepares agepi for later time-varying simulation work without adding
interpolation, reciprocity correction, population balancing, or simulator
integration. Contact schedules are not yet consumed directly by
`simulate_deterministic()`.

## Age-Vector Transformations

`aggregate_age_vector()` aggregates a numeric vector from one age structure to a
coarser age structure when each target age bin is an exact union of complete
source age bins.

`segregate_age_vector()` splits a coarse age vector into a finer age structure
when each target age bin is nested inside one source age bin. Splits use
explicit target weights within each source age bin.

`transform_age_vector()` is a convenience wrapper for exact age
transformations, including identity, aggregation, segregation, and mixed exact
transformations when the source and target bins align to a common set of
boundaries. Its `split_method` controls how source bins split across multiple
target bins.

## Demography Tables And Schedule Lookup

`validate_demography_table()` checks a tidy table with `time`, `age_group`, and
`population` columns against an age structure.

`Demography()` stores a validated demography table sorted by time and age-group
order. `demography_times()`, `demography_population_at()`,
`demography_population_vector()`, and `demography_population_table()` provide
exact-time accessors. Population interpolation is deliberately not implemented:
population tables may represent initial conditions, observed trajectories, or
projection targets, and those meanings need different interpolation policies.

Demographic-only helpers cover age grids, ageing operators, fertility,
mortality, migration, demographic process assembly, simulation, comparison, and
residual diagnostics. `simulate_demography()` uses exact-time schedule lookup
with Euler by default and offers opt-in interval-start stepwise lookup via
`time_policy = "step"` or bounded linear rate interpolation via
`time_policy = "linear"`.

Stepwise schedules are left-continuous: a schedule row at `t_i` applies from
`t_i` up to the next schedule time. Linear interpolation applies only to
rate-like fertility, mortality, and migration schedules between available
schedule times; it does not extrapolate and does not interpolate `Demography()`
population access. `migration_count` remains a per-time additive flow rather
than an interval total.

For adaptive `deSolve` runs with demographic schedules, `time_policy = "linear"`
is generally recommended. `time_policy = "step"` gives piecewise-constant
demographic rates, while `time_policy = "exact"` is generally unsuitable unless
all solver evaluation times are schedule times.

## SIR/SEIR-Demography Coupling

`simulate_deterministic()` can optionally use a `DemographicProcess()` for
first-pass SIR, SEIR, or generic compartment-demography coupling with Euler or
`deSolve`. This coupling is not WPP projection matching and does not support
time-varying contact matrices.

For SIR-demography coupling, fertility uses the current total age-specific
infection-state population `S + I + R` as its exposure, births enter only the
youngest susceptible compartment, mortality applies independently to `S`, `I`,
and `R`, and ageing moves `S`, `I`, and `R` independently through the ageing
operator. Net migration, including residual-derived migration schedules,
defaults to allocation entirely to `S`.

This `S`-only migration rule is an allocation convention for age-total net
migration inputs, not a mechanistic model of who moves while infected or
recovered. `simulate_deterministic(migration_policy = "proportional")` instead
allocates age-total net migration across compartments by current age-specific
compartment shares. `migration_policy = "error"` allows zero migration but
errors when non-zero migration would require an allocation choice.

The SEIR-demography policy follows the same convention: fertility exposure uses
`S + E + I + R`; births enter only the youngest `S`; mortality and ageing apply
independently to `S`, `E`, `I`, and `R`; migration is allocated by the same
`migration_policy`; and `E -> I` progression plus `I -> R` recovery remain
disease-model transitions. The force of infection depends on `I`, not `E`.

This policy does not add disease-induced mortality, vaccination, waning
immunity, compartment-specific demographic rates, or WPP projection matching.

## WPP-Style Data Semantics

`population_from_wpp()` and `demography_from_wpp()` convert WPP-style tidy
population tables into a `Demography()` object. Population values are
interpreted as caller-supplied counts with no scaling.

`standardise_wpp_fertility()` accepts already-computed age-specific fertility
rates, while `fertility_from_wpp_percent_asfr()` converts WPP 2024
`percentASFR1dt`-style fertility weights to agepi fertility rates using:

```text
fertility_rate = TFR * fraction / age_bin_width
```

Percent weights are divided by 100 first, and TFR is interpreted as births per
woman over the reproductive lifetime. The result is a `FertilitySchedule()`
whose values are annual births per female person-year.

`mortality_from_wpp_mx()` converts WPP-style central death rates (`mx`) into a
`MortalitySchedule()` using agepi's `annual_hazard` convention. Death
probabilities (`qx`) and survival probabilities are not converted because their
period and interval conventions require extra metadata.

These are adapter-layer helpers, not a complete WPP projection system. Inputs
should be pre-filtered to one country or location, and maternal age groups must
have finite age-bin widths, so open-ended maternal bins are rejected.
