# External Data Adapters Design Note

## Purpose

This note records the intended boundary for future external data adapters. The
adapter layer can connect agepi to WPP-style demographic data, socialmixr
contact matrices, and conmat contact predictions without changing the
deterministic SIR/SEIR simulation workflow.

External package integrations should remain optional initially. The core package
should continue to accept agepi-native objects and dependency-free inputs.

## Adapter Scope

### WPP Population Adapter

`population_from_wpp()`/`demography_from_wpp()` convert external WPP-like age,
year, and population tables into a validated `Demography()` object aligned to an
`AgeStructure()`.

The adapter should handle data-shape translation only:

- map external age labels or bounds to agepi age groups;
- map external year or period fields to agepi `time`;
- map population values to the required `population` column;
- call existing demography validation and construction helpers.

Population inputs are expected in tidy long format with one row per time and age
group. Values are interpreted as population counts in the caller's supplied
units; no scaling is applied. Inputs should be pre-filtered to one country,
location, or entity, or selected with the adapter's `location` and
`location_col` arguments. Multi-location inputs are not silently mixed.

The adapter does not implement demographic projection dynamics, interpolation,
or WPP projection matching. Separate WPP-style standardisers can adapt already
supplied fertility, mortality, and migration tables into agepi schedule objects
without requiring `wpp2024`.

WPP 2024 `percentASFR1dt`-style fertility schedules are weights over maternal
age groups, not direct per-capita fertility rates. `fertility_from_wpp_percent_asfr()`
converts those weights plus a time-specific TFR table into a `FertilitySchedule()`
using:

```r
fertility_rate = TFR * fraction / age_bin_width
```

Percent inputs are divided by 100 before conversion. TFR is interpreted as
births per woman over the reproductive lifetime, and the output schedule uses
agepi's `births_per_female_person_year` convention. The adapter validates that
the maternal-age weights sum to 1 for fractions or 100 for percentages at each
time point. Inputs should be pre-filtered to one country or location, and
open-ended maternal bins are rejected because conversion requires finite
age-bin widths. This remains a data-shape and unit-conversion adapter, not a
WPP projection engine.

WPP-style mortality inputs should use `mortality_from_wpp_mx()` when the input
contains age-specific central death rates (`mx`). These values are interpreted
as annual per-capita rates and stored in `MortalitySchedule()` using agepi's
`annual_hazard` convention, which is what `demographic_derivative()` expects.
`mortality_from_wpp()` currently supports only `quantity = "mx"`; age-specific
death probabilities (`qx`) and survival probabilities are rejected because their
period and interval conventions are ambiguous without additional metadata.
Mortality inputs should also be pre-filtered to one country, location, or entity.

WPP-derived schedule objects can be composed directly with the existing
demographic process framework. `build_demographic_process()` is the narrow
assembly helper for already-standardised schedules; it validates that supplied
fertility, mortality, and migration schedules share the requested
`AgeStructure()`, and then returns a `DemographicProcess()` usable by
`demographic_derivative()` and `simulate_demography()`.

```r
ages <- wpp_age_structure_5year()

population <- population_from_wpp(
  population_table,
  age_structure = ages,
  time_col = "year",
  age_group_col = "age",
  population_col = "population",
  location = "Exampleland",
  location_col = "location"
)

fertility <- fertility_from_wpp_percent_asfr(
  percent_asfr_table,
  age_structure = ages,
  time_col = "year",
  age_col = "age",
  weight_col = "percent_asfr",
  tfr_data = tfr_table,
  tfr_time_col = "year",
  tfr_col = "tfr"
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
demographic_derivative(initial_state, time = 2020, process = process)
simulate_demography(
  process = process,
  initial_state = initial_state,
  times = c(2020, 2021),
  time_policy = "step"
)
```

Schedule time grids do not have to be identical at construction time. When all
supplied schedules use the same grid, the process records it in `process$times`;
otherwise `process$times` is `NULL`. Evaluation follows the existing
demographic time-policy rules: exact lookup by default, or left-continuous
interval-start lookup with `time_policy = "step"`, or bounded linear
interpolation with `time_policy = "linear"`. Linear interpolation applies only
to rate-like fertility, mortality, and migration schedules, and only inside the
supported schedule range. WPP-derived annual or five-year rate schedules can be
used with either stepwise or linear rate interpolation, but nearest-year lookup,
qx conversion, survival-probability conversion, residual migration fitting, and
WPP projection dynamics are not implemented by these adapters.

Population `Demography()` accessors remain exact-time only. Population tables
can mean externally supplied initial states, observed trajectories, or
projection targets, so population interpolation needs a separate explicit
state-trajectory policy rather than the rate-schedule policy above.

### socialmixr Contact Adapter

A socialmixr adapter should convert contact-matrix outputs into agepi contact
matrices using agepi's recipient-row, source-column convention.

The first adapter should focus on ordinary socialmixr-style outputs that contain
a numeric `matrix` component. It should validate matrix dimensions and age-group
labels before returning a matrix suitable for existing agepi contact-matrix
helpers.

The adapter should preserve the matrix supplied by socialmixr. It should not add
or redo reciprocity correction, population balancing, age weighting, symmetry,
or per-capita transformations inside agepi.

Published or bundled source matrices should be loaded separately from any model
age-grid adaptation. `load_contact_matrix_source()` returns a source matrix, its
native `AgeStructure()`, and metadata. `adapt_contact_matrix_to_age_structure()`
then performs explicit exact adaptation: fine-to-coarse aggregation delegates to
`transform_contact_matrix()` and requires source-grid population weights, while
coarse-to-fine expansion copies source-band contacts to nested target ages. The
deprecated `contact_matrix_for_age_structure()` wrapper remains only for
compatibility with older examples.

### conmat Contact Adapter

A conmat adapter should convert generated or predicted contact matrices into
agepi contact matrices using the same recipient-row, source-column convention.

The first adapter can support conmat-style square matrices and long prediction
tables with `age_group_from`, `age_group_to`, and `contacts` columns. Long
tables should be reshaped so `age_group_to` becomes recipient rows and
`age_group_from` becomes source columns, matching the existing agepi convention.

The adapter should treat conmat modelling choices as upstream choices. It should
not refit models, apply household-size adjustments, add reciprocity correction,
or alter generated predictions.

## Core Workflow Boundary

External data adapters should prepare inputs for the existing model workflow.
They should not change:

- `force_of_infection()`;
- `simulate_deterministic()`;
- the deterministic SIR/SEIR state-vector convention;
- the current static-contact-matrix simulation assumptions.

Adapters should return objects that existing validation and simulation code can
already consume.

## Remaining Out Of Scope

The following remain out of scope for the initial adapter phase:

- mandatory WPP, socialmixr, or conmat dependencies;
- WPP projection matching;
- population interpolation or nearest-time population lookup;
- qx conversion or survival-probability conversion;
- residual migration fitting;
- additional infection-demography coupling beyond the current first-pass SIR/SEIR
  coupling;
- automatic demographic residual forcing;
- contact-matrix reciprocity correction;
- contact-matrix population balancing;
- source-bin splitting or general contact-matrix rebinning;
- stochastic simulation, vaccination, waning immunity, or event handling;
- changes to `force_of_infection()` semantics.

Residual-derived migration schedules remain age-total demographic inputs. In
the current SIR-demography coupling, those net migration values are allocated to
`S` only. The first SEIR-demography policy keeps the same `S`-only
allocation convention; proportional allocation across infection compartments is
out of scope for the adapter layer unless a future simulation option requests it
explicitly.
