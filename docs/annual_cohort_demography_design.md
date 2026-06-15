# Annual Cohort Demography Design

This note sketches an explicit annual cohort demography mode for `agepi`. It is
a design document only. It does not propose changing the current continuous
ageing implementation in place, and it does not treat the existing
`AgeingOperator()` derivative path as a WPP-style cohort projection mechanism.

## 1. Motivation

The current demographic process represents ageing as an ODE-compatible
compartmental outflow. `AgeingOperator()` computes finite-bin departure rates as

```text
departure_rate = 1 / age_bin_width
```

where `age_bin_width` is the inclusive width of the age group. The final
open-ended age group has zero departure and no destination. In
`demographic_derivative()` and `compartment_demographic_derivative()`, that
departure rate is applied continuously: population leaves one age bin and enters
the next at each derivative evaluation.

This rate convention gives a mean residence time equal to the bin width. It does
not make a cohort remain in the bin for exactly that width and then exit as a
block. Under a continuous solver, residence times are exponential. After one
bin width, a fraction `exp(-1)` remains in the starting bin. Under explicit
Euler, the result depends on the step size: annual Euler steps with 5-year bins
leave `(1 - 1 / 5)^5` in the starting bin after five years, while 1-year bins
with 1-year Euler steps empty the first bin only as a numerical artefact of the
step formula.

That distinction matters for WPP-style cohort projection. A population total can
match well while the age distribution remains wrong, because births and deaths
can balance in aggregate even as cohorts leak too slowly through broad age
groups. For WPP Kiribati total population, the fertility exposure correction
now improves the total trajectory, but age-structure diagnostics still show
material errors:

- 2030 total absolute age-distribution error: `5.988668`; max relative
  age-group error: `0.1378459`
- 2040 total absolute age-distribution error: `9.828091`; max relative
  age-group error: `0.2183497`
- 2050 total absolute age-distribution error: `10.632588`; max relative
  age-group error: `0.2374219`

The same diagnostic confirms the ageing mechanism directly:

- `deSolve` with one bin width leaves approximately `exp(-1) = 0.3678799` of a
  cohort in the starting bin.
- Euler with 5-year bins leaves `(1 - 1 / 5)^5 = 0.32768`.
- Euler with 1-year bins empties the bin only because the numerical step is
  exactly equal to the reciprocal departure rate.

The design implication is that WPP-style demographic projection should use an
explicit cohort update, not a patched continuous ageing rate.

## 2. Design Principle

Annual cohort demography should use operator splitting:

```text
state at start of year
  -> continuous epidemic ODE over [year, year + 1)
  -> annual demographic operator
  -> state at start of next year
```

The separation is:

```text
continuous epidemic process within annual intervals
+
discrete demographic cohort update at annual boundaries
```

In this design, demography is a distinct annual operator. It is not another
term added to the ODE derivative. The existing derivative-based demographic
path remains a continuous exponential ageing policy. The proposed cohort path
uses discrete survival, ageing, births, and optional migration at known annual
boundaries.

This framing also makes the coupled epidemic-demography semantics clearer:
infection, recovery, progression, and other epidemic transitions can remain
continuous between demographic updates, while demographic projection occurs at
calendar-year boundaries with deterministic cohort ageing.

## 3. Internal Age Grid

Annual cohort mode should use an internal demographic grid with:

- 1-year age bins, such as `0`, `1`, `2`, ..., up to a configured maximum finite
  age.
- A final open-ended age group, such as `100+`.
- Optional broader user-facing reporting groups.

The internal grid is the grid on which the demographic operator acts. User
inputs and outputs may still use 5-year or arbitrary age groups when a valid
mapping exists. For example, a WPP-style reporting grid can remain `0-4`,
`5-9`, ..., `100+`, while the internal state uses `0`, `1`, ..., `99`, `100+`.
Simulation output can include the internal 1-year trajectory and an aggregated
trajectory on the requested reporting grid.

Using 1-year internal bins avoids the exponential residence-time problem because
ageing is not represented by a continuous outflow. At the annual boundary,
survivors in age `a` are moved deterministically to age `a + 1`. Broader
reporting groups are then sums over complete internal ages, so 5-year WPP-style
reporting remains available without using 5-year bins as residence-time
compartments.

## 4. Annual Demographic Operator

The annual operator maps the state at the start of one year to the state at the
start of the next year. For a standalone population state, this is a vector by
single-year age. For a coupled epidemic state, it is a compartment by age state
using the existing compartment-major, age-group-minor convention internally.

### Deterministic Ageing Shift

For finite single-year ages:

```text
survivors at age a -> age a + 1
```

For the final open-ended age group:

```text
next final group =
  survivors already in final group
  + survivors from the previous finite age
```

Ageing alone should conserve population exactly, apart from numerical
round-off.

### Births

Births are computed from fertility rates and exposure. The existing WPP adapter
semantics should be preserved: `fertility_from_wpp_percent_asfr()` returns
annual births per female person-year, while a total-population simulation can
use `fertility_exposure_fraction = 0.5` as an approximate female exposure
correction.

Births are inserted into the youngest age group. In epidemic models, the birth
destination compartment must be explicit. For SIR and SEIR models this is
normally the susceptible compartment, but generic compartment models should use
the configured birth compartment rather than assuming the name `S` always
exists.

### Mortality

Background mortality should be applied as annual survival probabilities on the
internal 1-year grid. If source inputs are central death rates or hazards, the
operator should convert them to annual survival using a documented convention.
For example, an annual hazard `m_a` can map to survival

```text
S_a = exp(-m_a)
```

when rates are interpreted as constant hazards over the year. If source data
are death probabilities or survival probabilities, the period and interval
conventions should be explicit before conversion.

Recommended default ordering:

```text
apply survival over the year
-> age survivors
-> add births to age 0
-> add net migration if supplied
```

This convention treats mortality as survival during the interval lived by the
cohort before it enters the next age. Births are added at the end of the annual
interval as the age 0 population at the start of the next year. It is simple,
cohort-oriented, and avoids mixing newborns into the same year's mortality and
fertility exposure unless that is explicitly requested.

Two alternatives should remain documented:

- Add births before survival when input fertility represents births at the
  beginning of the interval and infant mortality should apply in the same
  annual operator.
- Apply half-year or midpoint exposure conventions when matching a specific
  demographic projection system. This would be a later refinement, not the
  first implementation default.

WPP-style infant and age 1 mortality require special care. Current examples use
WPP age `0` mortality as an approximation for `0-4` and omit age `1` on a
5-year grid. A 1-year internal grid should allow age `0` and age `1` mortality
to be represented separately when those inputs are available.

### Migration

Migration should be optional. The initial annual cohort design should treat
migration as annual net migration by age when supplied:

```text
next_state = state_after_births + net_migration_by_age
```

The representation must document whether values are interval totals, rates, or
residual-like adjustments. For coupled epidemic models, age-total migration
must also specify an allocation rule across compartments, such as susceptible
only, proportional, or error on non-zero migration, following the existing
`migration_policy` vocabulary.

## 5. Standalone Demography Mode

`simulate_demography()` could grow an explicit ageing policy option:

```text
ageing_policy = "exponential"
ageing_policy = "annual_cohort"
```

The initial default should remain `"exponential"` for backwards compatibility.
That default would preserve the current derivative pathway through
`demographic_derivative()`, `AgeingOperator()`, Euler, and optional `deSolve`.

When `ageing_policy = "annual_cohort"`, `simulate_demography()` should use the
discrete annual operator rather than the continuous derivative. This mode should
probably require annual update times, or at least require that demographic
operator times are an integer-year sequence. If users request non-annual output,
the design should decide whether to return only annual demographic states or to
interpolate/report between annual boundaries with no demographic change.

Annual cohort mode should require or internally construct a 1-year age grid. A
reasonable first milestone is:

- Accept a 1-year process grid directly.
- Add helper utilities to expand 5-year WPP-style inputs to the internal grid.
- Return output on the internal grid, with optional aggregation to the
  user-facing reporting grid.

WPP examples should eventually use `"annual_cohort"` once the helper utilities
and validation tests are in place. Until then, the examples should continue to
state that the current exponential pathway is diagnostic and not WPP projection
matching.

## 6. Coupled Epidemic-Demography Mode

Coupled deterministic simulations are a larger milestone than standalone
demography because they require two state representations to agree: the
epidemic ODE state and the annual demographic operator state.

The proposed coupled flow is:

```text
for each annual interval:
  solve epidemic ODE over [year, year + 1) with deSolve
  reshape final epidemic state by compartment and 1-year age
  apply annual demographic operator separately by compartment and age
  continue from the updated state at the next year
```

Epidemic transitions remain continuous between demographic update times.
`deSolve` should remain the continuous epidemic solver. At each annual boundary,
the demographic operator is applied outside the derivative path.

The coupled design needs clear conventions:

- Background mortality and disease mortality must be separate. Background
  mortality belongs to the annual demographic operator. Disease mortality, if
  later supported, belongs to the epidemic model or a clearly named disease
  transition.
- Births enter the configured birth compartment in the youngest internal age
  group.
- Ageing shifts each epidemic compartment independently unless a model defines
  a special cross-compartment demographic transition.
- Fertility exposure should use age-specific total population across epidemic
  compartments, scaled by `fertility_exposure_fraction` where appropriate.
- Contact matrices, susceptibility vectors, infectiousness vectors, and other
  age-specific epidemic rates must either be supplied on the internal 1-year
  grid or expanded to that grid using explicit assumptions.
- Outputs should be aggregatable back to the user-facing age groups after each
  requested output time.

This milestone also needs careful treatment of cumulative flow state variables.
If cumulative disease flows are tracked continuously, the annual demographic
operator should not age or kill cumulative counters as though they were living
population compartments.

## 7. Relationship To Euler

Euler should not be relied on as a demographic cohort mechanism.

With 1-year bins and annual steps, Euler can mimic deterministic cohort ageing
because the update subtracts `1 * state` from the starting bin and adds it to
the next bin. That is an accident of the numerical step size and the current
`1 / width` rate, not a demographic design.

With 5-year bins and annual steps, Euler still gives partial leakage:

```text
(1 - 1 / 5)^5 = 0.32768
```

remains in the starting bin after five years. With `deSolve`, the same
continuous rate gives the exponential result:

```text
exp(-1) = 0.3678799
```

`deSolve` should remain the continuous epidemic solver. Cohort demographic
updating should be an explicit discrete operator, not an artefact of Euler
stepping. If annual cohort mode becomes the recommended demographic path,
Euler support for demographic projection may eventually be narrowed or
documented as a legacy continuous-ageing approximation rather than a cohort
projection method.

## 8. Age-Grid Expansion And Aggregation Utilities

Annual cohort mode will need explicit mapping utilities. Likely helpers are:

`expand_age_grid_to_single_year()`

: Construct a 1-year internal age structure and mapping from a broader source
  or reporting age structure. It should preserve the final open-ended age group
  and reject ambiguous non-contiguous grids.

`expand_population_to_single_year()`

: Split broad age-group population counts onto the internal grid. This requires
  explicit weights within each source age group. Uniform splitting may be a
  convenience option, but WPP validation should prefer data-derived weights
  where available.

`expand_rates_to_single_year()`

: Expand fertility, mortality, susceptibility, infectiousness, or other
  age-specific vectors from broader bins to single-year ages. The helper must
  distinguish rates from counts. It should document whether broad-bin rates are
  copied to each internal age, converted by interval width, or rejected without
  metadata.

`expand_contact_matrix_to_single_year()`

: Expand a reporting-grid contact matrix to the internal 1-year grid. This is
  inherently assumption-heavy. Possible first policies include block-constant
  contacts within broad source and recipient bins, or requiring the caller to
  supply an already expanded matrix.

`aggregate_single_year_population()`

: Aggregate a single-year population vector or compartment by age matrix to a
  reporting age structure. Reporting bins must be exact unions of internal ages.

`aggregate_single_year_trajectory()`

: Aggregate a trajectory data frame from the internal grid to reporting age
  groups across all output times and, for epidemic outputs, all compartments.

`validate_age_grid_mapping()`

: Validate that source, internal, and reporting grids are contiguous, nested as
  required, and compatible with the final open-ended age group.

WPP-specific complications include:

- WPP mortality distinguishes age `0` and age `1` conventions in some tables.
  A 1-year internal grid should preserve those where available.
- WPP population and reporting tables are often supplied in 5-year age groups,
  so expansion needs within-bin weights or a documented approximation.
- The final open-ended age group must retain survivors and receive survivors
  from the previous finite age without needing a finite width.
- WPP fertility rates from `fertility_from_wpp_percent_asfr()` are annual
  births per female person-year. Total-population states still need
  `fertility_exposure_fraction`, or a later sex-structured model.

## 9. Testing And Validation Plan

Targeted tests should cover:

- Isolated single-cohort shift with no births, deaths, or migration.
- Final open age group retention.
- Population conservation under ageing alone.
- Birth-only updates into the youngest age group.
- Mortality-only updates using annual survival probabilities.
- Fertility exposure tests, including WPP female ASFR with total-population
  exposure and `fertility_exposure_fraction = 0.5`.
- WPP total population benchmarks.
- WPP age-structure benchmarks, especially the Kiribati diagnostics that show
  remaining age-distribution errors under the current continuous path.
- Coupled epidemic-demography regression tests once operator splitting is
  introduced.
- Aggregation consistency tests showing that internal 1-year totals match
  reporting-grid totals after aggregation.

Tests should include exact conservation cases before using WPP benchmarks, so
failures can be localized to the operator, expansion assumptions, or source
data conventions.

## 10. Implementation Milestones

1. Design and helper utilities only.
2. Standalone annual cohort `simulate_demography()` prototype.
3. WPP benchmark updated to annual cohort mode.
4. Tests for cohort demography and age aggregation.
5. Coupled deterministic epidemic-demography operator splitting.
6. Documentation, examples, and possible Euler deprecation discussion.

## 11. Open Questions

- Should mortality be applied before ageing, after ageing, or via annual
  survival during the interval? This note recommends annual survival before
  ageing as the first default, but WPP matching may motivate more detailed
  infant and midpoint conventions.
- Should annual cohort mode require 1-year input data, or automatically expand
  broader age groups? Automatic expansion is ergonomic but can hide strong
  assumptions.
- Should the package store both internal and reporting age grids in the
  demographic process object, or keep reporting aggregation as a simulator
  output option?
- How should contact matrices be expanded to 1-year age bins? Block-constant
  expansion is simple but may be too crude for some epidemic models.
- How should migration be represented: annual net counts, rates, residual
  adjustments, or a richer in/out migration structure?
- Should full sex-structured demography be a later extension? It would remove
  the need for approximate female exposure fractions in WPP fertility, but it
  would also expand the state dimension and coupling rules substantially.
- How should non-annual epidemic output times be handled in coupled mode when
  demography updates only at annual boundaries?
- How should cumulative flow states interact with annual demographic updates in
  coupled simulations?

