# Coupled Epidemic-Demography Operator Splitting Design

## Motivation

`agepi` has two demographic implementations:

- the original continuous derivative path, where ageing, fertility, mortality,
  and migration contribute to `dN/dt`;
- the annual-cohort path, where a full 1-year age grid is advanced by one
  discrete calendar year with survival, ageing, births, and migration.

The derivative path is useful for smooth ODE coupling and remains the default.
The annual-cohort path is closer to cohort projection inputs such as WPP-style
one-year age schedules. Coupling it to deterministic epidemic models needs an
operator split rather than another derivative term, because cohort ageing is a
discrete annual boundary operation.

## Current Architecture

The relevant pieces already exist:

- `simulate_deterministic()` maps long-form compartment-age state into the
  package's compartment-major vector order, evaluates transition rates, and
  integrates either by Euler or by the suggested `deSolve` backend.
- `compartment_demographic_derivative()` adds continuous demographic derivatives
  compartment-wise for the existing coupled path.
- `simulate_demography(ageing_policy = "annual_cohort")` applies
  `annual_cohort_demographic_step()` once per annual interval for a
  one-compartment population.
- `annual_cohort_demographic_step()` validates a complete 1-year age grid,
  applies annual-hazard survival, ages survivors, adds births, and applies
  count- or rate-based net migration.
- `AgeGridMapping()` and the output aggregation wrappers aggregate internal
  1-year outputs back to a reporting grid without changing totals.

## Proposed API

`simulate_deterministic()` gains two optional arguments:

```r
ageing_policy = c("exponential", "annual_cohort")
output_age_structure = NULL
```

`ageing_policy = "exponential"` is the default and preserves all existing
continuous derivative behavior. `ageing_policy = "annual_cohort"` opts into
operator splitting. If `demographic_process` is supplied, it must use the same
complete 1-year internal age grid as `age_structure`, starting at age 0 and
ending in a terminal open-ended age group. Requested output times must be spaced
exactly one year apart.

`output_age_structure` is optional. When supplied, trajectory and cumulative-flow
outputs are aggregated from the simulation grid to that reporting grid using the
existing age-grid mapping helpers. Reporting bins must be exact unions of the
internal bins.

## Algorithm

For annual times `t_0, t_1, ..., t_n`, with `t_{i+1} - t_i = 1`:

1. Start with the compartment-age state at `t_i`.
2. Integrate epidemic transition ODEs only over `[t_i, t_{i+1})`.
3. At `t_{i+1}`, apply the annual-cohort demographic operator to the ordinary
   compartment state.
4. Emit the post-demography state at `t_{i+1}`.
5. Repeat.

Cumulative disease flows are integrated during the epidemic interval and are not
changed by the demographic boundary operator.

## State Representation

The ordinary state remains the existing deterministic vector order:

```text
S_age1, S_age2, ..., I_age1, I_age2, ...
```

Auxiliary cumulative-flow states, when requested, are appended after the ordinary
state exactly as in the existing deterministic implementation. The annual
demographic operator only mutates the ordinary compartment state.

## Age-Grid Rules

Annual-cohort demographic coupling with a `demographic_process` requires:

- contiguous age groups;
- finite 1-year bins before the terminal group;
- a first age group starting at 0;
- a final open-ended age group;
- matching `age_structure` and `demographic_process$age_structure`;
- annual output times.

Broader reporting groups are supported after simulation with
`output_age_structure`, not as the internal residence-time grid.

## Births, Deaths, Ageing, and Migration

The compartment-wise helper applies the annual cohort operation in three parts:

- mortality and ageing are applied independently to every disease compartment;
- births are computed from total age-specific population exposure and added only
  to the configured birth compartment, usually `S`;
- net migration is computed as an age-total annual-cohort quantity and allocated
  by the existing `migration_policy`.

The policies are:

- `susceptible`: put all net migration in the model migration compartment,
  usually `S`;
- `proportional`: allocate net migration by post-survival/post-ageing/post-birth
  age-specific compartment shares;
- `error`: reject non-zero net migration.

Mortality uses the annual-cohort convention `survival = exp(-mortality_rate)`.
Disease-induced mortality is not part of this demographic operator.

## Contact and Rate Handling

During each annual interval, infection and disease-transition rates are evaluated
through the existing deterministic transition-rate pathway. Contact matrices,
`beta`, susceptibility, infectiousness, and generic transition rates keep their
current semantics.

The first implementation assumes static contact and disease-rate inputs over
the simulated horizon. Schedule lookup for demographic fertility, mortality,
and migration occurs at the interval start using `time_policy`.

## Cumulative Flows

Cumulative flows track disease transitions only. They are integrated over the
continuous epidemic interval and then carried unchanged through the demographic
boundary. Births, background deaths, ageing movements, and migration are not
currently exposed as cumulative-flow counters.

## Outputs

Without `output_age_structure`, output is on the internal simulation age grid.
With `output_age_structure`, ordinary deterministic trajectories and cumulative
flows are aggregated with:

- `aggregate_epidemic_trajectory_age_grid()`;
- `aggregate_cumulative_flows_age_grid()`.

Population summaries made from the returned trajectory, such as
`age_group_totals()` and `total_population()`, therefore inherit the reporting
grid when output aggregation is requested. The aggregation is count-preserving
and preserves non-age metadata columns.

## Validation

The test strategy covers:

- default `simulate_deterministic()` behavior unchanged;
- targeted compartment-wise annual cohort births, mortality, and ageing;
- zero-epidemic annual-coupled outputs matching
  `simulate_demography(ageing_policy = "annual_cohort")`;
- no-demography split runs matching the existing deterministic solver;
- SIR/SEIR, age-structured infection, and cumulative-flow equivalence where
  relevant;
- toy SIR with annual demography sanity checks;
- output aggregation conservation and metadata preservation.

## Implementation Milestones

1. Add this design document.
2. Add an internal compartment-wise annual cohort demographic helper.
3. Add an opt-in annual-interval deterministic split loop.
4. Preserve the existing derivative path as default.
5. Validate zero-epidemic annual-cohort coupling against standalone demography.
6. Validate split infection dynamics with demography absent.
7. Add a toy SIR plus annual demography sanity example.
8. Aggregate deterministic outputs back to a reporting grid.
9. Update user-facing docs, examples, and roadmap.

## Deferred Features

The first implementation intentionally does not include:

- stochastic annual-cohort coupling;
- subannual output inside an annual split interval;
- time-varying contact matrices;
- demographic event cumulative counters;
- sex-structured fertility exposure;
- disease-induced mortality in the demographic operator;
- residual forcing or WPP projection matching;
- automatic expansion of epidemic inputs from a reporting grid to the internal
  1-year grid.
