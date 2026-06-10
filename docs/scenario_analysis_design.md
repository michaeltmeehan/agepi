# Scenario Analysis Design

## Executive Summary

Scenario analysis should be implemented as a thin orchestration layer around
the existing simulators. The core solvers should continue to know only about
states, times, models, demographic processes, contact matrices, and cumulative
flow requests. Scenario utilities should prepare common arguments, extract a
branch state from an already completed trajectory, apply named overrides, call
`simulate_deterministic()` or `simulate_stochastic()`, and add scenario
metadata to the returned tables.

The first milestone should target deterministic simulations. It should support
the common workflow:

1. run a historical, calibration, or burn-in trajectory;
2. extract the ordinary compartment state at an exact branch time;
3. run named projection scenarios from that shared state;
4. combine trajectory and cumulative-flow tables with `scenario` metadata;
5. compute simple differences from a named baseline scenario.

The design should avoid disease-specific concepts. TB outcomes such as cases
averted, treatment initiations, deaths averted, prevalence reduction, and
cumulative incidence should be expressible as user-defined metrics over
trajectory and cumulative tables, not hard-coded into the scenario engine.

## Scenario-Analysis Use Cases

The scenario layer should support these general use cases:

- A common burn-in followed by several alternative projections.
- Parameter sweeps where each scenario changes one or more scalar or vector
  arguments such as `beta`, susceptibility, infectiousness, or transition-rate
  values embedded in a model object.
- Alternative model structures, such as adding or removing a transition, when
  the resulting model is compatible with the branch state.
- Alternative demographic schedules or contact matrices for the projection
  period.
- Different cumulative-flow requests for reporting, while keeping ordinary
  trajectory output compatible with the existing summary helpers.
- Scenario comparisons against a named baseline using final values, changes
  over an interval, cumulative totals, or user-supplied summary functions.
- Future stochastic scenario sets with multiple replicates per scenario and
  replicate-level summaries.

The package should not try to encode every policy-analysis convention. It
should make the common orchestration reliable and leave disease-specific
outcome definitions to caller-supplied metric functions.

## Current Package Capabilities Relevant To Scenarios

The current simulation surface is already scenario-friendly:

- `simulate_deterministic()` accepts long-form or numeric initial states and
  returns the initial state at the first requested time.
- Ordinary deterministic output is a data frame with `time`, `compartment`,
  `age_group`, and `value`.
- With `cumulative_flows`, deterministic output is a list with `trajectory` and
  `cumulative`.
- `simulate_stochastic()` has the same ordinary trajectory schema and can
  return `events` and `cumulative` tables.
- Cumulative-flow tables use `time`, `cumulative_name`, `transition_id`,
  `from`, `to`, `age_group`, and `value`.
- `state_long_to_vector()` and `state_vector_to_long()` define the
  compartment-major, age-group-minor ordering.
- `compartment_totals()`, `age_group_totals()`, and `total_population()` work
  on ordinary trajectory tables and should continue to ignore cumulative
  outputs unless users pass them explicitly.
- Deterministic cumulative flows are auxiliary outputs, not ordinary model
  compartments. Stochastic cumulative flows are derived from realised event
  logs.

These contracts imply that scenario utilities can be written without changing
the solver internals.

## Proposed Minimal API

The minimal public API should be:

```r
Scenario(
  name,
  description = NULL,
  overrides = list(),
  modifier = NULL,
  metadata = list()
)

ScenarioSet(
  scenarios,
  baseline = NULL,
  common_args = list(),
  branch_time = NULL,
  cumulative_policy = c("reset", "continue")
)

initial_state_from_simulation(
  simulation,
  time,
  model = NULL,
  age_structure = NULL,
  tolerance = 0
)

run_scenarios(
  scenario_set,
  simulator = c("deterministic", "stochastic"),
  historical = NULL,
  branch_time = scenario_set$branch_time,
  times = NULL,
  ...
)

combine_scenario_outputs(outputs)

summarise_scenarios(output, metrics)

compare_scenarios(summary, baseline = NULL)
```

For the first implementation, only `Scenario()`, `ScenarioSet()`,
`initial_state_from_simulation()`, `run_scenarios()`, and
`combine_scenario_outputs()` are essential. `summarise_scenarios()` and
`compare_scenarios()` can start small or remain documented patterns until the
output shape is proven by examples.

## State Extraction At Branch Time

`initial_state_from_simulation()` should extract ordinary compartment rows from
a completed simulation at an exact time. It should accept either:

- a trajectory data frame returned by `simulate_deterministic()` or
  `simulate_stochastic()`;
- a simulation list with a `trajectory` element.

The returned value should be a long-form state data frame with columns
`compartment`, `age_group`, and `value`, suitable for passing directly as
`initial_state`.

Recommended first-milestone rules:

- Require `time` to match exactly by default.
- Allow a numeric `tolerance` only for floating-point defensive matching.
- Reject branch times that are absent, duplicated unexpectedly, or incomplete.
- Validate one row per compartment-age cell when `model` and `age_structure`
  are supplied.
- Do not interpolate compartment states in the first milestone.

Exact extraction matches current package conventions: demographic population
accessors and contact schedules also prefer exact-time semantics unless a
specific interpolation policy exists.

## Scenario Representation

`Scenario()` should return a plain list with class `agepi_scenario`. A scenario
should contain:

- `name`: unique, non-empty character scalar;
- `description`: optional character scalar;
- `overrides`: named list of simulation argument replacements;
- `modifier`: optional function for advanced changes;
- `metadata`: optional named list copied to output metadata when practical.

`overrides` should be shallow replacements of common simulation arguments. For
example:

```r
Scenario("reduced_transmission", overrides = list(beta = 0.04))
Scenario("new_contacts", overrides = list(contact_matrix = projected_contacts))
Scenario("different_reporting", overrides = list(cumulative_flows = flows))
```

The optional `modifier` should be called with the resolved argument list and
the scenario object, and should return a modified argument list:

```r
modifier <- function(args, scenario) {
  args$beta <- args$beta * 0.8
  args
}
```

This gives advanced users a general hook without forcing the package to
invent override semantics for every object type.

## Scenario Set Representation

`ScenarioSet()` should return a plain list with class `agepi_scenario_set`.
It should contain:

- `scenarios`: named list of `Scenario()` objects;
- `baseline`: optional scenario name used by comparison helpers;
- `common_args`: named list of arguments shared across scenarios;
- `branch_time`: optional numeric scalar;
- `cumulative_policy`: default handling for cumulative counters after branch;
- `metadata`: optional list for project-level annotations.

The scenario names should be unique and should be used as stable output
identifiers. `baseline` should be a scenario in `scenarios`; it is the
projected baseline after the branch, not necessarily the historical trajectory.

## Scenario Execution

`run_scenarios()` should build one simulator call per scenario:

1. Resolve the branch state. If `historical` is supplied, call
   `initial_state_from_simulation(historical, branch_time, ...)`; otherwise
   require `common_args$initial_state` or a top-level `initial_state`.
2. Merge `common_args`, the extracted `initial_state`, and the projection
   `times`.
3. Apply scenario `overrides`.
4. Apply scenario `modifier`, if supplied.
5. Call `simulate_deterministic()` for the first milestone.
6. Normalize the output to a list with `trajectory`, optional `cumulative`,
   optional `events`, and metadata.
7. Add `scenario` and `scenario_description` columns to each returned table.

The scenario runner should not mutate the objects in `common_args`. It should
copy the resolved list and let each scenario return its own arguments and
outputs.

Projection `times` need an explicit convention. The simplest first milestone
is to require users to pass absolute simulation times beginning at
`branch_time`. The extracted branch state becomes the initial state and the
first projection time must equal `branch_time`. This preserves the existing
meaning of time in demographic schedules, contact schedules, and cumulative
flow output.

Relative projection time could be added later through an explicit
`time_origin = c("absolute", "relative")` argument. It should not be implicit,
because demographic and contact schedules are time-indexed.

## Common Arguments And Overrides

Common simulation arguments should be a named list whose names match the
target simulator:

```r
common_args <- list(
  model = model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.06,
  susceptibility = susceptibility,
  infectiousness = infectiousness,
  demographic_process = demographic_process,
  time_policy = "linear",
  method = "euler",
  cumulative_flows = cumulative_flows
)
```

Scenario-specific changes should be applied in this order:

1. common arguments;
2. branch `initial_state`;
3. projection `times`;
4. scenario `overrides`;
5. scenario `modifier`.

Later arguments win. This keeps the rules simple and gives the modifier final
control.

## What Scenarios May Modify

Scenarios should be allowed to modify any simulator argument by replacement:

- scalar or vector parameters such as `beta`, susceptibility, infectiousness,
  and stochastic `seed`;
- `model`, including transition rates or compartment definitions;
- `contact_matrix`;
- `demographic_process`;
- `time_policy`, `migration_policy`, or `method`;
- `initial_state`, for scenarios that intentionally start from a different
  state;
- `cumulative_flows`;
- stochastic-only arguments such as `population`, `return_events`, and later
  replicate controls.

The scenario layer should validate compatibility at the same boundary as the
simulator. For example, if a scenario supplies a model whose compartments do
not match the branch state, `simulate_deterministic()` should raise the same
state-mapping error it already raises. Optional preflight validation can make
errors clearer, but should not duplicate all solver validation.

The design should not add special intervention objects in the first
milestone. Intervention start and end times can be represented by supplying
time-varying objects only after the simulator supports them, or by custom
modifier functions that build static projection-period arguments. The package
can later add a generic intervention layer once time-varying inputs are part
of the solver contract.

## Output Combination

`combine_scenario_outputs()` should accept a named list of normalized scenario
outputs and return a list:

```r
list(
  trajectory = combined_trajectory,
  cumulative = combined_cumulative_or_NULL,
  events = combined_events_or_NULL,
  metadata = metadata
)
```

Combined trajectory rows should preserve existing columns and add:

- `scenario`;
- `scenario_description`, if available;
- later `replicate`, for stochastic runs.

Combined cumulative rows should preserve existing cumulative columns and add
the same metadata columns. Combined event rows should preserve event-log
columns and add scenario and replicate metadata.

The helper should not combine trajectory and cumulative rows into one table.
They have different meanings and schemas, and existing summary helpers are
safe only on ordinary trajectory tables.

When historical output is included, it should be marked as a phase rather than
a scenario projection. A useful combined trajectory schema is:

```text
phase, scenario, time, compartment, age_group, value
```

where `phase` is `"historical"` for the common trajectory and `"projection"`
for scenario-specific output. To avoid duplicating branch-time rows, the
combined full trajectory can either:

- include historical rows through `branch_time` and projection rows strictly
  after `branch_time`; or
- include all projection rows and omit the historical branch row.

The first option is easier to plot as continuous lines and avoids duplicate
time-scenario rows. The projection-only output should still retain the
scenario branch row because it is useful for checking common initial values.

## Scenario Comparison Summaries

`summarise_scenarios()` should be metric-driven rather than disease-specific.
A metric can be a named function that receives the normalized combined output
and returns a data frame containing `scenario`, optional grouping columns, and
one or more numeric values.

Examples of generic metric helpers that agepi can provide later:

- `final_compartment_value(compartment, by = NULL)`;
- `final_cumulative_value(cumulative_name, by = NULL)`;
- `interval_cumulative_increment(cumulative_name, start, end, by = NULL)`;
- `compartment_prevalence(numerator, denominator = "population", by = NULL)`;
- `cumulative_incidence(cumulative_name, population = NULL, by = NULL)`.

`compare_scenarios()` should take a summary table and compute differences from
a baseline scenario by matching on all columns except `scenario` and numeric
measure columns. It should support at least:

- absolute difference: `scenario - baseline`;
- averted: `baseline - scenario`;
- ratio: `scenario / baseline`;
- percent reduction: `(baseline - scenario) / baseline * 100`.

Baseline comparison should be explicit about sign. For example, "cases
averted" is not a different data source; it is a baseline-minus-scenario
comparison of a cumulative flow metric.

## Cumulative-Flow Handling

Cumulative counters should support both reset and continue semantics.

`cumulative_policy = "reset"` should be the default for projected scenario
outputs. Each scenario starts cumulative counters at zero at the branch time.
This makes scenario-period outcomes such as projected cases, treatments, or
deaths easy to compare.

`cumulative_policy = "continue"` should add the branch-time historical
cumulative value to each projected cumulative output for matching
`cumulative_name`, `transition_id`, and `age_group`. This produces continuous
lifetime or whole-horizon counters. It should require that historical and
projection cumulative specifications are compatible, or else error clearly.

The solver does not currently accept non-zero cumulative initial values. The
scenario layer can implement `continue` as a post-processing offset:

```text
continued_value = projected_reset_value + historical_branch_value
```

This avoids changing `simulate_deterministic()`. It also works for stochastic
event-log-derived cumulative outputs.

For combined full trajectories, historical cumulative rows through branch time
can be included with `phase = "historical"`, while projection rows can either
be reset or continued according to policy. Comparison helpers should default
to projection-period values; whole-horizon values should be requested
explicitly.

## Deterministic Workflow

A deterministic scenario workflow should look like:

```r
historical <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 20, by = 0.25),
  model = model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = beta,
  demographic_process = demographic_process,
  time_policy = "linear",
  cumulative_flows = flows
)

scenarios <- ScenarioSet(
  baseline = "status_quo",
  branch_time = 20,
  common_args = list(
    model = model,
    age_structure = ages,
    contact_matrix = contact_matrix,
    beta = beta,
    demographic_process = demographic_process,
    time_policy = "linear",
    method = "euler",
    cumulative_flows = flows
  ),
  scenarios = list(
    Scenario("status_quo"),
    Scenario("lower_transmission", overrides = list(beta = beta * 0.8))
  )
)

output <- run_scenarios(
  scenarios,
  simulator = "deterministic",
  historical = historical,
  times = seq(20, 30, by = 0.25)
)
```

This keeps the branch state common and lets scenarios differ only where their
definitions say they differ.

## Stochastic Workflow

Stochastic support should come after deterministic orchestration. The output
shape should reserve `replicate` from the start.

A future stochastic API can add:

```r
run_scenarios(
  scenario_set,
  simulator = "stochastic",
  replicates = 100,
  seeds = NULL
)
```

Recommended stochastic representation:

- one row per scenario-replicate-time-compartment-age in `trajectory`;
- one row per scenario-replicate-time-flow-age in `cumulative`;
- one row per scenario-replicate-event in `events`, when returned;
- `replicate` as an integer column;
- `seed` as metadata or a column if per-replicate seeds are generated.

Stochastic branch-state extraction from a deterministic historical trajectory
is risky because stochastic initial states must be whole-number counts. The
first stochastic milestone should require branch states to satisfy existing
`simulate_stochastic()` validation. Rounding or resampling deterministic
states should remain user code unless agepi later defines an explicit policy.

Stochastic summaries should distinguish:

- replicate-level summaries;
- scenario means, medians, quantiles, or probabilities over replicates;
- baseline comparisons either paired by replicate ID or unpaired by scenario
  distribution. Pairing should be explicit.

## Testing Strategy

Recommended tests for the first deterministic milestone:

- `Scenario()` validates name, description, overrides, modifier, and metadata.
- `ScenarioSet()` validates unique scenario names and baseline membership.
- `initial_state_from_simulation()` extracts exact branch rows from both a
  data-frame trajectory and a list output with `trajectory`.
- Missing, duplicated, or incomplete branch states error clearly.
- `run_scenarios()` runs two deterministic scenarios from the same branch state
  and preserves common branch-time values.
- Scenario overrides change only the intended simulator argument.
- Custom `modifier` functions can alter arguments.
- `combine_scenario_outputs()` preserves existing trajectory and cumulative
  columns plus scenario metadata.
- Cumulative reset starts projection counters at zero at branch time.
- Cumulative continue offsets projection counters by historical branch values.
- Historical plus projection combination avoids duplicate branch-time rows in
  full combined output.
- Summary helpers can compute final cumulative totals and compare them to a
  named baseline.

Later stochastic tests should cover replicate IDs, seed reproducibility,
event-log preservation, cumulative event counts by scenario-replicate, and
paired versus unpaired baseline comparison semantics.

## Implementation Roadmap

Milestone 1: deterministic orchestration document and constructors

- Add `R/scenarios.R`.
- Implement `Scenario()` and `ScenarioSet()` as lightweight validated lists.
- Add docs and examples without changing solvers.

Milestone 2: branch-state extraction

- Implement `initial_state_from_simulation()`.
- Test exact extraction and validation.

Milestone 3: deterministic runner and combiner

- Implement `run_scenarios(..., simulator = "deterministic")`.
- Implement internal argument resolution and output normalization.
- Implement `combine_scenario_outputs()`.
- Add a deterministic example using an existing generic model.

Milestone 4: cumulative policy

- Add reset and continue handling.
- Keep continue as post-processing offsets rather than solver changes.

Milestone 5: metric and comparison helpers

- Add small generic metric helpers or document metric-function patterns.
- Implement `compare_scenarios()` for summary tables.

Milestone 6: stochastic extension

- Add replicate handling around `simulate_stochastic()`.
- Preserve deterministic API shape by adding `replicate = NA_integer_` only
  when needed or by documenting that deterministic outputs omit it.

## Open Design Questions

- Should projection times be required to be absolute forever, or should
  relative projection times be supported with explicit `time_origin`?
- Should `Scenario()` offer special fields such as `parameter_overrides`,
  `model_overrides`, and `contact_matrix_overrides`, or is one `overrides`
  list plus `modifier` enough?
- Should scenario metadata become table columns automatically, and if so how
  should list-valued metadata be represented?
- Should `run_scenarios()` include historical rows by default, or should it
  return projection-only output and provide a separate helper to bind the
  historical trajectory?
- Should `compare_scenarios()` work only on summary tables, or should it also
  compare full time series directly?
- What is the right long-term representation for time-varying interventions
  once core simulators support time-varying inputs?
- Should stochastic scenario comparisons default to paired replicate
  comparisons when seeds are shared?

## Recommendation

Implement deterministic scenario orchestration first. Keep the scenario layer
as a package-managed wrapper over existing simulation calls, require exact
branch-time extraction, add `scenario` metadata to ordinary trajectory and
cumulative tables, and make cumulative reset versus continue an explicit
policy. Leave disease-specific outcome definitions to user-supplied metrics,
with small generic comparison helpers added only after the combined output
shape is stable.
