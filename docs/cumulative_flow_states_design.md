# Auxiliary Cumulative Flow States Design

## Scope

This document records a design inspection for package-managed cumulative flow
states. The narrow deterministic auxiliary-state milestone is implemented in
`simulate_deterministic(cumulative_flows = ...)`. Stochastic cumulative outputs
are implemented in `simulate_stochastic(cumulative_flows = ...)` by summarising
realised event logs; stochastic cumulative counters as state variables and
cumulative flows with demographic coupling remain out of scope.

The target use case is to let users request age-specific cumulative counters
for selected model transitions, for example infections, progressions, or
recoveries, while keeping those counters out of compartment totals,
force-of-infection denominators, demographic flows, and stochastic event
propensities.

## Current Architecture Summary

The current disease simulation stack is built around compartment-age state
cells and transition-rate rows:

```text
initial state
  -> state_vector_to_long()/state_long_to_vector()
  -> transition_rates()
  -> rates_to_derivative()
  -> deterministic_derivative()
  -> integrate_state_trajectory()
```

For stochastic simulation, the same `transition_rates()` rows are used as
Gillespie propensities:

```text
state
  -> transition_rates()
  -> stochastic_event_table()
  -> stochastic_apply_event()
```

The implementation currently assumes that the simulation state vector contains
only population compartments:

- `DiseaseModel` objects expose `model$compartments`; all state conversion,
  validation, derivative construction, stochastic indexing, and output use that
  vector.
- `state_vector_to_long()` and `state_long_to_vector()` require exactly one row
  per compartment-age pair.
- `transition_rates()` converts any numeric state through
  `state_vector_to_long(state, age_structure, model$compartments)`.
- `transition_population_by_age()` sums across `model$compartments`.
- `force_of_infection()` receives population denominators from those sums.
- `rates_to_derivative()` creates derivatives only for supplied compartments
  and subtracts/adds every transition row between those compartments.
- `simulate_deterministic()` validates the ordinary initial state length as
  `length(model$compartments) * n_age_groups`. When `cumulative_flows` is
  supplied, it augments the solver state internally and returns ordinary
  compartment trajectories separately from cumulative output.
- `simulate_stochastic()` validates integer compartment counts, checks fixed
  population totals from `model$compartments`, and turns each transition-rate
  row into a competing event.
- deterministic demographic coupling loops over `model$compartments` for
  births, deaths, ageing, migration allocation, and population exposure.
- summary helpers such as `compartment_totals()`, `age_group_totals()`, and
  `total_population()` currently sum all rows in simulation output.

Population totals are computed in three main places:

- `transition_population_by_age()`, used by `transition_rates()` and the
  stochastic population validator;
- `compartment_demographic_derivative()`, which sums all model compartments by
  age before fertility, mortality, and migration calculations;
- simulation-summary helpers, which sum the returned tidy rows without a state
  type distinction.

Fixed-population checks are performed by `validate_stochastic_population()`.
It derives initial age-specific totals from the model compartments and requires
any supplied `population` vector to match those totals. Stochastic conservation
tests also aggregate over all output compartments.

State variables are assumed to be non-negative compartment counts in
`transition_rates()`, `simulate_deterministic()`, Euler non-negativity checks,
`simulate_stochastic()`, stochastic integer checks, demographic derivative
coupling, and simulation-summary validation.

Naively adding cumulative states to `model$compartments` or to the ordinary
state vector would therefore break important assumptions. Counters would be
included in population denominators, total-population summaries, demographic
birth/death/ageing/migration flows, and stochastic propensities. In stochastic
simulation, a counter transition such as `S -> C_infections` would become a
separate competing event rather than an update attached to the real `S -> E`
event.

## State Mapping

The current state-mapping infrastructure does not distinguish compartments
from non-compartment states. It supports a compartment-major,
age-group-minor vector layout and long-form rows with columns `compartment`,
`age_group`, and `value`. State vector names are generated as
`<compartment>_<age_group>` but are intentionally ignored on input.

Age-specific cumulative states could be added to a deterministic state vector,
but they should not be represented as ordinary entries in `model$compartments`.
The safest design is a separate internal mapping for the augmented state:

- a base disease-compartment mapping using the existing `model$compartments`;
- an auxiliary cumulative-flow mapping with counter names and age groups;
- helper functions to split an augmented vector into `compartment_state` and
  `cumulative_state`, and to combine their derivatives and outputs.

The existing mapping could be extended only if it gains an explicit state type
or namespace, for example `state_type = "compartment"` versus
`state_type = "cumulative_flow"`. Without that distinction, downstream code
cannot know whether a row contributes to population, force of infection,
demography, or propensities.

Recommended output should also carry an explicit state kind. A compatible
option is to preserve existing rows and add auxiliary rows only when requested,
with columns such as:

```text
time, state_type, compartment, flow, age_group, value
```

However, adding columns to the existing trajectory output may affect callers
that expect exact column names. A lower-risk first implementation could return
a list when cumulative output is requested:

```r
list(
  trajectory = <existing time/compartment/age_group/value table>,
  cumulative = <time/flow/age_group/value table>
)
```

This preserves existing output shape when `cumulative_flows` is absent and
keeps counters out of population-summary helpers unless those helpers are
explicitly extended.

## Transition-Rate Table

`transition_rates()` now returns `from`, `to`, `age_group`, `rate`, and
`transition_id`.
Infection transitions and non-infection per-capita transitions both flow
through this table, so cumulative counters can be derived from one shared
source of transition logic.

Today, transition rows are stable enough for simple `from`/`to`/`age_group`
matching:

- SIR and SEIR have fixed unique transition pairs.
- `CompartmentModel()` validates unique `from`/`to` pairs separately for
  `infection_transitions` and `transitions`.
- `rates_to_derivative()` rejects exact duplicate `from`/`to`/`age_group`
  transition-rate rows as ambiguous.

There is still a future ambiguity risk. A model may eventually need two
distinct flows with the same `from` and `to` but different meanings, parameters,
or intervention schedules. Current validation rejects that shape, but adding
cumulative counters is a good point to introduce stable transition identifiers
or at least reserve room for them.

Recommended transition metadata:

```text
from, to, age_group, rate, transition_id
```

For current models, `transition_id` is generated deterministically from the
declared logical transition, for example `infection:S->E` or
`transition:E->IP`. Age-specific rows for the same logical transition share the
same identifier. The identifier distinguishes infection transitions from
ordinary per-capita transitions and distinguishes different `from`/`to` pairs.

The current architecture still does not support duplicate same-`from`/same-`to`
logical transitions with different meanings. `CompartmentModel()` rejects
duplicate `from`/`to` rows within `infection_transitions` and within
`transitions`, and `rates_to_derivative()` rejects duplicate
`from`/`to`/`age_group` transition-rate rows. Future duplicate-transition
support would need explicit declared identifiers rather than deriving IDs only
from transition type and endpoints.

Internal cumulative-flow validation helpers now normalize named-list and
data-frame `from`/`to` specifications to:

```text
cumulative_name, transition_id, from, to
```

Those helpers select existing logical transitions. Deterministic simulation now
uses them to augment the internal ODE state with auxiliary cumulative counters.
They do not add counters to `model$compartments`, create stochastic
propensities, or apply demographic processes to counters.

## Recommended Design

Cumulative flow states should be package-managed auxiliary states, not disease
compartments.

For deterministic simulation:

1. Validate requested cumulative flows against model transition metadata.
2. Build an augmented initial vector:
   `c(compartment_state, cumulative_state)`.
3. In the shared derivative function, split the incoming vector and pass only
   the compartment slice to `transition_rates()`.
4. Use the resulting transition-rate rows both for ordinary compartment
   derivatives and for cumulative derivatives.
5. For each requested flow and age group, set
   `dC_flow,age/dt = sum(matching transition rates for that flow and age)`.
6. Append cumulative derivatives to the compartment derivatives before
   returning to Euler or deSolve.
7. Output compartment trajectory and cumulative trajectory separately.

Implemented milestone: `simulate_deterministic()` now follows this design for
infection-only deterministic runs. With `cumulative_flows = NULL`, the return
value remains the ordinary tidy trajectory data frame. With cumulative flows
requested, the return value is:

```r
list(
  trajectory = <existing time/compartment/age_group/value table>,
  cumulative = <time/cumulative_name/transition_id/from/to/age_group/value table>
)
```

The cumulative derivative for each requested transition-age row is the same
`rate` value from `transition_rates()` that is passed to
`rates_to_derivative()` for ordinary compartment derivatives.

Runnable example: see
[`examples/deterministic_cumulative_flows.R`](../examples/deterministic_cumulative_flows.R)
for a SIR simulation that tracks cumulative infections and recoveries while
leaving the ordinary compartment trajectory as `time`, `compartment`,
`age_group`, and `value`.

For stochastic simulation, cumulative counters must not generate
propensities. Two safe approaches are available:

- derive cumulative counts from the event log by counting fired events matching
  each requested flow up to each output time;
- update an auxiliary counter vector immediately after the selected real event
  is applied.

Implemented stochastic milestone: `simulate_stochastic()` uses the event-log
approach. When cumulative flows are requested, it records realised event
metadata internally, counts matching events with `event_time <= output_time`,
and returns a cumulative table across all requested output times and age
groups. The event sampler, stochastic state vector, transition-rate rows, and
propensities are unchanged by a cumulative-output request.

Runnable example: see
[`examples/stochastic_cumulative_flows.R`](../examples/stochastic_cumulative_flows.R)
for a fixed-population generic stochastic simulation that tracks exposures and
clinical/subclinical onset counts from realised Gillespie events.

## API Recommendation

Use a simulation-time argument first:

```r
simulate_deterministic(
  ...,
  cumulative_flows = list(
    infections = list(from = "S", to = "E"),
    clinical = list(from = "E", to = "IP"),
    subclinical = list(from = "E", to = "IS"),
    recoveries = list(from = "IC", to = "R")
  )
)
```

Reasons:

- cumulative outputs are an observation request, not necessarily part of the
  biological model;
- the same model may be simulated with different reporting requirements;
- leaving `CompartmentModel()` unchanged avoids changing model validation and
  existing constructor outputs.

Later, `CompartmentModel(cumulative_flows = ...)` could be added as stored
default metadata, with a simulation-time override. `from`/`to` matching should
be accepted only when it resolves to one declared logical transition ID. A
data-frame representation is also supported internally and may be preferable
publicly for programmatic workflows:

```text
name, from, to
```

or, after transition IDs exist:

```text
name, transition_id
```

Recommended naming for cumulative output:

- use user-supplied names as stable `flow` labels;
- require unique non-empty names;
- avoid disease-specific core names such as `clinical_incidence`;
- if names are omitted in a data-frame API, generate names such as
  `S_to_E`, but reject collisions.

Allowing users to track all flows is reasonable as a later convenience:

```r
cumulative_flows = "all"
```

It should expand to all declared transition metadata rows, with generated
names that are unique and documented.

## Deterministic Strategy

The best insertion point is `deterministic_derivative()`, or a small wrapper
around it, because both Euler and deSolve use that path through
`integrate_state_trajectory()`.

The derivative path should be refactored into a shared core that accepts a
compartment-only state vector:

```text
compartment_state -> transition_rates() -> rates_to_derivative()
```

An augmented deterministic derivative can then:

```text
split augmented state
compute transition rates from compartment state only
compute compartment derivative
compute cumulative derivatives from selected transition-rate rows
append derivatives
```

This works naturally with deSolve because `integrate_state_trajectory_desolve()`
already passes a numeric vector of arbitrary length to the derivative callback
and extracts the same number of solved state columns. The output callback is
the part that must understand the augmented layout.

Initial values for cumulative counters should default to zero for every
requested flow and age group. A future option could accept non-zero initial
cumulative values for continuing simulations, but the first milestone should
keep that explicit and validated if supported at all.

Euler should not implement post-hoc flow accounting. Because Euler and deSolve
already share a derivative callback, Euler can receive the same augmented
derivative as deSolve. This means the Euler update for counters is
`C_next = C_current + dt * flow(t)`, matching the solver primitive rather than
an output-only reconstruction.

## Stochastic Strategy

Stochastic cumulative counters should be derived from realised events or
updated alongside the realised event. They must not be encoded as ordinary
transition-rate rows.

The implemented first stochastic design is:

1. Validate `cumulative_flows` against stochastic transition metadata.
2. Run the existing Gillespie simulation unchanged, with event logging enabled
   internally when cumulative output is requested.
3. At requested output times, derive cumulative counts by age group and flow
   from events with `time <= output_time`.
4. Return the usual trajectory plus a cumulative table, and the event log when
   `return_events = TRUE`.

This keeps cumulative counters out of `model$compartments`, state mapping,
fixed-population checks, `transition_rates()`, and `stochastic_apply_event()`.
If counters are instead held in state in a future design,
`stochastic_apply_event()` should return both the updated compartment state and
updated counters after a sampled event. The propensities must still be computed
from compartment state only.

The deterministic and stochastic APIs can share `cumulative_flows`, but their
interpretation differs:

- deterministic counters integrate transition rates over continuous time;
- stochastic counters count realised events.

That difference should be documented explicitly.

## Demography Strategy

Cumulative states should be excluded from births, deaths, ageing, migration,
population conservation checks, fertility exposure, migration allocation, and
force-of-infection denominators.

Deterministic infection-demography coupling can support cumulative counters if
the augmented derivative splits the state and sends only the compartment slice
to `compartment_demographic_derivative()`. The cumulative derivative should
come only from selected disease transition-rate rows unless future APIs add
explicit cumulative demographic flows.

For a first implementation, it is reasonable to reject
`cumulative_flows != NULL` when `demographic_process != NULL` if that keeps the
milestone small. A more useful first deterministic milestone can still support
demography safely by:

- computing disease transition rates from compartment state only;
- computing demographic derivatives from compartment state only;
- appending cumulative disease-flow derivatives;
- leaving demographic flows untracked.

The rejection route is safer if output shape and summary semantics are not yet
settled.

Stochastic demography is already unsupported, so no additional demographic
stochastic interaction is needed in the current fixed-population Gillespie
surface.

## Validation Strategy

Validation should happen at simulation time because `CompartmentModel()` does
not know the `AgeStructure()` and cumulative flows are primarily an output
request. If model-level default cumulative flows are later added, validate
their basic schema at construction time and validate their resolution against
expanded transition metadata at simulation time.

Required validation:

- `cumulative_flows` is `NULL`, `"all"`, a named list, or a supported data
  frame;
- flow names are unique, non-empty, non-missing, and do not collide after any
  generated naming;
- every flow has either a valid `transition_id` or a valid `from`/`to` pair;
- `from` and `to` name model compartments, not auxiliary states;
- each requested flow resolves to at least one existing transition-rate row
  across age groups;
- each requested `from`/`to` pair resolves to exactly one declared transition
  unless the API explicitly asks to sum multiple transition IDs;
- counters are age-specific over `age_structure$age_groups`;
- initial cumulative values, if supported, are finite, non-missing,
  non-negative, and length-compatible;
- stochastic cumulative outputs use integer event counts.

Invalid `from`/`to` pairs should error rather than silently produce zero
counter derivatives. Silent zero counters would hide misspelled compartment
names and model/API mismatches.

## Backward Compatibility

This feature can be added without changing existing behaviour when
`cumulative_flows` is absent:

- keep default arguments as `NULL`;
- keep current `simulate_deterministic()` and `simulate_stochastic()` return
  values unchanged when no cumulative flows are requested;
- keep SIR, SEIR, and `CompartmentModel()` constructors unchanged for the first
  milestone;
- keep current state mapping helpers compartment-only;
- keep existing summary helpers operating on ordinary trajectory output only.

The main compatibility risk is output shape when cumulative flows are
requested. Returning a list only when cumulative output is requested mirrors
the existing `simulate_stochastic(return_events = TRUE)` pattern and avoids
adding auxiliary rows to the ordinary compartment trajectory.

## Minimal Implementation Plan

Milestone 1: transition metadata and validation

- Add internal helpers to build declared transition metadata for SIR, SEIR, and
  `CompartmentModel()`.
- Include stable transition IDs in `transition_rates()` metadata while
  preserving the semantics of `from`, `to`, `age_group`, and `rate`.
- Add `validate_cumulative_flows()` with `from`/`to` matching and room for
  future transition IDs.
- Add focused tests for valid, invalid, and ambiguous requests. `"all"` remains
  deferred.

Milestone 2: deterministic auxiliary states without demography

- Add `cumulative_flows = NULL` to `simulate_deterministic()`.
- Build augmented initial vectors only when cumulative flows are requested.
- Split augmented state inside the derivative callback.
- Compute cumulative derivatives from transition-rate rows.
- Return `list(trajectory, cumulative)` when requested.
- Test deSolve and Euler paths where deSolve is available.

Milestone 3: deterministic demography policy

- Either reject cumulative flows with `demographic_process != NULL`, or support
  them by splitting compartment and cumulative state before both disease and
  demographic derivatives.
- Document that cumulative counters track disease transition flows only.
- Add tests ensuring counters are excluded from demographic population totals
  and demographic derivatives.

Milestone 4: stochastic cumulative outputs

- Add `cumulative_flows = NULL` to `simulate_stochastic()`.
- Derive cumulative tables from realised event logs at requested output times.
- Ensure internal event logging does not force event-log return unless
  `return_events = TRUE`.
- Status: implemented for fixed-population stochastic disease simulations.
- Test that counters equal event counts and that trajectory/population
  conservation is unchanged.

Milestone 5: output summaries and documentation

- Add helper functions for cumulative-output summaries if needed.
- Document the deterministic/stochastic interpretation difference.
- Decide whether to add model-level default cumulative-flow metadata.

## Open Questions and Risks

- Should transition IDs become part of the public `transition_rates()` output,
  or remain internal metadata until users need them?
- Should `from`/`to` matching be allowed long term, or only as a convenience
  while duplicate `from`/`to` transitions are disallowed?
- Should cumulative-output requests return a list or add auxiliary rows with a
  state-type column? The list design is safer for compatibility.
- Should deterministic cumulative counters be supported with demography in the
  first implementation, or initially rejected to reduce risk?
- Should users be able to provide non-zero initial counter values for resumed
  simulations?
- How should future demographic cumulative flows, such as births or deaths, be
  represented without mixing them with disease transition metadata?
- Existing summary helpers sum all rows they receive. They should not be used
  directly on combined compartment-plus-counter output unless output gains
  explicit state-type filtering.

## Recommendation

Proceed next with a small implementation milestone only after adding transition
metadata and cumulative-flow validation. The deterministic deSolve/Euler
strategy is straightforward if cumulative states are treated as an augmented
state slice and the existing compartment state mapping remains unchanged.
Stochastic support should wait until deterministic validation and output shape
are settled, then derive counters from realised event logs or update them
alongside sampled events.
