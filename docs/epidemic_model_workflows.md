# Epidemic And Compartmental Model Workflows

`agepi` treats epidemic modelling as a bookkeeping problem over compartments,
age groups, and time. This guide explains how the main pieces fit together:
how models are declared, how infection pressure is computed, how deterministic
and stochastic simulations differ, how cumulative flows are tracked, how
demography and contact matrices plug in, and which workflows are the right
entry point for a given task.

The short version is:

- `SIRModel()` and `SEIRModel()` are compact convenience constructors for the
  built-in age-structured epidemic models;
- `CompartmentModel()` is the generic interface for custom deterministic and
  supported fixed-population stochastic compartment structures;
- `simulate_deterministic()` runs deterministic epidemic simulation and can
  optionally couple in demography and cumulative flows;
- `simulate_stochastic()` runs fixed-population Gillespie simulation and can
  optionally return the event log and cumulative-flow counts;
- `compartment_totals()`, `age_group_totals()`, and `total_population()` are
  summary helpers for epidemic trajectories that have `compartment` and
  `value` columns;
- `CalibrationTarget()` and `evaluate_calibration_*()` are first-pass
  calibration-scaffold helpers that compare observed data with already-run
  simulation output.

## 1. Overview

An epidemic compartmental model represents people as counts in named states
such as susceptible, exposed, infectious, recovered, latent, or treated. In
`agepi`, every count is also indexed by age group and time. The core state is
therefore a compartment-by-age table evolving over time.

There are two broad modelling paths:

1. built-in convenience models such as `SIRModel()` and `SEIRModel()`;
2. custom generic models built with `CompartmentModel()`.

Use the built-in constructors when your structure really is SIR or SEIR and
you want the clearest possible workflow. Use `CompartmentModel()` when you need
additional compartments, more than one infectious state, custom per-capita
transitions, or a more specialised structure such as MSIR or a TB-style
compartment set.

The package currently focuses on deterministic prototype modelling, fixed-
population stochastic simulation, and first-pass epidemic-demography coupling.
It is not a full transmission-fitting framework and it does not introduce new
epidemic model families beyond the exported constructors.

## 2. Core Mental Model

The mechanical view is simple:

- compartments store numbers of people;
- transitions move people between compartments;
- transition rates are per-capita hazards or flows;
- force of infection turns infectious people in source age groups into new
  infection pressure for recipient age groups;
- age groups allow susceptibility, infectiousness, contact rates, and other
  parameters to vary by age;
- simulation output is tidy and can be summarised after the run.

This is bookkeeping for population flows, not a separate epidemiological
engine. A model only does what its compartments and transition rules allow.
The main task is to define those rules clearly and keep the age-grid
conventions consistent.

## 3. Age-Structured State Layout

`AgeStructure()` defines the age grid used throughout the simulation. The
package assumes the model state is ordered by compartment and age group. In
long form, the state has one row for each compartment-age cell. In vector form,
the order is compartment-major and age-group-minor. For `c("S", "I", "R")`
over `c("0-4", "5-9")`, the vector order is:

```text
S_0-4, S_5-9, I_0-4, I_5-9, R_0-4, R_5-9
```

State mapping helpers make that explicit:

- `state_long_to_vector()` converts long-form rows to a numeric state vector;
- `state_vector_to_long()` converts the other way;
- `initialise_compartments_from_proportions()` allocates an age-specific
  population into compartments from proportions;
- `validate_state_long()` and `validate_state_vector()` enforce the expected
  layout and missing-row checks.

The `initial_state` supplied to simulation functions must match the model age
structure and compartment order. Names on numeric state vectors are ignored on
input, so the position of each entry matters more than the label.

Why this matters: the state layout is used by the deterministic solver, the
stochastic event sampler, the force-of-infection calculation, the cumulative
flow accounting, and the simulation summary helpers. If the layout is wrong,
the model may still run but it will be simulating the wrong state.

## 4. Built-In SIR And SEIR Models

`SIRModel()` and `SEIRModel()` are the two convenience constructors currently
exported for the common built-in workflows.

- `SIRModel(gamma = ...)` defines `S -> I -> R`.
- `SEIRModel(sigma = ..., gamma = ...)` defines `S -> E -> I -> R`.

These constructors are useful when the model structure is standard and the
transition set is small. They are less verbose than a generic compartment
specification and they make the most common examples easier to read.

Use them when:

- your model really is SIR or SEIR;
- all infectious pressure is driven by the single infectious compartment `I`;
- you do not need custom static transitions or extra compartments.

Move to `CompartmentModel()` when:

- you need extra compartments such as maternal immunity, latent states, or
  treatment states;
- you need more than one infectious compartment;
- you need custom age-specific per-capita transitions;
- you want to represent a model closer to the paper or protocol you are
  reproducing.

The smallest deterministic SIR example is in
[examples/mock_sir_deterministic.R](../examples/mock_sir_deterministic.R).

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9", "10-14"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, 14)
)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 3),
  value = c(c(995, 1197, 898), c(5, 3, 2), c(0, 0, 0)),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(
  c(4, 2, 1,
    2, 5, 2,
    1, 2, 4),
  nrow = ages$n_age_groups,
  byrow = TRUE
)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 1, by = 0.1),
  model = SIRModel(gamma = 0.25),
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  method = "euler"
)
```

## 5. Generic CompartmentModel Workflow

`CompartmentModel()` is the general-purpose constructor for deterministic and
supported stochastic compartment structures.

The model declaration has three main parts:

- `compartments`: the compartment names and their order;
- `infection_transitions`: one or more `from -> to` infection transitions;
- `transitions`: per-capita non-infection transitions with rates.

The model can also declare:

- `infectious_compartments`: which compartments contribute to infectious
  pressure;
- `infectiousness_weights`: relative infectiousness by compartment;
- `birth_compartment`: where demographic births enter when demography is
  coupled in;
- `migration_compartment`: which compartment receives migration when the
  migration policy allocates all net movement to a single compartment.

The key distinction is between infection transitions and ordinary
per-capita transitions:

- infection transitions use the force of infection;
- ordinary transitions use the named rate supplied in the `rate` column or
  list-column entry.

Rates can be scalar or age-specific. When age-specific, they must be named by
age group. The rates are validated and aligned to the simulation age structure
when `transition_rates()` is evaluated.

`transition_rates()` returns the per-age logical transition table used by both
deterministic and stochastic simulation. Its rows carry `from`, `to`,
`age_group`, `rate`, and `transition_id`, and that transition identifier is
what cumulative-flow accounting uses to match declared flows back to model
transitions.

This is a rate-based model, not a branching-probability interface. If a source
compartment has multiple outgoing paths, each path needs its own transition
definition and rate. That is how `generic_seir.R`, `generic_msir.R`, and the
TB-style examples express competing routes.

Example `CompartmentModel()` declarations are in:

- [examples/generic_sir.R](../examples/generic_sir.R)
- [examples/generic_seir.R](../examples/generic_seir.R)
- [examples/generic_msir.R](../examples/generic_msir.R)
- [examples/tb_age_structured_demography.R](../examples/tb_age_structured_demography.R)

The generic SEIR example is the smallest comparison point against the built-in
constructor:

```r
model <- CompartmentModel(
  compartments = c("S", "E", "I", "R"),
  infection_transitions = data.frame(from = "S", to = "E"),
  transitions = data.frame(
    from = c("E", "I"),
    to = c("I", "R"),
    rate = c(0.4, 0.25)
  ),
  infectious_compartments = "I"
)
```

## 6. Force Of Infection

`force_of_infection()` is the age-structured transmission kernel. It computes
the infection pressure for each recipient age group from the current number
infectious in each source age group, the source-age population denominators,
and the contact matrix.

The convention is:

- rows are recipient age groups;
- columns are source or infectious age groups;
- `contact_matrix[a, b]` gives contacts made by recipient group `a` with
  source group `b`.

The conceptual form is:

```text
lambda = beta * susceptibility * (contact_matrix %*% (infectiousness * infectious / population))
```

In words:

- `beta` scales the overall transmission intensity;
- `infectious / population` is the infectious fraction in each source age
  group;
- `infectiousness` scales the source side;
- `susceptibility` scales the recipient side;
- the contact matrix translates source infectiousness into recipient force of
  infection.

For built-in SIR and SEIR models, infection transitions are generated from this
force-of-infection calculation. For generic compartment models, every declared
infection transition uses the same age-specific force of infection unless the
model itself is structured differently by compartment and weight.

See [docs/contact_matrix_workflows.md](contact_matrix_workflows.md) for how to
load and adapt contact matrices before they are used here.

## 7. Deterministic Simulation

`simulate_deterministic()` is the main deterministic simulator. It takes:

- an `initial_state`;
- `times`;
- a `model`;
- an `AgeStructure()`;
- a simulation-ready `contact_matrix`;
- `beta`, and optionally `susceptibility` and `infectiousness`.

When `demographic_process` is supplied, it also accepts demographic coupling
controls such as `time_policy`, `migration_policy`, `ageing_policy`, and
`output_age_structure`.

The simulator returns a tidy long data frame with columns `time`,
`compartment`, `age_group`, and `value`. If `cumulative_flows` is requested,
the return value is a list with:

- `trajectory`: the ordinary compartment trajectory;
- `cumulative`: cumulative-flow counts by time, age group, and named flow.

Deterministic simulation is appropriate when:

- you want a smooth expected trajectory;
- you want to compare model variants or scenarios;
- you want a run that is easier to inspect and summarise than a stochastic
  realisation.

It assumes the state evolves according to the chosen solver and the supplied
transition and demographic rates. It does not add random event noise.

The default solver path is explicit Euler unless `deSolve` is requested or
selected automatically for infection-only runs. The annual-cohort operator
splitting path is opt-in and deterministic-only.

## 8. Stochastic Simulation

`simulate_stochastic()` is the fixed-population Gillespie simulator. It is
event-based rather than derivative-based.

Current scope:

- SIR, SEIR, and supported generic `CompartmentModel()` structures;
- fixed population;
- no demography;
- no ageing;
- no fertility, mortality, or migration;
- no tau-leaping;
- optional event log output;
- optional cumulative flows derived from realised events.

The event sampler reuses `transition_rates()`. That means the same contact
matrix convention, force-of-infection semantics, and generic transition
structures apply. Each transition-rate row becomes a possible event within an
age group.

`seed` makes runs reproducible. `return_events = TRUE` returns the realised
event log alongside the trajectory. If no event occurs before an output time,
the current state is carried forward to that time.

With `return_events = FALSE`, the return value is the trajectory data frame.
With `return_events = TRUE`, the return value is a list containing
`trajectory` and `events`, and `cumulative_flows` adds a `cumulative` table as
well.

Use stochastic simulation when:

- the population is small enough that event noise matters;
- you need sample paths instead of means;
- you want to inspect event-level histories.

See:

- [examples/stochastic_sir.R](../examples/stochastic_sir.R)
- [examples/stochastic_seir.R](../examples/stochastic_seir.R)
- [examples/stochastic_cumulative_flows.R](../examples/stochastic_cumulative_flows.R)

## 9. Cumulative Flows And Event Accounting

Many workflows need incidence-like or event-count outputs that are not
compartments. `agepi` handles those as cumulative flows.

Cumulative flows are derived from transitions, not state compartments. They are
useful for quantities such as:

- infections;
- progression to active disease;
- treatment starts;
- recoveries;
- relapses.

Deterministic cumulative outputs are integrated alongside the ordinary
compartment trajectory. Stochastic cumulative outputs are counted from the
realised event log. In both cases, the cumulative flow names map to declared
transition endpoints through `transition_id`.

Relevant helpers:

- `cumulative_flow_totals()`;
- `cumulative_flow_totals_wide()`;
- `cumulative_flow_increments()`.

The cumulative helpers operate on the cumulative output, not on a raw epidemic
trajectory. They are distinct from `compartment_totals()`,
`age_group_totals()`, and `total_population()`, which expect epidemic
trajectories with `compartment` and `value` columns.

See [examples/deterministic_cumulative_flows.R](../examples/deterministic_cumulative_flows.R)
and [examples/stochastic_cumulative_flows.R](../examples/stochastic_cumulative_flows.R).

## 10. Summarising Simulation Output

`agepi` exposes a small set of summary helpers for epidemic trajectories:

- `compartment_totals()`: totals by time and compartment;
- `age_group_totals()`: totals by time and age group;
- `total_population()`: totals by time.

These helpers expect epidemic-trajectory output with `time`, `compartment`,
`age_group`, and `value` columns. That includes deterministic trajectories and
stochastic trajectories, but not standalone demography output from
`simulate_demography()`, which uses `time`, `age_group`, and `population`
instead.

For cumulative outputs, use the cumulative-flow helpers instead of the
compartment helpers.

Observed case data can be mapped into age groups with `as_agepi_cases()`. That
is useful when building calibration targets or when summarising observed
records by model age bands.

## 11. Adding Demography To Epidemic Models

`simulate_deterministic()` can couple a `DemographicProcess()` into an epidemic
simulation. The coupling is deterministic and schedule-based.

Mechanically, it does four things:

- births enter the configured birth compartment, often `S`;
- background mortality removes people from the active compartments;
- ageing moves each compartment independently between age groups;
- migration is allocated according to the selected migration policy.

The current coupling is a first-pass convention, not a full demographic
projection engine. It uses the current compartment-age state to determine how
many people are exposed to fertility and migration calculations, and then
updates the epidemic state as the solver advances.

This means:

- demography and epidemic dynamics interact through the changing age-specific
  population;
- migration allocation is a policy choice, not a biological fact;
- disease-induced mortality is not automatically added;
- a WPP-backed demographic process is still not the same as exact WPP
  reproduction.

See [docs/demography_workflows.md](demography_workflows.md) for the
mechanics of the demography layer and the age-grid assumptions it uses.

## 12. Adding Contact Matrices To Epidemic Models

The epidemic simulator expects a simulation-ready numeric contact matrix, not a
source object. The source/adapt workflow is therefore:

```r
source_contacts <- load_contact_matrix_source(...)
contact_matrix <- adapt_contact_matrix_to_age_structure(
  source_contacts,
  age_structure,
  population = ...
)
```

The adapted matrix must match the model `AgeStructure()`. `beta`,
susceptibility, and infectiousness complete the transmission setup.

The contact-matrix guide in
[docs/contact_matrix_workflows.md](contact_matrix_workflows.md) explains the
source loaders, orientation convention, age-grid adaptation, and limitations in
more detail.

## 13. Calibration-Facing Outputs

`agepi` includes a thin calibration scaffold, not a full optimisation system.
The main APIs are:

- `CalibrationTarget()`;
- `evaluate_calibration_target()`;
- `evaluate_calibration_objective()`.

These helpers compare already-run simulation output against observed data. They
do not run simulations themselves, and they do not estimate parameters.

Supported output types are:

- `trajectory`: epidemic trajectories with `compartment` and `value` columns;
- `cumulative`: cumulative-flow outputs with `cumulative_name` and `value`
  columns.

That makes the calibration layer useful for simple likelihood checks or for
building small objectives around observed case counts, cumulative event counts,
or compartment totals. It is still a scaffold: the package does not yet provide
a full calibration workflow with parameter search, priors, or MCMC.

The case-data adapter `as_agepi_cases()` is useful when observed records need
to be mapped into the same age bins as the model before calibration.

## 14. Worked Example A: Simple Age-Structured SIR

This example is dependency-free and shows the built-in SIR workflow end to end.

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9", "10-14"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, 14)
)

initial_state <- data.frame(
  compartment = rep(c("S", "I", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 3),
  value = c(
    c(995, 1197, 898),
    c(5, 3, 2),
    c(0, 0, 0)
  ),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(
  c(4, 2, 1,
    2, 5, 2,
    1, 2, 4),
  nrow = ages$n_age_groups,
  byrow = TRUE
)

sir <- SIRModel(gamma = 0.25)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 2, by = 0.1),
  model = sir,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  method = "euler"
)

compartment_totals(simulation)
total_population(simulation)
```

This is the easiest place to check that the age-grid, contact matrix, and
state layout all line up.

## 15. Worked Example B: Generic Compartment Model

This example shows when to use `CompartmentModel()` instead of the built-in
constructors.

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

model <- CompartmentModel(
  compartments = c("S", "E", "I", "R"),
  infection_transitions = data.frame(from = "S", to = "E"),
  transitions = data.frame(
    from = c("E", "I"),
    to = c("I", "R"),
    rate = c(0.4, 0.25)
  ),
  infectious_compartments = "I"
)

initial_state <- data.frame(
  compartment = rep(c("S", "E", "I", "R"), each = ages$n_age_groups),
  age_group = rep(ages$age_groups, times = 4),
  value = c(
    c(990, 1188, 896),
    c(4, 3, 2),
    c(6, 4, 2),
    c(0, 0, 0)
  ),
  stringsAsFactors = FALSE
)

contact_matrix <- matrix(
  c(4, 2, 1,
    2, 5, 2,
    1, 2, 4),
  nrow = ages$n_age_groups,
  byrow = TRUE
)

simulation <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 2, by = 0.1),
  model = model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  method = "euler"
)

tail(compartment_totals(simulation))
```

The generic model path is the right one once the model is no longer a plain
SIR or SEIR structure.

## 16. Worked Example C: Cumulative Flows

The cumulative-flow workflow is useful when you care about incidence-like
counts rather than only the compartment trajectory.

```r
output <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 10, by = 1),
  model = SIRModel(gamma = 0.25),
  age_structure = ages,
  contact_matrix = contact_matrix,
  beta = 0.08,
  cumulative_flows = list(
    infections = list(from = "S", to = "I"),
    recoveries = list(from = "I", to = "R")
  )
)

head(output$trajectory)
head(output$cumulative)
cumulative_flow_totals(output$cumulative)
cumulative_flow_totals_wide(output$cumulative)
```

The same request shape works for stochastic simulation, but the stochastic
cumulative output is counted from realised events rather than integrated from
transition rates.

## 17. Workflow Comparison Table

| workflow | main functions | use case | strengths | limitations | example file |
|---|---|---|---|---|---|
| Built-in SIR/SEIR | `SIRModel()`, `SEIRModel()`, `simulate_deterministic()` | Small age-structured epidemic prototypes | Very compact, easy to inspect | Only the built-in structures | `examples/mock_sir_deterministic.R` |
| Generic compartment model | `CompartmentModel()`, `transition_rates()`, `simulate_deterministic()` | Custom deterministic compartment structures | Flexible compartments and per-capita transitions | Still static-rate and compartment-based | `examples/generic_sir.R`, `examples/generic_seir.R`, `examples/generic_msir.R` |
| Deterministic simulation | `simulate_deterministic()` | Smooth epidemic trajectories | Supports demography, cumulative flows, and age structure | No stochastic event noise | `examples/mock_seir_demography.R` |
| Stochastic simulation | `simulate_stochastic()` | Event-level sample paths | Reproducible Gillespie simulation | Fixed population only; no demography | `examples/stochastic_sir.R`, `examples/stochastic_seir.R` |
| Cumulative-flow simulation | `simulate_deterministic(..., cumulative_flows = ...)`, `simulate_stochastic(..., cumulative_flows = ...)` | Incidence-like counts and derived outputs | Keeps event accounting separate from compartments | Flows must match declared transitions | `examples/deterministic_cumulative_flows.R`, `examples/stochastic_cumulative_flows.R` |
| Epidemic-demography coupling | `simulate_deterministic(..., demographic_process = ...)` | Age-structured epidemic turnover | Births, mortality, ageing, and migration in one run | Deterministic coupling; policy assumptions matter | `examples/annual_cohort_sir_demography.R`, `examples/kiribati_tb_realistic_demography.R` |
| Calibration-facing workflow | `as_agepi_cases()`, `CalibrationTarget()`, `evaluate_calibration_target()` | Small likelihood checks against observed data | Clear separation of observed data and model output | Not a full calibration engine | `examples/observed_case_age_groups.R` |

## 18. Limitations And Boundaries

The current epidemic model surface is intentionally limited:

- deterministic and stochastic paths are separate;
- stochastic simulation is fixed-population only;
- epidemic and demographic coupling is deterministic, schedule-based, and
  still a first-pass convention;
- contact matrices are interpreted as age-structured contact weights, not as a
  calibrated transmission model by themselves;
- contact matrices, `beta`, susceptibility, and infectiousness together define
  transmission pressure;
- disease-specific mortality conventions may still need explicit model choices;
- solver choice and time stepping matter for numerical accuracy;
- calibration support is a scaffold, not a full workflow;
- the examples are teaching or prototype scaffolds unless a file explicitly
  says it is a calibrated model.

## 19. Existing Examples

These scripts are the quickest way to see the package in action:

- `examples/mock_sir_deterministic.R`: smallest deterministic age-structured
  SIR example.
- `examples/generic_sir.R`: `CompartmentModel()` version of SIR, compared with
  the built-in constructor.
- `examples/generic_seir.R`: generic SEIR example with custom transitions.
- `examples/generic_msir.R`: MSIR-style custom model with maternal immunity.
- `examples/mock_seir_demography.R`: deterministic SEIR with demographic
  turnover.
- `examples/annual_cohort_sir_demography.R`: annual-cohort epidemic-demography
  coupling.
- `examples/deterministic_cumulative_flows.R`: deterministic cumulative
  infection and recovery tracking.
- `examples/stochastic_sir.R`: fixed-population stochastic SIR.
- `examples/stochastic_seir.R`: fixed-population stochastic SEIR.
- `examples/stochastic_cumulative_flows.R`: stochastic cumulative event counts.
- `examples/tb_age_structured_demography.R`: toy TB-style compartment model
  with demography and cumulative flows.
- `examples/kiribati_tb_realistic_demography.R`: public-data Kiribati TB
  scaffold with WPP demography, contact matrices, and provisional TB
  assumptions.

## 20. README Integration

The README now links to this guide from the documentation section. That keeps
the high-level package overview short while giving users a deeper workflow
explanation when they need it.

See also:

- [docs/demography_workflows.md](demography_workflows.md)
- [docs/contact_matrix_workflows.md](contact_matrix_workflows.md)

## 21. Consistency Notes

This guide follows the current exported API and output shapes:

- deterministic trajectories use `time`, `compartment`, `age_group`, and
  `value`;
- stochastic trajectories use the same shape, plus optional event logs and
  cumulative outputs;
- standalone demography uses `time`, `age_group`, and `population`;
- cumulative outputs use `time`, `cumulative_name`, `transition_id`, `from`,
  `to`, `age_group`, and `value`.

If you are reading this alongside the demography and contact-matrix guides,
the division of responsibilities is deliberate:

- demography workflows explain how population structure changes;
- contact-matrix workflows explain how age-specific exposure is loaded and
  adapted;
- this guide explains how epidemic compartments, transitions, and outputs fit
  together.
