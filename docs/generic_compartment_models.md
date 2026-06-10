# Generic Compartment Models

`CompartmentModel()` defines deterministic age-structured compartment models
without adding a new specialised constructor for every compartment layout. It
uses the same state layout, force-of-infection calculation, transition-rate
tables, derivative conversion, and deterministic simulation machinery as the
existing `SIRModel()` and `SEIRModel()` workflows.

`SIRModel()` and `SEIRModel()` remain supported convenience constructors for
the built-in SIR and SEIR models. Use `CompartmentModel()` when the compartment
structure is custom, such as MSIR or a model with additional static
per-capita transitions.

## Compartments and State Order

Compartments are supplied as a unique character vector:

```r
compartments = c("M", "S", "I", "R")
```

The deterministic state uses compartment-major, age-group-minor ordering. In
long form, this means all rows for the first compartment across ages, then all
rows for the second compartment, and so on. Compartment names are matched
literally, so examples and tests usually use short names such as `S`, `E`, `I`,
`R`, or `M`.

## Infection Transitions

Infection transitions are supplied separately from fixed per-capita
transitions:

```r
infection_transitions = data.frame(from = "S", to = "I")
```

For each recipient age group, infection flow is:

```r
lambda[age] * state[from, age]
```

where `lambda` is computed by `force_of_infection()`. Multiple infection
transitions can be supplied when a custom model needs them, but each transition
uses the same age-specific force of infection.

## Fixed Per-Capita Transitions

Non-infection transitions are supplied with `from`, `to`, and `rate` columns:

```r
transitions = data.frame(
  from = c("M", "I"),
  to = c("S", "R"),
  rate = c(0.15, 0.25)
)
```

Rates are non-negative per-capita rates. A rate may be a scalar or, when
transition rates are evaluated, a named age-specific numeric vector with one
value per age group. Age-specific vectors are supplied through a list column
and are reordered to match `age_structure$age_groups`:

```r
transitions <- data.frame(
  from = c("E", "E", "IP", "IC", "IS"),
  to = c("IP", "IS", "IC", "R", "R")
)
transitions$rate <- I(list(
  c("0-4" = 0.20, "5-9" = 0.35),
  c("0-4" = 0.30, "5-9" = 0.15),
  0.25,
  0.10,
  c("5-9" = 0.40, "0-4" = 0.30)
))
```

Names must match the age groups used for simulation; missing, unknown,
unnamed, negative, or non-finite age-specific rates are rejected when
`transition_rates()` expands the model. This deferred validation is necessary
because `CompartmentModel()` does not itself take an `AgeStructure()`.

`transition_rates()` also includes a stable `transition_id` metadata column for
each logical transition. Age-specific rows for the same logical transition
share the same ID, such as `infection:S->E` or `transition:E->IP`.

Multiple outgoing transitions from the same source compartment are supported
when they have different destinations. For example, `E -> IP` and `E -> IS`
can represent competing per-capita routes using age-specific rates. This is a
rate-based representation, not an explicit branching-probability interface.
Duplicate same-source/same-destination transitions with different meanings are
not currently supported.
For example, `M -> S` can represent fixed loss of maternal immunity and
`I -> R` can represent recovery.

Vaccination schedules and time-varying interventions are not yet implemented.
Static transitions such as `S -> V` or `M -> S` can be represented if supplied
as fixed per-capita rates.

## Cumulative Flow Tracking

For infection-only deterministic simulations, `simulate_deterministic()` can
track selected transition flows as auxiliary cumulative states without adding
them to `model$compartments`:

```r
output <- simulate_deterministic(
  initial_state = initial_state,
  times = seq(0, 10, by = 0.25),
  model = model,
  age_structure = ages,
  contact_matrix = contact_matrix,
  cumulative_flows = list(
    exposures = list(from = "S", to = "E"),
    symptomatic = list(from = c("E", "IP"), to = c("IP", "IC")),
    subclinical = list(from = "E", to = "IS")
  )
)

head(output$trajectory)
head(output$cumulative)
```

When cumulative flows are requested, the return value is a list with
`trajectory` containing the ordinary compartment trajectory and `cumulative`
containing `time`, `cumulative_name`, `transition_id`, `from`, `to`,
`age_group`, and `value`. A named-list entry can provide equal-length `from`
and `to` vectors to sum several disease transitions into one cumulative counter.
The cumulative derivatives use the same transition-rate rows as the ordinary
deterministic derivative. When
`demographic_process` is supplied, the cumulative counters still track disease
transition flows only and remain separate from demographic ageing, fertility,
mortality, and migration.

For fixed-population stochastic simulations, `simulate_stochastic()` can use
the same `cumulative_flows` specification to return cumulative counts derived
from realised event logs. These stochastic cumulative outputs count events with
`event_time <= output_time`; they are not stochastic state variables and do not
change propensities or event sampling.

## Infectious Compartments

`infectious_compartments` names the compartments that contribute to infectious
pressure:

```r
infectious_compartments = "I"
```

If the compartment vector contains `I`, it is used by default. Models with
other infectious states can name those states explicitly. When there is more
than one infectious compartment, named `infectiousness_weights` gives fixed
relative infectiousness weights for combining them before the force of
infection is calculated:

```r
infectious_compartments = c("IP", "IC", "IS")
infectiousness_weights = c(IP = 1, IC = 1, IS = 0.5)
```

These compartment weights are separate from the age-specific `infectiousness`
argument passed to `transition_rates()` or `simulate_deterministic()`.

## Contact Matrix Convention

Contact matrices use recipient rows and source columns. Entry `[i, j]`
contributes infectious contacts from source age group `j` to recipient age
group `i`. Age-specific `susceptibility` is therefore indexed by recipient age,
and age-specific `infectiousness` is indexed by source age.

## Demographic Coupling

When `simulate_deterministic()` is called with a `demographic_process`,
generic models use model-level target compartments for births and susceptible
policy migration:

```r
CompartmentModel(
  compartments = c("M", "S", "I", "R"),
  infection_transitions = data.frame(from = "S", to = "I"),
  transitions = data.frame(
    from = c("M", "I"),
    to = c("S", "R"),
    rate = c(0.15, 0.25)
  ),
  infectious_compartments = "I",
  birth_compartment = "M",
  migration_compartment = "S"
)
```

If `birth_compartment` is omitted and a compartment named `S` exists, births
enter `S`. If `migration_compartment` is omitted, it defaults to the birth
compartment. Under `migration_policy = "proportional"`, migration is allocated
across all compartments by current age-specific shares instead of using
`migration_compartment`.

## Current Limitations

- The generic interface currently supports deterministic, static-rate
  compartment models.
- Contact matrices, `beta`, susceptibility, infectiousness, compartment
  infectiousness weights, and fixed per-capita transition rates are static
  during a simulation.
- Vaccination schedules and time-varying intervention schedules are not yet
  implemented.
- Infection transitions all use the same force of infection.
- Cumulative flow tracking is supported for deterministic simulations with or
  without `demographic_process`, and as event-log-derived counts for
  fixed-population stochastic simulations.
- The built-in convenience constructors still cover the standard specialised
  SIR and SEIR paths; generic SIR and SEIR are mainly useful for validation and
  as templates for custom structures.

See `examples/generic_sir.R`, `examples/generic_seir.R`,
`examples/generic_msir.R`, `examples/deterministic_cumulative_flows.R`, and
`examples/stochastic_cumulative_flows.R` for complete runnable examples.
`examples/tb_age_structured_demography.R` gives a toy TB-style chronic
infection example using only the public API: `CompartmentModel()`, static
age-assortative mixing, closed demographic turnover, and cumulative flows,
including aggregation of `Lr -> I` and `Ld -> I` into one disease-onset output.
The TB example uses illustrative parameters only; it is not calibrated to
Kiribati or any other setting and should not be interpreted as policy evidence.
