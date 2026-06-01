# Stochastic Simulation Design

## Purpose

The stochastic simulator should add individual event stochasticity to disease
transitions while reusing agepi's existing age-structured model conventions:
compartment-major state ordering, age groups from `AgeStructure()`, recipient-row
and source-column contact matrices, `force_of_infection()` semantics, and tidy
long-format simulation output.

This document records the current stochastic design. `simulate_stochastic()` is
implemented for fixed-population disease events using the existing
transition-rate interface where the model supplies enough transition metadata.
It should not be read as committing to a fully stochastic demographic simulator.

## Current implementation: fixed-population supported compartment models

The implemented stochastic surface is deliberately narrow:

- SIR, SEIR, and supported generic `CompartmentModel()` transition structures;
- fixed population;
- no demography;
- no ageing;
- no migration;
- no fertility;
- no background mortality;
- Gillespie/direct-method continuous-time simulation;
- output aligned to requested observation times;
- reproducible trajectories when a random seed is supplied;
- optional event log output.

Each Gillespie propensity is built from `transition_rates()`. A transition-rate
row with `from`, `to`, `age_group`, and `rate` becomes one possible
individual-level event in that age group. When the event fires, one individual
is removed from `from` and added to `to`. This keeps event semantics aligned
with deterministic derivatives from `rates_to_derivative()` while preserving
integer counts.

The SIR event semantics are:

```text
infection in age group a: S_a -> I_a
recovery in age group a: I_a -> R_a
```

Infection propensities use the same force-of-infection convention as the
deterministic simulator:

```text
lambda_a(t) =
  beta(t) * susceptibility_a *
  sum_b C_ab(t) * infectiousness_b * I_b(t) / N_b(t)
```

Rows of `C` are recipient age groups `a`; columns are infectious source age
groups `b`. The infection event propensity for recipient age group `a` is then:

```text
lambda_a(t) * S_a(t)
```

The recovery event propensity for age group `a` is:

```text
gamma * I_a(t)
```

The implementation preserves integer counts and conserves total population
exactly for fixed-population SIR and SEIR.

The SEIR event semantics are:

```text
infection in age group a: S_a -> E_a
progression in age group a: E_a -> I_a
recovery in age group a: I_a -> R_a
```

The infection propensity uses the existing force-of-infection convention, with
infectious pressure coming from `I` rather than `E`. Progression and recovery
propensities are:

```text
sigma * E_a(t)
gamma * I_a(t)
```

For generic `CompartmentModel()` objects, infection transitions declared in
`infection_transitions` use the same force-of-infection path as deterministic
simulation, including named relative infectiousness weights for multiple
infectious compartments. Other declared `transitions` are interpreted as
per-capita within-age transitions with scalar or age-specific rates. Generic
stochastic support is therefore limited to fixed-population transitions that
can be represented as one source compartment and one destination compartment
in the same age group.

## Demography policy

Ageing must remain deterministic. Stochastic ageing should not be a goal.

Mortality may later support a stochastic implementation. Fertility may later
support a stochastic implementation. Migration may later support a stochastic
implementation, but only with a carefully documented policy.

The package should not assume that demography is entirely deterministic or
entirely stochastic. A future API might make this explicit with a mode such as:

```text
demography_mode = "none"
demography_mode = "deterministic"
demography_mode = "stochastic"
```

This is only a possible future API idea. It is not part of the current
fixed-population stochastic implementation.

## Stochastic mortality

A later stochastic mortality extension could represent death as a removal event
from each disease compartment:

```text
mortality event for compartment X and age group a: X_a -> removed
propensity: mu_a(t) * X_a(t)
```

Mortality would need to apply consistently across disease compartments. For
SIR, that means `S`, `I`, and `R`; for SEIR, `S`, `E`, `I`, and `R`. Any future
implementation should document whether mortality rates are shared across
compartments, disease-state-specific, or supplied by a separate policy.

## Stochastic fertility

A later stochastic fertility extension could represent births as events into
the youngest susceptible compartment:

```text
birth event: removed/external source -> S_youngest
```

The birth propensity could be based on the fertility schedule and the relevant
reproductive-age population. This requires careful interpretation because
agepi's current fertility convention is births per female person-year, while a
simulation state may contain total population, female population, or an
externally scaled fertility exposure. The stochastic fertility policy should
state which exposure is being used before implementation.

## Stochastic migration

Migration is more subtle than infection, recovery, mortality, or fertility.
Gross immigration and gross emigration can be represented as stochastic event
rates:

```text
immigration into compartment X and age group a: external source -> X_a
emigration from compartment X and age group a: X_a -> removed
```

Net migration is not naturally a Poisson event rate. Positive net migration
could be approximated as immigration, and negative net migration could be
approximated as emigration, but that is a modelling assumption. It loses
turnover information when simultaneous gross inflow and gross outflow produce a
small net value.

Stochastic migration should therefore not be implemented without an explicit
policy for gross versus net flows, compartment allocation, and how negative
inputs are interpreted.

## Hybrid deterministic/stochastic simulation

A later hybrid simulator could combine:

- stochastic disease events;
- deterministic ageing;
- optional deterministic or stochastic mortality, fertility, and migration;
- force of infection recomputed from the current state after each update.

If deterministic demographic updates are applied to a stochastic state, they
may create non-integer compartment values. That may be acceptable for a hybrid
approximation, but it should be documented clearly because the current
fixed-population implementation keeps integer counts.

There is an open design choice between:

- a fixed-step hybrid approximation, where stochastic disease events and
  deterministic demographic updates are combined over small intervals; and
- a continuous-time piecewise deterministic simulation, where deterministic
  flows evolve between stochastic events.

This choice should not be resolved prematurely. The current fixed-population
implementation does not require it.

## Output and validation expectations

The current stochastic implementation is validated against small, inspectable
cases and deterministic expectations. Targets include:

- a fixed seed gives reproducible trajectories;
- all states remain non-negative;
- fixed-population supported models conserve total population exactly;
- zero initial infections produce no epidemic events;
- averages over many stochastic simulations approximately track the
  deterministic trajectory;
- event times are nondecreasing;
- outputs are aligned to requested observation times.

The output should follow the existing agepi long-format convention with columns
compatible with deterministic simulation output wherever practical:
`time`, `compartment`, `age_group`, and `value`.
