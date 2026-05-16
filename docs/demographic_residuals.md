# Demographic Residual Diagnostics

This note documents the current convention for demographic residual diagnostics
and residual-derived migration schedules in agepi. These helpers are intended
for diagnostic workflows around demographic-only simulations, especially when
comparing closed or partially specified demographic dynamics to externally
supplied WPP-style observed or projected population trajectories.

They do not implement WPP projection matching, interpolation, automatic residual
forcing, or infection-demography coupling.

## Implied residuals

`implied_demographic_residual()` compares an observed population trajectory to
the one-step demographic update implied by a supplied demographic process. For
each observed interval `[time_start, time_end]`, it:

1. takes the observed population at `time_start` as the starting state;
2. evaluates `demographic_derivative()` at `time_start`;
3. predicts the interval endpoint using the current explicit Euler step rule;
4. compares the predicted endpoint to the observed endpoint.

For one age group in one interval:

```text
dt = time_end - time_start
predicted_end = observed_start + dt * model_derivative
residual_count = observed_end - predicted_end
residual_rate = residual_count / (dt * observed_start)
```

`residual_count` is therefore the total interval gap, not a per-time flow. It
may be positive or negative. `residual_rate` is `NA` when `observed_start` is
zero.

The function uses observed exact times only. There is no interpolation, nearest
time lookup, or age transformation. If the supplied process has schedules, they
must be available at the exact interval start times used in the observed
trajectory.

## Residual-derived migration schedules

`residual_to_migration_schedule()` converts residual diagnostics into a
`MigrationSchedule()` using the current Euler accounting convention.

For `use = "count"`:

```text
migration_count = residual_count / dt
schedule time = time_start
```

This division matters because `MigrationSchedule()` treats `migration_count` as
a per-time additive flow in `demographic_derivative()`, while
`residual_count` is the total gap accumulated over the interval.

For `use = "rate"`:

```text
migration_rate = residual_rate
schedule time = time_start
```

`residual_rate` values are used directly. Missing residual rates are not
silently replaced; `use = "rate"` errors when any `residual_rate` is `NA`.

## Interpretation

Residual-derived migration schedules can be useful diagnostic or forcing
objects for checking what age-specific net flow would reconcile an observed
trajectory under agepi's current one-step Euler semantics.

They should not be treated as automatically predictive migration estimates.
They are tied to:

- the supplied observed intervals;
- the current process and schedules;
- the explicit Euler update convention;
- exact-time schedule lookup at interval starts.

Residual schedules should not be used outside their observed intervals unless
the user makes an explicit modelling assumption about how those residual flows
or rates should continue, repeat, or be interpolated. agepi does not make that
assumption for the user.

## Practical use

A typical diagnostic workflow is:

1. build a demographic process from ageing, fertility, mortality, and any known
   migration schedules;
2. store an observed or projected population trajectory as a `Demography()`
   object;
3. call `implied_demographic_residual(observed, process)`;
4. inspect the residual counts and rates;
5. optionally call `residual_to_migration_schedule()` to create a migration
   schedule for a subsequent explicit diagnostic simulation.

This is especially useful when asking how far closed demographic dynamics, or a
partially specified demographic process, are from a WPP-style population
trajectory under the current agepi step rule.
