# Agepi Performance Baseline

This directory holds the development-only benchmark and profiling harness for
`agepi`.

## What The Benchmarks Cover

- `deterministic SIR`: a simple age-structured SIR simulation using the package
  deterministic solver path.
- `deterministic SEIR`: the same setup with an exposed compartment so the
  progression path is exercised.
- `generic CompartmentModel`: a custom multi-compartment model with age-specific
  transition rates and weighted infectious compartments.
- `stochastic Gillespie`: a fixed-population stochastic SIR simulation using the
  direct Gillespie sampler.
- `deterministic SIR with demography`: coupled epidemic and demographic
  dynamics with fertility, mortality, and migration schedules.

The fixtures are defined in `fixtures.R` and are intended to be stable across
development runs. The benchmark suite uses explicit Euler for the deterministic
cases and a short time grid so it stays practical to run during development.

## How To Run The Benchmarks

```r
Rscript development/benchmarks/run_benchmarks.R
```

The script prints a compact summary with median elapsed time, iterations per
second, memory allocation, and garbage-collection rate. The benchmark suite is
kept deliberately small enough for routine development use.

## How To Run The Profiler

```r
Rscript development/benchmarks/profile_benchmarks.R
```

The profiler writes temporary `.rprof` files and prints the top functions by
total and self time. The default run profiles the deterministic generic
compartment path and the stochastic Gillespie path.

## Baseline Results

Recorded from the current development environment:

| Benchmark | Median elapsed | Iterations/sec | Memory allocated |
| --- | ---: | ---: | ---: |
| deterministic SIR | pending | pending | pending |
| deterministic SEIR | pending | pending | pending |
| generic CompartmentModel | pending | pending | pending |
| stochastic Gillespie | pending | pending | pending |
| deterministic SIR with demography | pending | pending | pending |

Profile summaries from the current environment:

### Benchmark Baseline

| Benchmark | Median elapsed | Iterations/sec | Memory allocated |
| --- | ---: | ---: | ---: |
| deterministic SIR | 11.2 ms | 89.32 | 1.03 MB |
| deterministic SEIR | 13.87 ms | 72.13 | 282.38 KB |
| generic CompartmentModel | 40.89 ms | 24.46 | 889.5 KB |
| stochastic Gillespie | 9.01 s | 0.11 | 49.51 MB |
| deterministic SIR with demography | 26.49 ms | 37.75 | 384.06 KB |

These values came from the current development environment and should be used
as the reference point for future performance work.

## Interpreting The Profile

The deterministic profile showed runtime concentrated in:

- `transition_rates()`
- `rates_to_derivative()`
- `deterministic_derivative()`
- `generic_transition_rates()`
- validation and schedule lookup helpers
- `paste()`, `unique()`, and `data.frame()` inside validation and rate-table
  assembly

In this baseline, `force_of_infection()` is still part of the deterministic call
chain, but the surrounding validation and rate-table construction are more
visible in the sampled trace.

The stochastic profile showed runtime concentrated in:

- `stochastic_propensities()`
- `transition_rates()`
- `force_of_infection()`
- validation helpers
- `stochastic_event_table()`
- `duplicated()`, `%in%`, `unique()`, `paste()`, and `data.frame()` during
  event-table assembly and matching
- `validate_age_structure()` and related validation helpers

In other words, the dominant work is still in the rate/propensity machinery and
the data-frame and matching operations used to assemble and order those tables.

## Notes

- Do not commit generated `.rprof` files or large benchmark outputs.
- Use these measurements as a baseline for future optimisation work.
