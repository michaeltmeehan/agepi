# Demography Workflows

`agepi` supports demographic change as a deterministic schedule-based input.
This guide explains how the pieces fit together, how to choose between the
available workflows, and what each workflow assumes mechanically.

The short version is:

- `AgeStructure()` defines the age grid shared by demography, disease states,
  and contact matrices;
- `DemographicProcess()` stores fertility, mortality, migration, and ageing
  inputs on that grid;
- `simulate_demography()` evolves population counts only;
- `simulate_deterministic()` can consume a `demographic_process` alongside an
  epidemic model;
- `demographic_process_from_wpp()` is a WPP-backed convenience builder, not a
  projection engine of its own.

## 1. Overview

In `agepi`, demography means the non-disease processes that change population
structure over time:

- births;
- background mortality;
- ageing between adjacent age groups;
- migration, if supplied.

These flows interact with epidemic compartments. In a disease model, the state
is usually `compartment x age group`, so demographic changes can move people
between age groups while disease transitions move them between compartments.

The package currently treats demography as deterministic and schedule-based.
You provide age-specific fertility, mortality, and migration schedules, then
choose the workflow that matches the question:

- demographic-only projection or diagnostics;
- coupled epidemic-demography simulation;
- WPP-backed input assembly;
- WPP projection replay for comparison.

## 2. Core Mental Model

The central model is simple:

- the model state is a count for each compartment in each age group;
- births add people to a designated entry compartment or the youngest age
  group, depending on the workflow;
- mortality removes people from the active population;
- ageing moves survivors to the next age group;
- migration adds or removes people according to a declared policy.

This is different from holding the population fixed. A fixed-population model
keeps the age distribution static unless disease transitions change it. A
demographic model allows the age distribution itself to evolve, which then
changes the denominator for force of infection, the number of births, and the
relative weight of each age group over time.

## 3. Age Structures And Age Grids

`AgeStructure()` is the shared coordinate system for almost everything in the
package. Demography, epidemic compartments, and contact matrices all assume
the same age groups unless you explicitly aggregate or adapt them.

Common grids include:

- a model age grid, such as `0-4`, `5-9`, `10+`;
- a reporting age grid, which may be broader than the model grid;
- WPP single-year age grids, which are useful for annual cohort workflows;
- five-year or broader age groups for more compact models;
- open-ended terminal age groups such as `65+` or `100+`.

The package has three related ideas here:

1. `AgeStructure()` defines a grid.
2. `AgeGridMapping()` defines an exact nested mapping between two grids when
   the boundaries line up.
3. `aggregate_*_age_grid()` helpers collapse model output to a reporting grid.

This matters because demography and epidemic outputs are usually generated on a
model grid, but analysis may be done on a broader reporting grid. If the grids
do not line up exactly, use `AgeGridMapping()` and the relevant aggregation
helper instead of assuming implicit averaging.

WPP workflows often use `wpp_age_structure_1year()` or `wpp_age_structure_5year()`
to create grids compatible with the source tables they are derived from.

## 4. DemographicProcess Objects

`DemographicProcess()` is the package's container for demography inputs. It
stores:

- `age_structure`;
- `ageing_operator`;
- `fertility_schedule`;
- `fertility_exposure_fraction`;
- `mortality_schedule`;
- `migration_schedule`;
- `mode`;
- `times`.

The object does not simulate anything by itself. It is a validated data
container that downstream functions consume.

The schedule objects are time-indexed:

- fertility schedules hold annual fertility rates;
- mortality schedules hold annual hazards;
- migration schedules hold either counts or rates, depending on how they were
  standardised;
- `times` records a common schedule grid when all supplied schedules share the
  same times, otherwise it is `NULL`.

`time_policy` controls how scheduled values are read:

- `"exact"` requires the requested time to be present in the schedule;
- `"step"` uses left-continuous interval-start lookup;
- `"linear"` interpolates rate-like schedules between surrounding schedule
  times without extrapolation.

Current implementation boundaries matter here:

- there is no stochastic demographic simulation;
- there is no implicit interpolation of population tables;
- `linear` applies only to rate-like fertility, mortality, and migration
  schedules;
- `exact` is the strictest mode and is useful when you want schedule lookup to
  fail rather than silently infer a value.

`fertility_exposure_fraction` is a scalar multiplier on fertility exposure. The
default convention uses the full age-specific population as exposure. A value
such as `0.5` is an approximate female-exposure correction for total-population
models, not a substitute for sex-structured demography.

## 5. Standalone Demography Simulation

Use `simulate_demography()` when you want the population process only.
Typical uses include:

- testing fertility, mortality, and migration assumptions before coupling them
  to disease;
- comparing demographic projections against observed or reference population
  tables;
- validating age-grid logic on its own;
- building a background demographic trajectory for later reporting or
  diagnostics.

The function takes:

- a validated `DemographicProcess()`;
- an initial age-specific population vector;
- a sequence of times;
- optional `method` and `time_policy` settings;
- an optional `ageing_policy`.

The returned data frame has columns `time`, `age_group`, and `population`.
Each row is a model age group at a requested time.

The two ageing policies are:

- `"exponential"`: the standard derivative-based continuous-time path;
- `"annual_cohort"`: a discrete annual-cohort path that updates once per year
  and requires a complete 1-year age grid ending in an open-ended terminal
  group.

Annual-cohort demography is useful when you want a WPP-like cohort progression
story. It is not the same as the general continuous-time demographic path.

## 6. Coupled Epidemic-Demography Simulation

`simulate_deterministic()` accepts a `demographic_process` and applies it while
the epidemic model runs.

Mechanically, the coupling does four things:

- births enter the configured birth compartment, which is `S` for the built-in
  SIR and SEIR models;
- mortality removes people from the active compartments;
- ageing moves each compartment independently across age groups;
- migration is allocated according to `migration_policy`.

The available migration policies are:

- `"susceptible"`: put age-total migration into the configured susceptible or
  migration compartment;
- `"proportional"`: allocate age-total migration across compartments by their
  current shares in each age group;
- `"error"`: refuse non-zero migration unless an explicit allocation is not
  needed.

This is still a first-pass coupling. It does not add disease-induced mortality,
vaccination, waning immunity, or compartment-specific demographic rates.

The important practical point is that demography and transmission interact
through the evolving age-specific population. The force of infection sees the
current population in each source age group, so demographic change can alter
both denominators and the age profile of exposure over time.

## 7. WPP-Derived Demographic Workflow

`demographic_process_from_wpp()` is the main convenience wrapper for WPP-based
demography. It loads the relevant WPP inputs, standardises them to agepi
objects, and returns a structured object with these top-level fields:

- `country`;
- `years`;
- `age_structure`;
- `population`;
- `demographic_process`;
- `schedules`;
- `inputs`.

The components mean:

- `population` is the WPP-derived `Demography()` object, usually used as the
  reference population trajectory or as a source for initial values;
- `demographic_process` is the `DemographicProcess()` object you pass to
  simulation functions;
- `schedules` contains the fertility, mortality, and optional migration
  schedules used to build the process;
- `inputs` contains the standardised WPP tables used to create those schedules
  and the population trajectory.

The practical workflow is:

1. choose a country;
2. choose the year range;
3. choose the model age structure;
4. decide whether migration should be included;
5. choose a fertility exposure fraction when total-population fertility should
   be down-weighted toward female exposure;
6. extract the initial population at the first year;
7. simulate or compare the result against WPP reference values.

This workflow is useful for benchmarking and sanity checks. It is not an
automatic guarantee of exact WPP reproduction.

## 8. Projection-Backed WPP Workflow

There are two distinct WPP-style workflows in the package:

- simulate demography from fertility, mortality, and migration schedules;
- replay a precomputed WPP population projection trajectory.

The second path is what `examples/wpp_projection_backed_demography.R` shows.
It uses `population_trajectory_from_wpp()` to load the projected age-specific
population table and then `projection_population_vector()` to pull out one
exact year. That is useful when you want to inspect or reuse WPP's published
population trajectory directly.

This is not the same thing as simulating a mechanistic demographic process. It
reads back a trajectory that was already produced by another model.

## 9. Validation And Benchmarking

The validation workflow is about comparison, not automatic calibration.
`examples/wpp_demography_validation.R` shows how to interpret a benchmark run.

Useful comparisons are:

- total population by year;
- age-structure by year;
- component balance, such as births minus deaths and net change;
- the effect of migration on total and age-specific residuals;
- residuals converted back into migration-style flows for calibration
  experiments.

Two common mistakes are worth calling out:

- matching total population does not imply matching the age structure;
- very old age groups can show large relative errors because the denominator is
  small, so absolute errors are often more informative there.

## 10. Choosing The Right Workflow

| Workflow | Main functions | Use case | Strengths | Limitations | Example file |
| --- | --- | --- | --- | --- | --- |
| Static population epidemic model | `simulate_deterministic()` | Disease dynamics with fixed population | Simple, direct, easy to reason about | No demographic change | `examples/mock_sir_deterministic.R` |
| Mock demographic process | `DemographicProcess()`, `simulate_demography()` | Teaching and unit tests | Small, dependency-free, easy to inspect | Invented rates only | `examples/mock_demographic_workflow.R` |
| Standalone demographic simulation | `build_demographic_process()`, `simulate_demography()` | Population-only projection or diagnostics | Clear demographic mechanics | No disease coupling | `examples/mock_demographic_workflow.R` |
| Coupled epidemic-demography simulation | `simulate_deterministic()` with `demographic_process` | Disease and population change together | Shows births, deaths, ageing, and migration in one run | Still deterministic and first-pass | `examples/mock_seir_demography.R` |
| WPP-derived demographic process | `demographic_process_from_wpp()` | Use WPP as a demographic input source | Convenient, structured, benchmark-friendly | Not exact WPP reproduction | `examples/wpp_demography_validation.R` |
| WPP projection-backed replay | `population_trajectory_from_wpp()`, `projection_population_vector()` | Inspect or reuse published projections | Direct access to WPP trajectories | Not mechanistic demography | `examples/wpp_projection_backed_demography.R` |
| Kiribati TB scaffold | `demographic_process_from_wpp()`, `load_contact_matrix_source()`, `simulate_deterministic()` | Public-data model scaffold | Combines WPP demography and a contact-matrix source | Provisional and uncalibrated | `examples/kiribati_tb_realistic_demography.R` |

## 11. Limitations And Current Boundaries

The current boundaries are deliberate:

- demography is deterministic, not stochastic;
- population interpolation is not automatic;
- age-grid compatibility matters everywhere;
- annual-cohort demography assumes a complete 1-year grid and an open-ended
  terminal group;
- migration allocation is simplified and policy-driven;
- WPP-derived workflows are benchmarking tools, not exact-reproduction
  guarantees;
- disease-specific mortality conventions still need care when modelling
  diseases where mortality is not purely background mortality;
- sex structure is simplified or absent in current workflows;
- contact matrices are handled separately from demography.

## 12. Worked Mini-Example

The following example is intentionally small and dependency-free. It shows the
mechanics of a demographic-only workflow on a simple three-group grid. Because
the schedule values are meant to stay constant over the whole simulation
window, the example repeats the same rates at the final time point and uses
`time_policy = "step"` to carry them forward between listed schedule times.
Exact-time lookup would otherwise require schedule values at every simulated
time.

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9", "10+"),
  lower_bounds = c(0, 5, 10),
  upper_bounds = c(4, 9, Inf)
)

fertility <- FertilitySchedule(
  data.frame(
    time = c(0, 2),
    age_group = "10+",
    fertility_rate = c(0.02, 0.02),
    stringsAsFactors = FALSE
  ),
  ages
)

mortality <- MortalitySchedule(
  data.frame(
    time = c(0, 0, 0, 2, 2, 2),
    age_group = rep(ages$age_groups, times = 2),
    mortality_rate = c(0.01, 0.005, 0.02, 0.01, 0.005, 0.02),
    stringsAsFactors = FALSE
  ),
  ages
)

process <- build_demographic_process(
  age_structure = ages,
  fertility_schedule = fertility,
  mortality_schedule = mortality
)

initial_population <- c(500, 450, 800)
simulated <- simulate_demography(
  process = process,
  initial_state = initial_population,
  times = c(0, 1, 2),
  time_policy = "step"
)

total_population(simulated)
age_group_totals(simulated)
```

If you want a reporting grid, add a matching `AgeGridMapping()` and aggregate
the output:

```r
reporting <- AgeStructure(
  age_groups = c("0-9", "10+"),
  lower_bounds = c(0, 10),
  upper_bounds = c(9, Inf)
)

mapping <- AgeGridMapping(ages, reporting, open_ended = "include")
aggregate_demography_trajectory_age_grid(simulated, mapping)
```

## 13. Existing Examples

These scripts illustrate the main workflows:

- `examples/mock_demographic_workflow.R`: demographic-only workflow with
  diagnostics and residuals.
- `examples/mock_seir_demography.R`: epidemic model with demographic turnover.
- `examples/annual_cohort_sir_demography.R`: annual-cohort coupling path.
- `examples/wpp_demography_validation.R`: WPP-connected benchmark and residual
  analysis.
- `examples/wpp_projection_backed_demography.R`: replay of a published WPP
  projection trajectory.
- `examples/kiribati_tb_realistic_demography.R`: public-data Kiribati TB
  scaffold with WPP demography and a contact-matrix source.

## 14. Related Documentation

- [README.md](../README.md)
- [docs/model_conventions.md](model_conventions.md)
- [docs/external_data_adapters.md](external_data_adapters.md)
- [docs/demographic_residuals.md](demographic_residuals.md)
- [examples/mock_demographic_workflow.R](../examples/mock_demographic_workflow.R)
- [examples/wpp_demography_validation.R](../examples/wpp_demography_validation.R)
