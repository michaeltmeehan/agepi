# General Calibration Framework Design

## 1. Executive Summary

`agepi` should add calibration support as a small set of dependency-free data
contracts and evaluators, not as a full inference engine. The central idea is a
`CalibrationProblem` that combines:

- observed data targets;
- parameter specifications;
- a simulation wrapper;
- mappings from simulation outputs to comparable predictions;
- observation models or user-supplied objective functions.

The first implementation should support deterministic and stochastic outputs,
ordinary compartment trajectories, and cumulative-flow outputs. It should not
implement optimisers, MCMC, SMC, ABC, posterior diagnostics, plotting, or
disease-specific likelihoods. Those should remain in external packages or user
code that call a stable `agepi` objective evaluator.

The minimal public entry point should be:

```r
evaluate_calibration_objective(problem, values, ...)
```

where `values` are unconstrained or native-scale parameter values, depending on
the declared `ParameterSpec`s. The evaluator should inject parameters, run the
simulation, map outputs to each target, evaluate target-level log likelihoods
or losses, and return a structured result with total objective value and
per-target diagnostics.

## 2. Current Package Capabilities Relevant To Calibration

The current public simulation surface already provides useful calibration
building blocks:

- `simulate_deterministic()` returns ordinary compartment trajectories as a long
  data frame with `time`, `compartment`, `age_group`, and `value`.
- With `cumulative_flows`, deterministic simulation returns a list containing
  `trajectory` and `cumulative`. The cumulative table has `time`,
  `cumulative_name`, `transition_id`, `from`, `to`, `age_group`, and `value`.
- `simulate_stochastic()` returns the same trajectory shape for fixed-population
  stochastic runs. With `return_events = TRUE`, it also returns an event log.
  With `cumulative_flows`, it returns event-log-derived cumulative counts.
- `CompartmentModel()`, `SIRModel()`, `SEIRModel()`, `transition_rates()`, and
  `rates_to_derivative()` expose disease-model structure without TB-specific
  assumptions.
- `state_long_to_vector()` and `state_vector_to_long()` define the
  compartment-major, age-group-minor state convention.
- `as_agepi_cases()` maps observed event records into model age groups.
- `as_cfr_data()` converts selected cumulative flows into case/death interval
  tables for downstream severity workflows.
- `compartment_totals()`, `age_group_totals()`, and `total_population()` provide
  simple trajectory summaries.
- `compare_demography_to_observed()` and `summarise_demography_comparison()`
  show the existing pattern for diagnostic comparison helpers.
- Validation helpers consistently use plain data frames, lists, and explicit
  `stop(..., call. = FALSE)` errors.

The calibration layer should follow these conventions: small list/data-frame
constructors, explicit validation, no required new dependencies, and no
assumption that a model is tuberculosis.

## 3. Proposed Calibration Use Cases

Initial use cases should include:

- fitting a deterministic model to age-specific prevalence-like compartment
  targets;
- fitting cumulative incidence or transition-flow targets from deterministic
  cumulative outputs;
- fitting stochastic simulations to count data through repeated simulations or
  through a user-supplied summary/objective;
- combining several target types, such as total incidence, age distribution of
  incidence, and age-specific prevalence;
- calibrating scalar and age-specific parameters in model inputs, initial
  states, contact matrices, susceptibility, infectiousness, transition rates,
  or demographic schedule values;
- exposing an objective function that external optimisers or samplers can call.

Non-goals for the first calibration layer:

- TB-specific targets, priors, model constructors, or reporting conventions;
- internal optimisation or sampling algorithms;
- automatic identifiability analysis;
- posterior summaries, convergence diagnostics, or plotting;
- complex censoring, reporting-delay, nowcasting, or survey-design likelihoods.

## 4. Minimal Viable Calibration API

The minimal API should be:

```r
ParameterSpec(name, lower = -Inf, upper = Inf, transform = "identity",
              initial = NULL, fixed = FALSE, inject = NULL)

ParameterSet(values, specs)

CalibrationTarget(name, data, output_mapping, observation_model,
                  weight = 1, required = TRUE)

OutputMapping(source, filter = NULL, aggregate = NULL,
              time = NULL, age = NULL, value = "value",
              transform = "identity")

ObservationModel(family, ..., log = TRUE)

SimulationSpec(simulate, base_args, parameter_injector = NULL,
               stochastic = FALSE, replicates = 1, seed = NULL)

CalibrationProblem(parameters, simulation, targets)

evaluate_calibration_objective(problem, values, return_details = TRUE)
```

The constructors can return simple S3 lists. The first implementation can also
allow plain lists with validation, but exported constructors make user code less
fragile.

`evaluate_calibration_objective()` should be deterministic for deterministic
models. For stochastic models it should either use controlled seeds or accept
that the objective is noisy, with this behaviour explicit in `SimulationSpec`.

## 5. Object And Data Structures

### Observed Data

Observed data should be plain data frames. A minimal target data frame should
contain:

```text
time, value
```

Optional standard columns should include:

```text
age_group, target, denominator, sd, se, weight, lower, upper
```

The `target` column is useful when one table contains several measurement
series. `denominator` supports binomial and multinomial observations. `sd` or
`se` supports Gaussian observations with known measurement uncertainty.

Age-stratified data should use the existing `age_group` labels from an
`AgeStructure()`. Exact matching should be the first implementation. Existing
age transformation helpers can support explicit pre-processing outside
calibration; automatic rebinning during objective evaluation should be deferred.

### CalibrationTarget

A `CalibrationTarget` should bind one observed table to one mapping and one
observation model:

```r
CalibrationTarget(
  name = "annual_incidence_by_age",
  data = observed_incidence,
  output_mapping = OutputMapping(
    source = "cumulative",
    filter = list(cumulative_name = "infections"),
    aggregate = list(time = "interval", age_group = "keep")
  ),
  observation_model = ObservationModel("negative_binomial", size = 20)
)
```

Targets should evaluate independently, then combine by summing weighted
negative log likelihoods or losses.

### CalibrationProblem

`CalibrationProblem` should store validated specs, not fitted results:

```text
parameters: list of ParameterSpec
simulation: SimulationSpec
targets: list of CalibrationTarget
```

A later `CalibrationResult` could store optimiser output, accepted samples, or
posterior summaries, but that should not be part of the first milestone.

## 6. Parameter Handling

`ParameterSpec` should declare:

- `name`: stable unique parameter name;
- `lower`, `upper`: finite or infinite native-scale bounds;
- `transform`: `"identity"`, `"log"`, `"logit"`, or a user function pair;
- `initial`: optional starting value for external algorithms;
- `fixed`: whether the parameter is held constant;
- `inject`: optional injection rule or function.

The evaluator should support two scales:

- native scale: values are what the simulator expects;
- working scale: values are unconstrained or transformed for optimisers.

The first implementation should include helper functions such as:

```r
parameters_to_working(values, specs)
parameters_from_working(values, specs)
validate_parameter_values(values, specs, scale = c("native", "working"))
```

Bounds should be checked after back-transformation. For `log`, native values
must be positive. For `logit`, finite lower and upper bounds are required.

Parameter injection should be general. Good first options are:

- a named `parameter_injector(args, parameters)` function in `SimulationSpec`;
- optional per-parameter `inject` functions for small changes;
- no string-path mutation in the first milestone unless a robust recursive
  setter is designed and tested.

The recommended first pattern is:

```r
SimulationSpec(
  simulate = simulate_deterministic,
  base_args = list(...),
  parameter_injector = function(args, p) {
    args$beta <- p[["beta"]]
    args$model <- SIRModel(gamma = p[["gamma"]])
    args
  }
)
```

This keeps injection explicit and allows users to rebuild model objects,
contact matrices, age-specific vectors, or initial states without the package
guessing object internals.

## 7. Output Mapping

`OutputMapping` should convert simulation output into rows comparable with one
target data frame. It should support:

- selecting output source: `"trajectory"`, `"cumulative"`, `"events"`, or a
  user function;
- filtering by columns such as `compartment`, `age_group`, `cumulative_name`,
  `transition_id`, `from`, or `to`;
- aggregating across compartments, age groups, cumulative names, or time
  intervals;
- converting cumulative values to interval increments;
- computing proportions or rates when a denominator mapping is supplied;
- joining predictions to observations by exact keys.

The first implementation should avoid hidden interpolation. It should join on
exact `time`, `age_group`, and optional `target` keys unless the mapping
explicitly requests an aggregation policy.

Example mapping types:

```text
trajectory compartment value:
  source = "trajectory", filter = list(compartment = "I")

trajectory prevalence:
  numerator = compartment "I"; denominator = all compartments by time-age

cumulative incidence:
  source = "cumulative"; filter cumulative_name = "infections";
  time aggregation = interval increments

age distribution:
  source values by age; normalise within time to proportions
```

The mapping result should be a data frame with at least:

```text
time, observed, predicted
```

and optional keys such as `age_group`, `target`, `denominator`, and target-level
metadata.

## 8. Observation Models

The first observation models should cover common aggregate calibration targets:

- Gaussian for continuous quantities with additive error;
- lognormal for positive continuous quantities with multiplicative error;
- Poisson for counts;
- negative binomial for overdispersed counts;
- binomial for proportions or positive counts out of denominators;
- multinomial for age distributions or other compositions;
- user-supplied objective functions.

Each `ObservationModel` should evaluate one mapped target table and return a
data frame of pointwise contributions plus a scalar objective contribution.

Suggested required columns by family:

```text
gaussian: observed, predicted, sd or sigma
lognormal: observed > 0, predicted > 0, sdlog or sigma
poisson: observed count, predicted mean
negative_binomial: observed count, predicted mean, size or dispersion
binomial: observed successes or proportion, denominator, predicted probability
multinomial: observed counts, predicted probabilities, group key
user: mapped table passed to function
```

The package should use base R `stats` densities only. Negative binomial can use
`stats::dnbinom()`. No additional likelihood packages are needed.

## 9. Objective-Function Evaluation

`evaluate_calibration_objective()` should perform:

1. validate and transform parameter values;
2. inject parameters into simulation arguments;
3. run the simulation;
4. map simulation output to each target;
5. evaluate each observation model;
6. combine weighted target contributions.

The return value should be structured:

```r
list(
  value = <scalar objective, usually negative log likelihood>,
  log_likelihood = <scalar, if available>,
  parameters = <native-scale named vector/list>,
  target_results = <data frame with target, value, log_likelihood, n>,
  mapped = <optional list of mapped target tables>,
  simulation = <optional simulation output>
)
```

The objective convention should be minimisation by default:

```text
objective = -sum(weighted log_likelihood)
```

For user-supplied loss functions, the function should declare whether it
returns a log likelihood or a loss. Mixing both is possible but should be
explicit in diagnostics.

Failures during simulation or mapping should return an infinite objective only
when requested by `on_error = "infinite"`. The default for direct use should be
to error clearly, because silent infinite objectives can hide invalid model
configuration.

## 10. Deterministic Workflow

For deterministic models, the objective is a deterministic function of the
parameter values, solver method, and requested times. This is the simplest first
workflow:

```r
problem <- CalibrationProblem(
  parameters = list(
    ParameterSpec("beta", lower = 0, transform = "log"),
    ParameterSpec("gamma", lower = 0, transform = "log")
  ),
  simulation = SimulationSpec(
    simulate = simulate_deterministic,
    base_args = list(...),
    parameter_injector = function(args, p) {
      args$beta <- p[["beta"]]
      args$model <- SIRModel(gamma = p[["gamma"]])
      args
    }
  ),
  targets = list(...)
)

objective <- function(x) {
  evaluate_calibration_objective(problem, x, scale = "working")$value
}
```

External callers can pass `objective` to `stats::optim()`, `nloptr`,
`optimx`, `rstan`, `cmdstanr`, `TMB`, `pomp`, or other tools. `agepi` should
not wrap those tools initially.

## 11. Stochastic Workflow

Stochastic calibration needs a different policy because one simulation is a
random draw. The calibration framework should support, but not over-prescribe:

- single-replicate noisy objective evaluation;
- fixed-seed common-random-number evaluation for optimisation diagnostics;
- multiple replicates with mapped outputs summarised by mean, median, quantile,
  or empirical likelihood supplied by the user;
- user-supplied objective functions for simulation-based methods.

`SimulationSpec` should include:

```text
stochastic = TRUE
replicates = 1
seed = NULL or integer
replicate_seeds = optional integer vector
replicate_summary = "mean" or user function
```

The first implementation should not try to produce exact stochastic likelihoods
for latent Markov jump processes. It can evaluate observation likelihoods on
replicate summaries, or delegate to a user function that receives all replicate
outputs.

For stochastic cumulative outputs, the mapping should remember that cumulative
values are realised event counts. This makes Poisson, negative binomial,
binomial, and multinomial observation models natural, but the process model
uncertainty and observation uncertainty remain conceptually distinct.

## 12. Testing Strategy

Initial tests should focus on contracts, not inference performance:

- constructors validate required fields and reject malformed specs;
- parameter transformations round-trip and enforce bounds;
- parameter injection updates scalar arguments, model constructors, age-specific
  vectors, and initial states in small examples;
- trajectory mappings select compartments, aggregate over age, and join to
  observations by exact keys;
- cumulative mappings select `cumulative_name` or `transition_id`, convert
  cumulative values to increments, and preserve age groups;
- observation models return hand-checkable log likelihoods for Gaussian,
  lognormal, Poisson, negative binomial, binomial, and multinomial examples;
- deterministic objective evaluation is reproducible and decomposes by target;
- stochastic objective evaluation is reproducible when seeds are fixed;
- invalid simulations or invalid mapped predictions fail clearly;
- examples run through `test-public-api-examples.R` once user-facing examples
  are added.

Testing should use tiny SIR/SEIR/generic models already present in the test
suite. No TB fixture is needed for core calibration tests.

## 13. Implementation Roadmap

### Milestone 1: Contracts And Deterministic Evaluation

Add a new `R/calibration.R` or a small set of focused files:

- `R/calibration_parameters.R`
- `R/calibration_targets.R`
- `R/calibration_mapping.R`
- `R/calibration_observation_models.R`
- `R/calibration_evaluate.R`

Implement constructors, validators, parameter transformation, explicit
simulation injection, basic trajectory/cumulative mappings, Gaussian/Poisson
negative log likelihoods, and `evaluate_calibration_objective()`.

Recommended first observation families: Gaussian and Poisson, plus
user-supplied objective. These cover continuous trajectories and simple counts
while keeping implementation small. Negative binomial, binomial, lognormal, and
multinomial can follow once the mapping contract is stable.

### Milestone 2: Broader Observation Models

Add lognormal, negative binomial, binomial, and multinomial models. Add tests
for zero handling, denominators, invalid probabilities, dispersion parameters,
and grouped multinomial evaluation.

### Milestone 3: Stochastic Replicate Support

Add `replicates`, seed handling, replicate summaries, and stochastic-specific
diagnostics. Keep exact inference algorithms out of scope.

### Milestone 4: Convenience Adapters

Add small helpers for common target tables, such as:

- `as_calibration_target_cases()`;
- `as_calibration_target_prevalence()`;
- `as_calibration_target_cumulative_flow()`.

These should remain disease-neutral and should mostly validate or reshape data.

### Milestone 5: Vignettes And External Optimiser Examples

Add examples showing how to call `stats::optim()` or another optional external
tool from user code. Do not add dependencies or imports.

## 14. Open Design Questions

- Should `values` passed to `evaluate_calibration_objective()` default to
  native scale or working scale?
- Should `ParameterSpec` include priors now, or should priors be left entirely
  to external Bayesian packages until a later milestone?
- Should output mappings be list specifications only, or should formula-like
  helpers be added later for ergonomics?
- Should time interpolation ever happen inside calibration, or should users
  always request simulation times matching observation times?
- Should age aggregation be automatic when observation age bins are coarser
  than model bins, using existing `aggregate_age_vector()` semantics, or should
  observed data always be pre-aligned?
- How should missing observations be represented: omitted rows, `NA` values
  with explicit handling policy, or both?
- How much stochastic replicate summarisation belongs in `agepi` before users
  should move to `pomp`, particle MCMC, ABC-SMC, or custom simulation-based
  inference tools?
- Should mapped predictions be cached when external optimisers evaluate several
  target functions at the same parameter values?

## Recommended Boundary

Implement inside `agepi`:

- parameter specs, transforms, bounds, and injection hooks;
- observed-target validation;
- deterministic and stochastic simulation wrappers;
- output mapping for trajectory and cumulative tables;
- simple observation-model log likelihoods;
- structured objective evaluation and diagnostics.

Leave to external packages:

- optimisers and samplers;
- priors and posterior computation;
- SMC, ABC, particle filters, PMCMC, HMC, variational inference;
- convergence diagnostics and posterior plotting;
- complex survey/censoring/reporting-delay models;
- disease-specific calibration workflows.

