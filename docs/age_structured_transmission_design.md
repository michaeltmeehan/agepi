# Age-Structured Transmission Prototype: Design Document

## 1. Project aim

Build a small but extensible R prototype for age-structured infectious disease transmission models. The current package implementation supports deterministic age-structured SIR and SEIR simulation, including first-pass demographic coupling, with small utility layers for simulation summaries, age-vector transformation, demography-table access, demographic schedules, WPP-style adapters, and contact-matrix validation/coercion/aggregation.

Checkpoint note: this design document preserves much of the initial prototype
scope. Demographic-only ODE components, WPP-style demographic adapters, residual
diagnostics, SEIR infection-only simulation, linear schedule interpolation, and
optional deSolve backends have since been added; see `README.md`,
`docs/external_data_adapters.md`, `docs/demographic_residuals.md`, and
`examples/mock_demographic_workflow.R` for the current public surface.

The first implementation should be deliberately simple: a deterministic age-structured SIR model with one static contact matrix. However, the code should be designed so that later extensions can support:

- different disease compartment structures;
- different age-bin definitions;
- time-varying demography;
- age-specific susceptibility, infectiousness, morbidity, and mortality;
- different contact/mixing models;
- stochastic simulation;
- eventual porting of the simulation core to Julia.

The prototype should prioritise clean abstractions and testable model components over broad feature coverage.

---

## 2. Initial scope

### Included in first prototype

| Component | Current implementation |
|---|---|
| Language | R |
| Disease model | Age-structured SIR and SEIR |
| Simulation | Deterministic explicit Euler; optional deSolve for infection-only SIR/SEIR and SIR/SEIR-demography |
| Age bins | User-defined model age bins |
| Mixing | Static age-specific contact matrix |
| Susceptibility | Age-specific vector |
| Infectiousness | Age-specific vector, default all ones |

### Not currently implemented

| Component | Reason deferred |
|---|---|
| WPP projection matching, population interpolation, and automatic residual forcing | Schedule-level step and linear rate lookup exist, but population interpolation and projection matching remain out of scope |
| `socialmixr` or `conmat` package dependencies | Current contact support is dependency-free coercion and exact aggregation |
| Contact-matrix splitting or general rebinning | Current `transform_contact_matrix()` supports exact aggregation only |
| Other compartment structures | Current SIR and SEIR support use fixed compartment definitions |
| Additional adaptive or external ODE solver features | Optional deSolve support exists for the documented narrow combinations |
| Stochastic simulation | Design for it, but do not implement yet |
| Vaccination | Requires additional compartment conventions |
| Waning immunity | Can be added once transition framework is stable |
| Calibration/inference | Depends on stable model/simulation interface |
| Multi-country batching | Should come after one-country workflow is correct |
| Julia implementation | Design should make later port straightforward |

---

## 3. Core design principle

The prototype should be built around a **transition-rate interface**, not around hard-coded differential equations.

The model should expose a function that returns transition rates between compartments and age groups at time `t`.

For deterministic simulation, these transition rates are converted into derivatives.

For future stochastic simulation, the same transition rates can be interpreted as event intensities.

This design should support the following pathway:

```text
model specification
      ↓
transition rates
      ↓
deterministic derivatives
      ↓
Euler or deSolve simulation
```

Later:

```text
model specification
      ↓
transition rates
      ↓
event intensities
      ↓
Gillespie / tau-leaping / discrete-time stochastic simulation
```

---

## 4. Core epidemiological equation

The reusable force of infection should be:

$$
\lambda_a(t) =
\beta(t) s_a(t)
\sum_b C_{ab}(t)
\frac{\iota_b(t) I_b(t)}{N_b(t)}
$$

where:

| Symbol | Meaning |
|---|---|
| $a$ | recipient age group |
| $b$ | source age group |
| $\lambda_a(t)$ | force of infection for age group $a$ |
| $\beta(t)$ | transmission scaling parameter |
| $s_a(t)$ | age-specific susceptibility |
| $C_{ab}(t)$ | contact rate from age group $a$ to age group $b$ |
| $\iota_b(t)$ | age-specific infectiousness |
| $I_b(t)$ | infectious individuals in age group $b$ |
| $N_b(t)$ | total population in age group $b$ |

The function implementing this should not be tied to SIR specifically. It should accept the relevant state, population, contact matrix, susceptibility vector, infectiousness vector, and transmission parameter.

---

## 5. Main objects

Use simple R lists initially, but keep names and fields stable so they can later become formal classes or Julia structs.

### 5.1 `AgeStructure`

Defines model age bins.

Suggested fields:

```r
age_structure <- list(
  age_groups = c("0-4", "5-9", "10-14", "15-19", "20-29",
                 "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
  n_age_groups = 11,
  lower_bounds = c(0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80),
  upper_bounds = c(4, 9, 14, 19, 29, 39, 49, 59, 69, 79, Inf)
)
```

Required checks:

- `lower_bounds`, `upper_bounds`, and `age_groups` must have the same length.
- `n_age_groups` must equal `length(age_groups)`.
- Age bins must be non-overlapping.
- Age bins must be sorted.
- The final bin may be open-ended.
- Age groups must be unique.

---

### 5.2 `Demography`

Stores validated age-specific population tables aligned to the model age bins.
Current `Demography()` support is deliberately minimal: validation, sorted
storage, and exact-time population accessors. Demographic-only ODE components
are implemented separately and do not affect infection simulation dynamics.

Current constructor:

```r
Demography(demography, age_structure)
```

where `demography` is a data frame with:

| Field | Shape |
|---|---|
| `time` | numeric time point |
| `age_group` | one of `age_structure$age_groups` |
| `population` | non-negative finite population value |

Implemented helpers:

```r
validate_demography_table(demography, age_structure)
demography_times(demography)
demography_population_vector(demography, time)
demography_population_table(demography, time = NULL)
```

The current accessors require exact available time points. Interpolation,
nearest-year lookup, WPP integration, fertility, mortality, births, deaths,
ageing, migration, and demographic projection dynamics remain future work.

---

### 5.3 Contact matrices and future `MixingModel`

A `MixingModel` object is not currently implemented; current functions accept a
contact matrix directly. The package does include contact-matrix validation,
coercion, and exact aggregation helpers.

Current helpers:

```r
validate_contact_matrix(contact_matrix, age_structure = NULL)
as_agepi_contact_matrix(x, age_structure = NULL, orientation = ..., transpose = FALSE)
transform_contact_matrix(contact_matrix, from_age_structure, to_age_structure, population)
```

Expected dimensions:

```text
number_of_age_groups × number_of_age_groups
```

The entry `C[a, b]` represents contacts made by individuals in recipient age group `a` with individuals in source age group `b`. This matches the implemented force-of-infection convention: rows are recipients and columns are sources.

Implemented validation checks:

- Contact matrix must be square.
- Matrix dimensions must match the number of model age groups.
- Entries must be non-negative.
- Missing values are not allowed.

`as_agepi_contact_matrix()` currently supports numeric matrices, numeric data
frames, socialmixr-like lists with a numeric `matrix` element, and conmat-style
long data frames with `age_group_from`, `age_group_to`, and `contacts` columns.
It does not add `socialmixr` or `conmat` dependencies.

`transform_contact_matrix()` currently supports exact fine-to-coarse
aggregation only, where every target age bin is an exact union of complete
source age bins. Source-bin splitting and general rebinning remain future work.

Later versions may allow:

```r
get_contact_matrix <- function(t, mixing_model, demography) {
  ...
}
```

---

### 5.4 Future `RiskModel`

Would store age-specific susceptibility, infectiousness, and severity parameters. This object is not currently implemented; current functions accept susceptibility and infectiousness vectors directly.

Initial version:

```r
risk_model <- list(
  age_structure = age_structure,
  susceptibility = susceptibility,
  infectiousness = infectiousness,
  morbidity_risk = morbidity_risk,
  mortality_risk = mortality_risk
)
```

Required checks:

- All age-specific vectors must have length equal to the number of model age groups.
- Values must be non-negative.
- Susceptibility and infectiousness default to one for all age groups if not supplied.
- Morbidity and mortality risks may initially be used only for derived outputs.

---

### 5.5 `DiseaseModel`

Defines compartments and disease-specific transitions.

Initial SIR version:

```r
disease_model <- list(
  name = "SIR",
  compartments = c("S", "I", "R"),
  infectious_compartments = "I",
  transitions = list(
    infection = list(from = "S", to = "I"),
    recovery = list(from = "I", to = "R")
  ),
  parameters = list(
    beta = 0.03,
    gamma = 1 / 5
  )
)
```

Design requirements:

- Do not hard-code the assumption that there are exactly three compartments throughout the codebase.
- It is acceptable for the first derivative function to support only SIR, but state handling should be generic over compartments.
- Compartment ordering should be explicit and consistent.

---

### 5.6 Future `SimulationProblem`

Would combine all components needed to run a simulation. This object is not currently implemented; `simulate_deterministic()` currently takes explicit arguments.

```r
problem <- list(
  age_structure = age_structure,
  demography = demography,
  mixing_model = mixing_model,
  risk_model = risk_model,
  disease_model = disease_model,
  initial_state = initial_state,
  t_start = 0,
  t_end = 365,
  output_times = seq(0, 365, by = 1)
)
```

Required checks:

- All subcomponents must use the same `AgeStructure`.
- Initial state must match the model compartments and age groups.
- Time range must be valid.
- Output times must lie within the simulation interval.

---

## 6. State representation

The model should support a human-readable long-format state table and a solver-compatible numeric vector.

### Long-format state

```r
initial_state_long <- data.frame(
  age_group = c("0-4", "0-4", "0-4", "5-9", "5-9", "5-9"),
  compartment = c("S", "I", "R", "S", "I", "R"),
  value = c(100000, 10, 0, 110000, 5, 0)
)
```

### Solver vector

Use compartment-major ordering:

```text
S_1, S_2, ..., S_A,
I_1, I_2, ..., I_A,
R_1, R_2, ..., R_A
```

where `A` is the number of age groups.

Implemented helper functions:

```r
state_long_to_vector(state_long, age_structure, compartments)
state_vector_to_long(state_vector, age_structure, compartments)
aggregate_age_vector(values, from_age_structure, to_age_structure)
segregate_age_vector(values, from_age_structure, to_age_structure, weights)
transform_age_vector(values, from_age_structure, to_age_structure, weights = NULL, split_method = ...)
```

No matrix state helper is currently implemented. Numeric state vectors are interpreted by position only; names on numeric vectors are ignored when converting back to long form or simulating.

`aggregate_age_vector()` supports aggregation only when every target age bin is
an exact union of complete source age bins.

`segregate_age_vector()` supports coarse-to-fine splitting when every target age
bin is fully nested inside exactly one source age bin and explicit target
weights are supplied.

`transform_age_vector()` supports identity, aggregation, segregation, and mixed
exact transformations when source and target bins align to a common set of
boundaries. Its `split_method` controls source-bin splitting for mixed exact
transformations.

---

## 7. Transition-rate interface

The central model function should return a transition table.

Current function signature:

```r
transition_rates <- function(
  state,
  model,
  age_structure,
  contact_matrix,
  beta = 1,
  susceptibility = NULL,
  infectiousness = NULL
) {
  ...
}
```

Current output:

```r
data.frame(
  from = character(),
  to = character(),
  age_group = character(),
  rate = numeric()
)
```

For the current SIR model, include:

| Transition | From | To | Rate |
|---|---|---|---|
| Infection | `S[a]` | `I[a]` | `lambda[a] * S[a]` |
| Recovery | `I[a]` | `R[a]` | `gamma * I[a]` |

Notes:

- Ageing, births, background deaths, and infection-induced morbidity/mortality are not currently implemented.

---

## 8. Derivative construction

Implement:

```r
rates_to_derivative <- function(transition_rate_table, compartments, age_structure) {
  ...
}
```

This function should:

1. start with a zero derivative for every compartment-age cell;
2. subtract rates from origin cells;
3. add rates to destination cells;
4. return derivatives in the same compartment-major, age-group-minor ordering as the solver state vector.

This keeps the deterministic solver separated from the transition-rate logic.

---

## 9. Simulation function

Initial deterministic simulation function:

```r
simulate_deterministic <- function(
  initial_state,
  times,
  model,
  age_structure,
  contact_matrix,
  beta = 1,
  susceptibility = NULL,
  infectiousness = NULL,
  method = "euler"
) {
  ...
}
```

Expected behaviour:

- Validate inputs.
- Convert initial state to solver vector if needed.
- Use explicit Euler time steps by default. Optional `method = "deSolve"` is
  supported for the documented deterministic SIR/SEIR combinations when the
  suggested `deSolve` package is installed.
- Return output in long format, with columns:

```text
time, age_group, compartment, value
```

Optional derived outputs:

```text
time, age_group, incidence, morbidity, infection_mortality
```

Derived outputs are future work.

Implemented summary helpers for deterministic simulation output:

```r
compartment_totals(simulation_output)
age_group_totals(simulation_output)
total_population(simulation_output)
```

These helpers summarise the `time`, `compartment`, `age_group`, and `value` output produced by `simulate_deterministic()`. They do not add new model dynamics.

---

## 10. File structure

Use a lightweight R project structure.

```text
age-transmission-prototype/
  R/
    age_structure.R
    disease_model.R
    state_mapping.R
    age_transform.R
    demography.R
    contact_matrix.R
    force_of_infection.R
    transition_rates.R
    derivative.R
    simulate_deterministic.R
    simulation_summaries.R
    ...

  examples/
    mock_sir_deterministic.R
    mock_seir_demography.R

  tests/
    testthat/

  docs/
    age_structured_transmission_design.md
```

The current examples use invented mock data only. WPP-style adapters can prepare
population, fertility, mortality, and migration inputs, but WPP projection
matching remains future work.

---

## 11. Testing requirements

### 11.1 Age-structure tests

Check that:

- valid age structures pass;
- overlapping bins fail;
- unsorted bins fail;
- mismatched label lengths fail;
- duplicate age groups fail.

### 11.2 State mapping tests

Check that:

- long-format state converts to vector correctly;
- vector converts back to long format without changing values;
- total population by age is preserved;
- total population by compartment is preserved.

### 11.3 Force-of-infection tests

Use a small two-age-group example where expected values can be calculated manually.

Check that:

- output length equals number of age groups;
- values are non-negative;
- zero infectious individuals gives zero force of infection;
- increasing susceptibility increases force of infection;
- increasing infectiousness increases force of infection.

### 11.4 Transition-rate tests

For a small SIR model:

- infection rates are created for each age group;
- recovery rates are created for each age group;
- unsupported disease models are rejected.

### 11.5 Derivative tests

Check conservation properties in simplified cases:

| Case | Expected result |
|---|---|
| Infection + recovery only | Total population conserved |
| No infected individuals | No new infections |

### 11.6 Simulation smoke test

Run a small deterministic simulation and check that:

- the solver completes;
- outputs contain expected columns;
- all state values are non-negative;
- output includes all requested times;
- unsupported methods are rejected, and optional `deSolve` combinations run when
  the suggested package is available.

---

## 12. Implementation milestones for Codex

### Milestone 1: Skeleton and validation

Implemented:

- `AgeStructure`;
- validation helpers;
- state conversion helpers;
- exact age-vector aggregation helper;
- exact age-vector segregation and mixed transformation helper;
- basic tests.

Acceptance criteria:

- tests for age structure and state mapping pass;
- no disease simulation yet required.

---

### Milestone 2: Force of infection and SIR transition rates

Implemented:

- `force_of_infection()`;
- simple `DiseaseModel` constructor for SIR;
- `transition_rates()` for SIR;
- tests for force of infection and transition table.

Acceptance criteria:

- rates are correct in small manually checkable examples.

---

### Milestone 3: Deterministic derivative and simulation

Implemented:

- `rates_to_derivative()`;
- `simulate_deterministic()`;
- simulation summary helpers;
- one complete SIR example using mock demographic data.

Acceptance criteria:

- deterministic simulation runs end-to-end;
- smoke tests pass;
- output is returned in tidy long format.

---

### Milestone 4: Demographic input layer

Implemented:

- `validate_demography_table()`;
- `Demography()` for already-cleaned age-specific population tables;
- sorted storage by time and age-group order;
- exact-time population accessors.

Acceptance criteria:

- demography table validation catches missing, duplicate, extra, negative, or non-finite population rows;
- accessors return populations in age-structure order.

Still future work:

- population interpolation or nearest-year population lookup;
- WPP projection matching;
- automatic residual forcing;
- richer infection-demography policies beyond the current first-pass SIR/SEIR
  coupling.

### SEIR-demography policy checkpoint

Current SIR-demography coupling uses one explicit allocation policy: fertility
exposure is the total age-specific infection-state population `S + I + R`;
births enter the youngest `S` compartment only; background mortality removes
from `S`, `I`, and `R` independently; ageing moves `S`, `I`, and `R`
independently; and net migration is allocated entirely to `S`. Infection
transitions are added separately through the disease transition-rate pathway.

The first SEIR-demography extension preserves that policy shape. Fertility
exposure uses `S + E + I + R`; births enter the youngest `S` only; background
mortality and ageing apply independently to `S`, `E`, `I`, and `R`; and net
migration is allocated entirely to `S`.
Disease progression `E -> I` and recovery `I -> R` remain disease-model
transitions, and force of infection continues to depend on `I`, not `E`.

The `S`-only migration rule is an allocation convention for age-total net
migration inputs and residual-derived migration schedules. It is not a
mechanistic model of movement by infection status. Proportional migration across
`S`/`E`/`I`/`R` may be added later, but only behind an explicit option. This
policy checkpoint does not add WPP projection matching, disease-induced
mortality, vaccination, waning immunity, or compartment-specific demographic
rates.

---

### Milestone 5: Contact matrix integration

Implemented:

- `validate_contact_matrix()`;
- `as_agepi_contact_matrix()` for dependency-free coercion of supported inputs;
- `transform_contact_matrix()` for exact fine-to-coarse aggregation;
- compatibility with externally prepared age-aligned contact matrices.

Acceptance criteria:

- model accepts any valid age-aligned contact matrix;
- contact matrix construction remains separate from disease simulation.

Still future work:

- `MixingModel`;
- socialmixr/conmat package-specific adapters or dependencies;
- contact-matrix splitting and general rebinning;
- reciprocity correction and population balancing.

---

## 13. Coding style

Use simple, explicit R code.

Preferred style:

- small functions;
- clear argument names;
- explicit validation;
- minimal hidden global state;
- no premature S3/S4/R6 complexity;
- avoid clever metaprogramming;
- keep disease logic separate from data loading;
- keep plotting separate from simulation;
- write tests for each model component.

Keep the package structure lightweight and avoid adding dependencies until a milestone requires them.

---

## 14. Later Julia translation

The R design should map naturally onto Julia types later.

| R object/function | Julia equivalent |
|---|---|
| `age_structure` list | `struct AgeStructure` |
| `demography` list | `struct Demography` |
| `mixing_model` list | `abstract type AbstractMixingModel end` |
| `risk_model` list | `struct RiskModel` |
| `disease_model` list | `abstract type AbstractDiseaseModel end` |
| `transition_rates()` | multiple dispatch on disease model |
| future solver-backed `simulate_deterministic()` | `ODEProblem` using DifferentialEquations.jl |
| future stochastic simulation | Gillespie or tau-leaping using same transition rates |

The R prototype should therefore avoid design choices that would make later Julia translation awkward, such as relying on loosely structured global data frames inside the simulation engine.

---

## 15. First target example

The first complete runnable example should:

1. define model age bins;
2. create mock population and infection inputs;
3. create a simple age-specific contact matrix;
4. define age-specific susceptibility and infectiousness vectors;
5. define an SIR disease model;
6. seed infections in one or more age groups;
7. run the deterministic simulation;
8. return tidy long-format outputs.

The example should be small enough to run quickly and simple enough to use in tests.

---

## 16. Non-negotiable design constraints

1. Age groups must be explicit.
2. State-vector ordering must be documented.
3. Disease transitions must be represented through transition rates.
4. Force of infection must be implemented as a reusable function.
5. Simulation must be separated from data loading.
6. Mixing matrices must be separated from disease models.
7. Derived morbidity/mortality outputs must initially be separated from compartment transitions.
8. The code should be easy to port to Julia later.
