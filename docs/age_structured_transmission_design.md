# Age-Structured Transmission Prototype: Design Document

## 1. Project aim

Build a small but extensible R prototype for age-structured infectious disease transmission models using projected demographic data, initially from the World Population Prospects 2024 data ecosystem.

The first implementation should be deliberately simple: a deterministic age-structured SIR or SEIR model for one country, one demographic scenario, and one contact matrix. However, the code should be designed so that later extensions can support:

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

| Component | Initial implementation |
|---|---|
| Language | R |
| Disease model | Age-structured SIR, with design compatible with SEIR |
| Simulation | Deterministic ODE |
| Geography | One country |
| Demography | WPP 2024 population projections, initially one projection scenario |
| Age bins | User-defined model age bins |
| Mixing | Static age-specific contact matrix |
| Susceptibility | Age-specific vector |
| Infectiousness | Age-specific vector, default all ones |
| Infection morbidity/mortality | Derived outputs, not yet disease-state transitions |
| Births | Births enter youngest susceptible age group |
| Background mortality | Age-specific demographic mortality |
| Ageing | Movement between age groups |

### Not included initially

| Component | Reason deferred |
|---|---|
| Stochastic simulation | Design for it, but do not implement initially |
| Vaccination | Requires additional compartment conventions |
| Waning immunity | Can be added once transition framework is stable |
| Calibration/inference | Depends on stable model/simulation interface |
| Multi-country batching | Should come after one-country workflow is correct |
| Full R package structure | Start as a prototype R project; convert later if needed |
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
ODE simulation
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

Stores time-varying demographic information aligned to the model age bins.

Suggested fields:

```r
demography <- list(
  age_structure = age_structure,
  times = times,
  population = population_matrix,
  mortality = mortality_matrix,
  births = births_vector,
  ageing_rates = ageing_rate_vector
)
```

Expected dimensions:

| Field | Shape |
|---|---|
| `population` | `length(times) × number_of_age_groups` |
| `mortality` | `length(times) × number_of_age_groups` |
| `births` | `length(times)` |
| `ageing_rates` | `number_of_age_groups` |

Notes:

- `population[t, a]` is the projected population in age group `a` at time `t`.
- `mortality[t, a]` is the background mortality rate for age group `a` at time `t`.
- `births[t]` is the number/rate of births entering the youngest age group.
- `ageing_rates[a]` controls ageing out of age group `a`.

The first version may use interpolation or nearest-year lookup for time-varying demographic quantities.

---

### 5.3 `MixingModel`

Stores or generates the age-specific contact matrix.

Initial version:

```r
mixing_model <- list(
  age_structure = age_structure,
  contact_matrix = C,
  time_varying = FALSE
)
```

Expected dimensions:

```text
number_of_age_groups × number_of_age_groups
```

The entry `C[a, b]` should represent contacts made by individuals in recipient age group `a` with individuals in source age group `b`.

Required checks:

- Contact matrix must be square.
- Matrix dimensions must match the number of model age groups.
- Entries must be non-negative.
- Missing values are not allowed.

Later versions may allow:

```r
get_contact_matrix <- function(t, mixing_model, demography) {
  ...
}
```

---

### 5.4 `RiskModel`

Stores age-specific susceptibility, infectiousness, and severity parameters.

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

### 5.6 `SimulationProblem`

Combines all components needed to run a simulation.

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

Required helper functions:

```r
state_long_to_vector(state_long, age_structure, compartments)
state_vector_to_matrix(state_vector, age_structure, compartments)
state_vector_to_long(state_vector, age_structure, compartments)
```

The matrix representation should have:

```text
rows = compartments
columns = age groups
```

or vice versa, but the choice must be documented and used consistently.

Recommended:

```r
state_matrix[compartment, age_group]
```

---

## 7. Transition-rate interface

The central model function should return a transition table.

Suggested function signature:

```r
transition_rates <- function(t, state, problem) {
  ...
}
```

Suggested output:

```r
data.frame(
  transition = character(),
  age_from = character(),
  age_to = character(),
  compartment_from = character(),
  compartment_to = character(),
  rate = numeric()
)
```

For the initial SIR model, include:

| Transition | From | To | Rate |
|---|---|---|---|
| Infection | `S[a]` | `I[a]` | `lambda[a] * S[a]` |
| Recovery | `I[a]` | `R[a]` | `gamma * I[a]` |
| Ageing | `X[a]` | `X[a+1]` | `ageing_rate[a] * X[a]` |
| Background death | `X[a]` | outside system | `mortality[a] * X[a]` |
| Birth | outside system | `S[1]` | `births[t]` |

Notes:

- For births and deaths, `compartment_from` or `compartment_to` may be `NA`.
- For ageing, transitions should apply to all compartments.
- The final age group should not age into a further group.
- Infection-induced morbidity/mortality should initially be calculated as derived outputs from infections, not as state transitions.

---

## 8. Derivative construction

Implement:

```r
rates_to_derivative <- function(rates, state, problem) {
  ...
}
```

This function should:

1. start with a zero derivative for every compartment-age cell;
2. subtract rates from origin cells;
3. add rates to destination cells;
4. ignore missing origin or destination cells for births/deaths;
5. return a numeric vector in the same ordering as the solver state vector.

Then implement:

```r
sir_derivatives <- function(t, state, parameters) {
  problem <- parameters$problem
  rates <- transition_rates(t, state, problem)
  rates_to_derivative(rates, state, problem)
}
```

This keeps the deterministic solver separated from the transition-rate logic.

---

## 9. Simulation function

Initial deterministic simulation function:

```r
simulate_deterministic <- function(problem) {
  ...
}
```

Expected behaviour:

- Validate the problem.
- Convert initial state to solver vector if needed.
- Use `deSolve::ode()` or equivalent.
- Return output in long format, with columns:

```text
time, age_group, compartment, value
```

Optional derived outputs:

```text
time, age_group, incidence, morbidity, infection_mortality
```

---

## 10. File structure

Use a lightweight R project structure.

```text
age-transmission-prototype/
  R/
    age_structure.R
    demography.R
    mixing.R
    risk_model.R
    disease_model.R
    state.R
    force_of_infection.R
    transition_rates.R
    derivatives.R
    simulate_deterministic.R
    outputs.R
    validate.R

  scripts/
    01_define_age_structure.R
    02_load_or_mock_demography.R
    03_define_mixing_model.R
    04_run_sir_example.R

  tests/
    test_age_structure.R
    test_state_mapping.R
    test_force_of_infection.R
    test_transition_rates.R
    test_derivatives.R
    test_simulation_smoke.R

  docs/
    model_design.md
    age_mapping.md

  outputs/
```

For the very first implementation, mock demographic data may be used before connecting to WPP 2024. However, the code should be written so that WPP-derived demographic inputs can replace mock inputs without changing the simulation code.

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
- ageing transitions are created for all but the final age group;
- births enter the youngest susceptible group;
- background deaths remove people from each compartment and age group.

### 11.5 Derivative tests

Check conservation properties in simplified cases:

| Case | Expected result |
|---|---|
| Infection + recovery only | Total population conserved |
| Add background mortality | Total population decreases |
| Add births only | Total population increases |
| Ageing only | Total population conserved |
| No infected individuals | No new infections |

### 11.6 Simulation smoke test

Run a small deterministic simulation and check that:

- the solver completes;
- outputs contain expected columns;
- all state values are non-negative or close to non-negative within numerical tolerance;
- output includes all requested times;
- total population behaves sensibly under births/deaths.

---

## 12. Implementation milestones for Codex

### Milestone 1: Skeleton and validation

Implement:

- `AgeStructure`;
- validation helpers;
- state conversion helpers;
- basic tests.

Acceptance criteria:

- tests for age structure and state mapping pass;
- no disease simulation yet required.

---

### Milestone 2: Force of infection and SIR transition rates

Implement:

- `force_of_infection()`;
- simple `DiseaseModel` constructor for SIR;
- `transition_rates()` for SIR;
- tests for force of infection and transition table.

Acceptance criteria:

- rates are correct in small manually checkable examples.

---

### Milestone 3: Deterministic derivative and simulation

Implement:

- `rates_to_derivative()`;
- `simulate_deterministic()`;
- one complete SIR example using mock demographic data.

Acceptance criteria:

- deterministic simulation runs end-to-end;
- smoke tests pass;
- output is returned in tidy long format.

---

### Milestone 4: Demographic input layer

Implement:

- functions to construct a `Demography` object from already-cleaned population, mortality, and births data;
- interpolation or lookup of demographic quantities at time `t`;
- basic mock WPP-style data loader.

Acceptance criteria:

- simulation code does not care whether demography comes from mock data or WPP-derived data.

---

### Milestone 5: Contact matrix integration

Implement:

- `MixingModel`;
- contact matrix validation;
- placeholder support for matrices created externally by `socialmixr` or `conmat`.

Acceptance criteria:

- model accepts any valid age-aligned contact matrix;
- contact matrix construction remains separate from disease simulation.

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

Do not build a full R package yet unless requested. A structured R project with source files and tests is sufficient.

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
| `simulate_deterministic()` | `ODEProblem` using DifferentialEquations.jl |
| future stochastic simulation | Gillespie or tau-leaping using same transition rates |

The R prototype should therefore avoid design choices that would make later Julia translation awkward, such as relying on loosely structured global data frames inside the simulation engine.

---

## 15. First target example

The first complete runnable example should:

1. define model age bins;
2. create mock demographic projections for one country;
3. create a simple age-specific contact matrix;
4. define age-specific susceptibility and infectiousness vectors;
5. define an SIR disease model;
6. seed infections in one or more age groups;
7. run the deterministic simulation;
8. return tidy long-format outputs;
9. calculate derived incidence and infection mortality outputs.

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
