# agepi

`agepi` is an early-stage R prototype for age-structured epidemic models using projected demographic data.

The initial goal is to build a deterministic age-structured SIR model with explicit age groups, demographic turnover, age-specific mixing, susceptibility, infectiousness, morbidity, and mortality outputs.

The prototype should be designed around reusable model components rather than hard-coded disease equations, so that it can later support alternative compartment structures, stochastic simulation, and eventual porting of core simulation components to Julia.

## Initial implementation priorities

1. Define and validate age structures.
2. Define state-vector conventions and state mapping helpers.
3. Implement a reusable force-of-infection function.
4. Implement SIR transition rates.
5. Convert transition rates into deterministic derivatives.
6. Run a small deterministic simulation using mock demographic data.

## Design notes

See `docs/age_structured_transmission_design.md`.

## Current status

Project scaffold only. No production model implementation yet.
