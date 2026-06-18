# EpiDynamics validation examples

This folder contains small validation scripts that compare the public `agepi`
interface with deterministic examples from `EpiDynamics` 0.3.1. The aim is to
validate current behaviour, not to add new model structures to `agepi`.

Run scripts from the package root with `Rscript`, for example:

```sh
Rscript examples/validation/epidynamics_sir2_age_classes.R
Rscript examples/validation/epidynamics_sir.R
Rscript examples/validation/epidynamics_sir_birth_death.R
Rscript examples/validation/epidynamics_seir.R
Rscript examples/validation/finalsize_sir_final_size.R
```

The scripts load the local source tree using the same convention as the existing
examples: direct sourcing of `R/*.R` when running from the package root, or
`library(agepi)` as a fallback. They use
`EpiDynamics` only at runtime; the package is not added as a hard dependency.

## Exact validations

- `epidynamics_sir.R` compares `EpiDynamics::SIR()` with a one-age
  `agepi::SIRModel()`. With one age group and total population equal to one,
  the `agepi` frequency-dependent force of infection is the same as the
  `EpiDynamics` density-dependent proportion equations.
- `epidynamics_sir_birth_death.R` compares `EpiDynamics::SIRBirthDeath()` with
  `agepi::SIRModel()` plus one-age fertility and mortality schedules. Births
  enter susceptible people and mortality applies to all compartments in both
  models.
- `epidynamics_seir.R` compares `EpiDynamics::SEIR()` with `agepi::SEIRModel()`
  plus one-age fertility and mortality schedules.

These are non-age-structured reference validations. They are included as
baseline checks of disease dynamics and demographic coupling.

## Partial age-structured reproduction

- `epidynamics_sir2_age_classes.R` compares `EpiDynamics::sir2AgeClasses()`
  with a two-age `agepi::SIRModel()` plus demographic ageing, fertility, and
  mortality. The `EpiDynamics` object is a list with `model`, `pars`, `init`,
  `time`, and `results`; the comparable trajectory table is
  `reference$results`.

`EpiDynamics::sir2AgeClasses()` returns only `SC`, `IC`, `SA`, and `IA`.
Recovered states are omitted from the returned trajectory and from the source
equations; recovery leaves the reported four-state system. The current public
`agepi` SIR interface carries recovered children and adults explicitly, and the
demographic coupling applies ageing and mortality across disease compartments.
That difference is visible in the numerical comparison, so this script is a
worked partial reproduction rather than an exact validation.

There is also a transmission-convention mismatch: `EpiDynamics` uses terms such
as `SC * (betaCC * IC + betaCA * IA)`, while `agepi` computes force of infection
from infectious fractions. The script scales the contact-matrix columns by the
equilibrium child and adult population totals, so the `agepi` force of infection
matches the `EpiDynamics` density-dependent equations while those age totals
remain at equilibrium. This makes the comparison a transparent partial
reproduction rather than a claim that the two public interfaces have identical
semantics.

## Optional final-size comparison

- `finalsize_sir_final_size.R` compares a closed-population age-structured
  `agepi::SIRModel()` simulation with `finalsize::final_size()` when the
  optional `finalsize` package is installed. It simulates to a long horizon,
  tracks `S -> I` cumulative infections with `cumulative_flows`, and compares
  age-specific attack rates.

The script documents the compatibility assumptions in place: no births, deaths,
migration, waning, or reinfection; recipient-row/source-column contact-matrix
orientation; `finalsize` demography input normalised to population proportions;
`finalsize` contact-matrix normalisation so `contact_matrix * demography_vector`
has dominant real eigenvalue 1; and the translation from `agepi`'s `beta`/`gamma`
parameters to `finalsize`'s `r0` using the dominant eigenvalue of the fully
susceptible next-generation matrix.
`finalsize` is used only behind
`requireNamespace("finalsize", quietly = TRUE)` and is not required to run the
package.

## Investigated but skipped

- Lowercase `sir()`, `seir()`, `sirs()`, `sirDemog()`, and `seirDemog()` were
  requested for investigation but are not exported or documented under those
  names in the installed `EpiDynamics` 0.3.1 package. The corresponding
  available examples are `SIR()`, `SEIR()`, and `SIRBirthDeath()`.
- `SIR2Stages()` duplicates the same source equations as
  `sir2AgeClasses()` in this installed package, so adding another script would
  not provide an independent validation.
- `SEIR4AgeClasses()` uses four age groups with explicit annual ageing updates
  performed outside the ODE integration loop. That discrete ageing/event
  structure is not the same as the current continuous demographic coupling in
  `agepi`.
- `SIRSinusoidalForcing()`, `SIRTermTimeForcing()`, and related seasonal
  examples require time-varying or piecewise transmission forcing. Current
  public `agepi` deterministic simulation uses a static contact matrix and
  static `beta`.
- `SIRDemogStoch()`, `SISDemogStoch()`, `SIRTauLeap()`, and additive-noise
  examples are stochastic, while this validation folder is limited to
  deterministic comparisons.
- SIS and risk-group examples such as `SIS()`, `SIS2RiskGroups()`, and
  `SISnRiskGroups()` are not SIR/SEIR/SIRS-style age-structured validations and
  require disease structures outside the current built-in `agepi` SIR/SEIR
  focus. They may be candidates for future generic-compartment validation if a
  clean public-interface mapping is desired.

## Notes

All scripts print maximum absolute differences for shared state variables and
avoid writing output files or plots.
