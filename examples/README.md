# Examples

This folder is a grab bag of runnable scripts. The sections below group them by the kind of workflow they demonstrate.

## Basic Epidemic Workflows

- [generic_sir.R](generic_sir.R) - compares the built-in `SIRModel()` with a generic `CompartmentModel()` SIR specification; optional dependencies: none.
- [generic_seir.R](generic_seir.R) - compares the built-in `SEIRModel()` with a generic `CompartmentModel()` SEIR specification; optional dependencies: none.
- [mock_sir_deterministic.R](mock_sir_deterministic.R) - the smallest deterministic age-structured SIR example; optional dependencies: none.

## Custom Compartment Models

- [generic_msir.R](generic_msir.R) - a maternal-immunity SIR-style model built with `CompartmentModel()`; optional dependencies: none.
- [generic_outflows.R](generic_outflows.R) - a generic compartment model with explicit external outflows and matching cumulative removals; optional dependencies: none.

## Cumulative Outputs

- [deterministic_cumulative_flows.R](deterministic_cumulative_flows.R) - tracks deterministic cumulative infections and recoveries alongside the ordinary trajectory; optional dependencies: none.
- [generic_outflows.R](generic_outflows.R) - demonstrates cumulative removals from an external sink declared with `outflows`; optional dependencies: none.

## Stochastic Simulation

- [stochastic_sir.R](stochastic_sir.R) - small stochastic SIR example with a reproducible Gillespie trajectory and event log; optional dependencies: none.
- [stochastic_seir.R](stochastic_seir.R) - small stochastic SEIR example with a reproducible Gillespie trajectory and event log; optional dependencies: none.
- [stochastic_cumulative_flows.R](stochastic_cumulative_flows.R) - stochastic compartment model with event-log-derived cumulative flows; optional dependencies: none.

## Demography-Only Workflows

- [mock_demographic_workflow.R](mock_demographic_workflow.R) - invented WPP-like demographic schedules, simulation, residual diagnostics, and comparison output; optional dependencies: none.
- [demographic_ageing_diagnostics.R](demographic_ageing_diagnostics.R) - checks demographic ageing behaviour against theoretical expectations and optional WPP inputs; optional dependencies: `deSolve` for the ODE comparison, `wpp2024` for the WPP benchmark.
- [demography_plots.R](demography_plots.R) - exploratory plotting for demographic tables; optional dependencies: `ggplot2`.

## Coupled Epidemic-Demography

- [mock_seir_demography.R](mock_seir_demography.R) - deterministic SEIR with fertility, mortality, and migration coupling; optional dependencies: none.
- [annual_cohort_sir_demography.R](annual_cohort_sir_demography.R) - annual-cohort demographic coupling with a deterministic SIR model; optional dependencies: none.

## External-Data Adapters

- [epiparameter_seir.R](epiparameter_seir.R) - uses an `epiparameter` incubation-period object to parameterise SEIR progression; optional dependencies: `epiparameter`.
- [wpp_demography_validation.R](wpp_demography_validation.R) - benchmarks the WPP-adjacent demographic pipeline against WPP-style inputs; optional dependencies: `wpp2024`.
- [wpp_projection_backed_demography.R](wpp_projection_backed_demography.R) - replays a WPP population projection directly as a demographic trajectory; optional dependencies: `wpp2024`.

## Case Studies / Scaffolds

- [tb_age_structured_demography.R](tb_age_structured_demography.R) - a toy TB-style age-structured scaffold with demography and cumulative flows; optional dependencies: none.
- [kiribati_tb_realistic_demography.R](kiribati_tb_realistic_demography.R) - a public-data Kiribati TB scaffold that uses WPP-backed demography and contact-matrix adapters; optional dependencies: `wpp2024`, `contactdata`, and `socialmixr` for the contact-matrix fallback.

## Validation / Comparison Examples

- [validation/README.md](validation/README.md) - explains the validation folder and the comparison scripts it contains; optional dependencies: none.
- [validation/epidynamics_sir.R](validation/epidynamics_sir.R) - compares `agepi::SIRModel()` with `EpiDynamics::SIR()`; optional dependencies: `EpiDynamics`.
- [validation/epidynamics_seir.R](validation/epidynamics_seir.R) - compares `agepi::SEIRModel()` with `EpiDynamics::SEIR()`; optional dependencies: `EpiDynamics`.
- [validation/epidynamics_sir_birth_death.R](validation/epidynamics_sir_birth_death.R) - compares `agepi::SIRModel()` plus demography with `EpiDynamics::SIRBirthDeath()`; optional dependencies: `EpiDynamics`.
- [validation/epidynamics_sir2_age_classes.R](validation/epidynamics_sir2_age_classes.R) - partial age-structured reproduction of `EpiDynamics::sir2AgeClasses()`; optional dependencies: `EpiDynamics`.
- [validation/finalsize_sir_final_size.R](validation/finalsize_sir_final_size.R) - optional final-size comparison against `finalsize::final_size()`; optional dependencies: `finalsize`.

## Notes

- The `validation/validation_helpers.R` file is a support helper for the comparison scripts rather than a standalone example.
- Most scripts are written to run from the package root and will source the local `R/` files when `agepi` is not already loaded.
