# Docs

This folder mixes user-facing workflow notes, maintainer design notes, disease-specific notes, and status updates.

## Recommended Workflow / Reference Docs

- [model_conventions.md](model_conventions.md) - core implementation conventions for state ordering, force of infection, and contact-matrix handling.
- [epidemic_model_workflows.md](epidemic_model_workflows.md) - the main guide to SIR/SEIR, `CompartmentModel()`, deterministic and stochastic simulation, and workflow selection.
- [demography_workflows.md](demography_workflows.md) - the main guide to demographic-only and coupled demographic workflows.
- [contact_matrix_workflows.md](contact_matrix_workflows.md) - how to load, orient, validate, and adapt contact matrices.
- [generic_compartment_models.md](generic_compartment_models.md) - reference for building custom compartment models with `CompartmentModel()`.
- [external_data_adapters.md](external_data_adapters.md) - reference for WPP-style and other external input adapters.
- [demographic_residuals.md](demographic_residuals.md) - diagnostics for comparing demographic outputs and deriving residual-based migration schedules.
- [age_grid_helper_inventory.md](age_grid_helper_inventory.md) - implementation inventory for age-grid helpers and mapping layers.

## Design Notes

Files ending in `_design.md` are mainly design or maintainer notes. They are useful when extending the package, but they are not the first place new users should start.

- [age_structured_transmission_design.md](age_structured_transmission_design.md) - original age-structured transmission design note and prototype scope.
- [annual_cohort_demography_design.md](annual_cohort_demography_design.md) - design note for annual-cohort demographic updates.
- [calibration_framework_design.md](calibration_framework_design.md) - proposed calibration API and implementation plan.
- [contact_matrix_integration_design.md](contact_matrix_integration_design.md) - reconnaissance note for contact-matrix adapter work.
- [coupled_epidemic_demography_operator_splitting_design.md](coupled_epidemic_demography_operator_splitting_design.md) - design note for operator-splitting demographic coupling.
- [cumulative_flow_states_design.md](cumulative_flow_states_design.md) - design note for auxiliary cumulative flow state handling.
- [scenario_analysis_design.md](scenario_analysis_design.md) - design note for scenario orchestration around the simulators.
- [stochastic_simulation_design.md](stochastic_simulation_design.md) - design note for the stochastic simulation surface and limits.

## Disease Notes

- [disease_notes/influenza_parameterisation.md](disease_notes/influenza_parameterisation.md) - influenza parameterisation notes.
- [disease_notes/measles_parameterisation.md](disease_notes/measles_parameterisation.md) - measles parameterisation notes.
- [disease_notes/rsv_parameterisation.md](disease_notes/rsv_parameterisation.md) - RSV parameterisation notes.
- [disease_notes/tuberculosis_parameterisation.md](disease_notes/tuberculosis_parameterisation.md) - tuberculosis parameterisation notes.
- [disease_notes/kiribati_tb_modelling_inputs.md](disease_notes/kiribati_tb_modelling_inputs.md) - Kiribati TB modelling input inventory and assumptions.
- [disease_notes/kiribati_tb_data_inventory.md](disease_notes/kiribati_tb_data_inventory.md) - supporting Kiribati TB data inventory.
- [disease_notes/TEMPLATE_parameterisation.md](disease_notes/TEMPLATE_parameterisation.md) - template for adding a new disease note.

## Development / Status Notes

- [development_status.md](development_status.md) - current support, scope, limitations, and roadmap notes.
