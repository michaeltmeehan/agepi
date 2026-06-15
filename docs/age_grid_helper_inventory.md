# Age Grid Helper Inventory

This note inventories existing age-grid expansion, aggregation, and mapping
helpers relevant to a future annual-cohort demography mode. It is an
implementation gap analysis only.

## Scope Inspected

- `R/demographic_age_grid.R`
- `R/age_structure.R`
- `R/age_transform.R`
- `R/state_mapping.R`
- `R/contact_matrix.R`
- `R/force_of_infection.R`, limited to contact-matrix validation
- `R/demography.R`
- `R/demographic_rates.R`, limited to schedule validation/storage
- `R/demographic_process_builder.R`
- `R/simulation_summaries.R`
- `R/wpp_demographic_inputs.R`
- `docs/annual_cohort_demography_design.md`
- `docs/model_conventions.md`
- Related tests: `test-demographic-age-grid.R`, `test-age-transform.R`,
  `test-contact-matrix.R`, `test-contact-matrix-transform.R`,
  `test-state-mapping.R`, `test-simulation-summaries.R`,
  `test-demography.R`, `test-wpp-demographic-inputs.R`,
  `test-demographic-rates.R`, `test-demographic-process-builder.R`,
  and narrowly relevant demography validation/process tests.

`R/age_structured_transmission.R` was requested but is not present in this
checkout.

## Existing Helpers

| Helper | Location | Exported | Tests found | What it does | Annual-cohort suitability |
| --- | --- | --- | --- | --- | --- |
| `AgeStructure()` | `R/age_structure.R` | Yes | Yes | Builds an age grid from labels and inclusive lower/upper bounds. Allows a final open-ended bin. | Reuse as the core representation for internal and reporting grids. |
| `validate_age_structure()` | `R/age_structure.R` | Yes | Yes | Validates required fields, unique labels, sorted non-overlapping bounds, finite lower bounds, and only-final open-ended upper bound. | Reuse; add cross-grid mapping validation on top. |
| `wpp_age_structure_1year()` | `R/demographic_age_grid.R` | Yes | Yes | Creates `0`, `1`, ..., `max_age+` WPP-style single-year grid. | Directly reusable for the internal annual-cohort grid. |
| `wpp_age_structure_5year()` | `R/demographic_age_grid.R` | Yes | Yes | Creates `0-4`, `5-9`, ..., `max_age+` WPP-style five-year grid. | Directly reusable for WPP/reporting grid. |
| `parse_wpp_age_labels()` and lower-bound parsers | `R/wpp_demographic_inputs.R` | No | Indirect | Maps labels like `0`, `5`, `10+`, `0-4` to an existing `AgeStructure` by lower bound. | Reuse internally, but likely promote or generalize for mapping diagnostics. |
| `age_bin_widths()` | `R/wpp_demographic_inputs.R` | No | Indirect | Computes finite inclusive bin widths; rejects open-ended widths. | Useful for fertility conversion; not enough for general open-ended expansion. |
| `aggregate_age_vector()` | `R/age_transform.R` | Yes | Yes | Sums a vector from fine to coarse when every target bin is an exact union of source bins. | Reuse for 1-year to 5-year population/rate-total aggregation. |
| `segregate_age_vector()` | `R/age_transform.R` | Yes | Yes | Splits coarse vector values to nested finer bins using explicit target weights. | Reuse for count-like WPP population expansion if weights are supplied. |
| `transform_age_vector()` | `R/age_transform.R` | Yes | Yes | Wrapper for identity, aggregation, segregation, and mixed exact transforms; supports width/equal/weight splits, except width splits of open-ended bins. | Reuse for generic age-vector transforms; needs rate-specific wrappers and policy naming. |
| `AgeGridMapping()`, `aggregate_age_counts()`, `expand_age_counts()` | `R/age_grid_mapping.R` | No | Yes | Internal mapping object plus count/population table helpers for exact nested age grids. Aggregation sums finer bins into coarser bins; expansion currently splits uniformly across nested target bins. | First minimal mapping/count layer for future 1-year internal grids and aggregated reporting groups. Not a rate helper and not annual-cohort simulation. |
| `expand_age_rates()`, `expand_age_hazards()`, `expand_age_interval_probabilities()` | `R/age_grid_mapping.R` | No | Yes | Internal rate-expansion helpers for exact nested age grids. Annual per-person rates and hazards are copied to finer bins; interval probabilities are converted by survival composition over finite split widths. | Reusable for future fertility, mortality, and migration schedule expansion. Not annual-cohort simulation. |
| `aggregate_demography_trajectory_age_grid()`, `aggregate_epidemic_trajectory_age_grid()`, `aggregate_population_summary_age_grid()`, `aggregate_cumulative_flows_age_grid()` | `R/age_grid_mapping.R` | No | Yes | Internal output wrappers around `aggregate_age_counts()` for count-like demographic, epidemic, age-summary, and cumulative-flow tables. They validate each output shape, preserve non-age grouping columns such as `time`, `compartment`, `scenario`, `transition_id`, `cumulative_name`, `from`, and `to`, and sum internal age rows to reporting age groups. | Reusable after future 1-year internal simulations to produce reporting-grid outputs. Not for rates, hazards, proportions, percentages, or simulator behavior changes. |
| `expand_fertility_schedule_age_grid()`, `expand_mortality_schedule_age_grid()`, `expand_migration_schedule_age_grid()` | `R/demographic_rates.R` | No | Yes | Internal schedule wrappers that expand stored schedule data to a finer nested age grid, then rebuild validated schedule objects. Fertility rates are copied as annual births per female person-year; mortality rates are copied as annual hazards; migration rates are copied and migration counts are uniformly split across nested target ages, allowing negative net counts while conserving each source row total. | Reusable as a preparation step for future annual-cohort inputs. Not wired into any simulator path. |
| Internal transform validators (`can_aggregate_age_structures()`, `can_segregate_age_structures()`, overlap helpers) | `R/age_transform.R` | No | Indirect | Check exact unions, nested bins, overlap coverage, and split proportions. | Strong base for mapping validation, but not exported or packaged as a mapping object. |
| `transform_contact_matrix()` | `R/age_transform.R` | Yes | Yes | Aggregates contact matrices from fine to coarse exact-union grids using recipient-population weighting. | Reuse for internal-to-reporting contact aggregation; does not expand 5-year to 1-year. |
| `validate_contact_matrix()` | `R/force_of_infection.R` | Yes | Yes | Validates numeric finite non-negative square matrices and optional age-grid dimensions. | Reuse after any contact expansion/aggregation. |
| `as_agepi_contact_matrix()` | `R/contact_matrix.R` | Yes | Yes | Coerces matrix/data-frame/socialmixr-like/conmat-like inputs into recipient-source orientation and validates labels. | Reuse for input standardisation before grid transforms. |
| `ContactSchedule()` / `contact_matrix_at()` | `R/contact_matrix.R` | Yes | Yes | Stores exact-time contact matrices on one age grid; no interpolation. | Reuse after expansion to internal grid; no multi-grid schedule support. |
| `validate_demography_table()` | `R/demography.R` | Yes | Yes | Validates complete `time` x `age_group` population table on one grid. | Reuse for internal and reporting trajectories separately. |
| `Demography()` and accessors | `R/demography.R` | Yes | Yes | Stores sorted age-specific population trajectories and exact-time accessors. | Reuse for trajectory storage, but it does not transform between grids. |
| `population_from_wpp()` / `demography_from_wpp()` | `R/demography.R` | Yes | Yes | WPP-like population adapter; reshapes labels into a target age grid. No expansion, interpolation, or scaling. | Reuse as ingestion only; needs expansion step before 1-year mode. |
| `FertilitySchedule()`, `MortalitySchedule()`, `MigrationSchedule()` and validators | `R/demographic_rates.R` | Yes | Yes | Store age-specific schedules on one grid; mortality/migration require full coverage, fertility can be partial. | Reuse after rates are expanded to internal ages. |
| `standardise_wpp_fertility()`, `fertility_from_wpp_percent_asfr()` | `R/wpp_demographic_inputs.R` | Yes | Yes | Convert WPP-like fertility inputs into schedule objects; percent ASFR uses finite bin widths. | Reuse for 5-year schedules; not sufficient for 1-year age-specific rates without expansion. |
| `standardise_wpp_mortality()`, `mortality_from_wpp_mx()`, `mortality_from_wpp()` | `R/wpp_demographic_inputs.R` | Yes | Yes | Convert WPP-like mortality rates to `MortalitySchedule`; only `mx` supported. | Reuse for ingestion; needs 5-year to 1-year mortality policy. |
| `standardise_wpp_migration()` | `R/wpp_demographic_inputs.R` | Yes | Yes | Converts WPP-like migration counts/rates into `MigrationSchedule`. | Reuse for ingestion; expansion/allocation policy still needed. |
| `state_long_to_vector()` / `state_vector_to_long()` | `R/state_mapping.R` | Yes | Yes | Convert between long compartment-age rows and compartment-major vectors. | Reuse for coupled annual operator state reshaping. |
| `compartment_totals()`, `age_group_totals()`, `total_population()` | `R/simulation_summaries.R` | Yes | Yes | Sum deterministic output over existing age groups and/or compartments. | Reuse for totals only; does not aggregate internal ages to reporting groups. |
| `build_demographic_process()` | `R/demographic_process_builder.R` | Yes | Yes | Assembles age grid, ageing operator, and schedules into a derivative process. | Mostly not reusable for annual cohort except as a pattern for schedule consistency checks. |

## Coverage Summary

Existing tests cover the key current helpers: WPP age-grid constructors,
`AgeStructure` validation through many callers, age-vector transforms,
contact-matrix coercion and aggregation, state vector mapping, demography table
validation/accessors, WPP population/rate adapters, and simulation summary
totals. The tests are focused on current exact-grid behavior. They do not cover
annual-cohort requirements such as 5-year-to-1-year WPP expansion, rate
disaggregation semantics, contact-matrix expansion, or reporting-grid
aggregation of compartment outputs.

## Major Gaps

- A minimal internal `AgeGridMapping()` now records source/target grids, nested
  index relationships, supported aggregation/expansion directions, and
  open-ended policy for count/population tables and rate-like quantities. It
  does not yet cover contact matrices or annual-cohort operators.
- No public validator for "reporting bins are exact unions of internal bins"
  beyond internal transform helpers.
- No helper to expand 5-year WPP population counts to 1-year internal ages with
  documented weights, constraints, and open-ended handling.
- Internal helpers now expand annual per-person rates, annual hazards,
  interval probabilities, and fertility, mortality, and migration schedule
  objects to finer nested age grids. These wrappers only prepare schedule
  objects; they do not alter simulator behavior.
- Internal wrappers now aggregate 1-year internal count-like output tables back
  to reporting bins. They cover standalone demography trajectories,
  deterministic epidemic trajectories, age-stratified population summaries, and
  cumulative-flow output tables, but do not yet build complete `Demography`
  objects.
- No helper to expand coarse contact matrices to a 1-year contact matrix.
  Current `transform_contact_matrix()` only aggregates fine to coarse exact
  unions.
- Internal compartment-output aggregation now preserves compartments and other
  metadata while summing nested internal ages to reporting bins. Public summary
  helpers still only sum over labels already present in the output.
- Open-ended final groups are validated and aggregation-compatible, but
  width-based splitting of open-ended bins is intentionally rejected. Annual
  cohort mode needs an explicit terminal-age policy rather than implicit width
  splitting.

## Naming And Design Inconsistencies

- Age transforms use both "aggregate" and "segregate"; demography design text
  uses "expand". Consider standardising user-facing names around
  `aggregate_*()` and `expand_*()` while keeping `segregate_age_vector()` as a
  lower-level or backwards-compatible term.
- WPP label parsing is internal to `wpp_demographic_inputs.R`, but population
  ingestion also depends on it from `R/demography.R`. A general age-label parser
  or documented internal location would reduce coupling.
- `transform_age_vector()` can split count-like vectors by width/equal/weights,
  but rate expansion needs different semantics. Avoid reusing it directly for
  hazards or fertility rates without a rate-specific wrapper.
- Contact-matrix transformation is named generically but currently supports
  only fine-to-coarse aggregation. A future `aggregate_contact_matrix()` /
  `expand_contact_matrix()` split would make assumptions clearer.

## Recommended Minimal Helper Roadmap

1. Extend or specialise the internal mapping/count helpers if non-uniform
   expansion weights are needed for WPP population inputs.
2. Add contact-matrix expansion from coarse to internal ages with explicit
   assumptions, then reuse existing `transform_contact_matrix()` for
   internal-to-reporting aggregation.
3. Add focused tests for mapping validation, open-ended terminal handling,
   population/rate expansion, contact expansion, and compartment-output
   reporting aggregation before wiring annual-cohort simulation behavior.
