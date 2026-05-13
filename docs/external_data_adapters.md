# External Data Adapters Design Note

## Purpose

This note records the intended boundary for future external data adapters. The
next phase can connect agepi to WPP demographic data, socialmixr contact
matrices, and conmat contact predictions without changing the deterministic SIR
simulation workflow.

External package integrations should remain optional initially. The core package
should continue to accept agepi-native objects and dependency-free inputs.

## Adapter Scope

### WPP Demography Adapter

A WPP adapter should convert external WPP-like age, year, and population tables
into a validated `Demography()` object aligned to an `AgeStructure()`.

The adapter should handle data-shape translation only:

- map external age labels or bounds to agepi age groups;
- map external year or period fields to agepi `time`;
- map population values to the required `population` column;
- call existing demography validation and construction helpers.

It should not implement demographic projection dynamics, interpolation,
fertility, mortality, births, deaths, ageing, or migration.

### socialmixr Contact Adapter

A socialmixr adapter should convert contact-matrix outputs into agepi contact
matrices using agepi's recipient-row, source-column convention.

The first adapter should focus on ordinary socialmixr-style outputs that contain
a numeric `matrix` component. It should validate matrix dimensions and age-group
labels before returning a matrix suitable for existing agepi contact-matrix
helpers.

The adapter should preserve the matrix supplied by socialmixr. It should not add
or redo reciprocity correction, population balancing, age weighting, symmetry,
or per-capita transformations inside agepi.

### conmat Contact Adapter

A conmat adapter should convert generated or predicted contact matrices into
agepi contact matrices using the same recipient-row, source-column convention.

The first adapter can support conmat-style square matrices and long prediction
tables with `age_group_from`, `age_group_to`, and `contacts` columns. Long
tables should be reshaped so `age_group_to` becomes recipient rows and
`age_group_from` becomes source columns, matching the existing agepi convention.

The adapter should treat conmat modelling choices as upstream choices. It should
not refit models, apply household-size adjustments, add reciprocity correction,
or alter generated predictions.

## Core Workflow Boundary

External data adapters should prepare inputs for the existing model workflow.
They should not change:

- `force_of_infection()`;
- `simulate_deterministic()`;
- the deterministic SIR state-vector convention;
- the current static-contact-matrix simulation assumptions.

Adapters should return objects that existing validation and simulation code can
already consume.

## Not Yet Implemented

The following remain out of scope for the initial adapter phase:

- mandatory WPP, socialmixr, or conmat dependencies;
- demographic projection dynamics;
- WPP fertility, mortality, births, deaths, ageing, or migration dynamics;
- contact-matrix reciprocity correction;
- contact-matrix population balancing;
- source-bin splitting or general contact-matrix rebinning;
- changes to `force_of_infection()` or `simulate_deterministic()`.
