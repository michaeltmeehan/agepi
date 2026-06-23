# Contact Matrix Workflows

`agepi` treats contact matrices as one of the core inputs to age-structured
transmission. This guide explains what the matrix means, how the package
orients it, how source matrices are loaded and adapted, and how to decide which
workflow fits a given model.

The short version is:

- user-supplied matrices can be coerced and validated with
  `as_agepi_contact_matrix()`;
- published or synthetic sources can be loaded with
  `load_contact_matrix_source()`;
- source matrices are adapted to the model age grid with
  `adapt_contact_matrix_to_age_structure()`;
- `ContactSchedule()` stores externally supplied time-indexed matrices, but is
  not yet consumed directly by `simulate_deterministic()`;
- `contact_matrix_for_age_structure()` remains only as a deprecated
  compatibility wrapper.

## 1. Overview

A contact matrix represents how often people in one age group contact people in
another age group. In age-structured transmission models, it is the mechanism
that turns infectious people in source age groups into exposure risk for
recipient age groups.

`agepi` uses a fixed convention:

- rows are recipient age groups;
- columns are source or infectious age groups;
- `contact_matrix[a, b]` gives contacts made by recipient group `a` with source
  group `b`.

That orientation matters. If rows and columns are swapped, the force of
infection is attached to the wrong age groups and the simulation can still run
while giving the wrong answer.

## 2. Core Mental Model

The contact matrix is only one part of the transmission calculation. In the
current force-of-infection convention:

- recipient age group `a` is the group whose risk we are calculating;
- source age group `b` is the age group contributing infectious contacts;
- susceptibility scales the recipient side;
- infectiousness scales the source side;
- `I_b / N_b` is the source-age infectious fraction used in the denominator;
- `beta` scales the overall transmission intensity.

The practical consequence is that a contact matrix is not a free-floating
table. It is coupled to:

- the model age structure;
- the infectious counts by source age;
- the population by source age;
- susceptibility and infectiousness vectors when those are used.

If the age grid or the orientation is wrong, the force of infection is wrong
even if the matrix looks numerically plausible.

## 3. Age Structures And Contact Matrices

Contact matrices must line up with the model `AgeStructure()`. The model age
grid is the one used in `simulate_deterministic()` and the disease state
vector.

There are three common grids to keep distinct:

- the model age grid used in simulation;
- the source or native age grid used by a published contact matrix;
- the reporting age grid used for presentation or aggregation.

WPP-based demography often uses a one-year or five-year age grid, while a
contact source may live on five-year bands, broad survey bands, or a different
synthetic grid. The current package design adapts the contact source to the
model age grid, not the other way around. That keeps provenance explicit and
avoids silently re-creating a source matrix from a model grid.

If you need to aggregate epidemic or demographic outputs after simulation, use
the age-grid tools described in `docs/demography_workflows.md`. The contact
matrix workflow is separate from those output-aggregation workflows.

Open-ended age groups are allowed, but they must still be compatible with the
source and target age grids used in the adaptation step.

## 4. Contact-Matrix Layers In agepi

The package has a layered design for contact matrices:

1. `as_agepi_contact_matrix()` coerces and validates user-supplied matrix-like
   inputs into the agepi recipient-source convention.
2. `contact_matrix_from_socialmixr()` is a dependency-free adapter for
   `socialmixr`-shaped outputs.
3. `contact_matrix_from_conmat()` is a dependency-free adapter for conmat-style
   long tables.
4. `load_contact_matrix_source()` loads a native/source contact matrix with
   provenance metadata.
5. `validate_contact_matrix_source()` checks that a source object has the
   expected structure.
6. `adapt_contact_matrix_to_age_structure()` adapts a loaded source matrix to
   a target `AgeStructure()`.
7. `ContactSchedule()` and `contact_matrix_at()` store and retrieve
   time-indexed externally supplied matrices.
8. `contact_matrix_for_age_structure()` remains only for backward
   compatibility.

This separation is deliberate:

- user-supplied matrices are about coercion and validation;
- source loaders are about provenance and optional external packages;
- adaptation is about changing age grids.

## 5. Source Objects

`load_contact_matrix_source()` returns an object of class
`agepi_contact_matrix_source`. The intended structure is:

- `matrix`
- `age_structure`
- `source`
- `country`, if applicable
- `setting`, if applicable
- `orientation`
- `convention`
- `source_reference`
- `notes`
- `limitations`
- `metadata`
- class `agepi_contact_matrix_source`

The provenance fields are there because the same numeric matrix can mean very
different things depending on how it was generated:

- empirical survey data versus synthetic estimates;
- country-specific source versus proxy country;
- a broad `"all"` matrix versus a setting-specific matrix;
- a matrix that has already been filtered or symmetrised upstream.

When a matrix is used as a proxy, the limitation belongs in the source object
and in the narrative around the model. That is especially important when the
matrix is not country-specific or not directly observed in the place being
modelled.

`validate_contact_matrix_source()` is useful for custom source objects or for
asserting structure in tests. Loader functions already return validated source
objects.

## 6. Loading Published Or Synthetic Sources

`load_contact_matrix_source()` currently supports four source values:

- `"polymod_uk"`: compatibility alias for the United Kingdom POLYMOD source;
- `"polymod"`: empirical POLYMOD survey matrices via optional `socialmixr`;
- `"prem"`: Prem et al. synthetic country matrices via optional `contactdata`;
- `"conmat"`: matrices generated by optional `conmat` from caller-supplied
  population data.

### POLYMOD

POLYMOD is empirical European social-contact survey data. It is useful as a
survey-based contact source, but it is not country-specific to the country you
are modelling unless that country is one of the surveyed POLYMOD settings.
`"polymod_uk"` is only a compatibility alias for the United Kingdom case.

Optional dependency:

- `socialmixr`

Key arguments:

- `country`
- `setting`
- `age_limits`
- `filter`
- `symmetric`

Limitations:

- European and pre-pandemic survey context;
- may require explicit filter arguments for some settings;
- not automatically calibrated to non-survey populations.

### Prem / contactdata

Prem matrices are synthetic country-level contact matrices distributed through
`contactdata`. They are the preferred source when the requested country is
available in the installed dataset.

Optional dependency:

- `contactdata`

Key arguments:

- `country`
- `setting`
- `geographic_setting`
- `data_source`

Limitations:

- synthetic rather than directly observed;
- availability depends on the installed `contactdata` dataset;
- Kiribati may be absent from the installed dataset;
- nearby Pacific matrices should only be used as explicit proxy assumptions.

### conmat

`conmat` support is a generation workflow, not a country lookup. It uses
caller-supplied population data to generate a matrix.

Optional dependency:

- `conmat`

Key arguments:

- `population`
- `setting`
- `age_limits`
- `per_capita_household_size`

Limitations:

- depends on the supplied population and modelling defaults;
- not a simple published-country fetch;
- agepi does not refit or reinterpret the upstream conmat model.

When a source is unavailable because the optional package is missing, the
loader fails with a clear install/use message.

## 7. Adapting Source Matrices To A Model Age Grid

`adapt_contact_matrix_to_age_structure()` is the function that bridges a
native/source age grid and the model age grid.

The preferred method is `method = "source_band"`. The older
`method = "exact"` spelling is retained only as a deprecated alias.

The current source-band assumption is:

- the source matrix is loaded on its native age grid;
- the target model age grid is either nested inside or nests the source grid;
- fine-to-coarse aggregation uses recipient-population weighting;
- coarse-to-fine expansion treats contact rates as constant within each source
  age band.

This is a practical modelling assumption, not a claim that fine-age contact
structure is literally known. It lets agepi keep provenance and adaptation
separate while still producing a matrix on the age grid used by the disease
model.

### Population requirement

Population is required whenever fine-to-coarse aggregation is performed. The
adapter errors if `population` is missing. There is no silent equal-weight
fallback, because that would hide an assumption about how contacts should be
averaged.

If the source matrix already matches the target age grid exactly, the matrix is
returned unchanged apart from validation and metadata attachment.

### What the adapter does not do

The adapter does not:

- interpolate contacts across age bands;
- calibrate contacts to transmission data;
- balance reciprocity after import;
- re-fit source models;
- split source bins into arbitrary finer bins.

`transform_contact_matrix()` is still available for exact fine-to-coarse
aggregation, but it has narrower behaviour than the source-layer workflow. Use
it when you already have a plain matrix and a compatible exact aggregation
problem.

## 8. User-Supplied Matrices

If you already have a matrix, use `as_agepi_contact_matrix()` to coerce it to
the agepi convention and validate it against the chosen age structure.

Supported inputs include:

- a numeric matrix;
- a data frame coercible to a numeric matrix;
- a socialmixr-like list with a numeric `matrix` element;
- a conmat-style long data frame with `age_group_from`, `age_group_to`, and
  `contacts`.

Orientation controls matter:

- `orientation = "recipient_source"` means the input is already in agepi
  orientation;
- `orientation = "source_recipient"` means rows and columns should be
  transposed before validation;
- `transpose = TRUE` applies one additional transposition after orientation
  handling.

If an `age_structure` is supplied, the matrix is also checked against the age
labels and dimensions in that structure.

## 9. Time-Indexed Contact Matrices

`ContactSchedule()` stores externally supplied matrices indexed by time. The
input can be:

- a named list of matrices;
- a long data frame with `time`, `age_group_from`, `age_group_to`, and
  `contacts`.

`contact_matrix_at()` then retrieves the matrix at an exact available time.
There is no interpolation, smoothing, or automatic time bridging.

Important limitation: `ContactSchedule()` is a storage/accessor layer, not a
simulation input layer. It is not yet consumed directly by
`simulate_deterministic()`. If you want to use a schedule in simulation, you
currently need to extract the matrix for the relevant time and pass that
numeric matrix to the simulator.

## 10. Using Contact Matrices In Simulations

`simulate_deterministic()` expects a numeric `contact_matrix` aligned to the
model age structure. It does not accept a source object directly.

The practical simulation workflow is:

1. define the model `AgeStructure()`;
2. load or construct a contact matrix;
3. adapt it to the model age grid if needed;
4. pass the resulting numeric matrix to `simulate_deterministic()`.

The contact matrix then interacts with:

- `beta`, which scales overall transmission;
- susceptibility, which adjusts the recipient side;
- infectiousness, which adjusts the source side;
- the current age-specific population, which enters the force-of-infection
  denominator.

This is the same machinery used in simple SIR/SEIR examples and in the Kiribati
TB scaffold. The disease model may be simple, but the contact matrix still has
to match the model age grid and orientation exactly.

## 11. Choosing The Right Workflow

| Workflow | Main functions | Input type | Use case | Strengths | Limitations | Example file |
| --- | --- | --- | --- | --- | --- | --- |
| Simple manual matrix | `as_agepi_contact_matrix()` | Numeric matrix or data frame | Fastest path for toy models | No optional dependencies | User must manage provenance | `examples/mock_sir_deterministic.R` |
| socialmixr adapter | `contact_matrix_from_socialmixr()` | socialmixr-like list or matrix | Reuse a socialmixr output | Dependency-free coercion | Does not load source provenance | `examples/mock_sir_deterministic.R` |
| conmat-style long table adapter | `contact_matrix_from_conmat()` | Long table | Reuse a prediction table | Explicit orientation handling | No source provenance beyond the table | `examples/mock_sir_deterministic.R` |
| POLYMOD source loading | `load_contact_matrix_source()` | Optional `socialmixr` | Use POLYMOD as a source layer | Provenance attached | European survey context | `examples/generic_sir.R` |
| Prem/contactdata source loading | `load_contact_matrix_source()` | Optional `contactdata` | Use a country-specific synthetic source | Preferred when country exists | Kiribati may be absent | `examples/kiribati_tb_realistic_demography.R` |
| conmat generation | `load_contact_matrix_source()` | Optional `conmat` + population | Generate a matrix from population data | Explicit upstream generation | Requires caller-supplied population | `examples/kiribati_tb_realistic_demography.R` |
| Source adaptation | `adapt_contact_matrix_to_age_structure()` | Source object | Map a source grid to the model grid | Keeps provenance and adaptation separate | Requires population for aggregation | `examples/kiribati_tb_realistic_demography.R` |
| Time-indexed schedule | `ContactSchedule()`, `contact_matrix_at()` | Named list or long table | Store matrices by time | Exact lookup only | Not simulator-integrated yet | `examples/mock_seir_demography.R` |

## 12. Worked Examples

### Example A: Dependency-Free User-Supplied Matrix

This is the simplest path. It uses a plain matrix and makes the orientation
explicit.

```r
library(agepi)

ages <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

contacts <- matrix(
  c(4, 2,
    1, 5),
  nrow = 2,
  byrow = TRUE
)

contact_matrix <- as_agepi_contact_matrix(
  contacts,
  age_structure = ages,
  orientation = "recipient_source"
)

contact_matrix
```

This matrix now has rows as recipient age groups and columns as source age
groups, which is the orientation used by `force_of_infection()`.

### Example B: Source/Adapt Workflow

This example shows the intended pattern for a published source matrix.

```r
if (requireNamespace("contactdata", quietly = TRUE)) {
  countries <- contactdata::list_countries()
  if (length(countries) > 0) {
    country <- countries[[1]]

    source_contacts <- load_contact_matrix_source(
      source = "prem",
      country = country,
      setting = "all"
    )

    target_ages <- AgeStructure(
      age_groups = c("0-4", "5-9", "10+"),
      lower_bounds = c(0, 5, 10),
      upper_bounds = c(4, 9, Inf)
    )

    contact_matrix <- adapt_contact_matrix_to_age_structure(
      source_contacts,
      target_ages,
      method = "source_band"
    )
  }
}
```

If the source age grid is finer than the target model grid, supply
`population = ...` so the aggregation step can use recipient-population
weighting.

### Example C: Time-Indexed Contacts

```r
ages <- AgeStructure(
  age_groups = c("0-4", "5-9"),
  lower_bounds = c(0, 5),
  upper_bounds = c(4, 9)
)

schedule <- ContactSchedule(
  list(
    "0" = matrix(c(4, 2, 1, 5), nrow = 2, byrow = TRUE),
    "1" = matrix(c(3, 2, 2, 4), nrow = 2, byrow = TRUE)
  ),
  ages
)

contact_matrix_at(schedule, time = 1)
```

This returns the matrix at time 1 exactly. There is no interpolation.

## 13. Kiribati TB Scaffold

`examples/kiribati_tb_realistic_demography.R` uses the contact-matrix workflow
as a scaffold, not as a calibrated contact study.

The script:

- prefers a Prem/contactdata Kiribati matrix when it is available;
- does not silently substitute another Pacific country;
- falls back to POLYMOD UK only with an explicit limitation message;
- adapts the chosen source matrix to the model age grid;
- treats the matrix as one input to a provisional model, not as Kiribati
  social-contact truth.

That distinction matters. The model may be epidemiologically interesting, but
the contact matrix itself is still a proxy unless the source is actually
Kiribati-specific and available in the installed dataset.

## 14. Limitations And Current Boundaries

The current contact-matrix boundary is intentionally narrow:

- contact matrices are not automatically calibrated to local transmission;
- empirical social contacts are not the same thing as TB-relevant prolonged
  indoor exposure;
- source datasets may not include all countries;
- fine-to-coarse aggregation needs source-grid population weights;
- coarse-to-fine expansion assumes constant contacts within source bands;
- `ContactSchedule()` is not yet simulator-integrated;
- orientation must remain consistent from source loading through simulation;
- demography realism and contact realism are separate concerns.

## 15. Existing Examples

These scripts illustrate the main contact-matrix workflows:

- `examples/mock_sir_deterministic.R`: dependency-free manual matrix used in a
  simple SIR model.
- `examples/generic_sir.R`: manual matrix with a generic compartment model.
- `examples/mock_seir_demography.R`: manual matrix in a coupled SEIR and
  demography example.
- `examples/annual_cohort_sir_demography.R`: manual matrix in annual-cohort
  demographic coupling.
- `examples/kiribati_tb_realistic_demography.R`: source loading, adaptation,
  and a proxy-aware public-data scaffold.

## 16. Related Documentation

- [README.md](../README.md)
- [docs/model_conventions.md](model_conventions.md)
- [docs/external_data_adapters.md](external_data_adapters.md)
- [docs/demography_workflows.md](demography_workflows.md)
- [docs/age_structured_transmission_design.md](age_structured_transmission_design.md)
