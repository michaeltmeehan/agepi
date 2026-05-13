# Contact Matrix Integration Design Note

## Purpose

This note records a narrow reconnaissance pass for future contact-matrix
adapters. It does not propose adding package dependencies or changing the
current `force_of_infection()` convention.

Current agepi convention:

- `contact_matrix[a, b]` is the contact rate from recipient age group `a` with
  source age group `b`;
- rows are recipients;
- columns are infectious/source age groups;
- matrices passed to `force_of_infection()` must be square numeric matrices.

## Sources Inspected

Local agepi files:

- `docs/age_structured_transmission_design.md`
- `R/force_of_infection.R`
- `R/age_transform.R`

External package documentation:

- socialmixr `contact_matrix()` reference:
  <https://epiforecasts.io/socialmixr/reference/contact_matrix.html>
- conmat package site:
  <https://idem-lab.github.io/conmat/>
- conmat reference manual:
  <https://idem-lab.r-universe.dev/conmat/doc/manual.html>

Neither `socialmixr` nor `conmat` was installed in the local R library during
this pass.

## socialmixr Behaviour

### Object Types Returned

`socialmixr::contact_matrix()` returns a list-like object whose default printed
example contains:

- `matrix`: a numeric contact matrix;
- `participants`: a data table/data frame summarising participant counts and
  proportions by participant age group.

Optional arguments can add or change returned components:

- `return_demography = TRUE` can return the demography used for age weighting,
  symmetry, per-capita matrices, or split matrices;
- `return_part_weights = TRUE` can return participant weights;
- `split = TRUE` returns a decomposed representation rather than a single plain
  matrix;
- `counts = TRUE` returns counts rather than mean contacts;
- `per_capita = TRUE` returns per-capita contact rates, but is incompatible
  with `counts = TRUE` and `split = TRUE`.

For a first agepi adapter, the only supported input should be the ordinary
`contact_matrix()` result with a numeric `x$matrix`.

### Age-Group Representation

Age groups are represented as interval labels on row and column names, for
example `[0,1)`, `[1,5)`, `[5,15)`, `[15,Inf)`.

The `age_limits` argument supplies lower limits for the age groups. If omitted,
socialmixr infers limits from participant and contact ages. Missing participant
ages can be kept as a row labelled `NA`; missing contact ages can be kept as a
column labelled `NA`. Those outputs should be rejected by a first agepi adapter.

### Orientation

The socialmixr documentation describes `mij` as contact rates between
participants of age `i` in rows and contacts of age `j` in columns. The example
prints rows as `age.group` and columns as `contact.age.group`.

This is naturally:

- rows: survey participants / contact makers;
- columns: reported contacts.

For agepi, if survey participants are interpreted as susceptible recipients and
reported contacts as potential infectious sources, this orientation matches the
current `recipient-by-source` convention and should not be transposed.

The adapter should still make the orientation explicit in its name or metadata,
because other epidemiological contact-matrix conventions often use the opposite
wording.

### Squareness

With ordinary age limits and no kept missing ages, socialmixr returns a square
matrix over the same age bins for participant and contact age groups. However,
options that keep missing ages can introduce `NA`-labelled row or column groups.
A first adapter should require:

- numeric matrix;
- square dimensions;
- identical row and column age labels;
- no `NA` age-group labels.

### Population Correction and Reciprocity

`symmetric = TRUE` makes the matrix reciprocal according to:

```text
c_ij * N_i = c_ji * N_j
```

This requires survey population data through `survey_pop` or an implicit lookup
where supported. `weigh_age = TRUE` weights participants against the population
age distribution. `per_capita = TRUE` divides columns by the corresponding
population group to produce per-capita contact rates.

agepi should not redo these corrections inside a minimal adapter. The adapter
should preserve the matrix supplied by socialmixr and record, where possible,
whether the caller requested symmetry, age weighting, counts, split, or
per-capita output.

## conmat Behaviour

### Object Types Returned

conmat has several relevant output shapes:

- `predict_contacts()` returns a data frame with columns
  `age_group_from`, `age_group_to`, and `contacts`;
- `predict_contacts_1y()` returns a data frame with `age_from`, `age_to`,
  `contacts`, and `se_contacts`;
- `predictions_to_matrix()` converts predicted contacts to a square matrix;
- `predict_setting_contacts()` returns a list of setting matrices;
- `setting_prediction_matrix()` creates the matrix-list class used by outputs
  such as `extrapolate_polymod()` and `predict_setting_contacts()`;
- bundled Prem matrices are a named list of 16 by 16 matrices for `home`,
  `work`, `school`, `other`, and `all`.

For a first agepi adapter, the most useful supported inputs are:

- a plain square matrix produced by `predictions_to_matrix()`;
- a `predict_contacts()`-style data frame with `age_group_from`,
  `age_group_to`, and `contacts`;
- optionally, a named list of setting matrices where the caller explicitly
  chooses one setting or asks to sum settings.

### Age-Group Representation

conmat commonly represents age groups as interval factors such as `[0,5)`,
`[5,10)`, with `age_breaks` supplied as numeric lower bounds ending in `Inf`.
Population inputs can be data frames or `conmat_population` objects, which are
data frames with attributes identifying the age and population columns.

The conmat reference also describes `new_age_matrix()` as a numeric matrix with
an `age_breaks` attribute, usually from row names.

### Orientation

The conmat package site explicitly notes that its contact matrices are
transposed compared with Prem and Mossong: rows are "age group to" and columns
are "age group from".

The `predictions_to_matrix()` documentation says it converts predictions to
matrix format with contact age groups as rows and survey participant age groups
as columns. Its output is therefore:

- rows: `age_group_to` / contacted group;
- columns: `age_group_from` / participant or source group.

This matches agepi's current `recipient-by-source` convention if:

- `age_group_to` is interpreted as recipient/contacted age group;
- `age_group_from` is interpreted as source/contact-making age group.

No transposition should be needed for conmat matrix outputs that follow this
documented convention.

For long-form `predict_contacts()` data, agepi should pivot using
`age_group_to` as rows and `age_group_from` as columns.

### Squareness

`predictions_to_matrix()` returns a square matrix using the unique age groups
from the `age_group_from` and `age_group_to` columns. Bundled Prem matrices are
documented as 16 by 16 matrices over 5-year groups from 0 to 80.

A first adapter should still validate squareness and identical row/column age
labels after conversion.

### Population Correction and Reciprocity

conmat uses population data directly in prediction. `predict_contacts()` first
predicts contacts at one-year resolution and then aggregates predictions back to
the requested age resolution while weighting the contact rate by population.

Model fitting defaults to `symmetrical = TRUE`, using symmetric terms to reflect
the reciprocity of contacts. conmat also includes optional household adjustment
using per-capita household size, especially through setting-level workflows.

agepi should treat these as upstream modelling choices. A first adapter should
not refit models, enforce reciprocity, or apply household-size corrections.

## Recommended First Adapter API

Do not add these yet; this is the suggested shape for a later implementation.

```r
as_agepi_contact_matrix <- function(
  x,
  age_structure = NULL,
  source = c("matrix", "socialmixr", "conmat"),
  setting = NULL,
  combine_settings = c("error", "sum"),
  orientation = c("auto", "recipient_by_source", "participant_by_contact"),
  transpose = FALSE
) {
  # returns a numeric square matrix in agepi recipient-by-source orientation
}
```

Minimal source-specific helpers could be clearer:

```r
as_agepi_contact_matrix.socialmixr <- function(x, age_structure = NULL, transpose = FALSE) {
  # extract x$matrix; validate labels and dimensions
}

as_agepi_contact_matrix.conmat_predictions <- function(x, age_structure = NULL) {
  # pivot age_group_to rows by age_group_from columns
}

as_agepi_contact_matrix.conmat_matrix <- function(x, age_structure = NULL) {
  # validate matrix; read row/column names or age_breaks attribute
}
```

The first implementation should return only the numeric matrix expected by
`force_of_infection()`. If metadata is needed, return it as attributes or a
small list only after the current modelling API has a place for mixing metadata.

## Minimal Validation Needed

The first adapters should:

- extract or construct a numeric matrix;
- require square dimensions;
- require identical row and column age labels;
- reject missing, infinite, negative, or `NA`-labelled age groups;
- optionally compare labels with `age_structure$age_groups`;
- document whether transposition was applied;
- leave age aggregation/segregation to the separate age transformation API.

## Out of Scope for First Implementation

Keep the first implementation deliberately small. Defer:

- adding `socialmixr` or `conmat` to `DESCRIPTION`;
- changing `force_of_infection()`;
- constructing contact matrices from raw survey data;
- fitting conmat GAMs;
- downloading or managing population data;
- enforcing reciprocity after import;
- applying socialmixr age weighting, per-capita conversion, or split-matrix
  decomposition inside agepi;
- applying conmat household-size adjustments inside agepi;
- summing setting-specific matrices unless the caller explicitly requests it;
- transforming matrices across incompatible age structures before the age
  transformation API is final;
- supporting non-square or rectangular cross-classification matrices;
- storing time-varying contact matrices.

## Recommendation

The first contact-matrix integration should be an optional, dependency-free
coercion layer. It should accept already-created objects from socialmixr or
conmat, extract a plain numeric square matrix, validate age labels, and return
the matrix in agepi's recipient-by-source orientation.

For both socialmixr's ordinary `x$matrix` and conmat's documented matrix
outputs, the expected default is no transposition. The adapter should expose an
explicit `transpose` escape hatch because matrix orientation is a common source
of silent modelling errors.
