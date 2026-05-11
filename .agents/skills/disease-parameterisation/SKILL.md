# Disease Parameterisation Skill

## Purpose

Use this skill to create a standardised disease-specific parameterisation note for the `agepi` repository.

`agepi` is an early-stage R package for age-structured epidemic model prototypes. Disease-specific research notes should provide realistic values for examples and demonstrations, but must not drive or alter the core package architecture.

## When to use this skill

Use this skill when asked to research or draft plausible parameter values for a disease, such as:

- "Use the disease-parameterisation skill for measles."
- "Create a disease parameterisation note for RSV."
- "Research influenza parameters for the agepi prototype."

## Required input

The user must provide a disease name.

Optional useful inputs include:

- geographical setting;
- time period;
- target model type, such as SIR or SEIR;
- preferred age groups;
- preferred source types;
- whether vaccination, maternal immunity, or waning immunity should be considered.

If these optional inputs are not provided, use the default `agepi` model age groups:

```text
0-4, 5-9, 10-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
```

## Output location

Write the disease-specific note to:

```text
docs/disease_notes/<disease>_parameterisation.md
```

Use a lowercase, underscore-separated disease name for the filename.

Examples:

```text
docs/disease_notes/measles_parameterisation.md
docs/disease_notes/influenza_parameterisation.md
docs/disease_notes/rsv_parameterisation.md
```

## Required template

Use this template if present:

```text
docs/disease_notes/TEMPLATE_parameterisation.md
```

If the template is missing, create the disease note using the fixed structure defined below.

## Hard constraints

Do not:

- edit package source code;
- modify files in `R/`;
- change the package architecture;
- implement disease-specific model code;
- modify `DESCRIPTION`, `NAMESPACE`, tests, or package infrastructure unless explicitly asked;
- make disease-specific assumptions part of the generic `agepi` design;
- present uncertain values as universal defaults.

Only create or edit files under:

```text
docs/disease_notes/
```

unless the user explicitly asks otherwise.

## Source standards

Prioritise authoritative and traceable sources:

1. WHO, CDC, ECDC, national public health agencies;
2. major epidemiology or infectious disease textbooks;
3. systematic reviews and meta-analyses;
4. well-cited peer-reviewed modelling papers;
5. primary epidemiological studies, when needed.

Avoid relying on informal websites unless no better source is available. If a value comes from a weak or indirect source, label confidence as low.

## Citation requirements

For every numerical parameter, provide a source.

For every source-dependent claim, include enough citation detail for the user to verify it later. Prefer DOI, PubMed, official webpage title, report title, or stable URL where available.

Do not invent citations. If a source cannot be verified, state that the value is provisional and needs verification.

## Required output structure

Every disease note must use this structure:

1. Purpose
2. Recommended minimal model structure
3. Core natural-history parameters
4. Age-specific parameters
   - susceptibility
   - infectiousness
   - morbidity/severity
   - infection-induced mortality
5. Transmission and mixing assumptions
6. Initial conditions for demonstration runs
7. Parameters as package inputs
8. Implementation implications for agepi
9. Cautions and limitations
10. Reference table
11. Machine-readable parameter draft

## Parameter table requirements

Where possible, include:

- parameter name;
- symbol;
- suggested value;
- plausible range;
- unit;
- uncertainty or distribution;
- source;
- confidence level;
- notes.

Use confidence levels:

- High: robust, widely accepted value from authoritative or repeated sources;
- Medium: plausible value with some support but context-dependent;
- Low: weak evidence, indirect evidence, or placeholder assumption.

## Age-specific parameter requirements

Use the default age groups unless instructed otherwise:

```text
0-4, 5-9, 10-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
```

For each of these, try to provide:

- relative susceptibility;
- relative infectiousness;
- morbidity or severity risk;
- infection-induced mortality risk.

If evidence is unavailable, do not fabricate precision. Use transparent placeholders only where reasonable, and clearly mark them as placeholders.

## Parameter classification

Classify each parameter using one of these categories:

- `core`: required to run the transmission model;
- `risk_output`: used only for derived burden estimates;
- `optional_extension`: useful later but not required initially;
- `context_specific`: setting-dependent and should not have a universal default.

## Machine-readable draft

End the note with an R-style list named `disease_parameters`.

The list should be syntactically close to valid R, but it may contain `NA_real_` where values are unavailable or deliberately deferred.

Use this structure:

```r
disease_parameters <- list(
  disease = "<disease>",
  model = "<SIR/SEIR/etc>",
  natural_history = list(
    latent_period_days = NA_real_,
    infectious_period_days = NA_real_,
    incubation_period_days = NA_real_,
    gamma = NA_real_,
    sigma = NA_real_,
    R0 = NA_real_
  ),
  age_specific = list(
    age_groups = c("0-4", "5-9", "10-14", "15-19", "20-29",
                   "30-39", "40-49", "50-59", "60-69",
                   "70-79", "80+"),
    susceptibility = c(),
    infectiousness = c(),
    morbidity_risk = c(),
    mortality_risk = c()
  ),
  parameter_classification = list(
    core = c(),
    risk_output = c(),
    optional_extension = c(),
    context_specific = c()
  ),
  notes = list(
    source_summary = "",
    cautions = ""
  )
)
```

## Quality checks before finishing

Before reporting completion, check that:

- the note is saved under `docs/disease_notes/`;
- the filename matches the disease name;
- all required sections are present;
- numerical values have sources;
- placeholders are clearly labelled;
- age-specific vectors match the required age groups;
- the note does not propose changes to `agepi` core architecture except as non-binding future implications;
- no source-code files were modified.

## Completion report

After creating the note, report:

- output file path;
- sources consulted;
- any parameters with weak or missing evidence;
- any assumptions made;
- confirmation that no package source code was changed.
