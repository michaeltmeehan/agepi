# Disease parameterisation note: <Disease>

## 1. Purpose

This note summarises plausible parameter values for using `<Disease>` as an example disease in the `agepi` age-structured epidemic modelling prototype.

The purpose is to provide realistic values for testing and demonstration. This note does not define the core `agepi` model architecture and should not be treated as a definitive disease model.

---

## 2. Recommended minimal model structure

| Feature | Recommendation | Rationale | Source | Confidence |
|---|---|---|---|---|
| Minimal compartmental structure |  |  |  |  |
| Optional extensions |  |  |  |  |
| Initial prototype suitability | High / Medium / Low |  |  |  |

Notes:

- Recommended structure should be minimal and suitable for plugging into the generic `agepi` prototype.
- Optional extensions may include vaccination, maternal immunity, waning immunity, asymptomatic infection, hospitalisation, or disease-induced mortality, but these should not be assumed to be part of the core package.

---

## 3. Core natural-history parameters

| Parameter | Symbol | Suggested value | Plausible range | Unit | Distribution / uncertainty | Source | Confidence | Notes |
|---|---:|---:|---:|---|---|---|---|---|
| Latent period |  |  |  | days |  |  |  |  |
| Infectious period |  |  |  | days |  |  |  |  |
| Incubation period |  |  |  | days |  |  |  |  |
| Recovery rate |  |  |  | per day |  |  |  | Use `1 / mean infectious period` where appropriate. |
| Progression rate |  |  |  | per day |  |  |  | Use `1 / mean latent period` where appropriate. |
| Basic reproduction number | \(R_0\) |  |  | dimensionless |  |  |  |  |
| Duration of immunity |  |  |  | days / years / lifelong |  |  |  |  |

Where possible, convert durations into rates:

```text
rate = 1 / mean duration
```

---

## 4. Age-specific parameters

Use these default `agepi` model age groups unless a different age structure is explicitly requested:

```text
0-4, 5-9, 10-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
```

### 4.1 Susceptibility

| Age group | Suggested relative susceptibility | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 |  |  |  |  |
| 5-9 |  |  |  |  |
| 10-14 |  |  |  |  |
| 15-19 |  |  |  |  |
| 20-29 |  |  |  |  |
| 30-39 |  |  |  |  |
| 40-49 |  |  |  |  |
| 50-59 |  |  |  |  |
| 60-69 |  |  |  |  |
| 70-79 |  |  |  |  |
| 80+ |  |  |  |  |

### 4.2 Infectiousness

| Age group | Suggested relative infectiousness | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 |  |  |  |  |
| 5-9 |  |  |  |  |
| 10-14 |  |  |  |  |
| 15-19 |  |  |  |  |
| 20-29 |  |  |  |  |
| 30-39 |  |  |  |  |
| 40-49 |  |  |  |  |
| 50-59 |  |  |  |  |
| 60-69 |  |  |  |  |
| 70-79 |  |  |  |  |
| 80+ |  |  |  |  |

### 4.3 Morbidity / severity

Define the severity outcome before filling this table, for example symptomatic infection, hospitalisation, severe disease, ICU admission, or another disease-relevant outcome.

| Age group | Suggested morbidity / severity risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 |  |  |  |  |  |
| 5-9 |  |  |  |  |  |
| 10-14 |  |  |  |  |  |
| 15-19 |  |  |  |  |  |
| 20-29 |  |  |  |  |  |
| 30-39 |  |  |  |  |  |
| 40-49 |  |  |  |  |  |
| 50-59 |  |  |  |  |  |
| 60-69 |  |  |  |  |  |
| 70-79 |  |  |  |  |  |
| 80+ |  |  |  |  |  |

### 4.4 Infection-induced mortality

State clearly whether the value is an infection fatality risk, case fatality risk, disease-induced mortality rate, or another measure.

| Age group | Suggested infection-induced mortality risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 |  |  |  |  |  |
| 5-9 |  |  |  |  |  |
| 10-14 |  |  |  |  |  |
| 15-19 |  |  |  |  |  |
| 20-29 |  |  |  |  |  |
| 30-39 |  |  |  |  |  |
| 40-49 |  |  |  |  |  |
| 50-59 |  |  |  |  |  |
| 60-69 |  |  |  |  |  |
| 70-79 |  |  |  |  |  |
| 80+ |  |  |  |  |  |

---

## 5. Transmission and mixing assumptions

| Component | Recommendation | Source | Confidence | Notes |
|---|---|---|---|---|
| Contact matrix | Use external age-specific contact matrix |  |  | Candidate sources include `socialmixr`, `conmat`, or user-supplied matrices. |
| Transmission scaling | Calibrate beta to target \(R_0\) |  |  |  |
| Age-specific mixing |  |  |  |  |
| Seasonality | Include / exclude initially |  |  |  |
| Intervention effects | Defer / optional beta multiplier |  |  |  |
| Vaccination relevance |  |  |  | Mention only as an optional extension unless needed for basic realism. |

---

## 6. Initial conditions for demonstration runs

| Quantity | Suggested default | Source / rationale | Confidence | Notes |
|---|---:|---|---|---|
| Initial infected proportion |  |  |  |  |
| Initial infected age groups |  |  |  |  |
| Initial recovered / immune proportion |  |  |  |  |
| Initial susceptible assumption |  |  |  |  |
| Seeding approach |  |  |  |  |

---

## 7. Parameters as package inputs

Classify parameters according to how `agepi` should eventually treat them.

Use the following categories:

- `core`: required to run the transmission model;
- `risk_output`: used only for derived burden estimates;
- `optional_extension`: useful later but not required initially;
- `context_specific`: setting-dependent and should not have a universal default.

| Parameter | Type | Category | Suggested handling | Notes |
|---|---|---|---|---|
| beta | scalar or time-varying function | core | user-supplied or calibrated |  |
| gamma | scalar | core | disease default, user-overridable |  |
| sigma | scalar | core / optional_extension | disease default, user-overridable | Relevant for SEIR-like models. |
| susceptibility | age-specific vector | core | disease default, user-overridable |  |
| infectiousness | age-specific vector | core | disease default, user-overridable |  |
| morbidity risk | age-specific vector | risk_output | disease default, user-overridable |  |
| mortality risk | age-specific vector | risk_output | disease default, user-overridable |  |
| vaccination parameters | age-specific / time-varying | optional_extension / context_specific | defer initially |  |

---

## 8. Implementation implications for agepi

List implications for future examples or configuration only. Do not propose changes to the core architecture unless absolutely necessary.

| Observation | Possible implication | Priority | Notes |
|---|---|---|---|
|  |  | Low / Medium / High |  |

---

## 9. Cautions and limitations

Discuss cautions about:

- transferability across countries;
- temporal changes;
- vaccination history;
- demographic context;
- uncertainty in age-specific estimates;
- distinction between biological susceptibility and contact-driven risk;
- case fatality versus infection fatality;
- morbidity outcome definitions;
- compatibility between source age groups and `agepi` age groups.

---

## 10. Reference table

| Claim / parameter | Source | Source type | Confidence | Notes |
|---|---|---|---|---|
|  |  | WHO / CDC / ECDC / textbook / review / modelling study / primary study | High / Medium / Low |  |

---

## 11. Machine-readable parameter draft

This draft is intended as a provisional R-style list for later copying into an example or disease-parameter file. Use `NA_real_` where values are unavailable, uncertain, or deliberately deferred.

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
