# Disease parameterisation note: Measles

## 1. Purpose

This note summarises plausible parameter values for using measles as an example disease in the `agepi` age-structured epidemic modelling prototype.

The purpose is to provide realistic values for testing and demonstration. This note does not define the core `agepi` model architecture and should not be treated as a definitive disease model.

---

## 2. Recommended minimal model structure

| Feature | Recommendation | Rationale | Source | Confidence |
|---|---|---|---|---|
| Minimal compartmental structure | SEIR | Measles has a clear exposed/latent phase before infectiousness, followed by a short highly infectious period and long-lasting post-infection immunity. | CDC Clinical Overview; CDC Behind the Model measles simulator | High |
| Optional extensions | Vaccination, maternal immunity, age-specific prior immunity, importation/seeding, seasonality, disease-induced mortality, and outbreak isolation | These features are important in realistic measles applications but should be kept outside the minimal generic prototype. | WHO Measles fact sheet; WHO surveillance standards | High |
| Initial prototype suitability | High | Measles is a strong stress-test for age-structured transmission because it combines high transmissibility, durable immunity, vaccine relevance, and pronounced age/contact dependence. | Guerra et al. 2017; WHO Measles fact sheet | High |

Notes:

- A minimal SEIR model is preferable to SIR because the delay from infection to infectiousness is epidemiologically important.
- Vaccination and pre-existing immunity dominate real measles risk in most settings, so examples should be explicit about whether the simulated population is hypothetical and fully susceptible or partially immune.

---

## 3. Core natural-history parameters

| Parameter | Symbol | Suggested value | Plausible range | Unit | Distribution / uncertainty | Source | Confidence | Notes |
|---|---:|---:|---:|---|---|---|---|---|
| Latent period | \(1 / \sigma\) | 8 | 7-11 | days | Fixed default or gamma-distributed waiting time in stochastic examples | CDC Behind the Model measles simulator; WHO surveillance standards | Medium | CDC's recent simulator uses about 11 days for latent period and about 9 days for infectious period; a simpler SEIR demonstration can use 8 days if infectiousness begins about 4 days before rash and rash occurs around day 14. |
| Infectious period | \(1 / \gamma\) | 8 | 8-9 | days | Fixed default or gamma-distributed waiting time in stochastic examples | CDC Clinical Overview; WHO surveillance standards; CDC Behind the Model measles simulator | High | Public-health infectious window is usually from 4 days before to 4 days after rash onset. |
| Incubation period to prodrome |  | 11.5 | 7-23 | days | Use 11-12 days as central value | CDC Clinical Diagnosis Fact Sheet; WHO surveillance standards | High | CDC states 11-12 days to first symptoms; WHO gives 10-14 days usually, range 7-23 days. |
| Incubation period to rash |  | 14 | 7-21 | days | Use 14 days as central value | CDC Yellow Book; CDC About Measles | High | Rash usually appears about 14 days after exposure. |
| Recovery rate | \(\gamma\) | 0.125 | 0.111-0.125 | per day | Derived as \(1 / infectious_period_days\) | Derived from infectious-period sources above | High | Use `1 / 8` for the suggested value. |
| Progression rate | \(\sigma\) | 0.125 | 0.091-0.143 | per day | Derived as \(1 / latent_period_days\) | Derived from latent-period sources above | Medium | Use `1 / 8` for a simple SEIR default; set to `1 / 11` if aligning with CDC simulator definitions. |
| Basic reproduction number | \(R_0\) | 15 | 12-18 | dimensionless | Treat as context-dependent; calibrate beta to target value | Guerra et al. 2017; Anderson and May value cited therein | Medium | The 12-18 range is widely used, but the systematic review emphasises that measles \(R_0\) is not a biological constant. |
| Duration of immunity after infection |  | Lifelong | Lifelong / very long | years | Not modelled as waning in the minimal model | CDC About Measles; WHO Measles fact sheet | High | Waning of infection-derived immunity is not needed for a basic example. |

Where possible, convert durations into rates:

```text
rate = 1 / mean duration
```

---

## 4. Age-specific parameters

Use these default `agepi` model age groups:

```text
0-4, 5-9, 10-14, 15-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
```

### 4.1 Susceptibility

These values are placeholders for biological susceptibility among non-immune people. In real populations, age-specific susceptibility is driven mainly by vaccination coverage, prior infection, maternal antibodies in infants, and contact structure.

| Age group | Suggested relative susceptibility | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Use 1.00 for unvaccinated/non-immune demonstration populations; infants may have temporary maternal protection depending on age and maternal immunity. |
| 5-9 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only; setting-specific vaccination history should override. |
| 10-14 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 15-19 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 20-29 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 30-39 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 40-49 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 50-59 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 60-69 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 70-79 | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |
| 80+ | 1.00 | WHO Measles fact sheet; WHO surveillance standards | Low | Biological default only. |

### 4.2 Infectiousness

No robust general-purpose age-specific infectiousness vector was identified for measles. Use equal infectiousness by age in the minimal prototype unless a setting-specific study provides otherwise.

| Age group | Suggested relative infectiousness | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 5-9 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 10-14 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 15-19 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 20-29 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 30-39 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 40-49 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 50-59 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 60-69 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 70-79 | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |
| 80+ | 1.00 | CDC Clinical Overview; WHO surveillance standards | Low | Placeholder. |

### 4.3 Morbidity / severity

Outcome definition: probability of any recognised measles complication among reported or clinically apparent measles cases. This is not an infection fatality risk and is not directly comparable across surveillance systems.

| Age group | Suggested morbidity / severity risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.30 | Any complication among cases | CDC About Measles; WHO position-paper summary | Medium | CDC describes about 1 in 5 unvaccinated US cases hospitalised and complications as common in children under 5; WHO summary gives complications in approximately 30% of reported cases. |
| 5-9 | 0.20 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder lower than under-5 risk. |
| 10-14 | 0.15 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |
| 15-19 | 0.15 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |
| 20-29 | 0.20 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Adults are often noted as higher-risk than school-age children, but exact age mapping is context-dependent. |
| 30-39 | 0.20 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |
| 40-49 | 0.20 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |
| 50-59 | 0.25 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |
| 60-69 | 0.25 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |
| 70-79 | 0.30 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |
| 80+ | 0.30 | Any complication among cases | CDC About Measles; WHO position-paper summary | Low | Placeholder. |

### 4.4 Infection-induced mortality

Outcome definition: case fatality risk among reported or clinically apparent measles cases, not infection fatality risk. Values below are demonstration placeholders for a high-income setting with access to care. They should be replaced for low-resource, conflict-affected, malnourished, or immunocompromised populations.

| Age group | Suggested infection-induced mortality risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.0030 | Case fatality risk | CDC MMWR Colorado 2025; WHO Measles fact sheet; Portnoy et al. 2022 preprint | Low | CDC cites death in 0.1%-0.3% of US cases; global and LMIC risks can be much higher, especially in children under 5. |
| 5-9 | 0.0010 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 10-14 | 0.0010 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 15-19 | 0.0010 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 20-29 | 0.0015 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 30-39 | 0.0015 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 40-49 | 0.0020 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 50-59 | 0.0020 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 60-69 | 0.0030 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 70-79 | 0.0030 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |
| 80+ | 0.0030 | Case fatality risk | CDC MMWR Colorado 2025; Portnoy et al. 2022 preprint | Low | Placeholder. |

---

## 5. Transmission and mixing assumptions

| Component | Recommendation | Source | Confidence | Notes |
|---|---|---|---|---|
| Contact matrix | Use external age-specific contact matrix | POLYMOD-style matrices; `socialmixr`; `conmat` | Medium | Contact structure is the main age-specific transmission driver in the minimal model. |
| Transmission scaling | Calibrate beta to target \(R_0\) | Guerra et al. 2017 | High | Do not hard-code beta; solve or scale it for the chosen contact matrix, susceptibility vector, infectiousness vector, and target \(R_0\). |
| Age-specific mixing | Strong assortative school-age mixing is plausible | General contact-matrix literature | Medium | School and household mixing are important for measles-like transmission. |
| Seasonality | Exclude initially | Modelling simplification | Medium | Add later as a beta multiplier if needed for school-term or seasonal demonstrations. |
| Intervention effects | Defer or use optional beta multiplier / isolation effect | CDC Clinical Overview; WHO surveillance standards | Medium | Public-health isolation and contact tracing matter operationally but are optional for the first prototype. |
| Vaccination relevance | Important optional extension; not part of minimal SEIR | WHO Measles fact sheet | High | For realistic examples, vaccination or pre-existing immunity should be represented in initial conditions or a vaccination compartment. |

---

## 6. Initial conditions for demonstration runs

| Quantity | Suggested default | Source / rationale | Confidence | Notes |
|---|---:|---|---|---|
| Initial infected proportion | 0.00001 | Demonstration seed | Low | Equivalent to 1 infectious person per 100,000 population. |
| Initial infected age groups | 5-9 or 10-14 | School-age mixing demonstration | Low | Use one seeded age group for clear examples; importation into any age group is possible. |
| Initial recovered / immune proportion | 0.00 for fully susceptible demo; setting-specific otherwise | Modelling simplification; WHO vaccination context | Medium | A fully susceptible population is useful for checking \(R_0\)-driven behaviour but unrealistic in most modern settings. |
| Initial susceptible assumption | `1 - I0 - E0 - R0` within each age group | Mass-balance requirement | High | If using vaccination, subtract immune/vaccinated people from susceptibility. |
| Seeding approach | Seed one or a few exposed/infectious individuals | Demonstration seed | Medium | For outbreak simulations, importation pulses may be more realistic. |

---

## 7. Parameters as package inputs

| Parameter | Type | Category | Suggested handling | Notes |
|---|---|---|---|---|
| beta | scalar or time-varying function | core | user-supplied or calibrated | Calibrate to target \(R_0\), not a universal disease default. |
| gamma | scalar | core | disease default, user-overridable | Suggested 0.125 per day. |
| sigma | scalar | core | disease default, user-overridable | Suggested 0.125 per day for simple SEIR; user may choose 0.091 per day for 11-day latent period. |
| R0 target | scalar | context_specific | user-supplied or example default | Suggested 15 with plausible range 12-18. |
| susceptibility | age-specific vector | core | disease default, user-overridable | Equal biological placeholder; setting-specific immunity should override. |
| infectiousness | age-specific vector | core | disease default, user-overridable | Equal placeholder. |
| morbidity risk | age-specific vector | risk_output | disease default, user-overridable | Placeholder complication risks, not core transmission. |
| mortality risk | age-specific vector | risk_output | disease default, user-overridable | Placeholder case fatality risks; highly setting-dependent. |
| vaccination parameters | age-specific / time-varying | optional_extension / context_specific | defer initially | Essential for realistic measles scenarios but optional for generic prototype. |
| maternal immunity | age-specific infant modifier | optional_extension / context_specific | defer initially | Relevant if splitting the 0-4 group into infant ages. |

---

## 8. Implementation implications for agepi

List implications for future examples or configuration only. These are not proposed changes to the core architecture.

| Observation | Possible implication | Priority | Notes |
|---|---|---|---|
| Measles is better represented as SEIR than SIR. | Use measles as an example for latent-period support when the prototype supports an exposed compartment. | Medium | Do not make measles-specific code part of the core package. |
| \(R_0\) is high and contact-matrix dependent. | Provide examples that calibrate beta from a target \(R_0\). | Medium | Avoid universal beta defaults. |
| Real measles dynamics depend strongly on vaccination and prior immunity. | Keep vaccination as an optional example layer or initial-condition setting. | Medium | This should be context-specific. |
| Age-specific severity and mortality evidence is weak in the default age groups. | Mark risk outputs as illustrative and user-overridable. | Low | Do not overstate precision. |

---

## 9. Cautions and limitations

- Transferability across countries is limited. Mortality and complications vary strongly with nutrition, vitamin A status, health-care access, surveillance completeness, conflict, and age distribution.
- Temporal changes matter. Vaccination coverage, outbreak response, and immunity profiles change over time.
- Vaccination history dominates susceptibility. The equal susceptibility vector is only a biological placeholder for non-immune people.
- Maternal immunity is obscured by the broad `0-4` group and should not be inferred from this note.
- Biological susceptibility should not be confused with realised infection risk, which also depends on contacts and immunity.
- The mortality table uses case fatality risk placeholders, not infection fatality risks.
- Morbidity is defined as any recognised complication among cases; other outcomes such as hospitalisation, pneumonia, encephalitis, or SSPE need separate parameterisation.
- The age groups used here are `agepi` defaults and do not align exactly with many published measles severity or CFR groupings.

---

## 10. Reference table

| Claim / parameter | Source | Source type | Confidence | Notes |
|---|---|---|---|---|
| Incubation to prodrome 11-12 days; infectious 4 days before through 4 days after rash | CDC. "Measles Clinical Diagnosis Fact Sheet." https://www.cdc.gov/measles/hcp/communication-resources/clinical-diagnosis-fact-sheet.html | Public health agency | High | Used for incubation and infectious-period defaults. |
| Patients contagious from 4 days before to 4 days after rash; isolate for 4 days after rash | CDC. "Clinical Overview of Measles." https://www.cdc.gov/measles/hcp/clinical-overview/index.html | Public health agency | High | Used for infectious period. |
| Rash usually appears about 14 days after exposure | CDC Yellow Book. "Measles (Rubeola)." https://www.cdc.gov/yellow-book/hcp/travel-associated-infections-diseases/measles-rubeola.html | Public health agency / travel medicine | High | Used for incubation to rash. |
| Incubation usually 10-14 days, range 7-23 days; contagious around 4 days before to 4 days after rash | WHO. "Vaccine Preventable Diseases Surveillance Standards: Measles." https://www.who.int/publications/m/item/vaccine-preventable-diseases-surveillance-standards-measles | WHO surveillance standard | High | Used for plausible incubation range. |
| Measles is highly contagious, serious, airborne, and can lead to severe complications and death; 2024 deaths mostly under 5 | WHO. "Measles" fact sheet. https://www.who.int/news-room/fact-sheets/detail/measles | WHO fact sheet | High | Used for model rationale, severity caveats, and vaccination relevance. |
| R0 commonly cited as 12-18 but context-dependent | Guerra FM, Bolotin S, Lim G, et al. "The basic reproduction number (R0) of measles: a systematic review." Lancet Infect Dis. 2017;17(12):e420-e428. doi:10.1016/S1473-3099(17)30307-9 | Systematic review | Medium | Used for \(R_0\). |
| SEIRV simulator uses latent and infectious-period assumptions for measles | CDC. "Behind the Model: Interactive Measles Outbreak Simulator." https://www.cdc.gov/cfa-behind-the-model/php/data-research/interactive-measles-outbreak-simulator.html | Public health agency modelling documentation | Medium | Used for latent/infectious period cross-check. |
| Complications and deaths in US cases: pneumonia, encephalitis, death 0.1%-0.3% | CDC MMWR. "Measles Outbreak Associated with an Infectious Traveler - Colorado, May-June 2025." https://www.cdc.gov/mmwr/volumes/75/wr/mm7504a1.htm | Public health agency report | Medium | Used for mortality placeholder in high-income settings. |
| National-level measles CFR varies substantially in LMIC settings | Portnoy A, et al. "Estimating national-level measles case fatality ratios: an updated systematic review and modelling study." medRxiv. doi:10.1101/2022.10.05.22280730 | Systematic review / modelling preprint | Low | Preprint; used only to flag context-dependence of CFR. |

---

## 11. Machine-readable parameter draft

This draft is intended as a provisional R-style list for later copying into an example or disease-parameter file. Use `NA_real_` where values are unavailable, uncertain, or deliberately deferred.

```r
disease_parameters <- list(
  disease = "measles",
  model = "SEIR",
  natural_history = list(
    latent_period_days = 8,
    infectious_period_days = 8,
    incubation_period_days = 11.5,
    incubation_to_rash_days = 14,
    gamma = 1 / 8,
    sigma = 1 / 8,
    R0 = 15
  ),
  age_specific = list(
    age_groups = c("0-4", "5-9", "10-14", "15-19", "20-29",
                   "30-39", "40-49", "50-59", "60-69",
                   "70-79", "80+"),
    susceptibility = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    infectiousness = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    morbidity_risk = c(0.30, 0.20, 0.15, 0.15, 0.20,
                       0.20, 0.20, 0.25, 0.25, 0.30, 0.30),
    mortality_risk = c(0.0030, 0.0010, 0.0010, 0.0010, 0.0015,
                       0.0015, 0.0020, 0.0020, 0.0030, 0.0030, 0.0030)
  ),
  parameter_classification = list(
    core = c("beta", "gamma", "sigma", "susceptibility", "infectiousness"),
    risk_output = c("morbidity_risk", "mortality_risk"),
    optional_extension = c("vaccination", "maternal_immunity", "seasonality",
                           "importation", "isolation"),
    context_specific = c("R0", "contact_matrix", "prior_immunity",
                         "vaccination_coverage", "case_fatality_risk")
  ),
  notes = list(
    source_summary = paste(
      "Natural history values use CDC and WHO measles clinical and",
      "surveillance guidance. R0 uses Guerra et al. 2017 systematic review.",
      "Age-specific susceptibility and infectiousness are equal placeholders",
      "for non-immune people."
    ),
    cautions = paste(
      "Vaccination history and prior infection should override the",
      "susceptibility vector in realistic settings. Morbidity and mortality",
      "risks are illustrative case-based placeholders and are highly",
      "setting-dependent."
    )
  )
)
```
