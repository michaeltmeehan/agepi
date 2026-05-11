# Disease parameterisation note: Tuberculosis

## 1. Purpose

This note summarises plausible parameter values for using tuberculosis (TB) as an example disease in the `agepi` age-structured epidemic modelling prototype.

The purpose is to provide realistic values for testing and demonstration. This note does not define the core `agepi` model architecture and should not be treated as a definitive TB transmission model.

---

## 2. Recommended minimal model structure

| Feature | Recommendation | Rationale | Source | Confidence |
|---|---|---|---|---|
| Minimal compartmental structure | Susceptible-infected-active/recovered style only for very simple demonstrations; otherwise TB-specific latent infection and active disease states are needed | TB has long, variable latency, progression from infection to disease, and disease-driven infectiousness rather than a short acute infectious period. A standard acute SEIR model is a poor biological fit. | WHO TB fact sheet; CDC Clinical Overview of Tuberculosis; Houben and Dodd 2016 | Medium |
| Optional extensions | Latent TB infection, fast and slow progression, endogenous reactivation, reinfection, treatment, HIV/comorbidity risk multipliers, BCG, drug resistance, and disease-induced mortality | These are central to realistic TB applications but are too disease-specific for a minimal generic age-structured prototype. | WHO TB fact sheet; CDC Clinical Overview of Tuberculosis; Houben and Dodd 2016 | High |
| Initial prototype suitability | Medium | TB is useful for testing chronic infection and age-risk examples, but its natural history is not well represented by simple SIR/SEIR defaults. | Ma et al. 2018; Houben and Dodd 2016 | Medium |

Notes:

- For a generic `agepi` prototype, TB should be treated as an illustrative chronic-infection stress test, not as the first example for acute epidemic dynamics.
- If only SIR/SEIR model forms are currently available, use TB parameter values cautiously and describe the example as a simplified active pulmonary TB transmission scenario.
- For realistic TB, the model needs at least susceptible, latent infection, active infectious disease, and treated/recovered states, with progression and reactivation processes.

---

## 3. Core natural-history parameters

| Parameter | Symbol | Suggested value | Plausible range | Unit | Distribution / uncertainty | Source | Confidence | Notes |
|---|---:|---:|---:|---|---|---|---|---|
| Latent period, acute SEIR analogue | \(1 / \sigma\) | `NA_real_` | months to years | days / years | Not recommended as a single acute-model duration | CDC Clinical Overview of Tuberculosis; WHO TB fact sheet; Behr et al. 2018 | High | TB latency is highly variable; progression is most common in the first 2 years after infection but can occur later. |
| Infectious period, untreated active pulmonary TB | \(1 / \gamma\) | 365 | 180-730 | days | Use wide uncertainty or scenario-specific duration | Tiemersma et al. 2011; CDC treatment module | Low | This is a modelling placeholder for duration of infectious disease before cure, death, or treatment. In real settings, diagnosis and treatment strongly shorten infectiousness. |
| Incubation period |  | `NA_real_` | months to years | days / years | Not applicable as an acute incubation-period default | CDC Clinical Overview of Tuberculosis; Behr et al. 2018 | High | TB disease can appear soon after infection or after long latent infection; a single incubation value is misleading. |
| Recovery/removal rate from active infectious TB | \(\gamma\) | 1 / 365 | 1 / 730-1 / 180 | per day | Derived from active infectious duration placeholder | Derived from untreated-duration assumption above | Low | Use only for simplified active-TB examples. Treated TB should use treatment-specific removal rates. |
| Progression rate from latent infection to active disease |  | 0.0015 | 0.0005-0.005 | per day during early high-risk period | Better represented as age/time-since-infection specific | CDC Clinical Overview of Tuberculosis; WHO TB fact sheet; Behr et al. 2018 | Low | Approximate early progression placeholder, not a lifetime constant. Equivalent to about 0.55 per year; use only in short illustrative simulations. |
| Lifetime risk of active TB after infection |  | 0.05-0.10 | 0.05-0.15 | probability | Context- and risk-factor-specific | WHO TB fact sheet; CDC Clinical Overview of Tuberculosis | Medium | WHO describes about 5-10% of infected people eventually developing TB disease; risk is much higher with HIV and other risk factors. |
| Basic reproduction number | \(R_0\) | 1.5 | 1.0-4.0 | dimensionless | Context-specific; calibrate beta to contact matrix and active disease duration | Ma et al. 2018 | Low | Published TB reproduction-number estimates are heterogeneous and setting-specific; do not hard-code beta. |
| Duration of immunity |  | `NA_real_` | partial / imperfect | years | Reinfection possible; immunity is not lifelong sterilising | WHO TB fact sheet; Houben and Dodd 2016 | Medium | Previous infection may reduce but does not eliminate future TB risk; avoid a simple lifelong-immune recovered state for realistic models. |

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

These values are placeholders for relative risk of infection after effective exposure, not realised TB incidence. Realised risk depends strongly on household exposure, congregate settings, prior infection, BCG, HIV, diabetes, undernutrition, tobacco use, alcohol use, and local TB prevalence.

| Age group | Suggested relative susceptibility | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 1.20 | WHO TB fact sheet; Seddon and Shingadia 2014 | Low | Young children have high risk of disease after infection; infection susceptibility is hard to separate from progression risk. |
| 5-9 | 0.80 | Seddon and Shingadia 2014 | Low | School-age children often have lower progression risk than younger children and adults. |
| 10-14 | 0.90 | Seddon and Shingadia 2014 | Low | Placeholder. |
| 15-19 | 1.00 | WHO TB fact sheet; Seddon and Shingadia 2014 | Low | Adolescent risk begins to resemble adult pulmonary TB patterns. |
| 20-29 | 1.00 | WHO TB fact sheet | Low | Adult baseline placeholder. |
| 30-39 | 1.00 | WHO TB fact sheet | Low | Adult baseline placeholder. |
| 40-49 | 1.05 | WHO TB fact sheet | Low | Slightly higher placeholder to reflect accumulating comorbidity risk. |
| 50-59 | 1.10 | WHO TB fact sheet | Low | Placeholder. |
| 60-69 | 1.15 | WHO TB fact sheet | Low | Placeholder. |
| 70-79 | 1.20 | WHO TB fact sheet | Low | Placeholder. |
| 80+ | 1.20 | WHO TB fact sheet | Low | Placeholder. |

### 4.2 Infectiousness

Outcome definition: relative infectiousness among people in active TB disease. Children are less likely to have smear-positive cavitary pulmonary disease and are generally less infectious than adults; most transmission is from pulmonary or laryngeal TB.

| Age group | Suggested relative infectiousness | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 0.10 | CDC Clinical Overview of Tuberculosis; Seddon and Shingadia 2014 | Low | Young children are often paucibacillary and less infectious. |
| 5-9 | 0.20 | CDC Clinical Overview of Tuberculosis; Seddon and Shingadia 2014 | Low | Placeholder. |
| 10-14 | 0.50 | CDC Clinical Overview of Tuberculosis; Seddon and Shingadia 2014 | Low | Transitional age group. |
| 15-19 | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult-like pulmonary TB infectiousness placeholder. |
| 20-29 | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult baseline. |
| 30-39 | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult baseline. |
| 40-49 | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult baseline. |
| 50-59 | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult baseline. |
| 60-69 | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult baseline. |
| 70-79 | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult baseline. |
| 80+ | 1.00 | CDC Clinical Overview of Tuberculosis | Low | Adult baseline. |

### 4.3 Morbidity / severity

Outcome definition: probability of severe TB disease among active TB disease episodes, interpreted broadly as disseminated, meningeal, severe pulmonary, hospitalised, or otherwise clinically severe TB. These are illustrative placeholders, not validated age-specific risks.

| Age group | Suggested morbidity / severity risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.20 | Severe TB among active TB disease | WHO TB fact sheet; Seddon and Shingadia 2014 | Low | Young children have higher risk of severe forms including TB meningitis and disseminated TB. |
| 5-9 | 0.08 | Severe TB among active TB disease | Seddon and Shingadia 2014 | Low | Placeholder lower than under-5. |
| 10-14 | 0.08 | Severe TB among active TB disease | Seddon and Shingadia 2014 | Low | Placeholder. |
| 15-19 | 0.10 | Severe TB among active TB disease | WHO TB fact sheet | Low | Placeholder. |
| 20-29 | 0.10 | Severe TB among active TB disease | WHO TB fact sheet | Low | Adult baseline placeholder. |
| 30-39 | 0.10 | Severe TB among active TB disease | WHO TB fact sheet | Low | Adult baseline placeholder. |
| 40-49 | 0.12 | Severe TB among active TB disease | WHO TB fact sheet | Low | Placeholder reflecting higher comorbidity burden. |
| 50-59 | 0.15 | Severe TB among active TB disease | WHO TB fact sheet | Low | Placeholder. |
| 60-69 | 0.18 | Severe TB among active TB disease | WHO TB fact sheet | Low | Placeholder. |
| 70-79 | 0.20 | Severe TB among active TB disease | WHO TB fact sheet | Low | Placeholder. |
| 80+ | 0.25 | Severe TB among active TB disease | WHO TB fact sheet | Low | Placeholder. |

### 4.4 Infection-induced mortality

Outcome definition: case fatality risk among active TB disease episodes. Values below are illustrative, diagnosis-and-treatment-dependent placeholders. They are not infection fatality risks and should not be interpreted as age-specific estimates for any country.

| Age group | Suggested infection-induced mortality risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.08 | Case fatality risk among active TB disease | WHO TB fact sheet; Dodd et al. 2017; Tiemersma et al. 2011 | Low | High severe-disease risk in young children; mortality strongly depends on treatment access and HIV. |
| 5-9 | 0.03 | Case fatality risk among active TB disease | Dodd et al. 2017; Tiemersma et al. 2011 | Low | Placeholder. |
| 10-14 | 0.03 | Case fatality risk among active TB disease | Dodd et al. 2017; Tiemersma et al. 2011 | Low | Placeholder. |
| 15-19 | 0.04 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Placeholder. |
| 20-29 | 0.05 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Adult placeholder assuming some treatment access. |
| 30-39 | 0.06 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Placeholder. |
| 40-49 | 0.08 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Placeholder. |
| 50-59 | 0.10 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Placeholder. |
| 60-69 | 0.15 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Placeholder. |
| 70-79 | 0.20 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Placeholder. |
| 80+ | 0.25 | Case fatality risk among active TB disease | Tiemersma et al. 2011 | Low | Placeholder. |

---

## 5. Transmission and mixing assumptions

| Component | Recommendation | Source | Confidence | Notes |
|---|---|---|---|---|
| Contact matrix | Use external age-specific contact matrix, ideally with household and close-contact structure | General contact-matrix literature; CDC Clinical Overview of Tuberculosis | Medium | TB transmission often occurs through prolonged indoor exposure; generic social-contact matrices may miss household, workplace, prison, shelter, and health-care risks. |
| Transmission scaling | Calibrate beta to target \(R_0\) or incidence / prevalence target | Ma et al. 2018; Houben and Dodd 2016 | High | Do not hard-code beta; calibrate to the chosen contact matrix, infectiousness vector, active disease duration, and local epidemiology. |
| Age-specific mixing | Include household and adult/adolescent mixing when possible | CDC Clinical Overview of Tuberculosis; Seddon and Shingadia 2014 | Medium | Children are more often infected by adults than the main source of onward transmission. |
| Seasonality | Exclude initially | Modelling simplification | Medium | TB seasonality exists in some settings but is not needed for a first demonstration. |
| Intervention effects | Defer or represent as diagnosis/treatment rate and transmission multiplier | WHO TB fact sheet; CDC treatment module | Medium | Treatment rapidly changes infectiousness and survival, so intervention effects are central in realistic TB models. |
| Vaccination relevance | BCG is optional and context-specific | WHO TB fact sheet | Medium | BCG mainly protects young children from severe forms of TB; it is not a simple infection-blocking vaccine in this prototype. |

---

## 6. Initial conditions for demonstration runs

| Quantity | Suggested default | Source / rationale | Confidence | Notes |
|---|---:|---|---|---|
| Initial infected proportion | 0.00001 active infectious | Demonstration seed | Low | Equivalent to 1 active infectious person per 100,000 population. |
| Initial infected age groups | 20-29 or 30-39 | Demonstration seed | Low | Use an adult age group to reflect pulmonary TB transmission in a simple example. |
| Initial recovered / immune proportion | `NA_real_` or 0.00 in simple SIR examples | Modelling simplification | Medium | TB infection does not map cleanly to lifelong recovered immunity. |
| Initial latent infection proportion | Setting-specific; use 0.00 for disease-free introduction demo or 0.25 for global-background stress test | WHO TB fact sheet | Low | WHO estimates about one-quarter of the global population has been infected with TB bacteria; this is not appropriate for every country or age group. |
| Initial susceptible assumption | `1 - active - latent - treated/recovered` within each age group | Mass-balance requirement | High | If using only SIR/SEIR, document what latent infection and prior infection are being ignored. |
| Seeding approach | Seed one or a few active pulmonary TB cases | Demonstration seed | Medium | For endemic TB examples, initialise from prevalence or incidence rather than a single importation. |

---

## 7. Parameters as package inputs

| Parameter | Type | Category | Suggested handling | Notes |
|---|---|---|---|---|
| beta | scalar or time-varying function | core | user-supplied or calibrated | Calibrate to target \(R_0\), incidence, or prevalence; no universal TB default. |
| gamma | scalar | core for simplified active-TB model | disease default only for toy examples; user-overridable | Suggested `1 / 365` per day if using a simplified untreated active-disease duration. |
| sigma | scalar | optional_extension / context_specific | avoid as a universal acute-SEIR default | TB progression should be represented by progression/reactivation hazards, preferably age- and time-since-infection-specific. |
| R0 target | scalar | context_specific | user-supplied or example default | Suggested toy default 1.5 with plausible range 1.0-4.0. |
| susceptibility | age-specific vector | core | disease default, user-overridable | Low-confidence placeholders; progression risk may be a better age-specific parameter than susceptibility. |
| infectiousness | age-specific vector | core | disease default, user-overridable | Children have lower infectiousness placeholders; adults set to 1. |
| morbidity risk | age-specific vector | risk_output | disease default, user-overridable | Broad severe-TB placeholder. |
| mortality risk | age-specific vector | risk_output | disease default, user-overridable | Case fatality placeholders; strongly affected by diagnosis, treatment, HIV, and setting. |
| latent infection prevalence | age-specific vector | context_specific | user-supplied | Crucial for endemic TB examples. |
| progression / reactivation hazards | age-specific and time-since-infection parameters | optional_extension / context_specific | defer initially | Needed for realistic TB natural history. |
| diagnosis and treatment rates | scalar or age-specific rates | optional_extension / context_specific | defer initially | Treatment changes infectious duration, mortality, and onward transmission. |
| HIV/comorbidity risk multipliers | age-specific or subgroup-specific multipliers | optional_extension / context_specific | defer initially | Important in many settings, especially where HIV burden is high. |
| BCG parameters | age-specific protection parameters | optional_extension / context_specific | defer initially | Mainly relevant for severe paediatric TB outcomes. |

---

## 8. Implementation implications for agepi

List implications for future examples or configuration only. These are not proposed changes to the core architecture.

| Observation | Possible implication | Priority | Notes |
|---|---|---|---|
| TB is poorly represented by a short acute SEIR model. | Use TB as a later example once chronic/latent state support exists, or label any current example as a toy active-TB model. | High | Avoid implying TB has a single incubation or latent period like measles. |
| Age-specific infectiousness matters because children are usually less infectious. | Allow age-specific infectiousness vectors to affect the next-generation matrix. | Medium | This fits the generic age-structured model. |
| Latent infection and prior infection are common in many settings. | Demonstration initial conditions may need latent infection as a separate state or at least clear caveats. | Medium | WHO global latent infection estimates are not suitable as country defaults. |
| Treatment is central to TB dynamics. | Future examples could include a diagnosis/treatment removal rate or active-to-treated transition. | Medium | Keep outside the first generic core unless chronic infections are in scope. |

---

## 9. Cautions and limitations

- Transferability across countries is limited. TB incidence, latent infection prevalence, HIV burden, diabetes, undernutrition, tobacco use, health-care access, incarceration, homelessness, migration, and crowding all affect parameter values.
- Temporal changes matter. TB incidence, diagnosis, treatment access, drug resistance, and preventive therapy coverage change over time.
- BCG history affects severe paediatric TB outcomes but is not a simple universal susceptibility modifier.
- Biological susceptibility, progression from infection to disease, infectiousness once diseased, and realised incidence are distinct quantities.
- TB infection and disease are better viewed as a spectrum than a strict acute exposed/infectious/recovered sequence.
- The mortality table uses case fatality risk placeholders among active TB disease episodes, not infection fatality risks.
- The morbidity outcome is broad severe TB; meningitis, miliary TB, hospitalisation, and pulmonary severity need separate parameterisation in detailed work.
- The age groups used here are `agepi` defaults and do not align exactly with most published TB natural-history and burden estimates.
- Drug-resistant TB should not be modelled by changing infectiousness alone; delayed effective treatment can extend infectious duration.

---

## 10. Reference table

| Claim / parameter | Source | Source type | Confidence | Notes |
|---|---|---|---|---|
| TB spreads through the air; about one-quarter of the global population has been infected; about 5-10% of infected people develop TB disease; babies and children are higher risk | WHO. "Tuberculosis." https://www.who.int/news-room/fact-sheets/detail/tuberculosis | WHO fact sheet | High | Used for overall natural history, latent infection, and risk caveats. |
| Latent TB is not contagious; TB disease may spread; progression is most common within the first 2 years after infection; drug-resistant TB is not intrinsically more infectious but delayed recognition can prolong infectiousness | CDC. "Clinical Overview of Tuberculosis." https://www.cdc.gov/tb/hcp/clinical-overview/index.html | Public health agency | High | Used for latent/active distinction, progression timing, infectiousness, and drug-resistance caveats. |
| TB disease treatment takes at least 6 months in many regimens and prevents further transmission | CDC. "Self-Study Modules on Tuberculosis, Module 4: Treatment of Latent Tuberculosis Infection and Tuberculosis Disease." https://www.cdc.gov/tb/media/pdfs/Self_Study_Module_4_Treatment_of_Latent_TB_Infection_and_TB_Disease.pdf | Public health agency training module | Medium | Used for treatment and infectious-duration cautions. |
| TB natural history includes infection, latency, progression, reactivation, reinfection, active disease, and treatment | Houben RMGJ, Dodd PJ. "The global burden of latent tuberculosis infection: a re-estimation using mathematical modelling." PLoS Med. 2016;13(10):e1002152. doi:10.1371/journal.pmed.1002152; Ma et al. 2018 | Modelling literature | Medium | Used to caution against simple SEIR interpretation. |
| TB reproduction-number and serial-interval estimates are heterogeneous across settings and methods | Ma Y, Horsburgh CR Jr, White LF, Jenkins HE. "Quantifying TB transmission: a systematic review of reproduction number and serial interval estimates for tuberculosis." Epidemiol Infect. 2018;146(12):1478-1494. doi:10.1017/S0950268818001760 | Systematic review | Medium | Used for \(R_0\) caution and broad plausible range. |
| Most progression to active TB after infection occurs earlier than traditionally assumed, with later disease often reflecting reinfection or missed earlier disease in some settings | Behr MA, Edelstein PH, Ramakrishnan L. "Revisiting the timetable of tuberculosis." BMJ. 2018;362:k2738. doi:10.1136/bmj.k2738 | Review / perspective | Medium | Used for latency/progression cautions. |
| Untreated pulmonary TB has high mortality and long, variable duration | Tiemersma EW, van der Werf MJ, Borgdorff MW, Williams BG, Nagelkerke NJD. "Natural history of tuberculosis: duration and fatality of untreated pulmonary tuberculosis in HIV negative patients: a systematic review." PLoS One. 2011;6(4):e17601. doi:10.1371/journal.pone.0017601 | Systematic review | Medium | Used for active infectious duration and mortality placeholders. |
| Children have distinctive TB natural history, including high risk of severe disease in young children and lower infectiousness than adults | Seddon JA, Shingadia D. "Epidemiology and disease burden of tuberculosis in children: a global perspective." Infect Drug Resist. 2014;7:153-165. doi:10.2147/IDR.S45090 | Review | Medium | Used for paediatric susceptibility, infectiousness, morbidity, and mortality placeholders. |
| Global latent TB burden is large but uncertain and highly setting-specific | Houben RMGJ, Dodd PJ. "The global burden of latent tuberculosis infection: a re-estimation using mathematical modelling." PLoS Med. 2016;13(10):e1002152. doi:10.1371/journal.pmed.1002152 | Modelling study | Medium | Used for latent infection prevalence caveats. |
| Children are at high risk of severe and fatal TB without treatment | Dodd PJ, Gardiner E, Coghlan R, Seddon JA. "Burden of childhood tuberculosis in 22 high-burden countries: a mathematical modelling study." Lancet Glob Health. 2014;2(8):e453-e459. doi:10.1016/S2214-109X(14)70245-1 | Modelling study | Low | Used for paediatric mortality caveats. |

---

## 11. Machine-readable parameter draft

This draft is intended as a provisional R-style list for later copying into an example or disease-parameter file. Use `NA_real_` where values are unavailable, uncertain, or deliberately deferred.

```r
disease_parameters <- list(
  disease = "tuberculosis",
  model = "TB_latent_active",
  natural_history = list(
    latent_period_days = NA_real_,
    infectious_period_days = 365,
    incubation_period_days = NA_real_,
    gamma = 1 / 365,
    sigma = NA_real_,
    R0 = 1.5,
    lifetime_progression_risk = c(0.05, 0.10),
    early_progression_rate = 0.0015
  ),
  age_specific = list(
    age_groups = c("0-4", "5-9", "10-14", "15-19", "20-29",
                   "30-39", "40-49", "50-59", "60-69",
                   "70-79", "80+"),
    susceptibility = c(1.20, 0.80, 0.90, 1.00, 1.00,
                       1.00, 1.05, 1.10, 1.15, 1.20, 1.20),
    infectiousness = c(0.10, 0.20, 0.50, 1.00, 1.00,
                       1.00, 1.00, 1.00, 1.00, 1.00, 1.00),
    morbidity_risk = c(0.20, 0.08, 0.08, 0.10, 0.10,
                       0.10, 0.12, 0.15, 0.18, 0.20, 0.25),
    mortality_risk = c(0.08, 0.03, 0.03, 0.04, 0.05,
                       0.06, 0.08, 0.10, 0.15, 0.20, 0.25)
  ),
  parameter_classification = list(
    core = c("beta", "gamma", "susceptibility", "infectiousness"),
    risk_output = c("morbidity_risk", "mortality_risk"),
    optional_extension = c("latent_infection", "progression",
                           "reactivation", "reinfection", "treatment",
                           "BCG", "drug_resistance", "HIV_comorbidity"),
    context_specific = c("R0", "contact_matrix",
                         "latent_infection_prevalence",
                         "progression_rate", "diagnosis_rate",
                         "treatment_rate", "case_fatality_risk")
  ),
  notes = list(
    source_summary = paste(
      "TB values use WHO and CDC clinical summaries for latent infection,",
      "active disease, transmission, and treatment context. R0 and duration",
      "assumptions are broad modelling placeholders because TB dynamics are",
      "highly setting-specific."
    ),
    cautions = paste(
      "TB is not well represented by acute SEIR defaults. Use this parameter",
      "draft for simplified demonstrations only unless latent infection,",
      "progression, reactivation, treatment, and context-specific risk factors",
      "are modelled explicitly."
    )
  )
)
```
