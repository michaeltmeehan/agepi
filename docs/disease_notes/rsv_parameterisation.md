# Disease parameterisation note: RSV

## 1. Purpose

This note summarises plausible parameter values for using respiratory syncytial virus (RSV) as an example disease in the `agepi` age-structured epidemic modelling prototype.

The purpose is to provide realistic values for testing and demonstration. This note does not define the core `agepi` model architecture and should not be treated as a definitive RSV transmission model.

---

## 2. Recommended minimal model structure

| Feature | Recommendation | Rationale | Source | Confidence |
|---|---|---|---|---|
| Minimal compartmental structure | SEIRS | RSV is an acute respiratory infection with a short incubation period, short infectious period, frequent reinfection, and incomplete or waning immunity. | CDC How RSV Spreads; CDC Pink Book RSV chapter; Hall et al. 2013 review | Medium |
| Optional extensions | Maternal antibody or infant protection, vaccination/nirsevimab, age-specific clinical risk, seasonal forcing, household/day-care mixing, health-care or long-term-care outbreaks | These features matter for realistic RSV burden, but should remain optional example layers rather than core package assumptions. | CDC RSV prevention guidance; CDC ACIP RSV adult evidence review; Li et al. 2022 | High |
| Initial prototype suitability | High | RSV is a useful age-structured example because transmission is common across ages while severe outcomes concentrate in infants, young children, older adults, and people with risk conditions. | CDC RSV infants and young children; CDC ACIP RSV adult evidence review; Li et al. 2022 | High |

Notes:

- A minimal SEIRS model is preferable to SIR because reinfection after RSV is common and immunity is not lifelong sterilising.
- If the current prototype only supports SIR or SEIR, use the values below as a short acute-infection example and document that waning immunity is omitted.
- Prevention products introduced in recent years make modern real-world severity and susceptibility patterns time- and setting-dependent.

---

## 3. Core natural-history parameters

| Parameter | Symbol | Suggested value | Plausible range | Unit | Distribution / uncertainty | Source | Confidence | Notes |
|---|---:|---:|---:|---|---|---|---|---|
| Latent period | \(1 / \sigma\) | 3 | 2-5 | days | Fixed default or gamma-distributed waiting time | CDC How RSV Spreads; Lessler et al. 2009; Usher Institute RSV transmission rapid review | Low | Direct latent-period estimates are limited; use a value shorter than or similar to incubation because contagiousness can begin before symptoms. |
| Infectious period | \(1 / \gamma\) | 6 | 3-8 | days | Fixed default or gamma-distributed waiting time | CDC How RSV Spreads; CDC Pink Book RSV chapter | Medium | CDC describes most infected people as contagious for 3-8 days; shedding can last longer in young infants and immunocompromised people. |
| Incubation period |  | 4.5 | 2-8 | days | Use 4-6 days as central range | CDC RSV symptoms; CDC Pink Book RSV chapter; Lessler et al. 2009 | High | CDC states symptoms usually appear 4-6 days after infection; public-health manuals commonly give a 2-8 day range. |
| Recovery rate | \(\gamma\) | 0.167 | 0.125-0.333 | per day | Derived as \(1 / infectious_period_days\) | Derived from infectious-period sources above | Medium | Use `1 / 6` for the suggested value. |
| Progression rate | \(\sigma\) | 0.333 | 0.200-0.500 | per day | Derived as \(1 / latent_period_days\) | Derived from latent-period assumption above | Low | Use `1 / 3`; treat as an illustrative SEIR timing parameter. |
| Basic reproduction number | \(R_0\) | 2.5 | 1.2-4.0 | dimensionless | Context-specific; calibrate beta to contact matrix and seasonal setting | Reis and Shaman 2018; Otomaru et al. 2019; Usher Institute RSV transmission rapid review | Low | Published RSV transmissibility estimates vary by model, setting, season, and age structure. |
| Serial interval / generation-time proxy |  | 3.2 | 2.5-4.1 | days | Gamma distribution in household transmission examples | Otomaru et al. 2019; Usher Institute RSV transmission rapid review | Medium | Useful for checking generation speed, not directly a compartment duration. |
| Duration of immunity after infection |  | 365 | 180-730 | days | Waning-immunity scenario parameter | Hall et al. 1991; Hall et al. 2013 review | Low | Reinfection can occur within months to years; use this only for SEIRS demonstration, not as a universal biological constant. |

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

These values are low-confidence placeholders for relative susceptibility to infection after effective exposure. Realised infection risk is driven heavily by contact patterns, household structure, prior RSV exposure, infant maternal antibody or nirsevimab protection, adult vaccination, and seasonal circulation.

| Age group | Suggested relative susceptibility | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 1.20 | Hall et al. 2013 review; CDC RSV infants and young children | Low | Young children have frequent first infections; split infants from 1-4 years if infant protection is modelled. |
| 5-9 | 1.00 | Hall et al. 2013 review | Low | Placeholder. |
| 10-14 | 0.95 | Hall et al. 2013 review | Low | Placeholder. |
| 15-19 | 0.95 | Hall et al. 2013 review | Low | Placeholder. |
| 20-29 | 1.00 | Hall et al. 2013 review | Low | Adult baseline placeholder. |
| 30-39 | 1.00 | Hall et al. 2013 review | Low | Adult baseline placeholder; parents and carers may have high exposure. |
| 40-49 | 1.00 | Hall et al. 2013 review | Low | Adult baseline placeholder. |
| 50-59 | 1.05 | CDC ACIP RSV adult evidence review | Low | Slight increase is a placeholder for accumulating comorbidity risk, not biological susceptibility. |
| 60-69 | 1.10 | CDC ACIP RSV adult evidence review | Low | Placeholder; vaccination status should override in current applied examples. |
| 70-79 | 1.15 | CDC ACIP RSV adult evidence review | Low | Placeholder. |
| 80+ | 1.20 | CDC ACIP RSV adult evidence review | Low | Placeholder. |

### 4.2 Infectiousness

No robust general-purpose age-specific infectiousness vector was identified for the default `agepi` age bands. Children are often important introducers and amplifiers in households, but that is partly contact behaviour rather than per-infection infectiousness. Use equal infectiousness by age for the minimal prototype unless fitting a specific setting.

| Age group | Suggested relative infectiousness | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 1.10 | Otomaru et al. 2019; CDC How RSV Spreads | Low | Placeholder reflecting high viral burden/prolonged shedding in some young children, but contact matrix should carry most age structure. |
| 5-9 | 1.05 | Otomaru et al. 2019 | Low | Placeholder. |
| 10-14 | 1.00 | Otomaru et al. 2019 | Low | Placeholder. |
| 15-19 | 1.00 | Otomaru et al. 2019 | Low | Placeholder. |
| 20-29 | 1.00 | Otomaru et al. 2019 | Low | Adult baseline. |
| 30-39 | 1.00 | Otomaru et al. 2019 | Low | Adult baseline. |
| 40-49 | 1.00 | Otomaru et al. 2019 | Low | Adult baseline. |
| 50-59 | 1.00 | Otomaru et al. 2019 | Low | Adult baseline. |
| 60-69 | 1.00 | Otomaru et al. 2019 | Low | Adult baseline. |
| 70-79 | 1.00 | Otomaru et al. 2019 | Low | Adult baseline. |
| 80+ | 1.00 | Otomaru et al. 2019 | Low | Adult baseline. |

### 4.3 Morbidity / severity

Outcome definition: probability of RSV-associated hospitalisation among infections. Values are illustrative age-pattern placeholders, not validated infection-hospitalisation ratios for a particular country or season.

| Age group | Suggested morbidity / severity risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.0200 | Hospitalisation among infections | Li et al. 2022; CDC RSV infants and young children | Low | Under-5 burden is dominated by infants; the broad 0-4 age group hides much higher infant risk and lower toddler/preschool risk. |
| 5-9 | 0.0010 | Hospitalisation among infections | CDC RSV burden methods; Hall et al. 2013 review | Low | Placeholder for generally lower severe RSV risk in school-age children. |
| 10-14 | 0.0005 | Hospitalisation among infections | CDC RSV burden methods; Hall et al. 2013 review | Low | Placeholder. |
| 15-19 | 0.0005 | Hospitalisation among infections | CDC RSV burden methods; Hall et al. 2013 review | Low | Placeholder. |
| 20-29 | 0.0010 | Hospitalisation among infections | CDC RSV burden methods; Hall et al. 2013 review | Low | Placeholder. |
| 30-39 | 0.0010 | Hospitalisation among infections | CDC RSV burden methods; Hall et al. 2013 review | Low | Placeholder. |
| 40-49 | 0.0015 | Hospitalisation among infections | CDC RSV burden methods; Hall et al. 2013 review | Low | Placeholder. |
| 50-59 | 0.0030 | Hospitalisation among infections | CDC ACIP RSV adult evidence review | Low | Placeholder; chronic heart/lung disease and immune compromise matter more than age alone. |
| 60-69 | 0.0060 | Hospitalisation among infections | CDC ACIP RSV adult evidence review | Low | Placeholder for higher risk in older adults. |
| 70-79 | 0.0100 | Hospitalisation among infections | CDC ACIP RSV adult evidence review | Low | Placeholder. |
| 80+ | 0.0150 | Hospitalisation among infections | CDC ACIP RSV adult evidence review | Low | Placeholder. |

### 4.4 Infection-induced mortality

Outcome definition: infection fatality risk among RSV infections. Values below are illustrative placeholders for demonstration burden outputs. They are not case fatality risks among hospitalised people and should not be used as country defaults.

| Age group | Suggested infection-induced mortality risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.00020 | Infection fatality risk | Li et al. 2022; CDC RSV infants and young children | Low | Global under-5 mortality risk varies strongly by setting; infant and low-resource risks can be much higher. |
| 5-9 | 0.000005 | Infection fatality risk | CDC RSV burden methods; Li et al. 2022 | Low | Placeholder. |
| 10-14 | 0.000003 | Infection fatality risk | CDC RSV burden methods | Low | Placeholder. |
| 15-19 | 0.000003 | Infection fatality risk | CDC RSV burden methods | Low | Placeholder. |
| 20-29 | 0.000005 | Infection fatality risk | CDC RSV burden methods | Low | Placeholder. |
| 30-39 | 0.000005 | Infection fatality risk | CDC RSV burden methods | Low | Placeholder. |
| 40-49 | 0.000010 | Infection fatality risk | CDC RSV burden methods | Low | Placeholder. |
| 50-59 | 0.000030 | Infection fatality risk | CDC ACIP RSV adult evidence review | Low | Placeholder. |
| 60-69 | 0.000080 | Infection fatality risk | CDC ACIP RSV adult evidence review | Low | Placeholder. |
| 70-79 | 0.000200 | Infection fatality risk | CDC ACIP RSV adult evidence review | Low | Placeholder. |
| 80+ | 0.000500 | Infection fatality risk | CDC ACIP RSV adult evidence review | Low | Placeholder; risk is much higher in frail, institutionalised, or medically complex adults. |

---

## 5. Transmission and mixing assumptions

| Component | Recommendation | Source | Confidence | Notes |
|---|---|---|---|---|
| Contact matrix | Use an external age-specific close-contact matrix with household and child-care/school contacts where possible | Otomaru et al. 2019; general contact-matrix literature | Medium | Household and child contacts are central for RSV spread; generic all-contact matrices may underrepresent infant-care contacts. |
| Transmission scaling | Calibrate beta to target \(R_0\) or seasonal attack rate | Reis and Shaman 2018; Otomaru et al. 2019 | High | Do not hard-code beta; solve or scale it for the chosen contact matrix, susceptibility vector, infectiousness vector, immunity assumptions, and target \(R_0\). |
| Age-specific mixing | Include strong household and young-child mixing if available | Otomaru et al. 2019; CDC RSV infants and young children | Medium | Young children often contribute to household transmission, but evidence is setting-specific. |
| Seasonality | Include as an optional sinusoidal or time-varying beta multiplier | CDC RSV surveillance; Reis and Shaman 2018 | High | RSV has marked seasonal epidemics in temperate settings; seasonality should be optional and setting-specific. |
| Intervention effects | Defer or use optional susceptibility/severity multipliers for vaccination, maternal vaccination, or nirsevimab | CDC RSV prevention guidance; CDC ACIP RSV adult evidence review | High | Prevention products primarily affect severe disease and may affect infection risk depending on product and endpoint. |
| Vaccination / passive immunisation relevance | Important optional extension; not part of minimal SEIRS | CDC RSV prevention guidance | High | Current examples should state whether infant nirsevimab, maternal vaccination, or adult RSV vaccination is included. |

---

## 6. Initial conditions for demonstration runs

| Quantity | Suggested default | Source / rationale | Confidence | Notes |
|---|---:|---|---|---|
| Initial infected proportion | 0.00001 | Demonstration seed | Low | Equivalent to 1 infectious person per 100,000 population. |
| Initial infected age groups | 0-4 or 5-9 for household/child-driven demo; 70-79 for severe-adult demo | Demonstration seed informed by RSV age-risk pattern | Low | Choose seed ages to match the example purpose. |
| Initial recovered / immune proportion | 0.50 for endemic seasonal demo; 0.00 for fully susceptible stress test | Hall et al. 2013 review; modelling simplification | Low | Most people have repeated RSV infections over life, but immunity is incomplete and hard to map to a simple recovered state. |
| Initial susceptible assumption | `1 - E0 - I0 - R0` within each age group | Mass-balance requirement | High | If using vaccination or nirsevimab, represent protection as susceptibility or severity modifiers rather than removing everyone from susceptibility. |
| Seeding approach | Seed one or a few exposed/infectious individuals, or use a seasonal importation pulse | Demonstration seed | Medium | For seasonal examples, repeated low-level seeding may be more realistic than a single importation. |

---

## 7. Parameters as package inputs

| Parameter | Type | Category | Suggested handling | Notes |
|---|---|---|---|---|
| beta | scalar or time-varying function | core | user-supplied or calibrated | Calibrate to target \(R_0\), seasonal attack rate, or observed incidence; no universal RSV beta. |
| gamma | scalar | core | disease default, user-overridable | Suggested 0.167 per day. |
| sigma | scalar | core for SEIR/SEIRS | disease default, user-overridable | Suggested 0.333 per day; lower confidence than incubation and infectious-period values. |
| R0 target | scalar | context_specific | user-supplied or example default | Suggested 2.5 with plausible range 1.2-4.0. |
| waning immunity rate | scalar | optional_extension | scenario default, user-overridable | Suggested `1 / 365` per day for toy SEIRS examples; low confidence. |
| susceptibility | age-specific vector | core | disease default, user-overridable | Low-confidence placeholders; prevention products and prior exposure should override. |
| infectiousness | age-specific vector | core | disease default, user-overridable | Low-confidence placeholders; contact matrix should carry most age structure. |
| morbidity risk | age-specific vector | risk_output | disease default, user-overridable | Hospitalisation-risk placeholders; strongly setting- and season-dependent. |
| mortality risk | age-specific vector | risk_output | disease default, user-overridable | Infection-fatality placeholders; replace for any applied burden analysis. |
| seasonal forcing amplitude and phase | scalar / time-varying function | optional_extension / context_specific | user-supplied | Important for realistic RSV, especially outside purely pedagogic examples. |
| maternal vaccination / nirsevimab / infant protection | age-specific modifier | optional_extension / context_specific | defer initially | Relevant mainly for infants; broad 0-4 age group is too coarse for detailed infant products. |
| adult RSV vaccination | age-specific modifier | optional_extension / context_specific | defer initially | Relevant for 60+, 75+, or risk-based groups depending on policy and year. |

---

## 8. Implementation implications for agepi

List implications for future examples or configuration only. These are not proposed changes to the core architecture.

| Observation | Possible implication | Priority | Notes |
|---|---|---|---|
| RSV reinfection is common and immunity wanes or is incomplete. | Use RSV as an example for optional SEIRS-style waning immunity once supported. | Medium | Avoid making waning immunity a universal core assumption. |
| Severe burden is concentrated in infants and older adults, but the default `0-4` age group is coarse. | Future examples may benefit from infant-specific age groups when modelling RSV severity or prevention products. | Medium | Keep the default age groups for generic prototypes unless the user opts into finer infant strata. |
| Seasonality is central to RSV epidemics. | Example scripts could allow beta seasonality as a generic time-varying transmission option. | Medium | This is a general feature, not RSV-specific code. |
| Prevention products are policy- and year-dependent. | Treat vaccination/nirsevimab as optional context-specific modifiers rather than disease defaults. | Low | Real-world recommendations have changed recently and should be checked before applied use. |

---

## 9. Cautions and limitations

- Transferability across countries is limited. RSV seasonality, testing, hospital admission thresholds, infant protection coverage, older-adult vaccine uptake, comorbidities, and health-care access vary widely.
- Temporal changes matter. Maternal RSV vaccination, infant nirsevimab, and adult RSV vaccination can change severity patterns across seasons.
- The broad `0-4` age group hides very high risk in infants, especially young infants, and lower risk in older preschool children.
- Biological susceptibility, realised infection risk, and severe-outcome risk are distinct. Contact matrices should carry much of the age-specific transmission pattern.
- The morbidity table uses illustrative infection-hospitalisation risks, not case severity among medically attended infections.
- The mortality table uses illustrative infection fatality risks. Death certification and attribution differ across surveillance systems.
- Waning immunity is represented as a toy SEIRS parameter; real RSV immunity is partial, endpoint-specific, and affected by repeated exposures.
- Published \(R_0\) estimates are highly model-dependent. Calibrate beta for the chosen setting instead of treating the suggested value as universal.

---

## 10. Reference table

| Claim / parameter | Source | Source type | Confidence | Notes |
|---|---|---|---|---|
| People with RSV are usually contagious for 3-8 days and may become contagious 1-2 days before symptoms | CDC. "How RSV Spreads." https://www.cdc.gov/rsv/causes/index.html | Public health agency | High | Used for infectious period and latent-period rationale. |
| Symptoms usually appear 4-6 days after infection; most RSV infections are mild but infants can develop severe disease | CDC. "RSV in Infants and Young Children." https://www.cdc.gov/rsv/infants-young-children/ | Public health agency | High | Used for incubation and infant severity rationale. |
| RSV chapter summarises incubation, communicability, seasonality, and prevention | CDC. "Epidemiology and Prevention of Vaccine-Preventable Diseases: RSV." 2025 Pink Book chapter. https://www2.cdc.gov/vaccines/ed/pinkbook/2025/PB_RSV/PB_RSV.pdf | Public health agency textbook | High | Used for incubation, communicability, and seasonality cross-checks. |
| RSV causes substantial hospitalisations and deaths in older adults; adults 65+ have higher rates than adults 60-64 | CDC ACIP. "Evidence to Recommendations for Use of GSK Adjuvanted RSVPreF3 Vaccine in Adults Ages 60 and Older." https://www.cdc.gov/acip/evidence-to-recommendations/gsk-adjuvanted-rsvpref3-adults-etr.html | Public health agency evidence review | High | Used for older-adult morbidity and mortality risk pattern. |
| CDC burden estimates are model-based and account for under-detection and out-of-hospital deaths | CDC. "How CDC Estimates the Burden of RSV in the US." https://www.cdc.gov/rsv/php/surveillance/about-burden-estimates.html | Public health agency methods | Medium | Used to caution that burden outputs depend on surveillance methods. |
| Global 2019 under-5 RSV burden includes tens of millions of RSV acute lower respiratory infection episodes and millions of hospital admissions | Li Y, Wang X, Blau DM, et al. "Global, regional, and national disease burden estimates of acute lower respiratory infections due to respiratory syncytial virus in children younger than 5 years in 2019: a systematic analysis." Lancet. 2022;399:2047-2064. doi:10.1016/S0140-6736(22)00478-0 | Systematic analysis | High | Used for under-5 severity and mortality pattern. |
| RSV incubation-period estimates around 4.4 days and uncertainty in respiratory-virus incubation periods | Lessler J, Reich NG, Brookmeyer R, Perl TM, Nelson KE, Cummings DAT. "Incubation periods of acute respiratory viral infections: a systematic review." Lancet Infect Dis. 2009;9(5):291-300. doi:10.1016/S1473-3099(09)70069-6 | Systematic review | Medium | Used for incubation cross-check. |
| Reinfection occurs and immunity is incomplete; RSV can recur throughout life | Hall CB, Walsh EE, Long CE, Schnabel KC. "Immunity to and frequency of reinfection with respiratory syncytial virus." J Infect Dis. 1991;163(4):693-698. doi:10.1093/infdis/163.4.693; Hall CB, et al. RSV review, N Engl J Med. 2013;368:588-598. doi:10.1056/NEJMra0804877 | Cohort study / review | Medium | Used for waning-immunity and SEIRS rationale. |
| RSV model inferred mean infection duration near 5.2 days and mean \(R_0\) near 2.8 in US regional simulations | Reis J, Shaman J. "Simulation of four respiratory viruses and inference of epidemiological parameters." Infect Dis Model. 2018;3:23-34. doi:10.1016/j.idm.2018.03.001 | Modelling study | Low | Used for \(R_0\) and infectious-duration plausibility only. |
| Household RSV study used a 3.2-day serial interval and estimated setting-specific transmission dynamics in children under 5 | Otomaru H, Kamigaki T, Tamaki R, et al. "Transmission of Respiratory Syncytial Virus Among Children Under 5 Years in Households of Rural Communities, the Philippines." Open Forum Infect Dis. 2019;6(3):ofz045. doi:10.1093/ofid/ofz045 | Household transmission study | Medium | Used for serial interval and household transmission rationale. |
| Rapid review summarises RSV transmission parameters including incubation, serial interval, infectiousness, and \(R_0\) | Usher Institute, University of Edinburgh. "What are the parameters (attack rates, incubation period, serial interval, R0, generation time, etc.) and modes of transmission of RSV?" 2025 rapid review. https://usher.ed.ac.uk/sites/default/files/atoms/files/rsv_transmission_parameters_rr_-_version_3_formatted.pdf | Rapid review | Medium | Used as a recent parameter cross-check. |

---

## 11. Machine-readable parameter draft

This draft is intended as a provisional R-style list for later copying into an example or disease-parameter file. Use `NA_real_` where values are unavailable, uncertain, or deliberately deferred.

```r
disease_parameters <- list(
  disease = "rsv",
  model = "SEIRS",
  natural_history = list(
    latent_period_days = 3,
    infectious_period_days = 6,
    incubation_period_days = 4.5,
    gamma = 1 / 6,
    sigma = 1 / 3,
    R0 = 2.5,
    serial_interval_days = 3.2,
    immunity_duration_days = 365,
    omega = 1 / 365
  ),
  age_specific = list(
    age_groups = c("0-4", "5-9", "10-14", "15-19", "20-29",
                   "30-39", "40-49", "50-59", "60-69",
                   "70-79", "80+"),
    susceptibility = c(1.20, 1.00, 0.95, 0.95, 1.00,
                       1.00, 1.00, 1.05, 1.10, 1.15, 1.20),
    infectiousness = c(1.10, 1.05, 1.00, 1.00, 1.00,
                       1.00, 1.00, 1.00, 1.00, 1.00, 1.00),
    morbidity_risk = c(0.0200, 0.0010, 0.0005, 0.0005, 0.0010,
                       0.0010, 0.0015, 0.0030, 0.0060, 0.0100, 0.0150),
    mortality_risk = c(0.000200, 0.000005, 0.000003, 0.000003, 0.000005,
                       0.000005, 0.000010, 0.000030, 0.000080, 0.000200,
                       0.000500)
  ),
  parameter_classification = list(
    core = c("beta", "gamma", "sigma", "susceptibility", "infectiousness"),
    risk_output = c("morbidity_risk", "mortality_risk"),
    optional_extension = c("waning_immunity_rate", "seasonal_forcing",
                           "maternal_vaccination", "nirsevimab",
                           "adult_rsv_vaccination"),
    context_specific = c("R0_target", "beta", "initial_immunity",
                         "seasonality_phase", "prevention_product_coverage")
  ),
  notes = list(
    source_summary = "CDC clinical and surveillance sources for timing, spread, and severe-risk groups; Li et al. 2022 for global under-5 burden; Lessler et al. 2009 for incubation; Hall et al. for reinfection and immunity; Reis and Shaman 2018 and Otomaru et al. 2019 for modelling and household transmission parameters.",
    cautions = "Age-specific susceptibility, infectiousness, hospitalisation risk, and infection fatality vectors are illustrative placeholders. Replace them for applied work, especially where infant protection, adult vaccination, health-care access, or local RSV seasonality differ."
  )
)
```
