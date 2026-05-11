# Disease parameterisation note: Influenza

## 1. Purpose

This note summarises plausible parameter values for using seasonal influenza as an example disease in the `agepi` age-structured epidemic modelling prototype.

The purpose is to provide realistic values for testing and demonstration. This note does not define the core `agepi` model architecture and should not be treated as a definitive influenza transmission model.

---

## 2. Recommended minimal model structure

| Feature | Recommendation | Rationale | Source | Confidence |
|---|---|---|---|---|
| Minimal compartmental structure | SEIR for a single-season outbreak; SIRS/SEIRS as an optional extension | Influenza has a short incubation period, infectiousness can begin before symptoms, and population immunity changes through infection, vaccination, and antigenic drift. | WHO seasonal influenza fact sheet; CDC Pink Book influenza chapter; CDC Yellow Book influenza chapter | High |
| Optional extensions | Vaccination, waning or strain-specific immunity, asymptomatic infection, antiviral treatment, age-specific severity, seasonal forcing, importation | These features matter for realistic seasonal influenza burden, but should remain optional example layers rather than core package assumptions. | WHO seasonal influenza fact sheet; CDC influenza burden methods; Biggerstaff et al. 2014 | High |
| Initial prototype suitability | High | Influenza is a useful age-structured example because infection rates are often high in children while severe outcomes concentrate in young children, older adults, pregnant people, and people with underlying conditions. | WHO influenza surveillance standards; CDC FluSurv-NET 2010-2023 report; Tokars et al. 2018 | High |

Notes:

- A minimal SEIR model is preferable to SIR if the prototype can represent an exposed class, because the latent/incubation delay is short but epidemiologically meaningful.
- For a single illustrative epidemic, immunity can be represented through initial recovered/immune proportions or susceptibility multipliers rather than a full multi-strain immune-history model.
- For applied seasonal influenza work, subtype, season, vaccination coverage, prior immunity, and surveillance ascertainment should be specified explicitly.

---

## 3. Core natural-history parameters

| Parameter | Symbol | Suggested value | Plausible range | Unit | Distribution / uncertainty | Source | Confidence | Notes |
|---|---:|---:|---:|---|---|---|---|---|
| Latent period | \(1 / \sigma\) | 1.5 | 1-2 | days | Fixed default or gamma-distributed waiting time | CDC How Flu Spreads; CDC surveillance manual; Cowling et al. 2010 | Medium | Infectiousness may begin about 1 day before symptoms, while incubation averages about 2 days. |
| Infectious period | \(1 / \gamma\) | 4 | 3-7 | days | Fixed default or gamma-distributed waiting time | CDC How Flu Spreads; CDC Yellow Book; CDC Pink Book influenza chapter | Medium | CDC describes infectiousness from about 1 day before symptoms to about 5-7 days after symptom onset; peak shedding is early. |
| Incubation period |  | 2 | 1-4 | days | Use 2 days as central value | WHO seasonal influenza fact sheet; CDC How Flu Spreads; CDC Pink Book influenza chapter | High | WHO and CDC describe a usual incubation period around 2 days, with a 1-4 day range. |
| Recovery rate | \(\gamma\) | 0.250 | 0.143-0.333 | per day | Derived as \(1 / infectious_period_days\) | Derived from infectious-period sources above | Medium | Use `1 / 4` for the suggested value. |
| Progression rate | \(\sigma\) | 0.667 | 0.500-1.000 | per day | Derived as \(1 / latent_period_days\) | Derived from latent-period assumption above | Medium | Use `1 / 1.5`; align with the short delay from infection to infectiousness. |
| Basic reproduction number | \(R_0\) | 1.5 | 1.2-2.0 | dimensionless | Context-specific; calibrate beta to contact matrix and immunity setting | Biggerstaff et al. 2014; WHO pandemic influenza risk management guidance | Medium | Seasonal and pandemic influenza estimates vary by subtype, setting, contact structure, and immunity. |
| Serial interval / generation-time proxy |  | 2.6 | 1.5-4.0 | days | Gamma or lognormal distribution in stochastic examples | Cowling et al. 2010; Biggerstaff et al. 2014 | Medium | Useful for checking epidemic speed, not directly a compartment duration. |
| Duration of immunity |  | NA | months to years, strain-dependent | days / years | Defer unless using SIRS or multi-strain examples | WHO seasonal influenza fact sheet; CDC vaccine effectiveness and burden methods | Low | Immunity depends on prior infection, vaccination, antigenic drift, subtype, and outcome. Do not use a universal duration for the minimal prototype. |

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

These values are low-confidence placeholders for relative susceptibility to infection after effective exposure in a seasonal influenza demonstration. Realised infection risk is driven heavily by contact patterns, prior immunity, vaccination, circulating subtype, and school or household mixing.

| Age group | Suggested relative susceptibility | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 1.40 | WHO influenza surveillance standards; Tokars et al. 2018; Cauchemez et al. 2009 | Low | Children have higher annual attack rates and household susceptibility than adults in several studies; vaccination and prior immunity should override. |
| 5-9 | 1.60 | WHO influenza surveillance standards; Tokars et al. 2018; Cauchemez et al. 2009 | Low | School-age children often have high symptomatic incidence and intense contact mixing. |
| 10-14 | 1.50 | WHO influenza surveillance standards; Tokars et al. 2018 | Low | Placeholder for elevated child/adolescent infection risk. |
| 15-19 | 1.25 | Tokars et al. 2018; Cauchemez et al. 2009 | Low | Placeholder between school-age children and adult baseline. |
| 20-29 | 1.00 | Tokars et al. 2018 | Low | Adult baseline. |
| 30-39 | 1.00 | Tokars et al. 2018 | Low | Adult baseline; parent/carer exposure may be higher in household models. |
| 40-49 | 0.95 | Tokars et al. 2018 | Low | Placeholder. |
| 50-59 | 0.90 | Cauchemez et al. 2009; Tokars et al. 2018 | Low | Older adults may have lower infection risk in some seasons because of prior immunity and lower contact rates. |
| 60-69 | 0.85 | Cauchemez et al. 2009; Tokars et al. 2018 | Low | Placeholder; do not confuse lower infection risk with lower severity. |
| 70-79 | 0.85 | Cauchemez et al. 2009; Tokars et al. 2018 | Low | Placeholder. |
| 80+ | 0.85 | Cauchemez et al. 2009; Tokars et al. 2018 | Low | Placeholder; long-term-care exposure can change this substantially. |

### 4.2 Infectiousness

No robust general-purpose age-specific infectiousness vector was identified for the default `agepi` age bands. Children can shed virus for longer and are often important in household and school spread, but much of this effect belongs in the contact matrix and duration assumptions.

| Age group | Suggested relative infectiousness | Source | Confidence | Notes |
|---|---:|---|---|---|
| 0-4 | 1.20 | CDC Pink Book influenza chapter; CDC How Flu Spreads; Cauchemez et al. 2009 | Low | Placeholder reflecting longer shedding in young children and high household contribution. |
| 5-9 | 1.15 | CDC Pink Book influenza chapter; Cauchemez et al. 2009 | Low | Placeholder. |
| 10-14 | 1.10 | CDC Pink Book influenza chapter; Cauchemez et al. 2009 | Low | Placeholder. |
| 15-19 | 1.05 | CDC Pink Book influenza chapter | Low | Placeholder. |
| 20-29 | 1.00 | CDC Pink Book influenza chapter | Low | Adult baseline. |
| 30-39 | 1.00 | CDC Pink Book influenza chapter | Low | Adult baseline. |
| 40-49 | 1.00 | CDC Pink Book influenza chapter | Low | Adult baseline. |
| 50-59 | 1.00 | CDC Pink Book influenza chapter | Low | Adult baseline. |
| 60-69 | 0.95 | CDC Pink Book influenza chapter | Low | Placeholder; lower contact rates, not necessarily lower biological infectiousness. |
| 70-79 | 0.95 | CDC Pink Book influenza chapter | Low | Placeholder. |
| 80+ | 0.95 | CDC Pink Book influenza chapter | Low | Placeholder; institutional settings may reverse this assumption. |

### 4.3 Morbidity / severity

Outcome definition: probability of influenza-associated hospitalisation among symptomatic influenza illnesses. Values are illustrative demonstration inputs derived from CDC 2022-2023 burden estimates by dividing estimated hospitalisations by estimated symptomatic illnesses within broad age groups, then mapping to the default `agepi` bands.

| Age group | Suggested morbidity / severity risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.0070 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates; CDC FluSurv-NET 2010-2023 report | Medium | CDC 2022-2023 estimate: about 20,772 hospitalisations / 2,979,491 symptomatic illnesses for ages 0-4. |
| 5-9 | 0.0015 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates; CDC FluSurv-NET 2010-2023 report | Low | Mapped from broader 5-17 age group; risk is generally low in school-age children. |
| 10-14 | 0.0015 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates; CDC FluSurv-NET 2010-2023 report | Low | Mapped from broader 5-17 age group. |
| 15-19 | 0.0020 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Transitional placeholder between adolescent and adult burden groups. |
| 20-29 | 0.0030 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broad 18-49 group. |
| 30-39 | 0.0030 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broad 18-49 group. |
| 40-49 | 0.0040 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Placeholder for increasing comorbidity with age. |
| 50-59 | 0.0080 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates; CDC FluSurv-NET 2010-2023 report | Low | Mapped from broad 50-64 group. |
| 60-69 | 0.0150 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates; CDC FluSurv-NET 2010-2023 report | Low | Transitional placeholder between 50-64 and 65+ estimates. |
| 70-79 | 0.0350 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates; CDC FluSurv-NET 2010-2023 report | Low | Mapped from broad 65+ group; risk rises strongly in older adults. |
| 80+ | 0.0600 | Hospitalisation among symptomatic illnesses | CDC 2022-2023 influenza burden estimates; CDC FluSurv-NET 2010-2023 report | Low | Placeholder for very old adults; replace with setting-specific FluSurv-NET or burden data where possible. |

### 4.4 Infection-induced mortality

Outcome definition: influenza-associated death risk among symptomatic influenza illnesses. These are not infection fatality risks for all infections, because asymptomatic and unrecognised infections are not included in the denominator. Values are illustrative high-income-setting placeholders based primarily on CDC 2022-2023 burden estimates.

| Age group | Suggested infection-induced mortality risk | Outcome definition | Source | Confidence | Notes |
|---|---:|---|---|---|---|
| 0-4 | 0.000066 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Medium | CDC 2022-2023 estimate: about 198 deaths / 2,979,491 symptomatic illnesses for ages 0-4. |
| 5-9 | 0.000005 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broader 5-17 age group; deaths are rare but not zero. |
| 10-14 | 0.000005 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broader 5-17 age group. |
| 15-19 | 0.000010 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Transitional placeholder. |
| 20-29 | 0.000020 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broad 18-49 group. |
| 30-39 | 0.000030 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broad 18-49 group. |
| 40-49 | 0.000060 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Placeholder for increasing risk with age and comorbidity. |
| 50-59 | 0.000200 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broad 50-64 group. |
| 60-69 | 0.000600 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Transitional placeholder between 50-64 and 65+ groups. |
| 70-79 | 0.002000 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Mapped from broad 65+ group. |
| 80+ | 0.006000 | Death among symptomatic illnesses | CDC 2022-2023 influenza burden estimates | Low | Placeholder for oldest adults; long-term-care residence and comorbidities matter greatly. |

---

## 5. Transmission and mixing assumptions

| Component | Recommendation | Source | Confidence | Notes |
|---|---|---|---|---|
| Contact matrix | Use an external age-specific close-contact matrix with household, school, workplace, and community layers where possible | WHO influenza surveillance standards; Cauchemez et al. 2009; general contact-matrix literature | Medium | Child and school contacts are central to many influenza epidemics. |
| Transmission scaling | Calibrate beta to target \(R_0\), observed growth rate, or seasonal attack rate | Biggerstaff et al. 2014; CDC burden methods | High | Do not hard-code beta; solve or scale it for the chosen contact matrix, susceptibility vector, infectiousness vector, immunity assumptions, and target \(R_0\). |
| Age-specific mixing | Include strong school-age and household mixing when available | Cauchemez et al. 2009; WHO influenza surveillance standards | Medium | Age-specific incidence partly reflects contact patterns rather than biological susceptibility. |
| Seasonality | Include as an optional sinusoidal, school-term, climate, or empirical beta multiplier | WHO seasonal influenza fact sheet; CDC Yellow Book | High | Influenza seasonality is strong in temperate settings and less regular in tropical settings. |
| Intervention effects | Defer or use optional beta, susceptibility, infectiousness, or severity multipliers for vaccination, antivirals, isolation, or school closure | WHO seasonal influenza fact sheet; CDC burden prevented methods | Medium | Interventions are policy-, season-, and coverage-dependent. |
| Vaccination relevance | Important optional extension; not part of minimal SEIR | WHO seasonal influenza fact sheet; CDC vaccine effectiveness and burden methods | High | Vaccination can affect infection, symptomatic disease, hospitalisation, and death endpoints differently. |

---

## 6. Initial conditions for demonstration runs

| Quantity | Suggested default | Source / rationale | Confidence | Notes |
|---|---:|---|---|---|
| Initial infected proportion | 0.00001 | Demonstration seed | Low | Equivalent to 1 infectious person per 100,000 population. |
| Initial infected age groups | 5-9 or 10-14 for school-driven demo; 20-49 for workplace/community demo | Demonstration seed informed by influenza contact patterns | Low | Choose seed ages to match the example purpose. |
| Initial recovered / immune proportion | 0.30 for a seasonal demo; 0.00 for fully susceptible stress test | Modelling simplification; WHO influenza epidemiology | Low | Real immunity is subtype- and history-dependent. A single recovered compartment is an approximation. |
| Initial susceptible assumption | `1 - E0 - I0 - R0` within each age group | Mass-balance requirement | High | If using vaccination or prior immunity, represent protection as age-specific initial immunity or susceptibility modifiers. |
| Seeding approach | Seed one or a few exposed/infectious individuals, or use low-level seasonal importation | Demonstration seed | Medium | Repeated importation can be more realistic for seasonal examples than a single seed. |

---

## 7. Parameters as package inputs

| Parameter | Type | Category | Suggested handling | Notes |
|---|---|---|---|---|
| beta | scalar or time-varying function | core | user-supplied or calibrated | Calibrate to target \(R_0\), observed growth, or attack rate; no universal influenza beta. |
| gamma | scalar | core | disease default, user-overridable | Suggested 0.250 per day. |
| sigma | scalar | core for SEIR/SEIRS | disease default, user-overridable | Suggested 0.667 per day. |
| R0 target | scalar | context_specific | user-supplied or example default | Suggested 1.5 with plausible range 1.2-2.0. |
| susceptibility | age-specific vector | core | disease default, user-overridable | Low-confidence placeholders; prior immunity and vaccination should override. |
| infectiousness | age-specific vector | core | disease default, user-overridable | Low-confidence placeholders; contact matrix should carry most age structure. |
| morbidity risk | age-specific vector | risk_output | disease default, user-overridable | Hospitalisation among symptomatic illnesses; use setting-specific data where possible. |
| mortality risk | age-specific vector | risk_output | disease default, user-overridable | Deaths among symptomatic illnesses; use setting-specific data where possible. |
| seasonal forcing amplitude and phase | scalar / time-varying function | optional_extension / context_specific | user-supplied | Important for realistic influenza seasonality. |
| vaccination coverage and effect | age-specific / time-varying modifiers | optional_extension / context_specific | defer initially | Product, season, subtype match, and endpoint-specific effectiveness matter. |
| waning or strain-specific immunity | scalar / matrix / history state | optional_extension / context_specific | defer initially | Useful later for SIRS or multi-strain examples; not required for a first prototype. |

---

## 8. Implementation implications for agepi

List implications for future examples or configuration only. These are not proposed changes to the core architecture.

| Observation | Possible implication | Priority | Notes |
|---|---|---|---|
| Influenza has a short exposed phase and rapid generation time. | Use influenza as a compact SEIR example for testing fast epidemic dynamics. | Medium | Time-step choices should be checked if using discrete-time simulations. |
| Severity rises steeply in older adults while infection incidence is often high in children. | Influenza is useful for demonstrating separation between transmission outputs and burden outputs. | Medium | Keep severity vectors as risk-output parameters rather than core transmission assumptions. |
| Seasonality and vaccination are central in real applications. | Future examples may include generic time-varying beta and age-specific protection modifiers. | Medium | These should remain generic features, not influenza-specific package architecture. |
| A single recovered compartment is a crude representation of influenza immunity. | Document whether examples are single-season, single-strain, or toy SIRS scenarios. | Low | Multi-strain immune history is outside the minimal prototype. |

---

## 9. Cautions and limitations

- Transferability across countries is limited. Influenza subtype mix, vaccine coverage, health-care access, testing practice, school calendars, demographics, and climate vary widely.
- Temporal changes matter. Annual vaccine composition, antigenic drift, prior immunity, and post-pandemic behavioural changes can shift attack rates and severity.
- The broad `0-4` and `80+` groups hide important heterogeneity, especially infants, people with pregnancy, long-term-care residents, and people with chronic conditions.
- Biological susceptibility, realised infection risk, symptomatic illness risk, and severe-outcome risk are distinct. Contact matrices should carry much of the age-specific transmission pattern.
- The morbidity table uses hospitalisation among symptomatic illnesses, not infection-hospitalisation risk among all infections.
- The mortality table uses death risk among symptomatic illnesses, not infection fatality risk among all infections.
- Age-specific susceptibility and infectiousness values are illustrative placeholders and should not be treated as measured universal influenza biology.
- Published \(R_0\) estimates are model-dependent. Calibrate beta for the chosen setting instead of treating the suggested value as universal.

---

## 10. Reference table

| Claim / parameter | Source | Source type | Confidence | Notes |
|---|---|---|---|---|
| Seasonal influenza incubation period is about 2 days, range 1-4 days; annual attack rates are higher in children than adults | WHO. "Influenza (seasonal)." https://www.who.int/news-room/fact-sheets/detail/influenza-(seasonal) | Public health agency | High | Used for incubation, epidemiology, and vaccination context. |
| Influenza incubation period is 1-4 days; annual attack rates are about 5-10% in adults and 20-30% in children | WHO. "Influenza: Vaccine Preventable Diseases Surveillance Standards." https://www.who.int/publications/m/item/vaccine-preventable-diseases-surveillance-standards-influenza | Public health agency surveillance standard | High | Used for age attack-rate pattern. |
| Influenza viruses can be detected from about 1 day before symptoms through 5-7 days after illness onset; children and immunocompromised people may be contagious longer | CDC. "How Flu Spreads." https://www.cdc.gov/flu/spread/index.html | Public health agency | High | Used for latent and infectious-period rationale. |
| Incubation period ranges 1-4 days; peak shedding usually occurs from 1 day before onset to 3 days after; acute symptoms often last 2-7 days | CDC. "Chapter 6: Influenza." Manual for the Surveillance of Vaccine-Preventable Diseases. https://www.cdc.gov/surv-manual/php/table-of-contents/chapter-6-influenza.html | Public health agency surveillance manual | High | Used for incubation, infectiousness timing, and symptoms. |
| Influenza virus shedding is 5-10 days, peak 1-3 days after illness onset; adults can transmit from day before symptoms to 5-7 days after | CDC. "Chapter 12: Influenza." Pink Book. https://www.cdc.gov/pinkbook/hcp/table-of-contents/chapter-12-influenza.html | Public health agency textbook | High | Used for infectious period and prolonged shedding cautions. |
| Adults with influenza are usually infectious from day before symptom onset to about 5-7 days after symptom onset | CDC. "Influenza." Yellow Book. https://www.cdc.gov/yellow-book/hcp/travel-associated-infections-diseases/influenza.html | Public health agency travel medicine manual | High | Used for infectious period and travel/seasonality context. |
| Reproduction-number estimates for seasonal and pandemic influenza are generally modest but vary by setting and method | Biggerstaff M, Cauchemez S, Reed C, Gambhir M, Finelli L. "Estimates of the reproduction number for seasonal, pandemic, and zoonotic influenza: a systematic review of the literature." BMC Infect Dis. 2014;14:480. doi:10.1186/1471-2334-14-480 | Systematic review | Medium | Used for \(R_0\) range and caution. |
| Serial interval and generation timing for influenza are short; infectiousness is concentrated around symptom onset | Cowling BJ, Fang VJ, Riley S, Peiris JSM, Leung GM. "Estimation of the serial interval of influenza." Epidemiology. 2009/2010; often cited estimate around 2.6 days. doi:10.1097/EDE.0b013e3181b91f06 | Peer-reviewed transmission study | Medium | Used for serial interval and latent-period plausibility. |
| Children have higher symptomatic influenza incidence than adults | Tokars JI, Olsen SJ, Reed C. "Seasonal Incidence of Symptomatic Influenza in the United States." Clin Infect Dis. 2018;66(10):1511-1518. doi:10.1093/cid/cix1060 | Peer-reviewed incidence analysis | Medium | Used for susceptibility placeholder pattern. |
| Household contacts aged 18 years or younger were estimated to be more susceptible than adults aged 19-50 in 2009 H1N1 household data | Cauchemez S, Donnelly CA, Reed C, et al. "Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States." N Engl J Med. 2009;361:2619-2627. doi:10.1056/NEJMoa0905498 | Household transmission study | Medium | Used for age susceptibility and child transmission rationale; pandemic-specific. |
| US influenza burden estimates provide age-specific symptomatic illnesses, medical visits, hospitalisations, and deaths | CDC. "Preliminary Estimated Flu Disease Burden 2022-2023 Flu Season." https://www.cdc.gov/flu-burden/php/data-vis/2022-2023.html | Public health agency burden estimates | Medium | Used for illustrative morbidity and mortality risks. |
| Adults aged 65+ consistently have the highest influenza-associated hospitalisation rates, followed in most seasons by children aged 0-4 | CDC. "Laboratory-Confirmed Influenza-Associated Hospitalizations Among Children and Adults - FluSurv-NET, United States, 2010-2023." MMWR Surveill Summ. 2024;73(6):1-17. https://www.cdc.gov/mmwr/volumes/73/ss/ss7306a1.htm | Public health agency surveillance report | High | Used for age severity pattern. |
| CDC burden methods adjust surveillance hospitalisations for under-detection and estimate symptomatic illnesses, hospitalisations, and deaths | CDC. "How CDC Estimates the Burden of Seasonal Flu in the United States." https://www.cdc.gov/flu-burden/php/about/index.html | Public health agency methods | Medium | Used for caution about burden-derived risks. |

---

## 11. Machine-readable parameter draft

This draft is intended as a provisional R-style list for later copying into an example or disease-parameter file. Use `NA_real_` where values are unavailable, uncertain, or deliberately deferred.

```r
disease_parameters <- list(
  disease = "influenza",
  model = "SEIR",
  natural_history = list(
    latent_period_days = 1.5,
    infectious_period_days = 4,
    incubation_period_days = 2,
    gamma = 0.25,
    sigma = 0.667,
    R0 = 1.5
  ),
  age_specific = list(
    age_groups = c("0-4", "5-9", "10-14", "15-19", "20-29",
                   "30-39", "40-49", "50-59", "60-69",
                   "70-79", "80+"),
    susceptibility = c(1.40, 1.60, 1.50, 1.25, 1.00,
                       1.00, 0.95, 0.90, 0.85, 0.85, 0.85),
    infectiousness = c(1.20, 1.15, 1.10, 1.05, 1.00,
                       1.00, 1.00, 1.00, 0.95, 0.95, 0.95),
    morbidity_risk = c(0.0070, 0.0015, 0.0015, 0.0020, 0.0030,
                       0.0030, 0.0040, 0.0080, 0.0150, 0.0350, 0.0600),
    mortality_risk = c(0.000066, 0.000005, 0.000005, 0.000010, 0.000020,
                       0.000030, 0.000060, 0.000200, 0.000600, 0.002000, 0.006000)
  ),
  parameter_classification = list(
    core = c("beta", "gamma", "sigma", "susceptibility", "infectiousness"),
    risk_output = c("morbidity_risk", "mortality_risk"),
    optional_extension = c("seasonal_forcing", "vaccination", "waning_immunity",
                           "asymptomatic_infection", "antiviral_treatment"),
    context_specific = c("R0", "beta", "initial_immunity", "vaccine_coverage",
                         "vaccine_effectiveness", "subtype_mix")
  ),
  notes = list(
    source_summary = "Core timing values are from WHO and CDC influenza guidance. R0 is from influenza reproduction-number reviews. Age-specific severity values are illustrative values derived mainly from CDC 2022-2023 burden estimates and broader FluSurv-NET age patterns.",
    cautions = "Age-specific susceptibility and infectiousness are low-confidence placeholders. Morbidity and mortality risks are based on symptomatic illness denominators and should be replaced for applied analyses."
  )
)
```
