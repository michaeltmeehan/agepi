# Kiribati TB modelling inputs and collaborator data requests

## 1. Executive summary

This note converts the Kiribati TB evidence inventory into actionable modelling
inputs for a future age-structured TB model in `agepi`. It is a documentation
and design note only. It does not implement model code, add dependencies, or
commit external datasets.

Main conclusion:

- A first public-data-only model is feasible as a national, age-structured,
  latent-active TB model calibrated to WHO national burden estimates,
  notifications, age/sex incidence estimates, treatment outcomes, MDR/RR-TB
  estimates, WPP demography, and WUENIC BCG coverage.
- A public-data-only model can support method development, broad sensitivity
  analysis, and transparent assumptions. It should not be treated as a
  definitive Kiribati policy model because public data do not identify
  subnational transmission, latent infection prevalence, contact-tracing
  cascades, PEARL intervention effects, or treatment-pathway timing.
- A collaborator-enhanced model becomes substantially more useful if Kiribati
  NTP, MHMS, PEARL, or other collaborators can provide de-identified
  notification, contact tracing, TPT, treatment, drug-susceptibility,
  subnational denominator, and intervention activity data.

## 2. Public-data-only modelling pathway

The public-data-only pathway should be deliberately modest. It can establish a
reproducible national modelling scaffold and identify which posterior
uncertainties are driven by missing local data.

Sufficient public inputs for a first model:

| Input | Public source | Use | Sufficient for first model? | Caveat |
|---|---|---|---|---|
| National population by age/sex | UN WPP 2024 | Initial population and demographic schedules | Yes | National only; subnational population needs census/NSO tables. |
| Fertility, mortality, migration | UN WPP 2024 | Background demographic process | Yes | Use WPP-style schedules; not projection matching. |
| TB incidence and mortality | WHO Global TB Database burden estimates | Core calibration targets | Yes | WHO estimates have uncertainty and depend on surveillance adjustment. |
| TB notifications | WHO Global TB Database case notifications / burden CSV | Diagnosis/reporting target | Yes | Detection changes may reflect active case finding rather than incidence. |
| Age/sex incidence | WHO age/sex/risk-factor incidence estimates | Age pattern target | Yes | Age bins need mapping to model age groups. |
| Treatment outcomes | WHO treatment outcome CSV | Treatment success/death/loss constraints | Partly | Aggregate data do not give treatment timing. |
| MDR/RR-TB burden | WHO MDR/RR-TB estimates | Scenario bound or optional resistance state | Partly | Small counts and wide uncertainty; do not overfit. |
| Contact/TPT aggregate indicators | WHO contact/TPT CSV | Prevention cascade constraints | Partly | Denominator definitions and missingness need checking. |
| BCG coverage | WUENIC / WHO Immunization Data | Birth cohort BCG history | Yes for simple paediatric severe-TB modifier | Not direct protection against adult pulmonary TB. |
| Contact matrix | Synthetic/projection literature, if Kiribati available | Age mixing | Partly | No public empirical Kiribati contact matrix identified. |
| South Tarawa evidence | Census report/atlas, DFAT/UNDP/PEARL public narrative | Sensitivity/scenario motivation | Partly | Not enough for calibrated two-patch transmission. |

Recommended public-data-only product:

- national deterministic age-structured TB model;
- age groups chosen to align with WHO incidence and WPP data as much as possible;
- compartments for susceptible, latent infection, active TB, on treatment, and
  post-treatment/recovered;
- demography from WPP;
- a static synthetic age contact matrix or sensitivity set;
- calibration to national incidence, notifications, mortality, age/sex incidence
  distribution, and treatment outcomes;
- no calibrated South Tarawa/outer-island split;
- no PEARL effect estimate except as an external scenario sensitivity;
- MDR/RR-TB handled as a reporting outcome or scenario modifier unless enough
  DST data are obtained.

## 3. Collaborator-enhanced modelling pathway

Collaborator data would allow the model to move from a national scaffold to a
policy-facing Kiribati TB model.

Key additions with NTP/PEARL/collaborator access:

| Data access | New modelling capability | Why it changes the model |
|---|---|---|
| De-identified notification line list | Age, sex, time, site-of-disease, and geography-specific calibration | Separates age/geographic burden from national WHO estimates. |
| Subnational denominators and migration/mobility | South Tarawa/rest-of-Kiribati or island-group model | Allows explicit geographic strata and movement assumptions. |
| Contact tracing and TPT register | Household-contact infection and prevention module | Links index cases, contacts, screening, infection testing, TPT start, and completion. |
| PEARL screening and infection data | Calibration to South Tarawa active disease and TB infection prevalence | Provides contemporary intervention-denominator and latent infection evidence. |
| Treatment timing and outcomes | Diagnosis/treatment pathway and infectious duration | Converts treatment from a crude removal rate to a realistic pathway. |
| DST/lab register | MDR/RR-TB transmission or treatment-pathway module | Distinguishes true low resistance from low testing coverage. |
| Household/crowding/census microdata | TB-relevant mixing and crowding modifiers | Better represents prolonged indoor exposure than generic contact matrices. |
| Intervention activity logs | Time-varying case-finding and programme scenario inputs | Separates detection intensity from underlying incidence. |

## 4. Calibration target table

| Observed data stream | Source | Model output to compare | Scale | Stratification | Calibration role | Public-only or collaborator? | Notes |
|---|---|---|---|---|---|---|---|
| Estimated TB incidence count/rate | WHO burden estimates | New active TB disease episodes | Annual | National | Primary target | Public | Fit transmission/progression scale; include uncertainty. |
| Estimated TB mortality count/rate | WHO burden estimates | TB deaths among active/on-treatment states | Annual | National | Primary target | Public | Requires disease-induced mortality or case-fatality mapping. |
| TB notifications count/rate | WHO case notifications; NTP | Diagnosed/reported active TB | Annual/monthly if NTP | National public; subnational if NTP | Primary target | Public national; collaborator subnational | Detection rate and active case finding both influence this target. |
| Case detection/treatment coverage | WHO burden estimates | Diagnosed active TB / incident active TB | Annual | National | Secondary target | Public | Derived target; avoid double-counting if also fitting incidence and notifications. |
| Age/sex incidence distribution | WHO age/sex/risk-factor estimates | Incident active TB by age/sex | Annual | Age/sex national | Primary age-pattern target | Public | Mapping to model age groups is required. |
| Treatment success | WHO treatment outcomes; NTP | Completed/cured among treatment starts | Annual cohort | National public; age/geography if NTP | Secondary target | Public/collaborator | Treatment duration and outcome definitions must align. |
| Death during treatment | WHO treatment outcomes; NTP | Death among treatment cohort | Annual cohort | National public; age/geography if NTP | Secondary target | Public/collaborator | May overlap with mortality target. |
| Lost to follow-up/failure | WHO treatment outcomes; NTP | Treatment interruption/failure states | Annual cohort | National public; age/geography if NTP | Secondary target | Public/collaborator | Useful only if model has treatment outcome states. |
| MDR/RR-TB incidence/proportion | WHO MDR/RR estimates; DST register | RR/MDR incident cases or proportion among tested/notified | Annual | National public; age/geography if DST register | Scenario/optional target | Public/collaborator | Small numbers; use weakly unless detailed DST data exist. |
| Contacts screened/evaluated | WHO contact/TPT CSV; NTP | Contact investigation coverage | Annual | National public; age/geography if NTP | Scenario target | Public/collaborator | Requires denominator definition. |
| TPT starts/completions | WHO contact/TPT CSV; NTP/PEARL | Preventive treatment initiation/completion | Annual | National public; age/geography if NTP/PEARL | Scenario target | Public/collaborator | Completion is more useful than starts for effect. |
| BCG coverage | WUENIC | Vaccinated birth cohort proportion | Annual birth cohort | National | Fixed input / weak target | Public | Use only for paediatric severe TB protection scenarios. |
| TB infection prevalence / ARTI | PEARL TST/ARTI outputs | Latent infection prevalence or annual infection risk | Survey/intervention rounds | Age, school, South Tarawa | Primary latent/transmission target | Collaborator | Critical for identifying latent infection pool. |
| South Tarawa/rest notifications | NTP/PEARL | Diagnosed TB by geographic stratum | Annual/monthly | Geography, age if available | Primary spatial target | Collaborator | Needed for a calibrated two-patch model. |

## 5. Parameter-source table

| Parameter | Model role | Initial source | Public-data-only handling | Collaborator-enhanced handling | Priority |
|---|---|---|---|---|---|
| Initial population by age | Denominators and state sizes | UN WPP 2024 | Use national WPP age-sex population | Use census/NSO subnational denominators for spatial model | Critical |
| Fertility rate by age/time | Births into youngest age group | UN WPP 2024 | Fixed demographic input | Same; optionally reconcile with census/NSO | High |
| Mortality rate by age/time | Background deaths | UN WPP 2024 | Fixed demographic input | Same; optionally subnational sensitivity | High |
| Net migration | Demographic process | UN WPP 2024 | National net migration only | Subnational mobility/migration if available | Medium |
| Contact matrix | Force of infection | Synthetic/projected contact matrix | Use synthetic matrix and sensitivity analysis | Add household/crowding/PEARL contact-informed modifiers | Critical |
| Transmission coefficient | Force of infection scale | Calibrated | Fit to incidence/notifications | Fit by geography/time with intervention history | Critical |
| Relative susceptibility by age | Infection risk | TB literature; age-specific incidence target | Use literature prior and calibrate weakly | Re-estimate with line list and infection prevalence | High |
| Relative infectiousness by age | Contribution to transmission | TB literature | Use adult-high/child-low literature assumptions | Refine with pulmonary/bacteriological confirmation by age | High |
| Fast progression hazard | Recent infection to active TB | TB natural-history literature | Literature prior; not identifiable alone | Constrain with TST/ARTI/contact data | High |
| Reactivation hazard | Latent infection to active TB | TB natural-history literature | Literature prior; wide uncertainty | Refine with latent prevalence and age-specific incidence | High |
| Initial latent infection prevalence | Initial condition | Literature/WHO household-contact estimates | Scenario or calibrated latent pool with broad prior | Use PEARL TST/IGRA/ARTI data | Critical |
| Diagnosis rate | Active TB to treatment/reporting | WHO notifications vs incidence | Calibrate annually or by period | Calibrate by geography/time using NTP and ACF logs | Critical |
| Treatment initiation delay | Infectious duration before treatment | Literature / assumed | Fixed or broad prior | Estimate from symptom/diagnosis/treatment dates | High |
| Treatment success probability | Treatment outcome | WHO treatment outcomes | Fit aggregate outcome probabilities | Stratify by age, geography, resistance, treatment history | High |
| TB mortality hazard/CFR | Mortality output | WHO mortality estimates; treatment outcomes | Fit aggregate mortality/CFR | Stratify by age, treatment status, delay, resistance | High |
| Relapse/retreatment rate | Recurrent TB pathway | WHO notifications/treatment history; literature | Optional broad prior | Estimate from line list and treatment history | Medium |
| TPT initiation/completion | Prevention scenario | WHO contact/TPT CSV | Use national annual aggregate | Estimate by contact age, geography, regimen, exposure | Critical for prevention scenarios |
| TPT efficacy | Reduced progression from infection | TB prevention literature | Literature prior | Same, with local completion/adherence data | High |
| BCG coverage | Paediatric vaccine modifier | WUENIC | Fixed annual birth-cohort input | Add subnational EPI/survey data if available | Medium |
| BCG protection against severe childhood TB | Severity/mortality modifier | TB vaccine literature | Literature prior | Same; local severe paediatric TB can validate | Medium |
| MDR/RR proportion | Resistance module | WHO MDR/RR estimates | Scenario modifier | Estimate from DST/lab register | Medium |

## 6. Scenario-input table

| Scenario | Model mechanism | Required inputs | Public-data-only feasibility | Collaborator-enhanced inputs | Priority |
|---|---|---|---|---|---|
| Status quo continuation | Baseline diagnosis, treatment, demography, transmission | WHO burden/notifications/treatment, WPP demography | Feasible | NTP time series improves fit | Critical |
| Increased active case finding | Higher diagnosis rate; possibly earlier treatment | ACF dates, screened population, yield, target area | Only crude sensitivity from public narrative | ACF activity log by location/date/age and cases found | Critical |
| PEARL-like South Tarawa screen-and-treat | One-off or phased screening, treatment, TPT | Population enumerated, screening coverage, diagnostic algorithm, TB/TBI yield, TPT uptake/completion | Not calibratable; scenario placeholder only | PEARL aggregate or individual outputs | Critical |
| Expanded contact tracing | More contacts evaluated and active/TBI cases found | Contacts per index case, evaluation coverage, positivity, TPT start/completion | Crude national aggregate from WHO contact/TPT CSV | Contact tracing register linked to index cases | Critical |
| Expanded TPT for eligible contacts | Reduced progression among infected contacts | Eligibility, starts, completion, regimen, efficacy | Feasible as broad national sensitivity | Age/geography/regimen-specific cascade | High |
| Treatment strengthening | Shorter infectious duration, higher success, lower death/LTFU | Treatment success/death/LTFU; treatment delay | Feasible at national aggregate | Treatment dates, outcome, adherence, DOTS support by geography | High |
| Rapid diagnostic expansion | Shorter diagnosis delay and improved confirmation | Lab coverage, Xpert/culture/smear testing, positivity | Weak sensitivity only | Lab register and diagnostic access by time/geography | Medium |
| South Tarawa targeted control | Geographic transmission/detection shift | South Tarawa population, cases, intervention coverage | Not calibrated with public data | Subnational NTP/PEARL/census data | Critical |
| Outer-island service strengthening | Higher diagnosis/treatment access outside Tarawa | Outer-island notifications, service coverage, transfer/referral data | Not calibrated with public data | Outer-island line list and programme activity | High |
| MDR/RR-TB response | Separate resistance treatment pathway | RR/MDR incidence, DST coverage, treatment outcomes | Broad scenario only | DST/lab and DR-TB treatment data | Medium |
| BCG coverage change | Paediatric severe TB modifier | WUENIC birth-cohort coverage | Feasible for severity output only | Subnational EPI/survey coverage | Low/medium |

## 7. Data request checklist for collaborators

| Data item | Preferred format | Time period | Age stratification | Geographic stratification | Individual-level or aggregate? | Model use | Priority |
|---|---|---|---|---|---|---|---|
| TB notification line list | CSV with data dictionary; de-identified stable IDs | Ideally 2000-present; minimum 2010-present | Exact age or age in years; sex | Island, village, health facility, South Tarawa/rest | Individual preferred; aggregate acceptable for first pass | Incidence/notification calibration, age and geography patterns | Critical |
| Case classification and disease site | Columns in notification line list | Same as notifications | Exact age/age group | Same as notifications | Individual preferred | Infectious pulmonary fraction and paediatric infectiousness | Critical |
| Bacteriological confirmation and diagnostic test | CSV lab-linked table or line-list columns | 2010-present; include Xpert/culture introduction dates | Age/sex | Facility/island | Individual or aggregate by year/age/geography | Diagnostic pathway and infectiousness | High |
| Treatment dates and outcomes | CSV linked to notification ID | 2010-present | Exact age/age group | Facility/island | Individual preferred | Treatment rate, duration, success, death, LTFU | Critical |
| Relapse/retreatment history | Line-list columns or linked episode table | 2010-present | Exact age/age group | Island/village | Individual preferred | Recurrent TB and relapse/retreatment pathway | High |
| TB deaths | Line-list outcome plus vital registration if available | 2010-present | Exact age/age group | Island/village | Individual preferred; aggregate acceptable | Mortality calibration | High |
| Active case finding activities | CSV by activity/event | 2000-present if available; minimum 2018-present | Age of screened people if available | Activity location, island/village/hotspot | Aggregate by event acceptable | Detection-rate changes and ACF scenarios | Critical |
| Screening denominator and yield | CSV by campaign/round | Same as ACF/PEARL periods | Age group and sex | Location/round | Aggregate acceptable; individual better | Screen-and-treat scenarios | Critical |
| Contact tracing register | CSV linked to index case if possible | 2010-present; minimum recent 5 years | Contact exact age/age group; index age | Household/village/island | Individual preferred | Household transmission and contact cascade | Critical |
| TPT start/completion register | CSV linked to contact or aggregate cascade | 2010-present; minimum recent 5 years | Exact age/age group | Household/village/island | Individual preferred | TPT effectiveness and prevention scenarios | Critical |
| PEARL enumeration/screening outputs | De-identified individual CSV or aggregate tables | Full PEARL baseline/intervention/follow-up | Exact age/age group and sex | South Tarawa village/cluster if shareable | Aggregate may be enough; individual ideal | TB/TBI prevalence, ARTI, PEARL scenario | Critical |
| TST/IGRA/ARTI data | CSV with test dates/results | PEARL and any earlier surveys | Exact age/age group | School/village/island | Individual preferred; school aggregate acceptable | Latent infection and transmission calibration | Critical |
| Drug susceptibility testing | CSV lab register | 2010-present; minimum 2017-present | Age/sex | Facility/island | Individual preferred | MDR/RR pathway and DST coverage | High |
| DR-TB treatment outcomes | CSV linked to DST/patient ID | 2010-present | Age/sex | Facility/island | Individual preferred | Resistance scenario and mortality/treatment duration | Medium |
| Subnational population denominators | CSV table | 2000-present annual if possible; census years minimum | Five-year age/sex or exact age | Island/village/South Tarawa/rest | Aggregate | Geographic rates and two-patch model | Critical |
| Mobility/migration between islands | CSV or report tables | 2010-present; census migration windows | Age/sex if available | Origin-destination island | Aggregate | Spatial coupling | Medium/high |
| Household composition/crowding | Census aggregate tables or controlled microdata | 2020 census; earlier census if available | Age composition if possible | Village/island | Aggregate by household type acceptable | Household/crowding mixing modifier | High |
| Comorbidity summaries | Aggregate CSV | 2010-present or survey years | Age/sex | National/subnational if available | Aggregate preferred | Diabetes, smoking, HIV sensitivity analyses | Medium |

## 8. Minimum viable model specification

Purpose:

- produce a transparent national Kiribati TB calibration using public data only;
- test `agepi` age-structured demography/contact workflows for chronic TB-like
  dynamics;
- quantify which conclusions are sensitive to missing latent infection,
  contact, and subnational data.

Recommended structure:

| Component | Minimum viable choice | Rationale |
|---|---|---|
| Geography | National only | Public TB data are strongest nationally. |
| Age structure | Align to WHO/WPP-compatible bins; document mapping to `agepi` defaults if used | Avoid unnecessary age splitting beyond targets. |
| Demography | WPP population, fertility, mortality, migration | Public and compatible with current adapter direction. |
| Compartments | Susceptible, latent infection, active TB, on treatment, post-treatment/recovered | Captures minimum TB biology better than acute SEIR. |
| Transmission | Age contact matrix plus calibrated beta | Public contact matrix is synthetic/assumed, so beta must be calibrated. |
| Progression | Fast/slow or single latent-to-active hazard with age modifier | Public data cannot identify rich latency structure. |
| Diagnosis | Active TB to treatment/notification rate | Needed to match incidence and notifications separately. |
| Treatment | Aggregate treatment success/death/loss probabilities | Public treatment outcomes support aggregate constraints. |
| Mortality | Background mortality plus TB mortality/CFR output | Required for WHO mortality target. |
| BCG | Optional paediatric severe TB output modifier | Keep out of transmission unless explicitly modelled. |
| TPT | Scenario-only aggregate coverage | Public data can support broad sensitivity, not individual targeting. |
| MDR/RR-TB | Scenario/reporting output only | Low counts and limited public detail. |
| Calibration targets | WHO incidence, notifications, mortality, age/sex incidence, treatment outcomes | All public. |

Parameters to calibrate in the minimum model:

- transmission coefficient;
- initial latent infection prevalence or latent pool scaling;
- diagnosis/reporting rate by broad period;
- progression/reactivation scale within literature-informed bounds;
- TB mortality or CFR scale;
- treatment outcome probabilities if not fixed from WHO outcomes.

Parameters to fix or use as priors:

- WPP demography;
- age contact matrix;
- age-specific infectiousness and susceptibility patterns;
- TPT efficacy;
- BCG protection against severe childhood TB;
- treatment duration assumptions.

## 9. Upgraded policy model specification

Purpose:

- evaluate Kiribati-relevant intervention scenarios, especially South Tarawa
  screen-and-treat, contact tracing, TPT, active case finding, and treatment
  strengthening;
- separate incidence changes from detection changes;
- represent geography, household exposure, and programme history with local
  evidence.

Recommended upgrades:

| Component | Upgraded choice | Data required |
|---|---|---|
| Geography | South Tarawa plus rest of Kiribati; optionally island groups | Subnational cases, denominators, mobility. |
| Time-varying diagnosis | Period/location-specific diagnosis rates | ACF logs, lab expansion dates, notifications by geography. |
| PEARL intervention | Explicit screening/TPT/treatment pulse in South Tarawa | PEARL enumeration, screening, TB/TBI results, TPT uptake/completion. |
| Household/contact module | Contact tracing cascade and household risk | Index-contact register, household composition, contact outcomes. |
| Latent infection calibration | Age/geography-specific latent prevalence or ARTI | PEARL TST/IGRA/ARTI or school survey data. |
| Treatment pathway | Delay, start, completion, failure, death, relapse | Treatment dates/outcomes linked to cases. |
| Resistance | Optional RR/MDR-TB state or treatment pathway | DST register and DR-TB outcomes. |
| Subnational demography | Island/village age-sex denominators and migration | Census/NSO annual estimates or agreed interpolation. |
| Crowding modifier | Transmission multiplier by household/crowding/geography | Census/PEARL household data. |
| Comorbidity/risk groups | Optional diabetes/smoking/HIV modifiers | Aggregate risk-factor prevalence by age/geography. |

Policy outputs:

- incidence, notifications, mortality, and treatment outcomes by age and
  geography;
- active TB averted by intervention scenario;
- TPT courses needed and TB disease averted;
- expected effect of PEARL-like screening on active TB and TB infection;
- South Tarawa versus outer-island burden and spillover;
- sensitivity of conclusions to latent prevalence, contact matrix, and
  detection assumptions.

## 10. Remaining uncertainties

- Public data do not directly identify latent TB infection prevalence in
  Kiribati. This is the largest uncertainty for a chronic TB model.
- No empirical Kiribati age-contact matrix was identified. Synthetic matrices
  should be treated as structural assumptions, especially because TB transmission
  depends on prolonged indoor and household exposure.
- National WHO estimates are essential but cannot resolve South Tarawa versus
  outer-island heterogeneity.
- Notifications are affected by active case finding, lab changes, programme
  access, and underlying incidence; calibration should not treat them as pure
  incidence.
- Public TPT and contact investigation data need careful denominator checks
  before being used as targets.
- MDR/RR-TB counts are small, so a resistance module should use wide uncertainty
  unless detailed DST data are available.
- BCG should not be used as a simple infection-blocking vaccine. Its most
  defensible initial use is as a paediatric severe TB modifier.
- A policy-facing model requires governance decisions about whether individual
  de-identified data can be shared, or whether aggregate tables are the
  acceptable unit of collaboration.

