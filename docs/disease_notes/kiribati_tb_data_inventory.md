# Kiribati tuberculosis evidence and data inventory

## 1. Executive summary

This note inventories candidate evidence and data sources for a future
age-structured tuberculosis (TB) model for Kiribati in `agepi`. It is a research
and scoping note only; it does not define a final model or imply that public
aggregate data are sufficient for policy analysis.

Key facts from public sources:

- WHO's Global TB Database is the critical public source for annual Kiribati TB
  burden estimates, notifications, treatment outcomes, MDR/RR-TB estimates,
  contact investigation, and TB preventive treatment (TPT) indicators. The
  2025 data release includes national time series and age/sex/risk-factor
  incidence estimates.
- The WHO 2025 burden CSV gives Kiribati 2024 estimated TB incidence of 945 per
  100,000 population, about 1300 incident cases, and estimated TB mortality
  excluding HIV of 37 per 100,000, with wide uncertainty.
- The WHO 2025 age/sex/risk-factor CSV includes Kiribati 2024 incidence by age
  group and sex. Public age groups are useful but do not exactly match the
  current `agepi` default age structure.
- Kiribati 2020 Population and Housing Census materials are the critical public
  source for subnational population, household, crowding, age-sex structure,
  and South Tarawa versus outer-island demography.
- UN World Population Prospects (WPP) 2024 is the critical public source for
  national age-sex population estimates/projections, fertility, mortality, and
  migration schedules compatible with `agepi`'s WPP-style adapter boundary.
- No Kiribati-specific empirical age contact matrix was identified. Candidate
  mixing inputs are synthetic country matrices, projected matrices, or a
  bespoke household/census-informed mixing structure.
- Public evidence consistently points to South Tarawa as the main high-burden
  setting, with lower notification rates in outer islands, but detailed
  subnational TB time series will likely require Kiribati National TB Programme
  (NTP), Ministry of Health and Medical Services (MHMS), PEARL study, or
  collaborator access.

## 2. Candidate model use case

Candidate first policy-relevant model:

- national age-structured TB model with demography, latent infection, active
  infectious disease, diagnosis/treatment, relapse/retreatment, and mortality;
- optional subnational extension with South Tarawa and outer-island strata;
- intervention scenarios for active case finding, contact tracing, TPT, and
  treatment strengthening;
- optional BCG component focused on paediatric severe TB rather than adult
  infection blocking;
- optional MDR/RR-TB compartment or scenario modifier for treatment pathway and
  effective duration of infectiousness.

Facts versus assumptions:

- Fact: Public sources identify high TB burden in Kiribati and a concentration
  of intervention activity in South Tarawa.
- Fact: Public WHO data include annual national estimates and some age/sex
  burden estimates.
- Assumption: A serious TB model should include latent infection and
  reactivation, not a simple acute SEIR interpretation.
- Assumption: South Tarawa versus outer-island heterogeneity is likely important
  enough to test, but the minimum useful geographic split depends on access to
  subnational case and population denominators.

## 3. Core epidemiological data

| Data item | Relevance to model | Likely source | Access | Temporal coverage | Age stratification | Geographic stratification | Priority | Notes/caveats |
|---|---|---|---|---|---|---|---|---|
| Annual TB incidence, mortality, case fatality, treatment coverage | Calibration targets for national incidence, mortality, and diagnosis/treatment gap | WHO Global TB Database burden estimates | Public | WHO annual time series; 2025 release includes estimates through 2024 | Aggregate national burden; separate age/sex incidence file available | National only in public CSV | Critical | WHO data page states estimates include mortality, incidence, age/sex/risk-factor disaggregation, HIV, rifampicin resistance, case fatality, diagnosis/treatment coverage, and household-contact infection estimates. |
| National notifications and notification rate | Calibration/validation target for detected TB; intervention-sensitive outcome | WHO Global TB Database case notifications; Kiribati NTP | Public aggregate; detailed NTP likely private | Annual; WHO CSV updated regularly | Some public age/sex notification outcomes; finer case line list likely private | National public; subnational likely private | Critical | WHO 2025 burden CSV gives Kiribati 2024 notification rate `c_newinc_100k = 444` and case detection/treatment coverage estimate around 47%. |
| Age/sex incidence estimates | Age-specific force of infection/progression calibration | WHO TB incidence estimates disaggregated by age group, sex, and risk factor | Public | Annual; 2025 release through 2024 | Yes, but in WHO age groups such as 0-4, 5-9, 5-14, 15-24, 25-34, etc. | National only | Critical | The 2024 WHO age/sex CSV gives total incidence about 1300, female about 580, male about 690; age groups need mapping to `agepi`. |
| Treatment outcomes | Treatment success, death, failure, loss-to-follow-up, relapse/retreatment calibration | WHO Global TB Database treatment outcomes; historical Kiribati paper | Public aggregate; NTP line-list private | WHO annual; historical paper reports 2000-2011/2012 | WHO has treatment outcomes by age group and sex CSV | National public; subnational likely private | High | Viney et al. report high treatment commencement and success in 2000-2011, but current programmatic time series should come from WHO/NTP. |
| TB/HIV indicators | Comorbidity pathway and mortality/progression modifier | WHO TB burden estimates and TB/HIV programme data | Public aggregate; details private | Annual | Limited public stratification | National | Medium | WHO 2024 estimates imply very low TB/HIV contribution in Kiribati, but confirm with NTP/HIV programme if model includes HIV. |
| Latent TB infection prevalence | Initial conditions and reactivation pool | PEARL study; TST/IGRA/school ARTI data; WHO household-contact infection estimates | Mostly private/unclear; WHO household-contact estimates public | PEARL contemporary; WHO annual household-contact estimates | PEARL likely individual age; public WHO aggregate | PEARL South Tarawa; WHO national | Critical | Public aggregate incidence alone cannot identify latent prevalence and reactivation without strong assumptions. |
| Smear/culture/Xpert confirmation, pulmonary/extrapulmonary split | Infectiousness and diagnostic pathway | WHO case notifications; NTP register/lab data | Public aggregate; detailed private | Annual | Limited | National public; facility/geography private | High | Needed to distinguish infectious pulmonary disease from all notified TB. |

## 4. Demographic data

| Data item | Relevance to model | Likely source | Access | Temporal coverage | Age stratification | Geographic stratification | Priority | Notes/caveats |
|---|---|---|---|---|---|---|---|---|
| National population by single year/five-year age and sex | Population denominators, initial state, age structure | UN WPP 2024 | Public | 1950-2100 estimates/projections | Yes, one-year age/time since WPP 2022/2024 | National | Critical | Best fit for `agepi` WPP-style population, fertility, mortality, and migration adapters. |
| Fertility schedules/TFR | Birth inflow in demographic process | UN WPP 2024 | Public | 1950-2100 | Maternal age schedules | National | Critical | Use with existing WPP-style fertility helpers after checking units. |
| Mortality schedules | Background mortality by age | UN WPP 2024; World Bank for summary indicators | Public | 1950-2100 WPP | Yes in WPP | National | Critical | Use WPP central death rates where available; avoid qx conversion unless metadata are explicit. |
| Net migration | Demographic residuals and population projection consistency | UN WPP 2024; census comparison | Public national; subnational private/unclear | 1950-2100 WPP | WPP age/sex schedules where available | National only | High | South Tarawa and Kiritimati movement may be epidemiologically important but is not captured by national net migration alone. |
| 2020 census population, households, age/sex, island/village, urban/rural | Subnational denominators and crowding proxies | Kiribati 2020 Population and Housing Census, Kiribati NSO/SPC/Pacific Data Hub | Public report/atlas; microdata access controlled | Census night 7 Nov 2020 | Yes in report/tables/microdata | Island, village, urban/rural, districts | Critical | Pacific Data Hub describes report and atlas with island, village, urban/rural groupings and CAPI enumeration. |
| Household size/crowding/housing | TB transmission risk proxy and household-contact model | 2020 census report/atlas; HIES; NTP contact data | Public aggregate; microdata controlled | Census 2020; HIES 2019 | Household-level if microdata accessed | Island/village possible | High | Important for South Tarawa high-density transmission; model may need household or close-contact multiplier. |
| World Bank demographic indicators | Cross-check denominators, fertility, mortality, age dependency | World Bank Data / HealthStats / WDI | Public | Annual historical | Summary age bands | National | Medium | Useful for quick QA, but WPP/census should be primary. |

## 5. Contact/mixing data

| Data item | Relevance to model | Likely source | Access | Temporal coverage | Age stratification | Geographic stratification | Priority | Notes/caveats |
|---|---|---|---|---|---|---|---|---|
| Kiribati empirical age contact matrix | Direct age-assortative mixing input | None identified publicly | Not found | Not found | Would be yes if collected | Unknown | Critical gap | No Kiribati-specific diary/contact survey found in this search. |
| Synthetic country contact matrices | First-pass age mixing | Prem et al. synthetic matrices; `socialmixr`; `conmat`-style predictions | Public depending on dataset/package | Usually baseline pre-pandemic; projection variants exist | Yes | National synthetic | High | Need verify whether Kiribati is included; if absent, use nearest Pacific proxy or generated synthetic matrix with Kiribati demography. |
| Projected contact matrices from demographic/household data | Alternative when no empirical matrix exists | Prem et al. 2017/2021 methods; Arregui et al. demographic projection methods | Public methods/data | Depends on chosen baseline/projection | Yes | National | Medium | Synthetic matrices may underrepresent TB-relevant prolonged indoor exposure. |
| Household/close-contact mixing matrix | TB-specific transmission structure | 2020 census household composition; NTP contact-tracing data; PEARL household enumeration | Census aggregate public; NTP/PEARL private | 2020 census; PEARL contemporary | Yes if microdata/PEARL accessed | South Tarawa and islands possible | High | Better biological fit for TB than all social contacts; likely requires controlled data access. |
| Location-specific contacts (home, school, work, other) | Scenario analysis for intervention or setting-specific transmission | Synthetic matrices; census school/work data | Public/unclear | Baseline | Yes | National synthetic; census subnational proxies | Medium | Current `agepi` accepts a static contact matrix; setting-specific decomposition would be upstream. |

## 6. Intervention and programme data

| Data item | Relevance to model | Likely source | Access | Temporal coverage | Age stratification | Geographic stratification | Priority | Notes/caveats |
|---|---|---|---|---|---|---|---|---|
| Active case finding history | Time-varying detection rate and scenario design | Kiribati NTP; UNDP/Global Fund reports; DFAT project design; PEARL | Public narrative; detailed operational data private | 2000s-present, uneven public detail | Usually no in public reports | South Tarawa/hotspots described; detailed geography private | Critical | UNDP reported a 2018 hotspot approach with active case finding, immediate treatment, contact tracing, and DOTS arrangements. |
| Contact tracing | Household transmission/intervention impact | WHO contact/TPT CSV; NTP; DFAT/UNDP reports | Public aggregate; detailed private | WHO annual; NTP historical | Limited public | National public; local private | Critical | WHO contact/TPT CSV has variables for screened/evaluated contacts and preventive treatment, but interpretation requires data dictionary and completeness checks. |
| TB preventive treatment (TPT) | Prevention scenario; reduces progression from infection | WHO contact/TPT CSV; PEARL; NTP | Public aggregate; PEARL/NTP detailed private | WHO annual; PEARL intervention period | Limited public; PEARL likely age | National public; PEARL South Tarawa | Critical | Public WHO CSV includes contact TPT and short regimen variables; PEARL is the most important detailed source if accessible. |
| Diagnosis tools and lab strengthening | Detection delay and bacteriological confirmation | WHO lab diagnostic services CSV; NTP; project reports | Public aggregate; detailed private | WHO annual | No/limited | National public | High | Historical paper notes enhanced case finding and liquid culture contributed to increased notifications. |
| DOTS/treatment support | Treatment success, duration infectious, relapse | NTP; DFAT QTBECP design; WHO treatment outcomes | Public narrative/aggregate; detailed private | 2000s-present | Limited public | South Tarawa/outer islands likely private | High | DFAT project design discusses DOTS expansion and outer-island service sustainability. |
| PEARL screen-and-treat intervention | Major contemporary scenario and impact evaluation | PEARL protocol, publications, investigators, MHMS | Protocol public; participant/results data private/unclear | Protocol published 2022; intervention in South Tarawa | Likely individual age | South Tarawa intervention; rest of Kiribati comparator | Critical | Protocol screens residents aged >=3 years and selected exposed children <3 years using TST, symptoms, CXR/CAD, and Xpert Ultra; offers TPT for TB infection. |

## 7. Drug-resistance evidence

| Data item | Relevance to model | Likely source | Access | Temporal coverage | Age stratification | Geographic stratification | Priority | Notes/caveats |
|---|---|---|---|---|---|---|---|---|
| MDR/RR-TB burden estimates | MDR/RR compartment or scenario modifier | WHO MDR/RR-TB burden estimates | Public | Annual; 2025 release through 2024 | No/limited | National | High | WHO 2024 MDR/RR CSV estimates RR/MDR among new TB about 1.1% and among retreatment about 12%, with wide uncertainty; estimated incident RR/MDR-TB about 18 cases. |
| Drug susceptibility testing coverage/results | Diagnostic pathway and resistance ascertainment | WHO drug resistance testing CSV; NTP lab data | Public aggregate; detailed private | Since 2017 in WHO CSV | No/limited | National public | High | Needed to interpret whether low resistance estimates reflect true burden or testing coverage. |
| Molecular epidemiology and resistance | Transmission clustering; lineage/resistance assumptions | First molecular epidemiology study of M. tuberculosis in Kiribati | Public peer-reviewed | Study period around 2009-2011 | Limited patient data | Mostly Tarawa/Tungaru Central Hospital sample | Medium | Study found no MDR strains among tested MTBC strains and little first-line resistance; it is older and not enough for current MDR assumptions. |
| Historical drug-resistant cases | Scenario bounds and surveillance caution | Viney et al. 2000-2012 paper; NTP | Public paper; detailed private | 2000-2012 | No | National | Medium | Viney et al. report several drug-resistant samples including one MDR-TB case who died; reconcile with molecular study and WHO recent estimates. |

## 8. Vaccination/prevention evidence

| Data item | Relevance to model | Likely source | Access | Temporal coverage | Age stratification | Geographic stratification | Priority | Notes/caveats |
|---|---|---|---|---|---|---|---|---|
| BCG national coverage | Paediatric severe TB protection and birth cohort history | WHO/UNICEF WUENIC; WHO Immunization Data portal | Public | WUENIC data from 1980-2024; country PDF available at least through 2023 | Birth cohort/infant coverage, not individual age | National | High | 2023 country profile reports BCG estimate 96% for 2021-2023 after lower estimates in 2014-2016; 2024 release should be pulled from portal for final modelling. |
| BCG subnational coverage | Geographic differences in paediatric protection | EPI programme, surveys, MICS, census-linked immunisation data | Public survey aggregate/unclear; programme private | Survey/programme dependent | Children | Island/region if survey design supports | Medium | Important if South Tarawa/outer-island paediatric burden differs; likely not central for adult pulmonary TB transmission. |
| TPT among contacts | Prevention intervention and calibration | WHO contact/TPT CSV; NTP; PEARL | Public aggregate; detailed private | WHO annual; PEARL intervention | Limited public | National public; PEARL South Tarawa | Critical | Need eligibility denominator, regimen, completion, and targeting by age/household exposure. |
| School annual risk of infection (ARTI) / TST | Infection transmission calibration | PEARL protocol and study data | Results private/unclear | Contemporary PEARL | School children | South Tarawa | High | PEARL protocol uses annual risk of TB infection in primary school children as a co-primary transmission outcome. |

## 9. Geography: Tarawa versus outer islands

| Data item | Relevance to model | Likely source | Access | Temporal coverage | Age stratification | Geographic stratification | Priority | Notes/caveats |
|---|---|---|---|---|---|---|---|---|
| South Tarawa versus rest of Kiribati TB notifications | Two-patch model calibration and intervention contrast | NTP; PEARL; DFAT project design | Mostly private; public narrative and figures | Historical 2000-2010 figures in DFAT design; PEARL contemporary | Unknown | Yes | Critical | DFAT design states outer-island notification rates are consistently lower than South Tarawa and lists several Gilbert islands with recurrent cases. |
| Island/village denominators | Geographic population denominators | 2020 census report/atlas; Pacific Data Hub/PopGIS | Public aggregate; microdata controlled | 2020 census | Yes in census tables/microdata | Island/village/urban-rural | Critical | Needed to calculate rates consistently rather than use raw cases. |
| Mobility between South Tarawa, outer islands, and Kiritimati | Spatial coupling and importation/reactivation risk | Census migration questions, transport data, NTP/collaborators | Public partial; detailed likely private | Census and programme periods | Possibly | Island-level | High | DFAT design notes travel between Tarawa and outer islands and Kiritimati resettlement as relevant risks. |
| Crowding and household composition in South Tarawa | Transmission multiplier | Census, HIES, PEARL enumeration | Public aggregate/private detail | 2020 census; PEARL | Household age composition if microdata | Village/island | High | Molecular paper notes crowded South Tarawa conditions as relevant for ongoing chains of transmission. |

## 10. Calibration targets

Candidate public calibration targets:

| Target | Public source | Use | Caveat |
|---|---|---|---|
| National TB incidence rate and count by year | WHO burden estimates CSV | Fit transmission/progression scale | WHO estimates depend on surveillance adjustments and are uncertain. |
| National notifications and notification rate by year | WHO case notification CSV / burden CSV | Fit diagnosis/reporting rate | Notifications respond to active case finding and lab changes. |
| Age/sex incidence distribution | WHO age/sex/risk-factor incidence CSV | Fit age-specific burden and mixing/progression assumptions | WHO age groups need mapping and may be modelled estimates, not raw cases. |
| TB mortality and CFR | WHO burden estimates CSV | Fit mortality/treatment outcome pathway | Death registration quality and TB attribution may be uncertain. |
| Treatment success/loss/death | WHO treatment outcome CSV; Viney et al. historical paper | Fit treatment pathway | Public aggregate may lack subnational and age detail. |
| MDR/RR-TB incidence/proportion | WHO MDR/RR estimates CSV | Fit resistance scenario | Small numbers yield wide uncertainty. |
| Contact investigation and TPT counts | WHO contact/TPT CSV | Fit prevention coverage | Must check denominator definitions and reporting completeness. |
| BCG birth cohort coverage | WUENIC | Paediatric severe TB modifier | Does not measure protection against adult pulmonary disease. |
| South Tarawa notification rate before/after PEARL | PEARL/NTP | Evaluate intervention scenario | Likely requires collaborator access. |
| School-child ARTI or TST conversion | PEARL | Fit transmission intensity | Likely requires PEARL collaborator access. |

## 11. Data gaps and collaborator requests

Likely collaborator/NTP requests:

| Request | Why it matters | Minimum useful fields | Priority | Notes |
|---|---|---|---|---|
| De-identified TB notification line list | Age/geography/time calibration | Diagnosis date/year, age, sex, island/village, pulmonary/site, bacteriological confirmation, new/retreatment, outcome | Critical | Essential for South Tarawa versus outer-island modelling. |
| Case finding activity log | Time-varying detection intervention | Dates, locations, population screened, screening algorithm, tests, cases found, costs if available | Critical | Needed to separate real incidence from detection changes. |
| Contact tracing and TPT register | Household transmission/prevention | Index case ID/linkage, contact age/sex, household/village, screening result, infection test, regimen, start/completion | Critical | Public WHO aggregates are not enough for age-specific household TPT effect. |
| PEARL individual or aggregate outputs | Intervention scenario and infection prevalence | Enumeration denominator, age/sex, TST/CXR/Xpert results, TB/TBI diagnosis, TPT uptake/completion, geography | Critical | Could transform model from generic national calibration to policy-relevant South Tarawa analysis. |
| Drug susceptibility/lab register | MDR/RR pathway | Test type/date, rifampicin/isoniazid and first-line results, previous treatment, geography | High | Small numbers require careful interpretation. |
| Subnational denominators and mobility | Two-patch model | Annual island population, age/sex, migration/mobility flows | High | Census gives 2020; annual denominators and flows may require NSO/MHMS collaboration. |
| Treatment pathway timing | Infectious duration and outcomes | Symptom onset, diagnosis date, treatment start, smear conversion, completion/outcome, relapse | High | Needed for realistic infectious-duration assumptions. |
| Comorbidity data linked or aggregated | Risk stratification | Diabetes, smoking, undernutrition, HIV where ethically/shareably available | Medium | Diabetes and smoking may be important in Kiribati but add complexity. |

## 12. Recommended next steps

1. Build a public-data extraction notebook outside the package source tree that
   pulls only filtered Kiribati rows from WHO TB CSVs, WPP, WUENIC, and census
   tables; do not commit raw datasets.
2. Decide whether the first model is national-only or two-patch
   South-Tarawa/rest-of-Kiribati. If two-patch, prioritise NTP/PEARL access.
3. Map WHO age groups and WPP/census age groups to an explicit `agepi` age
   structure before any calibration.
4. Choose an interim contact matrix strategy: synthetic Kiribati matrix if
   available, otherwise a demography-projected matrix plus a household/crowding
   sensitivity analysis.
5. Request NTP/PEARL aggregate tables first, then line-list access only if the
   governance path is feasible.
6. Treat MDR/RR-TB, BCG, and TPT as scenario modules until the core latent-active
   TB model is stable.

## Sources checked

- WHO Global Tuberculosis Programme data page and CSV catalogue:
  <https://www.who.int/teams/global-tuberculosis-programme/data>
- WHO TB dashboard:
  <https://data.who.int/dashboards/tuberculosis>
- WHO Global tuberculosis report 2025:
  <https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/tb-reports/global-tuberculosis-report-2025>
- Kiribati 2020 Population and Housing Census, Pacific Data Hub:
  <https://microdata.pacificdata.org/index.php/catalog/767/related-materials>
- UN World Population Prospects 2024:
  <https://population.un.org/wpp/>
- World Bank Kiribati data:
  <https://data.worldbank.org/country/kiribati>
- WHO/UNICEF WUENIC Q&A:
  <https://www.who.int/news-room/questions-and-answers/item/who-unicef-estimates-of-national-immunization-coverage>
- WHO Immunization Data portal:
  <https://immunizationdata.who.int/>
- Kiribati WUENIC 2023 country profile PDF:
  <https://cdn.who.int/media/docs/default-source/country-profiles/immunization/2024-country-profiles/immunization-2024-kir.pdf>
- PEARL study website:
  <https://www.thepearlstudy.org/>
- PEARL protocol, BMJ Open 2022:
  <https://bmjopen.bmj.com/content/12/4/e055295>
- UNDP hotspot case-finding story:
  <https://www.undp.org/pacific/news/innovative-hotspot-approach-aims-bring-tb-under-control-kiribati>
- DFAT Kiribati TB elimination project design/funding proposal:
  <https://www.dfat.gov.au/sites/default/files/kiribati-tb-elimination-project-design-funding-proposal.pdf>
- Viney et al. "Battling tuberculosis in an island context...", 2000-2012:
  <https://openresearch-repository.anu.edu.au/items/a4b1721d-df8a-4ed0-b6b8-b808c241b642>
- First molecular epidemiology study of Mycobacterium tuberculosis in Kiribati:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3561247/>
- Prem et al. synthetic contact matrix literature:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC5609774/>
- Arregui et al. contact matrix projection method:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC6300299/>

## Repository files inspected

- `README.md`
- `ROADMAP.md`
- `docs/model_conventions.md`
- `docs/external_data_adapters.md`
- `docs/disease_notes/TEMPLATE_parameterisation.md`
- `docs/disease_notes/tuberculosis_parameterisation.md`
- repository-wide search hits for TB, Kiribati, demography, WPP, contact
  matrices, calibration, vaccination, BCG, Tarawa, outer islands, and MDR terms

