# Roadmap

## Milestone 1: Project skeleton and age/state infrastructure

Implement:

- source folder structure;
- age-structure constructor and validator;
- state long-to-vector and vector-to-long mapping;
- tests for age validation and state mapping.

Do not implement disease transmission yet.

## Milestone 2: Force of infection and SIR transition rates

Implement:

- reusable force-of-infection function;
- SIR disease model constructor;
- transition-rate table for infection, recovery, ageing, births, and background deaths;
- tests using small manually checkable examples.

## Milestone 3: Deterministic simulation

Implement:

- conversion from transition rates to derivatives;
- deterministic ODE simulation;
- tidy long-format simulation outputs;
- smoke test using mock demographic data.

## Milestone 4: Demography layer

Implement:

- demography object constructor;
- demographic quantity lookup/interpolation;
- mock WPP-style input format;
- tests showing simulation code is independent of data source.

## Milestone 5: Contact matrix integration

Implement:

- mixing model constructor;
- contact matrix validation;
- compatibility with externally generated contact matrices.
