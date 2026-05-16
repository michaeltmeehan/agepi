#' Build a demographic process from schedules
#'
#' Convenience builder for combining an age structure and already-constructed
#' demographic schedule objects into a [DemographicProcess()]. The builder
#' constructs [AgeingOperator()] from `age_structure` and delegates validation
#' to [DemographicProcess()].
#'
#' This function does not parse WPP data, simulate demography, infer missing
#' fertility, mortality, or migration schedules, implement interpolation, or
#' implement residual forcing. It is only a small assembly helper for schedule
#' objects that have already been standardised or constructed.
#'
#' @param age_structure Age structure validated by [validate_age_structure()].
#' @param fertility_schedule Optional [FertilitySchedule()] object.
#' @param mortality_schedule Optional [MortalitySchedule()] object.
#' @param migration_schedule Optional [MigrationSchedule()] object.
#' @param mode Process mode, either `"closed"` or `"migration"`.
#'
#' @return An `agepi_demographic_process` object.
#' @examples
#' ages <- wpp_age_structure_5year(max_age = 10)
#' process <- build_demographic_process(ages)
#' process$mode
#'
#' mortality <- MortalitySchedule(
#'   data.frame(
#'     time = 2020,
#'     age_group = ages$age_groups,
#'     mortality_rate = rep(0.01, ages$n_age_groups)
#'   ),
#'   ages
#' )
#' build_demographic_process(
#'   age_structure = ages,
#'   mortality_schedule = mortality
#' )
#' @export
build_demographic_process <- function(age_structure,
                                      fertility_schedule = NULL,
                                      mortality_schedule = NULL,
                                      migration_schedule = NULL,
                                      mode = c("closed", "migration")) {
  mode <- match.arg(mode)
  ageing_operator <- AgeingOperator(age_structure)

  DemographicProcess(
    age_structure = age_structure,
    ageing_operator = ageing_operator,
    fertility_schedule = fertility_schedule,
    mortality_schedule = mortality_schedule,
    migration_schedule = migration_schedule,
    mode = mode
  )
}
