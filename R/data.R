#' Childhood vaccination rates in Norway (2020 borders)
#'
#' We conveniently package vaccine coverage data taken from "Kommunehelsa statistikkbank".
#' This data was last updated on 2019-04-09.
#'
#' This dataset contains national/county/municipality level (5 year average) vaccination coverage rates
#' for 16 year olds for the childhood vaccination program (diphtheria, hpv, measles,
#' mumps, poliomyelitis, pertussis, rubella, tetanus).
#'
#' Municipalities are updated for the 2020 borders.
#'
#' @format
#' \describe{
#' \item{calyear}{The middle year of a 5 year range (e.g. 2011 is the average of data from 2009-2013).}
#' \item{location_code}{The location code.}
#' \item{age}{The population age.}
#' \item{vaccine}{The vaccine.}
#' \item{vaccination_coverage_pr100}{Percentage of people who are vaccinated.}
#' \item{vaccination_coverage_pr1}{Proportion of people who are vaccinated.}
#' \item{imputed}{FALSE if real data. TRUE if it is the national average.}
#' }
#' @source \url{http://khs.fhi.no/webview/}
"nor_childhood_vax_b2020"

#' A fake commuter dataset for Norway as a single entity in 2017
#'
#' For use with single_entity_seiiar_2017
#'
#' @format
#' \describe{
#' \item{from}{Location code.}
#' \item{to}{Location name.}
#' \item{n}{Number of people.}
#' }
"single_entity_fake_commuters_2017"

#' SEIIaR data.frame for Norway as a single entity with no one infected and everyone susceptible (2020 borders)
#'
#' For use with single_entity_fake_commuters_2017
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E}{Number of exposed people.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' }
"single_entity_seiiar_2017"

#' Daily number of commuters from/to municipalities in Norway in 2017 (2020 borders)
#'
#' @format
#' \describe{
#' \item{from}{Location code.}
#' \item{to}{Location name.}
#' \item{n}{Number of people.}
#' }
#' @source \url{https://www.ssb.no/statbank/table/03321}
"nor_commuters_2017_b2020"

#' SEIIaR data.frame for Norway with no one infected and everyone susceptible (2020 borders)
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E}{Number of exposed people.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' }
"nor_seiiar_noinfected_2017_b2020"

#' SEIIaR data.frame for Norway with 10 people infected in Oslo and everyone susceptible (2020 borders)
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E}{Number of exposed people.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' }
"nor_seiiar_oslo_2017_b2020"

#' SEIIaR data.frame for Norway with no one infected and real measles susceptibility (2020 borders)
#'
#' Measles vaccination coverate rates for 16 year olds in the 5 year average
#' from 2014 to 2018 were used as the proportion of recovered people.
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E}{Number of exposed people.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' }
#' @source \url{http://khs.fhi.no/webview/}
"nor_seiiar_measles_noinfected_2017_b2020"

#' SEIIaR data.frame for Norway with 10 people infected in Oslo and real measles susceptibility (2020 borders)
#'
#' Measles vaccination coverate rates for 16 year olds in the 5 year average
#' from 2014 to 2018 were used as the proportion of recovered people.
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E}{Number of exposed people.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' }
#' @source \url{http://khs.fhi.no/webview/}
"nor_seiiar_measles_oslo_2017_b2020"

#' Convert blank seiiar with vax
#'
#' Takes a fully susceptible population and
#' proportion of people vaccinated per location code
#' and allocates an appropriate amount of people to
#' recovered.
#' For more information, look at \code{vignette("including_vax","spread")}.
#' @param seiiar SEIIAR data.table representing a fully susceptible population
#' @param vax data.table containing proportion of people vaccinated per location code ("vaccination_coverage_pr1")
#' @examples
#' vax_measles <- spread::nor_childhood_vax_b2020[
#'   calyear == 2016 &
#'   vaccine == "measles"
#' ]
#'
#' nor_seiiar_measles_noinfected_2017_b2020 <- spread::convert_blank_seiiar_with_vax(
#'   seiiar = spread::nor_seiiar_noinfected_2017_b2020,
#'   vax = vax_measles
#' )
#' @export
convert_blank_seiiar_with_vax <- function(seiiar, vax) {
  R <- NULL
  S <- NULL
  proportion <- NULL
  validate_seiiar(seiiar)
  validate_vax_prev(vax)

  retval <- copy(seiiar)
  setDT(retval)
  retval[vax, on = "location_code", R := round(S * vaccination_coverage_pr1)]
  retval[, S := S - R]

  return(retval)
}

#' A fake population for asymmetric mobility
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E}{Number of exposed people.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' }
"asymmetric_mobility_dummy_seiiar_pop"

#' Fake mobility matrixes for asymmetric mobility
#'
#' A list of 20 matrices (1 for each time period)
#'
#' @format
#' \describe{
#' \item{from}{Location code.}
#' \item{to}{Location code.}
#' \item{n}{Number of people.}
#' }
"asymmetric_mobility_dummy_mobility_matrix"

#' Fake betas for asymmetric mobility
#'
#' A vector of 20 betas (1 for each time period)
#'
"asymmetric_mobility_dummy_betas"

#' Fake dynamic seeds for asymmetric mobility
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{day}{Day seeding occurs.}
#' \item{n}{Number of people.}
#' }
"asymmetric_mobility_dummy_dynamic_seeds"

create_asymmetric_mobility_dummy_files <- function(base_loc) {
  seiiar_pop <- data.table::data.table(
    "location_code" = c("a", "b", "c"),
    "S" = c(1000, 1000, 2000),
    "E" = c(0, 0, 0),
    "I" = c(50, 0, 0),
    "Ia" = c(0, 0, 0),
    "R" = c(0, 0, 0)
  )

  temp <- data.table::data.table(
    from = c("a", "a", "b", "b", "c", "c"),
    to = c("b", "c", "a", "c", "a", "b"),
    n = c(50, 10, 10, 10, 10, 10)
  )
  mobility_matrix <- vector("list", length = 20)
  for (i in seq_along(mobility_matrix)) {
    mobility_matrix[[i]] <- data.table::copy(temp)
    data.table::setnames(mobility_matrix[[i]], c("from", "to", "n"))
  }

  betas <- rep(0.6, 20)

  asymmetric_mobility_dummy_seiiar_pop <- seiiar_pop
  save(asymmetric_mobility_dummy_seiiar_pop, file = file.path(base_loc, "asymmetric_mobility_dummy_seiiar_pop.rda"), compress = "xz")

  asymmetric_mobility_dummy_mobility_matrix <- mobility_matrix
  save(asymmetric_mobility_dummy_mobility_matrix, file = file.path(base_loc, "asymmetric_mobility_dummy_mobility_matrix.rda"), compress = "xz")

  asymmetric_mobility_dummy_betas <- betas
  save(asymmetric_mobility_dummy_betas, file = file.path(base_loc, "asymmetric_mobility_dummy_betas.rda"), compress = "xz")

  asymmetric_mobility_dummy_dynamic_seeds <- data.table::data.table(
    "location_code" = c("a"),
    "day" = c(3, 7),
    "n" = c(3, 5)
  )
  save(asymmetric_mobility_dummy_dynamic_seeds, file = file.path(base_loc, "asymmetric_mobility_dummy_dynamic_seeds.rda"), compress = "xz")
}

#' A fake population for asymmetric mobility se1e2iiar
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E1}{Number of exposed, asymptomatic people, not infectious.}
#' \item{E2}{Number of exposed, presymptomatic people, infectious.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' }
"asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_pop"

#' Fake mobility matrixes for asymmetric mobility se1e2iiar
#'
#' A list of 20 matrices (1 for each time period)
#'
#' @format
#' \describe{
#' \item{from}{Location code.}
#' \item{to}{Location code.}
#' \item{n}{Number of people.}
#' }
"asymmetric_mobility_se1e2iiar_dummy_mobility_matrix"

#' Fake betas for asymmetric mobility se1e2iiar
#'
#' A data table of 60 betas (1 for each time period for each of the three locations)
#'
"asymmetric_mobility_se1e2iiar_dummy_betas"

#' Fake dynamic seeds for asymmetric mobility se1e2iiar
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{day}{Day seeding occurs.}
#' \item{n}{Number of people.}
#' }
"asymmetric_mobility_se1e2iiar_dummy_dynamic_seeds"

#' A fake population for asymmetric mobility se1e2iiar 2 strains
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{S}{Number of susceptible people.}
#' \item{E1}{Number of exposed, asymptomatic people, not infectious.}
#' \item{E2}{Number of exposed, presymptomatic people, infectious.}
#' \item{I}{Number of infectious and symptomatic people.}
#' \item{Ia}{Number of infectious and asymptomatic people.}
#' \item{R}{Number of recovered people.}
#' \item{E1_b}{Number of exposed, asymptomatic people, not infectious of strain B.}
#' \item{E2_b}{Number of exposed, presymptomatic people, infectious of strain B.}
#' \item{I_b}{Number of infectious and symptomatic people of strain B.}
#' \item{Ia_b}{Number of infectious and asymptomatic people of strain B.}
#' }
"asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop"

#' Fake mobility matrixes for asymmetric mobility se1e2iiar
#'
#' A list of 20 matrices (1 for each time period)
#'
#' @format
#' \describe{
#' \item{from}{Location code.}
#' \item{to}{Location code.}
#' \item{n}{Number of people.}
#' }
"asymmetric_mobility_se1e2iiar_dummy_mobility_matrix"

#' Fake betas for asymmetric mobility se1e2iiar
#'
#' A data table of 60 betas (1 for each time period for each of the three locations)
#'
"asymmetric_mobility_se1e2iiar_dummy_betas"

#' Fake dynamic seeds for asymmetric mobility se1e2iiar
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{day}{Day seeding occurs.}
#' \item{n}{Number of people.}
#' \item{n_b}{Number of people strain B.}
#' }
"asymmetric_mobility_se1e2iiar_2strains_dummy_dynamic_seeds"


