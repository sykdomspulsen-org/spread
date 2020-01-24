#' Daily number of commuters from/to municipalities in Norway in 2017 (2020 borders)
#'
#' @format
#' \describe{
#' \item{from}{Location code.}
#' \item{to}{Location name.}
#' \item{n}{Number of people.}
#' }
#' @source \url{https://www.ssb.no/statbank/table/03321}
"norway_commuters_2017_b2020"

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
"norway_seiiar_noinfected_2017_b2020"

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
"norway_seiiar_oslo_2017_b2020"

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
"norway_seiiar_measles_noinfected_2017_b2020"

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
"norway_seiiar_measles_oslo_2017_b2020"

create_blank_norway_2017 <- function() {
  . <- NULL
  year <- NULL
  n <- NULL
  from <- NULL
  to <- NULL
  level <- NULL
  pop <- NULL
  location_code <- NULL
  municip_code_current <- NULL
  weighting <- NULL
  S <- NULL

  dirData <- system.file("extdata", package = "spread")

  di_edge_list <- data.table(readxl::read_excel(file.path(dirData, "di_edge_list_2017.xlsx"), skip = 3))
  setnames(di_edge_list, c("from", "to", "n"))
  di_edge_list[, year := 2017]
  di_edge_list <- di_edge_list[!is.na(n)]
  di_edge_list[, from := zoo::na.locf(from)]

  di_edge_list[, from := stringr::str_extract(from, "^[0-9][0-9][0-9][0-9]")]
  di_edge_list[, to := stringr::str_extract(to, "^[0-9][0-9][0-9][0-9]")]

  di_edge_list[, from := sprintf("municip%s", from)]
  di_edge_list[, to := sprintf("municip%s", to)]

  norwayMunicipMerging <- fhidata::norway_municip_merging_b2020
  sum(di_edge_list$n)
  di_edge_list <- merge(
    di_edge_list,
    norwayMunicipMerging,
    by.x = c("from", "year"),
    by.y = c("municip_code_original", "year"),
    all.x = T
  )
  sum(di_edge_list$n)
  di_edge_list[, n := n * weighting]
  di_edge_list <- di_edge_list[!is.na(n) & n > 0]
  sum(di_edge_list$n, na.rm = T)
  di_edge_list[, from := NULL]
  di_edge_list[, weighting := NULL]
  setnames(di_edge_list, "municip_code_current", "from")

  di_edge_list <- merge(
    di_edge_list,
    norwayMunicipMerging,
    by.x = c("to", "year"),
    by.y = c("municip_code_original", "year"),
    all.x = T
  )
  di_edge_list[, n := n * weighting]
  di_edge_list <- di_edge_list[!is.na(n) & n > 0]
  sum(di_edge_list$n, na.rm = T)
  di_edge_list[, to := NULL]
  di_edge_list[, weighting := NULL]
  setnames(di_edge_list, "municip_code_current", "to")

  di_edge_list <- di_edge_list[, .(n = sum(n)), keyby = .(
    from,
    to
  )]
  di_edge_list <- di_edge_list[from != to]
  di_edge_list <- di_edge_list[n > 0]
  sum(di_edge_list$n)
  seiiar <- fhidata::norway_population_b2020[year == 2017 & level == "municipality", .(
    S = sum(pop),
    E = 0,
    I = 0,
    Ia = 0,
    R = 0
  ), by = .(location_code)]

  return(list(
    seiiar = seiiar,
    commuters = di_edge_list
  ))
}

#' Convert blank seiiar with vax
#'
#' Takes a fully susceptible population and
#' proportion of people vaccinated per location code
#' and allocates an appropriate amount of people to
#' recovered.
#' For more information, look at \code{vignette("including_vax","spread")}.
#' @param seiiar SEIIAR data.table representing a fully susceptible population
#' @param vax data.table containing proportion of people vaccinated per location code
#' @examples
#' vax_measles <- fhidata::norway_childhood_vax_b2020[
#'   year == 2016 &
#'     stringr::str_detect(location_code, "^municip") &
#'     vax == "measles",
#'   c("location_code", "proportion")
#' ]
#'
#' norway_seiiar_measles_noinfected_2017_b2020 <- spread::convert_blank_seiiar_with_vax(
#'   seiiar = spread::norway_seiiar_noinfected_2017_b2020,
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
  retval[vax, on = "location_code", R := round(S * proportion)]
  retval[, S := S - R]

  return(retval)
}

create_data_files_norway_2017 <- function(base_loc) {
  . <- NULL
  year <- NULL
  n <- NULL
  from <- NULL
  to <- NULL
  level <- NULL
  pop <- NULL
  location_code <- NULL
  S <- NULL
  vax <- NULL

  x <- create_blank_norway_2017()

  seiiar <- x[["seiiar"]]
  norway_commuters_2017_b2020 <- x[["commuters"]]

  save(norway_commuters_2017_b2020, file = file.path(base_loc, "norway_commuters_2017_b2020.rda"), compress = "xz")

  norway_seiiar_noinfected_2017_b2020 <- seiiar
  save(norway_seiiar_noinfected_2017_b2020, file = file.path(base_loc, "norway_seiiar_noinfected_2017_b2020.rda"), compress = "xz")

  norway_seiiar_oslo_2017_b2020 <- copy(seiiar)
  norway_seiiar_oslo_2017_b2020[location_code == "municip0301", I := 10]
  norway_seiiar_oslo_2017_b2020[location_code == "municip0301", S := S - I]
  save(norway_seiiar_oslo_2017_b2020, file = file.path(base_loc, "norway_seiiar_oslo_2017_b2020.rda"), compress = "xz")

  # measles
  vax_prev <- fhidata::norway_childhood_vax_b2020[year == 2016 & stringr::str_detect(location_code, "^municip") & vax == "measles"]
  norway_seiiar_measles_noinfected_2017_b2020 <- convert_blank_seiiar_with_vax(seiiar, vax_prev)

  save(norway_seiiar_measles_noinfected_2017_b2020, file = file.path(base_loc, "norway_seiiar_measles_noinfected_2017_b2020.rda"), compress = "xz")

  norway_seiiar_measles_oslo_2017_b2020 <- copy(norway_seiiar_measles_noinfected_2017_b2020)
  norway_seiiar_measles_oslo_2017_b2020[location_code == "municip0301", I := 10]
  norway_seiiar_measles_oslo_2017_b2020[location_code == "municip0301", S := S - I]
  save(norway_seiiar_measles_oslo_2017_b2020, file = file.path(base_loc, "norway_seiiar_measles_oslo_2017_b2020.rda"), compress = "xz")
}
