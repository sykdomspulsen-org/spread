#' Daily number of commuters from/to municipalities in Norway in 2017
#'
#' @format
#' \describe{
#' \item{from}{Location code.}
#' \item{to}{Location name.}
#' \item{n}{Number of people.}
#' }
#' @source \url{https://www.ssb.no/statbank/table/03321}
"norway_commuters_2017"

#' SEIIaR data.frame for Norway with no one infected and everyone susceptible.
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
"norway_seiiar_noinfected_2017"

#' SEIIaR data.frame for Norway with 10 people infected in Oslo and everyone susceptible.
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
"norway_seiiar_oslo_2017"

#' SEIIaR data.frame for Norway with no one infected and real measles susceptibility.
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
"norway_seiiar_measles_noinfected_2017"

#' SEIIaR data.frame for Norway with 10 people infected in Oslo and real measles susceptibility.
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
"norway_seiiar_measles_oslo_2017"

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

  norwayMunicipMerging <- fhidata::norway_municip_merging

  di_edge_list[fhidata::norway_municip_merging, on = c("from==municip_code_original", "year==year"), municip_code_current := municip_code_current]
  di_edge_list[, from := NULL]
  setnames(di_edge_list, "municip_code_current", "from")

  di_edge_list[fhidata::norway_municip_merging, on = c("to==municip_code_original", "year==year"), municip_code_current := municip_code_current]
  di_edge_list[, to := NULL]
  setnames(di_edge_list, "municip_code_current", "to")

  di_edge_list <- di_edge_list[, .(n = sum(n)), keyby = .(
    from,
    to
  )]
  di_edge_list <- di_edge_list[from != to]
  di_edge_list <- di_edge_list[n > 0]

  seiiar <- fhidata::norway_population_current[year == 2017 & level == "municipality", .(
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
#' recovered
#' For more information, look at \code{vignette("including_vax","spread")}
#' @param seiiar SEIIAR data.table representing a fully susceptible population
#' @param vax data.table containing proportion of people vaccinated per location code
#' @export
convert_blank_seiiar_with_vax <- function(seiiar, vax){
  R <- NULL
  S <- NULL
  proportion <- NULL
  validate_seiiar(seiiar)
  validate_vax_prev(vax)

  retval <- copy(seiiar)
  setDT(retval)
  retval[vax,on="location_code",R:=round(S*proportion)]
  retval[,S:=S-R]

  return(retval)
}

create_data_files_norway_2017 <- function(base_loc = "data") {
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
  norway_commuters_2017 <- x[["commuters"]]

  save(norway_commuters_2017, file = file.path(base_loc, "norway_commuters_2017.rda"), compress = "xz")

  norway_seiiar_noinfected_2017 <- seiiar
  save(norway_seiiar_noinfected_2017, file = file.path(base_loc, "norway_seiiar_noinfected_2017.rda"), compress = "xz")

  norway_seiiar_oslo_2017 <- copy(seiiar)
  norway_seiiar_oslo_2017[location_code == "municip0301", I := 10]
  norway_seiiar_oslo_2017[location_code == "municip0301", S := S - I]
  save(norway_seiiar_oslo_2017, file = file.path(base_loc, "norway_seiiar_oslo_2017.rda"), compress = "xz")

  # measles
  vax_prev <- fhidata::norway_childhood_vax[year==2016 & stringr::str_detect(location_code,"^municip") & vax=="measles"]
  norway_seiiar_measles_noinfected_2017 <- convert_blank_seiiar_with_vax(seiiar, vax_prev)

  save(norway_seiiar_measles_noinfected_2017, file = file.path(base_loc, "norway_seiiar_measles_noinfected_2017.rda"), compress = "xz")

  norway_seiiar_measles_oslo_2017 <- copy(norway_seiiar_measles_noinfected_2017)
  norway_seiiar_measles_oslo_2017[location_code == "municip0301", I := 10]
  norway_seiiar_measles_oslo_2017[location_code == "municip0301", S := S - I]
  save(norway_seiiar_measles_oslo_2017, file = file.path(base_loc, "norway_seiiar_measles_oslo_2017.rda"), compress = "xz")
}


