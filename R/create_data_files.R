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
  E <- NULL
  Ia <- NULL
  R <- NULL

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

  # norway as a single entity
  single_entity_fake_commuters_2017 <- data.table(from = "norge", to = "x", n = 1)
  save(single_entity_fake_commuters_2017, file = file.path(base_loc, "single_entity_fake_commuters_2017.rda"), compress = "xz")

  single_entity_seiiar_2017 <- spread::norway_seiiar_noinfected_2017_b2020[, .(
    location_code = "norge",
    S = sum(S),
    E = sum(E),
    I = sum(I),
    Ia = sum(Ia),
    R = sum(R)
  )]
  single_entity_seiiar_2017 <- rbind(single_entity_seiiar_2017, single_entity_seiiar_2017)
  single_entity_seiiar_2017[2, location_code := "x"]
  single_entity_seiiar_2017[2, S := 1]
  save(single_entity_seiiar_2017, file = file.path(base_loc, "single_entity_seiiar_2017.rda"), compress = "xz")
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

create_asymmetric_mobility_se1e2iiar_dummy_files <- function(base_loc) {
  se1e2iiar_pop <- data.table::data.table(
    "location_code" = c("a", "b", "c"),
    "S" = c(1000, 1000, 2000),
    "E1" = c(0, 0, 0),
    "E2" = c(0, 0, 0),
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

  dayEach <- rep(1:5, each = 4)
  timeEach <- rep(c(0, 6, 12, 18), 5)
  location_codes <- c(rep("a", 5 * 4), rep("b", 5 * 4), rep("c", 5 * 4))
  betas <- c(rep(0.6, 5 * 4), rep(0.4, 5 * 2), rep(0.2, 5 * 2), rep(0.7, 5 * 3), rep(0.4, 5))
  days <- rep(dayEach, 3)
  times <- rep(timeEach, 3)
  betas <- data.table("location_code" = location_codes, "day" = days, "time" = times, "beta" = betas)

  asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_pop <- se1e2iiar_pop
  save(asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_pop, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_pop.rda"), compress = "xz")

  asymmetric_mobility_se1e2iiar_dummy_mobility_matrix <- mobility_matrix
  save(asymmetric_mobility_se1e2iiar_dummy_mobility_matrix, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_dummy_mobility_matrix.rda"), compress = "xz")

  asymmetric_mobility_se1e2iiar_dummy_betas <- betas
  save(asymmetric_mobility_se1e2iiar_dummy_betas, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_dummy_betas.rda"), compress = "xz")

  asymmetric_mobility_se1e2iiar_dummy_dynamic_seeds <- data.table::data.table(
    "location_code" = c("a"),
    "day" = c(3, 7),
    "n" = c(3, 5)
  )
  save(asymmetric_mobility_se1e2iiar_dummy_dynamic_seeds, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_dummy_dynamic_seeds.rda"), compress = "xz")
}






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

create_asymmetric_mobility_se1e2iiar_2strains_dummy_files <- function(base_loc) {
  se1e2iiar_2strains_pop <- data.table::data.table(
    "location_code" = c("a", "b", "c"),
    "S" = c(1000, 1000, 2000),
    "E1" = c(0, 0, 0),
    "E2" = c(0, 0, 0),
    "I" = c(50, 0, 0),
    "Ia" = c(0, 0, 0),
    "R" = c(0, 0, 0),
    "E1_b" = c(0, 0, 0),
    "E2_b" = c(0, 0, 0),
    "I_b" = c(50, 0, 0),
    "Ia_b" = c(0, 0, 0)
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

  dayEach <- rep(1:5, each = 4)
  timeEach <- rep(c(0, 6, 12, 18), 5)
  location_codes <- c(rep("a", 5 * 4), rep("b", 5 * 4), rep("c", 5 * 4))
  betas <- c(rep(0.6, 5 * 4), rep(0.4, 5 * 2), rep(0.2, 5 * 2), rep(0.7, 5 * 3), rep(0.4, 5))
  days <- rep(dayEach, 3)
  times <- rep(timeEach, 3)
  betas <- data.table("location_code" = location_codes, "day" = days, "time" = times, "beta" = betas)

  asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop <- se1e2iiar_2strains_pop
  save(asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop.rda"), compress = "xz")

  asymmetric_mobility_se1e2iiar_dummy_mobility_matrix <- mobility_matrix
  save(asymmetric_mobility_se1e2iiar_dummy_mobility_matrix, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_dummy_mobility_matrix.rda"), compress = "xz")

  asymmetric_mobility_se1e2iiar_dummy_betas <- betas
  save(asymmetric_mobility_se1e2iiar_dummy_betas, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_dummy_betas.rda"), compress = "xz")

  asymmetric_mobility_se1e2iiar_2strains_dummy_dynamic_seeds <- data.table::data.table(
    "location_code" = c("a"),
    "day" = c(3, 7, 8),
    "n" = c(3, 5, 0),
    "n_b" = c(1, 0, 5)
  )
  save(asymmetric_mobility_se1e2iiar_2strains_dummy_dynamic_seeds, file = file.path(base_loc, "asymmetric_mobility_se1e2iiar_2strains_dummy_dynamic_seeds.rda"), compress = "xz")
}
