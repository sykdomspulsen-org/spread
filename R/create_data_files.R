#' Names of areas in Norway that currently exist.
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{location_name}{Location name.}
#' }
#' @source \url{https://snl.no/kommunenummer}
"norway_di_edge_list_2017"

#' Names of areas in Norway that previously existed.
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{location_name}{Location name.}
#' }
#' @source \url{https://no.wikipedia.org/wiki/Liste_over_norske_kommunenummer}
"norway_pop_wo_com_2017"

#' Names of areas in Norway that currently exist.
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{location_name}{Location name.}
#' }
#' @source \url{https://snl.no/kommunenummer}
"start_points_blank"

#' Names of areas in Norway that previously existed.
#'
#' @format
#' \describe{
#' \item{location_code}{Location code.}
#' \item{location_name}{Location name.}
#' }
#' @source \url{https://no.wikipedia.org/wiki/Liste_over_norske_kommunenummer}
"start_points_oslo"

create_data_files_norway_2017 <- function(base_loc) {

  . <- NULL
  year <- NULL
  n <- NULL
  from <- NULL
  to <- NULL
  level <- NULL
  pop <- NULL
  location_code <- NULL

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

  nrow(di_edge_list)
  di_edge_list <- merge(di_edge_list, norwayMunicipMerging[, c("year", "municip_code_current", "municip_code_original")],
                        by.x = c("from", "year"),
                        by.y = c("municip_code_original", "year")
  )
  nrow(di_edge_list)
  di_edge_list[, from := NULL]
  setnames(di_edge_list, "municip_code_current", "from")

  nrow(di_edge_list)
  di_edge_list <- merge(di_edge_list, norwayMunicipMerging[, c("year", "municip_code_current", "municip_code_original")],
                        by.x = c("to", "year"),
                        by.y = c("municip_code_original", "year")
  )
  nrow(di_edge_list)
  di_edge_list[, to := NULL]
  setnames(di_edge_list, "municip_code_current", "to")

  di_edge_list <- di_edge_list[, .(n = sum(n)), keyby = .(
    from,
    to
  )]
  di_edge_list <- di_edge_list[from != to]
  di_edge_list <- di_edge_list[n > 0]

  commuters <- di_edge_list[, .(n = sum(n)), keyby = .(from)]

  pop_wo_com <- fhidata::norway_population_current[year == 2017 & level=="municipality", .(
    pop = sum(pop)
  ), by = .(location_code)]

  pop_wo_com <- merge(pop_wo_com, commuters, by.x = c("location_code"), by.y = c("from"), all.x = T)
  pop_wo_com[!is.na(n), pop := pop - n]
  pop_wo_com[, n := NULL]

  norway_di_edge_list_2017 <- di_edge_list
  save(norway_di_edge_list_2017, file = file.path(base_loc, "norway_di_edge_list_2017.rda"), compress = "xz")

  norway_pop_wo_com_2017 <- pop_wo_com
  save(norway_pop_wo_com_2017, file = file.path(base_loc, "norway_pop_wo_com_2017.rda"), compress = "xz")

  start_points_blank <- copy(norway_pop_wo_com_2017)
  start_points_blank[,pop:=NULL]
  start_points_blank[,I:=0]
  save(start_points_blank, file = file.path(base_loc, "start_points_blank.rda"), compress = "xz")

  start_points_oslo <- copy(start_points_blank)
  start_points_oslo[location_code=="municip0301",I:=10]
  save(start_points_oslo, file = file.path(base_loc, "start_points_oslo.rda"), compress = "xz")

}
