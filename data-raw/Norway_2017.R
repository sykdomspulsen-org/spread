devtools::load_all()

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

  di_edge_list[, from := sprintf("municip_nor%s", from)]
  di_edge_list[, to := sprintf("municip_nor%s", to)]

  norwayMunicipMerging <- csdata::nor_locations_redistricting()
  sum(di_edge_list$n)
  di_edge_list <- merge(
    di_edge_list,
    norwayMunicipMerging,
    by.x = c("from", "year"),
    by.y = c("location_code_original", "calyear"),
    all.x = T
  )
  sum(di_edge_list$n)
  di_edge_list[, n := n * weighting]
  di_edge_list <- di_edge_list[!is.na(n) & n > 0]
  sum(di_edge_list$n, na.rm = T)
  di_edge_list[, from := NULL]
  di_edge_list[, weighting := NULL]
  setnames(di_edge_list, "location_code_current", "from")

  di_edge_list <- merge(
    di_edge_list,
    norwayMunicipMerging,
    by.x = c("to", "year"),
    by.y = c("location_code_original", "calyear"),
    all.x = T
  )
  di_edge_list[, n := n * weighting]
  di_edge_list <- di_edge_list[!is.na(n) & n > 0]
  sum(di_edge_list$n, na.rm = T)
  di_edge_list[, to := NULL]
  di_edge_list[, weighting := NULL]
  setnames(di_edge_list, "location_code_current", "to")

  di_edge_list <- di_edge_list[, .(n = round(sum(n))), keyby = .(
    from,
    to
  )]
  di_edge_list <- di_edge_list[from != to]
  di_edge_list <- di_edge_list[n > 0]
  sum(di_edge_list$n)
  seiiar <- csdata::nor_population_by_age_cats()[calyear == 2017 & granularity_geo == "municip", .(
    S = sum(pop_jan1_n),
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

# create_data_files_norway_2017

x <- create_blank_norway_2017()

seiiar <- x[["seiiar"]]
nor_commuters_2017_b2020 <- x[["commuters"]]

usethis::use_data(nor_commuters_2017_b2020, overwrite = TRUE, version = 3, compress = "xz")

nor_seiiar_noinfected_2017_b2020 <- seiiar
usethis::use_data(nor_seiiar_noinfected_2017_b2020, overwrite = TRUE, version = 3, compress = "xz")

nor_seiiar_oslo_2017_b2020 <- copy(seiiar)
nor_seiiar_oslo_2017_b2020[location_code == "municip_nor0301", I := 10]
nor_seiiar_oslo_2017_b2020[location_code == "municip_nor0301", S := S - I]
usethis::use_data(nor_seiiar_oslo_2017_b2020, overwrite = TRUE, version = 3, compress = "xz")

# measles
vax_prev <- nor_childhood_vax_b2020[calyear == 2016 & stringr::str_detect(location_code, "^municip") & vaccine == "measles"]
nor_seiiar_measles_noinfected_2017_b2020 <- convert_blank_seiiar_with_vax(seiiar, vax_prev)
usethis::use_data(nor_seiiar_measles_noinfected_2017_b2020, overwrite = TRUE, version = 3, compress = "xz")

nor_seiiar_measles_oslo_2017_b2020 <- copy(nor_seiiar_measles_noinfected_2017_b2020)
nor_seiiar_measles_oslo_2017_b2020[location_code == "municip_nor0301", I := 10]
nor_seiiar_measles_oslo_2017_b2020[location_code == "municip_nor0301", S := S - I]
usethis::use_data(nor_seiiar_measles_oslo_2017_b2020, overwrite = TRUE, version = 3, compress = "xz")

# norway as a single entity
single_entity_fake_commuters_2017 <- data.table(from = "norge", to = "x", n = 1)
usethis::use_data(single_entity_fake_commuters_2017, overwrite = TRUE, version = 3, compress = "xz")

single_entity_seiiar_2017 <- nor_seiiar_noinfected_2017_b2020[, .(
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
usethis::use_data(single_entity_seiiar_2017, overwrite = TRUE, version = 3, compress = "xz")
