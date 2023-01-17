library(data.table)

d <- fread("data-raw/SYSVAK_2019-04-09-14-17.csv")
d[GEO == 0, location_code := "nation_nor"]
d[GEO > 0 & GEO < 100, location_code := glue::glue("county_nor{X}", X = formatC(GEO, width = 2, flag = "0"))]
d[GEO >= 100, location_code := glue::glue("municip_nor{X}", X = formatC(GEO, width = 4, flag = "0"))]
d <- d[!is.na(location_code)]

d[, calyear := as.numeric(stringr::str_extract(AAR, "^[0-9][0-9][0-9][0-9]")) + 2]
d <- d[SPVFLAGG == 0 & ALDER == "16_16"]
d[, year_merging := 2019]

dnorge <- d[stringr::str_detect(location_code, "nation_nor")]
dcounty <- d[stringr::str_detect(location_code, "county")]
dmunicip <- d[stringr::str_detect(location_code, "municip")]

d_new <- merge(
  d,
  csdata::nor_locations_redistricting(),
  by.x = c("location_code", "year_merging"),
  by.y = c("location_code_original", "calyear")
)

d <- d_new[,.(
  RATE = sum(RATE * weighting) / sum(weighting)
), keyby = .(
  location_code_current,
  KJONN,
  ALDER,
  VAKSINE,
  calyear
)]

setnames(d, "location_code_current", "location_code")

d[, age := "016"]
setnames(d, "RATE", "vaccination_coverage_pr100")
d[, vaccine := as.character(forcats::fct_recode(VAKSINE,
                                            "measles" = "Meslinger",
                                            "diphtheria" = "Difteri",
                                            "hpv" = "HPV",
                                            "pertussis" = "Kikhoste",
                                            "mumps" = "Kusma",
                                            "poliomyelitis" = "Poliomyelitt",
                                            "rubella" = "Rodehunder",
                                            "tetanus" = "Stivkrampe"
))]
d <- d[, c("location_code", "calyear", "age", "vaccine", "vaccination_coverage_pr100")]
d[, imputed := FALSE]
national_results <- d[location_code == "nation_nor", .(
  national_vaccination_coverage_pr100 = mean(vaccination_coverage_pr100)
),
keyby = .(calyear, vaccine)
]
skeleton <- data.table(expand.grid(
  location_code = csdata::nor_locations_names()[granularity_geo %in% c("county", "municip")]$location_code,
  calyear = unique(d$calyear),
  vaccine = unique(d$vaccine),
  stringsAsFactors = F
))
d <- merge(d, skeleton,
           by = c("location_code", "calyear", "vaccine"), all = T
)
d <- merge(d, national_results, by = c("calyear", "vaccine"))
d[is.na(age), age := "016"]
d[is.na(vaccination_coverage_pr100), imputed := TRUE]
d[is.na(vaccination_coverage_pr100), vaccination_coverage_pr100 := national_vaccination_coverage_pr100]
d[, vaccination_coverage_pr1 := vaccination_coverage_pr100 / 100]
d[, national_vaccination_coverage_pr100 := NULL]

setcolorder(d, c(
  "calyear",
  "location_code",
  "age",
  "vaccine",
  "vaccination_coverage_pr100",
  "vaccination_coverage_pr1",
  "imputed"
))

nor_childhood_vax_b2020 <- d

usethis::use_data(nor_childhood_vax_b2020, overwrite = TRUE, version = 3, compress = "xz")

