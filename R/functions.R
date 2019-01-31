#' CreateDataFiles
#' @export CreateDataFiles
CreateDataFiles <- function() {
  dirData <- system.file("extdata", package = "commuter")
  pop_wo_com <- data.table(readxl::read_excel(file.path(dirData, sprintf("%s.xlsx", "pop_wo_com"))))
  di_edge_list <- data.table(readxl::read_excel(file.path(dirData, sprintf("%s.xlsx", "di_edge_list"))))

  popFiles <- sum(pop_wo_com$pop)
  popReal <- 5000000
  popMultiplier <- 1 # popReal/popFiles

  loc <- pop_wo_com[, c("kommuneNameOld", "location")]
  setnames(loc, "kommuneNameOld", "from")
  loc[, from := factor(from, levels = from)]

  nrow(di_edge_list)
  di_edge_list <- merge(di_edge_list, loc, by = "from")
  nrow(di_edge_list)
  di_edge_list[, from := NULL]
  setnames(di_edge_list, "location", "from")

  setnames(loc, "from", "to")
  nrow(di_edge_list)
  di_edge_list <- merge(di_edge_list, loc, by = "to")
  nrow(di_edge_list)
  di_edge_list[, to := NULL]
  setnames(di_edge_list, "location", "to")
  setcolorder(di_edge_list, c("from", "to", "n"))

  pop_wo_com[, kommuneNameOld := NULL]
  setcolorder(pop_wo_com, c("location", "pop"))



  aMaster <- sykdomspuls::GenNorwayMunicipMerging()
  for (i in unique(aMaster$year)) {
    a <- unique(aMaster[year == i, c("municip", "municipEnd")])
    setnames(a, c("from", "fromNew"))
    nrow(di_edge_list)
    di_edge_list <- merge(di_edge_list, a, by = "from", all.x = T)
    nrow(di_edge_list)
    di_edge_list[!is.na(fromNew), from := fromNew]
    di_edge_list[, fromNew := NULL]

    setnames(a, c("to", "toNew"))
    nrow(di_edge_list)
    di_edge_list <- merge(di_edge_list, a, by = "to", all.x = T)
    nrow(di_edge_list)
    di_edge_list[!is.na(toNew), to := toNew]
    di_edge_list[, toNew := NULL]

    setnames(a, c("location", "locationNew"))
    nrow(pop_wo_com)
    pop_wo_com <- merge(pop_wo_com, a, by = "location", all.x = T)
    nrow(pop_wo_com)
    pop_wo_com[!is.na(locationNew), location := locationNew]
    pop_wo_com[, locationNew := NULL]
  }
  di_edge_list <- di_edge_list[from != to]

  return(list(
    "di_edge_list" = di_edge_list,
    "pop_wo_com" = pop_wo_com
  ))
}

#' SetupCPPAndStructure
#' @import data.table
#' @importFrom readxl read_excel
#' @importFrom sykdomspuls GenNorwayMunicipMerging
#' @importFrom withr with_dir
#' @importFrom processx run
#' @export SetupCPPAndStructure
SetupCPPAndStructure <- function() {
  dirTemp <- RAWmisc::TempDir()

  dirCode <- system.file("cpp", package = "commuter")
  file.copy(file.path(dirCode, list.files(dirCode)), dirTemp)

  x <- CreateDataFiles()
  di_edge_list <- x[["di_edge_list"]]
  pop_wo_com <- x[["pop_wo_com"]]

  for (i in c("di_edge_list", "pop_wo_com")) {
    fwrite(get(i),
      file = file.path(dirTemp, sprintf("%s.txt", i)),
      sep = " ",
      col.names = F
    )
  }

  res <- withr::with_dir(dirTemp, {
    processx::run(
      command = "g++",
      args = c("-std=c++11", "-oinfl_kommuner.exe", "infl_kommuner.cpp"), echo = T
    )
  })

  return(dirTemp)

  aMaster <- aMaster[year == max(year), c("municipEnd", "county", "region")]

  setnames(aMaster, "municipEnd", "location")
  pop_wo_com <- merge(pop_wo_com, aMaster[, c("location", "county")], by = "location")
  pop_wo_com[, location := NULL]
  setnames(pop_wo_com, "county", "location")

  setnames(aMaster, "location", "to")
  di_edge_list <- merge(di_edge_list, aMaster[, c("to", "county")], by = "to")
  di_edge_list[, to := NULL]
  setnames(di_edge_list, "county", "to")

  setnames(aMaster, "to", "from")
  di_edge_list <- merge(di_edge_list, aMaster[, c("from", "county")], by = "from")
  di_edge_list[, from := NULL]
  setnames(di_edge_list, "county", "from")

  pop_wo_com <- pop_wo_com[, .(pop = round(sum(pop) * popMultiplier)), keyby = .(location)]
  di_edge_list <- di_edge_list[, .(n = round(sum(n) * popMultiplier)), keyby = .(from, to)]

  length(unique(pop_wo_com$location))
  length(unique(d$location))

  setorder(di_edge_list, from, to)
  setorder(pop_wo_com, location)

  # regions
  regions <- unique(aMaster[, c("county", "region")])
  setnames(regions, "county", "location")

  # now we add in the age structure
  pop <- fread("https://data.ssb.no/api/v0/dataset/1076.csv?lang=en")
  pop[, county := sprintf("county%s", stringr::str_extract(region, "^[0-9][0-9]"))]
  pop[, region := NULL]
  pop[, agex := as.numeric(stringr::str_extract(age, "^[0-9][0-9][0-9]"))]
  setnames(pop, "Population 1 January. Calculated figures for newest municipality division, by region, age, year and contents", "pop")
  pop[, contents := NULL]
  pop[agex %in% 0:4, age := "0-4"]
  pop[agex %in% 5:19, age := "5-19"]
  pop[agex %in% 20:64, age := "20-64"]
  pop[agex %in% 65:200, age := "65+"]
  pop <- pop[, .(pop = sum(pop)), by = .(age, county)]

  d <- readRDS(fhi::DashboardFolder("data_raw", "resYearLineMunicip_influensa.RDS"))
  d[, pop := NULL]
  d[age %in% c("5-14", "15-19"), age := "5-19"]
  d[age %in% c("20-29", "30-64"), age := "20-64"]
  d <- d[, .(n = sum(n)), by = .(x, week, year, wkyr, county, age)]
  d <- merge(d, pop, by = c("county", "age"))
  setnames(d, "county", "location")
  d[week >= 30, season := sprintf("%s/%s", year, year + 1)]
  d[is.na(season), season := sprintf("%s/%s", year - 1, year)]
  d[week %in% c(25:35), offSeason := mean(n), by = .(season, location)]
  d[, offSeason := mean(offSeason, na.rm = T), by = .(season, location)]
  d[, nMinusOffSeason := floor(n - offSeason)]
  d[nMinusOffSeason < 0, nMinusOffSeason := 0]

  d <- merge(d, regions, by = c("location"))

  pop_wo_com <- merge(pop_wo_com, unique(d[, c("location", "age", "pop")]), by = "location")

  di_edge_list <- merge(di_edge_list, unique(d[, c("location", "age", "pop")]), by.x = "from", by.y = "location", allow.cartesian = T)
  setnames(di_edge_list, "pop", "pop_from")
  setnames(di_edge_list, "age", "age_from")
  di_edge_list <- merge(di_edge_list, unique(d[, c("location", "age", "pop")]), by.x = "to", by.y = "location", allow.cartesian = T)
  setnames(di_edge_list, "pop", "pop_to")
  setnames(di_edge_list, "age", "age_to")
  di_edge_list[!(age_from %in% c("20-64") & age_to %in% c("20-64")), n := 0]
  di_edge_list[from == to, n := 0]
  # di_edge_list[n>0,pop_from_total:=sum(pop_from),by=.(from,to)]
  # di_edge_list[n>0,n_total:=sum(n),by=.(from,to)]
  # di_edge_list[!is.na(pop_from_total),n:=ceiling(n*pop_from/pop_from_total)]

  di_edge_list[, n_commuting := sum(n), by = .(from)]

  for (i in colnames(mixingMatrix)) for (j in colnames(mixingMatrix)) {
      if (i == j) next
      di_edge_list[from == to & age_from == i & age_to == j & n == 0, n := floor(pop_from * 0.1 * mixingMatrix[i, j])]
    }

  di_edge_list[from == "county01" & to == "county02"]
  di_edge_list[from == "county01" & to == "county01"]

  pop_wo_com <- di_edge_list[, .(n_commuting = sum(n)), by = .(from, age_from)]

  di_edge_list[, to := sprintf("%s_%s", to, age_to)]
  di_edge_list[, from := sprintf("%s_%s", from, age_from)]
  di_edge_list <- di_edge_list[from != to, c("from", "to", "n")]
  setorder(di_edge_list, from, to)

  pop_wo_com <- merge(pop_wo_com, unique(d[year == max(year), c("location", "age", "pop")]),
    by.x = c("from", "age_from"), by.y = c("location", "age")
  )
  pop_wo_com[, pop := pop - n_commuting]
  pop_wo_com[, n_commuting := NULL]
  pop_wo_com[, from := sprintf("%s_%s", from, age_from)]
  pop_wo_com[, age_from := NULL]
  setorder(pop_wo_com, from)

  for (i in c("di_edge_list", "pop_wo_com")) {
    fwrite(get(i),
      file = file.path(CONFIG_DIR$DIR_TMP, sprintf("%s.txt", i)),
      sep = " ",
      col.names = F
    )
  }

  res <- withr::with_dir(CONFIG_DIR$DIR_TMP, {
    processx::run(
      command = "g++",
      args = c("-std=c++11", "-oinfl_kommuner.exe", "infl_kommuner_region_betas.cpp"), echo = T
    )
  })


  # retval <- aMaster[year==max(year) & municipEnd %in% pop_wo_com$location,c("region","municipEnd")]
  # setnames(retval,"municipEnd","location")
  setorder(regions, location)
  regions[, betas := 0.6]
  return(list("regions" = regions, "d" = d))
}


#' RunSim
#' @param id fileid
#' @param startVals a
#' @param R0 R0
#' @param a days latent
#' @param gammaTheoretical days infectious
#' @param gammaEffective days infectious
#' @param asymptomaticProb Probability someone is asymptomatic
#' @param asymptomaticRelativeInfectiousness How infectious is this person relative to symptomatic
#' @param M days simulated
#' @param verbose True/False
#' @export
RunSim <- function(
                   id = 1,
                   startVals,
                   R0 = 1.88,
                   a = 1.9,
                   gammaTheoretical = 3,
                   gammaEffective = gammaTheoretical,
                   asymptomaticProb = 0.33,
                   asymptomaticRelativeInfectiousness = 0.5,
                   M = 100,
                   verbose = F) {
  fwrite(startVals[, "value"],
    file = file.path(dirTemp, sprintf("start_infected_%s.txt", id)),
    sep = " ",
    col.names = F
  )

  symptomaticProb <- 1 - asymptomaticProb
  beta <- R0 / gammaTheoretical / (symptomaticProb + asymptomaticProb * asymptomaticRelativeInfectiousness)
  beta <- round(beta, 3)

  # beta;  // infection parameter, 0.6
  # a;  // 1/latent period, 1/1.9
  # gamma; // 1/infectious period, 1/3
  # asymptomaticProb
  # asymptomaticRelativeInfectiousness
  res <- withr::with_dir(dirTemp, {
    processx::run(
      command = file.path(dirTemp, "infl_kommuner.exe"),
      args = c(
        as.character(id),
        as.character(beta),
        as.character(a),
        as.character(gammaEffective),
        as.character(asymptomaticProb),
        as.character(asymptomaticRelativeInfectiousness),
        as.character(M)
      ),
      echo_cmd = verbose,
      echo = verbose
    )
  })

  loc <- fread(file.path(dirTemp, "pop_wo_com.txt"))
  setnames(loc, c("location", "pop"))
  loc[, kn := 1:.N - 1]

  m <- fread(file.path(dirTemp, sprintf("cpp_res_series_%s.txt", id)))
  setnames(m, c("kn", "S", "E", "I", "IA", "R", "INCIDENCE"))
  # m <- m[seq(1,nrow(d),2)]
  m[, day := rep(1:(.N / 2), each = 2), by = kn]

  m <- merge(m, loc, by = "kn")

  m <- m[, .(
    S = base::mean(S),
    E = base::mean(E),
    I = base::mean(I),
    IA = base::mean(IA),
    R = base::mean(R),
    INCIDENCE = sum(INCIDENCE)
  ),
  by = .(
    location, day
  )
  ]

  return(m)
}
