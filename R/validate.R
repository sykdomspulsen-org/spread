validate_seiiar <- function(seiiar) {
  wanted_names <- c(
    "location_code",
    "S",
    "E",
    "I",
    "Ia",
    "R"
  )

  if (!all.equal(names(seiiar), wanted_names)) {
    stop(glue::glue("names are not c({glue::glue_collapse(wanted_names, sep=', ')})"))
  }
}

validate_vax_prev <- function(vax) {
  needed_names <- c("location_code", "proportion")
  if (sum(!needed_names %in% names(vax)) > 0) {
    stop(glue::glue("names are not c({glue::glue_collapse(needed_names, sep=', ')})"))
  }
}
