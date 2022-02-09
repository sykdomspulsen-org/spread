.onAttach <- function(libname, pkgname) {
  version <- tryCatch(
    utils::packageDescription("spread", fields = "Version"),
    warning = function(w){
      1
    }
  )

  packageStartupMessage(paste0(
    "spread ",
    version,
    "\n",
    "https://docs.sykdomspulsen.no/spread"
  ))
}
