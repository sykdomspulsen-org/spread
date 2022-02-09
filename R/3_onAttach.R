.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "spread ",
    # utils::packageDescription("spread")$Version,
    "\n",
    "https://docs.sykdomspulsen.no/spread"
  ))
}
