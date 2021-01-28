.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste(
    "plnr",
    utils::packageDescription("plnr")$Version,
    "https://folkehelseinstituttet.github.io/plnr"
  ))
  packageStartupMessage("\nCommuter model and C++ code developed by Solveig Engebretsen and Andreas Nyg\u00E5rd Osnes")
  packageStartupMessage("Asymmetric mobility model and C++ code developed by Solveig Engebretsen")
  packageStartupMessage("Commuter ported to RCPP by Richard White")
  packageStartupMessage("Branching process developed by Gunnar Ro")
}
