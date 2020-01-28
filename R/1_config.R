set_config <- function() {
  if (!foreach::getDoParRegistered()) foreach::registerDoSEQ()
}
