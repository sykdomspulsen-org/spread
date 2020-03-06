#' asymmetric mobility
#'
#' For more information, look at \code{vignette("commuter_model","spread")}
#' @param seiiar_pop Data frame containing `location_code`, `S`, `E`, `I`, `Ia`, and `R` for the entire population
#' @param mobility_matrix List (1 entry for each time period) with each entry containing a data frame with `from`, `to`, `n` for the number of people who travel. Each data frame must be complete.
#' @param betas Vector (1 entry for each time period). Float, infection parameter, 0.6
#' @param latent_period Float, 1.9
#' @param infectious_period Float, 3
#' @param asymptomatic_prob Float, Proportion/probability of asymptomatic given infectious
#' @param asymptomatic_relative_infectiousness Float, Relative infectiousness of asymptomatic infectious
#' @param N Int = 1 int, Number of internal simulations (average taken). This is generally used for parameter fitting.
#' @examples
#' spread::asymmetric_mobility(
#'   seiiar_pop = spread::asymmetric_mobility_dummy_seiiar_pop,
#'   mobility_matrix = spread::asymmetric_mobility_dummy_mobility_matrix,
#'   betas = spread::asymmetric_mobility_dummy_betas,
#'   latent_period = 1.9,
#'   infectious_period = 3.0,
#'   asymptomatic_prob = 0,
#'   asymptomatic_relative_infectiousness = 0,
#'   N = 1
#' )
#' @import data.table
#' @export
asymmetric_mobility <- function(
  seiiar_pop = spread::asymmetric_mobility_dummy_seiiar_pop,
  mobility_matrix = spread::asymmetric_mobility_dummy_mobility_matrix,
  betas = spread::asymmetric_mobility_dummy_betas,
  latent_period = 1.9,
  infectious_period = 3.0,
  asymptomatic_prob = 0,
  asymptomatic_relative_infectiousness = 0,
  N = 1) {

  stopifnot(length(mobility_matrix) == length(betas))

  a <- 1 / latent_period
  gamma <- 1 / infectious_period
  days_simulation <- length(betas)/4

  retval <- asymmetric_mobility_cpp(
    seiiar_pop = seiiar_pop,
    mobility_matrix = mobility_matrix,
    betas = betas,
    a = a,
    gamma = gamma,
    asymptomaticProb = asymptomatic_prob,
    asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness,
    N = N,
    M = days_simulation
  )
  retval <- copy(retval)
  return(retval)
}
