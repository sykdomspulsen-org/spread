#' asymmetric mobility SE1E2IIaR
#'
#' For more information, look at \code{vignette("commuter_model","spread")}
#' @param seiiar_pop Data frame containing `location_code`, `S`, `E`, `I`, `Ia`, and `R` for the entire population
#' @param mobility_matrix List (1 entry for each time period) with each entry containing a data frame with `from`, `to`, `n` for the number of people who travel. Each data frame must be complete.
#' @param seed_matrix matrix of seeding cases per date
#' @param betas Vector (1 entry for each time period). Float, infection parameter, 0.6
#' @param latent_period Float, 2.0
#' @param presymptomatic_period Float 3.0
#' @param infectious_period Float, 5.0
#' @param presymptomatic_relative_infectiousnes Float, relative infectiousness of presymptomatic
#' @param asymptomatic_prob Float, Proportion/probability of asymptomatic given infectious
#' @param asymptomatic_relative_infectiousness Float, Relative infectiousness of asymptomatic infectious
#' @param N Int = 1 int, Number of internal simulations (average taken). This is generally used for parameter fitting.
#' @examples
#' spread::asymmetric_mobility_se1e2iiar(
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
asymmetric_mobility_se1e2iiar <- function(
  seiiar_pop = spread::asymmetric_mobility_dummy_se1e2iiar_pop,
  mobility_matrix = spread::asymmetric_mobility_dummy_se1e2iiar_mobility_matrix,
  seed_matrix= NULL,
  betas = spread::asymmetric_mobility_dummy_se1e2iiar_betas,
  latent_period = 2.0,
  presymptomatic_period = 3.0,
  infectious_period = 5.0,
  presymptomatic_relative_infectiousness = 1.25,
  asymptomatic_prob = 0.4,
  asymptomatic_relative_infectiousness = 0.5,
  N = 1) {
  stopifnot(length(mobility_matrix) == length(betas))

  a1 <- 1 / latent_period
  a2 <- 1 / presymptomatic_period
  gamma <- 1 / infectious_period
  days_simulation <- length(betas) / 4

  retval <- asymmetric_mobility_se1e2iiar_cpp(
    seiiar_pop = seiiar_pop,
    mobility_matrix = mobility_matrix,
    seed_matrix=seed_matrix,
    betas = betas,
    a1 = a1,
    a2 = a2,
    gamma = gamma,
    presymptomaticRelativeInfectiousness = presymptomatic_relative_infectiousness,
    asymptomaticProb = asymptomatic_prob,
    asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness,
    N = N,
    M = days_simulation
  )
  retval <- retval[
    time == 4,
    .(
      b_S = b_S,
      b_E1 = b_E1,
      b_E2 = b_E2,
      b_I = b_I,
      b_Ia = b_Ia,
      b_R = b_R,
      c_S = c_S,
      c_E1 = c_E1,
      c_E2 = c_E2,
      c_I = c_I,
      c_Ia = c_Ia,
      c_R = c_R,
      c_incidence = sum(c_incidence)
    ),
    keyby = .(
      location_code,
      week,
      day
    )
    ]
  retval <- copy(retval)
  retval[, time := "23:59"]
  setcolorder(retval, c("location_code", "week", "day", "time"))

  return(retval)
}


se1e2iiaR_calculate_beta_from_r0 <- function(
  r0,
  a2 = a2,
  gamma = gamma,
  presymptomaticRelativeInfectiousness = presymptomatic_relative_infectiousness,
  asymptomaticProb = asymptomatic_prob,
  asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness) {

  denominator <- presymptomaticRelativeInfectiousness*(1-asymptomaticProb)/a2 +
                 (1 - asymptomaticProb)/gamma +
                (asymptomaticRelativeInfectiousness*asymptomaticProb)/gamma

  beta <- r0 /denominator
  beta <- round(beta, 3)

  return(beta)
}

