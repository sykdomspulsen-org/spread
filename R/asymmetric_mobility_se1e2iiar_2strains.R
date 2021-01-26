#' asymmetric mobility SE1E2IIaR two strains
#'
#' For more information, look at \code{vignette("commuter_model","spread")}
#'
#' @param se1e2iiar_2strains_pop Data frame containing `location_code`, `S`, `E1`, `E2`, `I`, `Ia`, `R`, `E1_b`, `E2_b`, `I_b`, and`Ia_b` for the entire population
#' @param mobility_matrix List (1 entry for each time period) with each entry containing a data frame with `from`, `to`, `n` for the number of people who travel. Each data frame must be complete.
#' @param dynamic_seeds Data.table containing `location_code`, `day`, `n`, `n_b` for dynamic seeding (or `NULL`)
#' @param betas Data table with `locaton code`, `day`, `time`, `beta` (1 entry for each location for each time period). Float, infection parameter, 0.6
#' @param inputSeed Integer with seed for cpp code
#' @param latent_period Float, 3.0
#' @param presymptomatic_period Float 2.0
#' @param infectious_period Float, 5.0
#' @param b_relative_infectiousness Float, relative infectiousness of strain B
#' @param presymptomatic_relative_infectiousness Float, relative infectiousness of presymptomatic
#' @param asymptomatic_prob Float, Proportion/probability of asymptomatic given infectious
#' @param asymptomatic_relative_infectiousness Float, Relative infectiousness of asymptomatic infectious
#' @param N Int = 1 int, Number of internal simulations (average taken). This is generally used for parameter fitting.
#' @return A data.table containing the following variables:
#' \describe{
#' \item{location_code}{Location code}
#' \item{week}{Week number}
#' \item{day}{Day number}
#' \item{time}{Time of reporting (23:59)}
#' \item{b_S}{Susceptibles belonging to location code}
#' \item{b_E1}{E1s belonging to location code}
#' \item{b_E2}{E2s belonging to location code}
#' \item{b_I}{Infectious (symptomatic) belonging to location code}
#' \item{b_Ia}{Infectious (asymptomatic) belonging to location code}
#' \item{b_R}{Recovered belonging to location code}
#' \item{b_E1_b}{E1_bs belonging to location code}
#' \item{b_E2_b}{E2_bs belonging to location code}
#' \item{b_I_b}{Infectious (symptomatic) strain B belonging to location code}
#' \item{b_Ia_b}{Infectious (asymptomatic) strain B belonging to location code}
#' \item{c_S}{Susceptibles currently in this location code}
#' \item{c_E1}{E1s currently in this location code}
#' \item{c_E2}{E2s currently in this location code}
#' \item{c_I}{Infectious (symptomatic) currently in this location code}
#' \item{c_Ia}{Infectious (asymptomatic) currently in this location code}
#' \item{c_R}{Recovered currently in this location code}
#' \item{c_E1_b}{E1_bs currently in this location code}
#' \item{c_E2_b}{E2_bs currently in this location code}
#' \item{c_I_b}{Infectious (symptomatic) strain B currently in this location code}
#' \item{c_Ia_b}{Infectious (asymptomatic) strain B currently in this location code}
#' \item{c_symp_incidence_a}{Transition from E2 to I currently in this location code}
#' \item{c_asymp_incidence_a}{Transition from E1 to Ia currently in this location code}
#' \item{c_symp_incidence_b}{Transition from E2_b to I_b currently in this location code}
#' \item{c_asymp_incidence_b}{Transition from E1_b to Ia_b currently in this location code}
#' \item{c_symp_incidence}{Transition from E2 and E2_b to I and I_b currently in this location code}
#' \item{c_asymp_incidence}{Transition from E1 and E1_b to Ia and Ia_b currently in this location code}
#' }
#' @examples
#' spread::asymmetric_mobility_se1e2iiar_2strains(
#'   se1e2iiar_2strains_pop = spread::asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop,
#'   mobility_matrix = spread::asymmetric_mobility_se1e2iiar_dummy_mobility_matrix,
#'   dynamic_seeds = spread::asymmetric_mobility_se1e2iiar_2strains_dummy_dynamic_seeds,
#'   betas = spread::asymmetric_mobility_se1e2iiar_dummy_betas,
#'   inputSeed = 123,
#'   latent_period = 3.0,
#'   presymptomatic_period = 2.0,
#'   infectious_period = 5.0,
#'   b_relative_infectiousness = 1.5,
#'   presymptomatic_relative_infectiousness = 1,
#'   asymptomatic_prob = 0,
#'   asymptomatic_relative_infectiousness = 0,
#'   N = 1
#' )
#' @import data.table
#' @export
asymmetric_mobility_se1e2iiar_2strains <- function(
  se1e2iiar_2strains_pop = spread::asymmetric_mobility_se1e2iiar_dummy_se1e2iiar_2strains_pop,
  mobility_matrix = spread::asymmetric_mobility_se1e2iiar_dummy_mobility_matrix,
  dynamic_seeds = NULL,
  betas = spread::asymmetric_mobility_se1e2iiar_dummy_betas,
  inputSeed = as.numeric(Sys.time()),
  latent_period = 3.0,
  presymptomatic_period = 2.0,
  infectious_period = 5.0,
  b_relative_infectiousness = 1.5,
  presymptomatic_relative_infectiousness = 1.25,
  asymptomatic_prob = 0.4,
  asymptomatic_relative_infectiousness = 0.5,
  N = 1) {
  stopifnot(length(mobility_matrix) >= length(unique(betas$day)))

  a1 <- 1 / latent_period
  a2 <- 1 / presymptomatic_period
  gamma <- 1 / infectious_period
  days_simulation <- length(unique(betas$day))

  if (!inherits(betas, "data.table")) {
    betas <- data.table(betas)
  }
  stopifnot(identical(
    names(betas),
    c("location_code", "day", "time", "beta")
  ))

  if (!inherits(se1e2iiar_2strains_pop, "data.table")) {
    se1e2iiar_2strains_pop <- data.table(se1e2iiar_2strains_pop)
  }
  stopifnot(identical(
    names(se1e2iiar_2strains_pop),
    c("location_code", "S", "E1", "E2", "I", "Ia", "R", "E1_b", "E2_b", "I_b", "Ia_b")
  ))
  data.table::setorder(se1e2iiar_2strains_pop, location_code)

  stopifnot(inherits(mobility_matrix, "list"))
  for (i in seq_along(mobility_matrix)) {
    if (!inherits(mobility_matrix[[i]], "data.table")) {
      mobility_matrix[[i]] <- data.table(mobility_matrix[[i]])
    }
    stopifnot(identical(
      names(mobility_matrix[[i]]),
      c("from", "to", "n")
    ))
    data.table::setorder(mobility_matrix[[i]], from, to)
  }

  # create seed_matrix from dynamic_seeds
  location_codes <- se1e2iiar_2strains_pop$location_code
  if(!is.null(dynamic_seeds)){
    dynamic_seeds_a <- data.table::copy(dynamic_seeds)[,"n_b":=NULL]
    dynamic_seeds_b <- data.table::copy(dynamic_seeds)[,"n":=NULL]
    names(dynamic_seeds_b)[names(dynamic_seeds_b) == "n_b"] = "n"
  }
  else{
    dynamic_seeds_a = NULL
    dynamic_seeds_b = NULL
  }
  seed_matrix <- convert_dynamic_seeds_to_seed_matrix(
    dynamic_seeds = dynamic_seeds_a,
    location_codes = location_codes,
    days = 1:days_simulation
  )

  seed_matrix_b <- convert_dynamic_seeds_to_seed_matrix(
    dynamic_seeds = dynamic_seeds_b,
    location_codes = location_codes,
    days = 1:days_simulation
  )

  # create beta_matrix from betas data frame
  location_codes <- se1e2iiar_2strains_pop$location_code
  beta_matrix <- convert_beta_to_matrix(
    betas = betas,
    location_codes = location_codes,
    days = 1:days_simulation,
    times = c(0, 6, 12, 18)
  )


  retval <- asymmetric_mobility_se1e2iiar_2strains_cpp(
    se1e2iiar_2strains_pop = se1e2iiar_2strains_pop,
    mobility_matrix = mobility_matrix,
    seed_matrix = seed_matrix,
    seed_matrix_b = seed_matrix_b,
    betas = beta_matrix,
    inputSeed = inputSeed,
    a1 = a1,
    a2 = a2,
    gamma = gamma,
    relativeInfectiousnessB = b_relative_infectiousness,
    presymptomaticRelativeInfectiousness = presymptomatic_relative_infectiousness,
    asymptomaticProb = asymptomatic_prob,
    asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness,
    N = N,
    M = days_simulation
  )
  retval = cbind(retval[[1]], retval[[2]])
  retval_incidence <- retval[, .(
    c_symp_incidence_a = sum(c_symp_incidence_a),
    c_asymp_incidence_a = sum(c_asymp_incidence_a),
    c_symp_incidence_b = sum(c_symp_incidence_b),
    c_asymp_incidence_b = sum(c_asymp_incidence_b),
    c_symp_incidence = sum(c_symp_incidence),
    c_asymp_incidence = sum(c_asymp_incidence)
  ),
  keyby = .(
    location_code,
    week,
    day
  )
  ]

  retval <- retval[
    time == 4,
    .(
      b_S = b_S,
      b_E1 = b_E1,
      b_E2 = b_E2,
      b_I = b_I,
      b_Ia = b_Ia,
      b_R = b_R,
      b_E1_b = b_E1_b,
      b_E2_b = b_E2_b,
      b_I_b = b_I_b,
      b_Ia_b = b_Ia_b,
      c_S = c_S,
      c_E1 = c_E1,
      c_E2 = c_E2,
      c_I = c_I,
      c_Ia = c_Ia,
      c_R = c_R,
      c_E1_b = c_E1_b,
      c_E2_b = c_E2_b,
      c_I_b = c_I_b,
      c_Ia_b = c_Ia_b
    ),
    keyby = .(
      location_code,
      week,
      day
    )
    ]
  retval[retval_incidence, on = c(
    "location_code",
    "week",
    "day"
  ), c_symp_incidence_a := c_symp_incidence_a]
  retval[retval_incidence, on = c(
    "location_code",
    "week",
    "day"
  ), c_asymp_incidence_a := c_asymp_incidence_a]
  retval[retval_incidence, on = c(
    "location_code",
    "week",
    "day"
  ), c_symp_incidence_b := c_symp_incidence_b]
  retval[retval_incidence, on = c(
    "location_code",
    "week",
    "day"
  ), c_asymp_incidence_b := c_asymp_incidence_b]
  retval[retval_incidence, on = c(
    "location_code",
    "week",
    "day"
  ), c_symp_incidence := c_symp_incidence]
  retval[retval_incidence, on = c(
    "location_code",
    "week",
    "day"
  ), c_asymp_incidence := c_asymp_incidence]
  retval <- copy(retval)
  retval[, time := "23:59"]
  setcolorder(retval, c("location_code", "week", "day", "time"))

  return(retval)
}
