convert_dynamic_seeds_to_seed_matrix <- function(
                                                 dynamic_seeds,
                                                 location_codes,
                                                 days) {
  skeleton <- data.table(expand.grid(
    location_code = location_codes,
    day = days,
    stringsAsFactors = FALSE
  ))

  if (!is.null(dynamic_seeds)) {
    retval <- merge(
      skeleton,
      dynamic_seeds,
      by = c("location_code", "day"),
      all.x = T
    )
    retval[is.na(n), n := 0]
  } else {
    retval <- skeleton
    retval[, n := 0]
  }
  retval <- dcast.data.table(retval, day ~ location_code, value.var = "n")
  retval[, day := NULL]
  retval <- as.matrix(retval)
  retval
}

convert_beta_to_matrix <- function(
                                   betas,
                                   location_codes,
                                   days,
                                   times) {
  retval <- dcast.data.table(data.table(betas), day + time ~ location_code, value.var = "beta")
  retval[, day := NULL]
  retval[, time := NULL]
  retval <- as.matrix(retval)
  retval
}


#' asymmetric mobility
#'
#' For more information, look at \code{vignette("commuter_model","spread")}
#' @param seiiar_pop Data frame containing `location_code`, `S`, `E`, `I`, `Ia`, and `R` for the entire population
#' @param mobility_matrix List (1 entry for each time period) with each entry containing a data frame with `from`, `to`, `n` for the number of people who travel. Each data frame must be complete.
#' @param dynamic_seeds Data.table containing `location_code`, `day`, and `n` for dynamic seeding (or `NULL`)
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
#'   dynamic_seeds = spread::asymmetric_mobility_dummy_dynamic_seeds,
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
                                dynamic_seeds = NULL,
                                betas = spread::asymmetric_mobility_dummy_betas,
                                latent_period = 1.9,
                                infectious_period = 3.0,
                                asymptomatic_prob = 0,
                                asymptomatic_relative_infectiousness = 0,
                                N = 1) {
  stopifnot(length(mobility_matrix) == length(betas))

  a <- 1 / latent_period
  gamma <- 1 / infectious_period
  days_simulation <- length(betas) / 4

  if (!inherits(seiiar_pop, "data.table")) {
    seiiar_pop <- data.table(seiiar_pop)
  }
  stopifnot(identical(
    names(seiiar_pop),
    c("location_code", "S", "E", "I", "Ia", "R")
  ))
  data.table::setorder(seiiar_pop, location_code)

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
  location_codes <- seiiar_pop$location_code
  days <- seq_along(mobility_matrix)
  seed_matrix <- convert_dynamic_seeds_to_seed_matrix(
    dynamic_seeds = dynamic_seeds,
    location_codes = location_codes,
    days = days
  )

  retval <- asymmetric_mobility_cpp(
    seiiar_pop = seiiar_pop,
    mobility_matrix = mobility_matrix,
    seed_matrix = seed_matrix,
    betas = betas,
    a = a,
    gamma = gamma,
    asymptomaticProb = asymptomatic_prob,
    asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness,
    N = N,
    M = days_simulation
  )
  retval <- retval[
    time == 4,
    .(
      b_S = b_S,
      b_E = b_E,
      b_I = b_I,
      b_Ia = b_Ia,
      b_R = b_R,
      c_S = c_S,
      c_E = c_E,
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
