commuter_calculate_beta_from_r0 <- function(
                                            r0,
                                            gamma,
                                            symptomatic_prob,
                                            asymptomatic_prob,
                                            asymptomatic_relative_infectiousness) {
  beta <- r0 * gamma / (symptomatic_prob + asymptomatic_prob * asymptomatic_relative_infectiousness)
  beta <- round(beta, 3)

  return(beta)
}

commuter_convert_seiiar <- function(
                                    seiiar = spread::norway_seiiar_oslo_2017,
                                    commuters = spread::norway_commuters_2017) {
  . <- NULL
  S <- NULL
  E <- NULL
  I <- NULL
  Ia <- NULL
  R <- NULL
  from <- NULL
  pop <- NULL
  n <- NULL
  num_commuters <- NULL
  com_S <- NULL
  com_E <- NULL
  com_I <- NULL
  com_Ia <- NULL
  com_R <- NULL
  cum_R <- NULL
  R_over_limit <- NULL
  R_over_limit_1 <- NULL
  to <- NULL

  seiiar_with_com <- copy(seiiar)
  setDT(seiiar_with_com)

  commutersx <- copy(commuters)
  setDT(commutersx)

  seiiar_with_com[, pop := S + E + I + Ia + R]
  com <- commutersx[, .(num_commuters = sum(n)), keyby = .(from)]
  seiiar_with_com <- merge(seiiar_with_com, com, by.x = "location_code", by.y = "from", all = T)
  seiiar_with_com[is.na(num_commuters), num_commuters := 0]
  seiiar_with_com[, com_S := round(S * num_commuters / pop)]
  seiiar_with_com[, com_E := 0]
  seiiar_with_com[, com_I := 0]
  seiiar_with_com[, com_Ia := 0]
  seiiar_with_com[, com_R := round(R * num_commuters / pop)]

  seiiar_home <- copy(seiiar_with_com)
  seiiar_home[, S := S - com_S]
  seiiar_home[, E := E]
  seiiar_home[, I := I]
  seiiar_home[, Ia := Ia]
  seiiar_home[, R := R - com_R]
  seiiar_home[, pop := NULL]
  seiiar_home[, num_commuters := NULL]
  seiiar_home[, com_S := NULL]
  seiiar_home[, com_E := NULL]
  seiiar_home[, com_I := NULL]
  seiiar_home[, com_Ia := NULL]
  seiiar_home[, com_R := NULL]

  seiiar_commuters <- merge(
    commutersx,
    seiiar_with_com[, c("location_code", "num_commuters", "com_S", "com_E", "com_I", "com_Ia", "com_R")],
    by.x = "from",
    by.y = "location_code",
    all.x = T
  )
  seiiar_commuters[, S := 0]
  seiiar_commuters[, E := 0]
  seiiar_commuters[, I := 0]
  seiiar_commuters[, Ia := 0]
  seiiar_commuters[, R := com_R * n / num_commuters]
  setorder(seiiar_commuters, from, -R)
  seiiar_commuters[, R := ceiling(R)]
  seiiar_commuters[, cum_R := cumsum(R), by = .(from)]
  seiiar_commuters[, R_over_limit := cum_R > com_R]
  seiiar_commuters[, R_over_limit_1 := shift(R_over_limit), by = from]
  seiiar_commuters[is.na(R_over_limit_1), R_over_limit_1 := FALSE]
  # the row before it is wrong, and so is this one
  seiiar_commuters[R_over_limit_1 == TRUE & R_over_limit == TRUE, R := 0]
  # the row before it was right, and this one goes over
  seiiar_commuters[R_over_limit_1 == FALSE & R_over_limit == TRUE, R := R - (cum_R - com_R)]

  seiiar_commuters[, S := n - R]
  seiiar_commuters[, n := NULL]
  seiiar_commuters[, num_commuters := NULL]
  seiiar_commuters[, com_S := NULL]
  seiiar_commuters[, com_E := NULL]
  seiiar_commuters[, com_I := NULL]
  seiiar_commuters[, com_Ia := NULL]
  seiiar_commuters[, com_R := NULL]
  seiiar_commuters[, cum_R := NULL]
  seiiar_commuters[, R_over_limit := NULL]
  seiiar_commuters[, R_over_limit_1 := NULL]

  setorder(seiiar_commuters, from, to)

  return(list(
    "seiiar_home" = seiiar_home,
    "seiiar_commuters" = seiiar_commuters
  ))
}

check_seiiar <- function(seiiar) {
  if (!all.equal(names(seiiar), c("location_code", "S", "E", "I", "Ia", "R"))) {
    stop('names(seiiar) is not c("location_code","S","E","I","Ia","R")')
  }
}

check_commuters <- function(commuters) {
  if (!all.equal(names(commuters), c("from", "to", "n"))) {
    stop('names(commuters) is not c("from","to","n")')
  }
}

#' commuter
#'
#' This model is a stochastic SEIIaR (susceptible, exposed, infectious, infectious asymptomatic, recovered)
#' metapopulation model. Each location has a local infection system, while the locations are connected
#' by people who commute each day. The model differentiates between day and night. During the day you
#' can infect/be infected in the location where you work, while during the night you can infect/be infected
#' in the location where you live. It is the same commuters who travel back and forth each day. At the
#' start of a day, all commuters are sent to their work location, where they mix for 12 hours. The
#' commuters are then sent to their respective home locations, where they mix for 12 hours. The model
#' is based upon a published model.
#' @param seiiar Data frame
#' @param commuters Data frame
#' @param r0 Float, basic reproduction number
#' @param beta Float, infection parameter, 0.6
#' @param latent_period Float, 1.9
#' @param infectious_period Float, 3
#' @param asymptomatic_prob Float, Proportion/probability of asymptomatic given infectious
#' @param asymptomatic_relative_infectiousness Float, Relative infectiousness of asymptomatic infectious
#' @param N Int = 1 int, Number of repetitions
#' @param M Int, Number of days
#' @import data.table
#' @export
commuter <- function(
                     seiiar = spread::norway_seiiar_oslo_2017,
                     commuters = spread::norway_commuters_2017,
                     r0 = NULL,
                     beta = NULL,
                     latent_period = 1.9,
                     infectious_period = 3.0,
                     asymptomatic_prob = 0,
                     asymptomatic_relative_infectiousness = 0,
                     N = 1,
                     M = 7 * 8) {
  . <- NULL
  INCIDENCE <- NULL
  location_code <- NULL
  day <- NULL
  is_6pm <- NULL

  check_seiiar(seiiar)
  check_commuters(commuters)

  a <- 1 / latent_period
  gamma <- 1 / infectious_period

  if (!is.null(r0) & !is.null(beta)) {
    stop("You cannot specify both 'r0' and 'beta' simultaneously")
  } else if (is.null(r0) & is.null(beta)) {
    stop("You need to specify either 'r0' or 'beta'")
  } else if (!is.null(r0) & is.null(beta)) {
    beta <- commuter_calculate_beta_from_r0(
      r0 = r0,
      gamma = gamma,
      symptomatic_prob = 1 - asymptomatic_prob,
      asymptomatic_prob = asymptomatic_prob,
      asymptomatic_relative_infectiousness = asymptomatic_relative_infectiousness
    )
  }

  x <- commuter_convert_seiiar(
    seiiar = seiiar,
    commuters = commuters
  )

  d <- commuter_cpp(
    seiiar_home = x[["seiiar_home"]],
    seiiar_commuters = x[["seiiar_commuters"]],
    beta = beta,
    a = a,
    gamma = gamma,
    asymptomaticProb = asymptomatic_prob,
    asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness,
    N = N,
    M = M
  )

  d[, incidence := sum(incidence), by = .(location_code, day)]
  d <- d[is_6pm == 1]
  d[, pop := S + E + I + Ia + R]

  return(d)
}
