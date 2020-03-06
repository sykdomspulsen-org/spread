#'
#'
#' Samples from a branching process model with a constant reproduction number
#'
#' R0 - basic reproduction number
#' serial_interval a distcrete object
#'
#' @param initial_cases Number of initial cases
#' @param R0 R0
#' @param dispersion Dispersion parameter
#' @param serial_interval Serial interval (distribution in days)
#' @param days_simulation Days simulated
#' @param simulations Number of simulations
#' @export
branching_process <- function(
                              initial_cases = 1,
                              R0 = 3,
                              dispersion = 1,
                              serial_interval = distcrete::distcrete(
                                "gamma",
                                interval = 1,
                                shape = epitrix::gamma_mucv2shapescale(mu = 8.4, cv = 3.4 / 8.4)$shape,
                                scale = epitrix::gamma_mucv2shapescale(mu = 8.4, cv = 3.4 / 8.4)$scale,
                                w = 0
                              ),
                              days_simulation = 10,
                              simulations = 1000) {
  incidences <- matrix(nrow = days_simulation, ncol = simulations)
  incidences[1, ] <- initial_cases
  for (day in 2:days_simulation) {
    if (day == 2) {
      FI <- incidences[1:(day - 1), ] * rev(serial_interval$d(1:(day - 1))) * R0
    } else {
      FI <- colSums(incidences[1:(day - 1), ] * rev(serial_interval$d(1:(day - 1)))) * R0
    }
    new <- stats::rnbinom(simulations, mu = FI, size = dispersion)
    incidences[day, ] <- new
  }
  class(incidences) <- append(class(incidences), "bp-incidence")
  return(incidences)
}


#' Fit parameters on cumulative incidence using approximate Bayesian computation (ABC)
#'
#'
#' @param cases_min a
#' @param cases_max a
#' @param param_list a
#' @param simulations Number of simulations
#' @export
fit_params_bp <- function(cases_min, cases_max, param_list, simulations = 100) {
  cumulative <- NULL
  run_bp <- function(param) {
    force(param)
    incidences <- branching_process(
      initial_cases = param$initial_cases,
      R0 = param$R0,
      dispersion = param$dispersion,
      serial_interval = param$serial_interval,
      days_simulation = param$days_simulation,
      simulations = simulations
    )
    cum <- colSums(incidences)
    df <- data.table(
      cumulative = cum,
      R0 = param$R0,
      dispersion = param$dispersion,
      initial_cases = param$initial_cases,
      time = param$end_time,
      serial_interval = param$serial_interval$parameters$shape
    )

    return(df)
  }
  # results_list <- lapply(param_list,run_bp)
  results_list <- future.apply::future_lapply(param_list, run_bp)
  results <- rbindlist(results_list)

  fitting_results <- results[cumulative >= cases_min & cumulative < cases_max]

  class(fitting_results) <- append(class(fitting_results), "bp-fit")
  return(fitting_results)
}





#' plot_quantiles
#'
#' @param da a
#' @param x a
#' @param max_v a
#' @param min a
#' @import ggplot2
#' @export
plot_quantiles_bp <- function(da, x = NULL, max_v = NULL, min = 0) {
  stopifnot(inherits(da, c("bp-incidence", "bp-cumulative")))

  med <- matrixStats::rowQuantiles(da, p = c(0.05, 0.5, 0.95))
  if (is.null(x)) {
    x <- 1:nrow(med)
  }
  if (!is.null(max_v)) {
    med[med[, 3] > max_v, 3] <- max_v
  } else {
    max_v <- max(med[, 3])
  }
  ggplot() +
    geom_line(aes(x = x, y = med[, 2], group = 1)) +
    geom_ribbon(aes(x = x, ymin = med[, 1], ymax = med[, 3], group = 1), alpha = 0.3) +
    ylim(min, max_v) +
    theme_minimal()
}


#' summarize_bp
#'
#' @param incidences A `bp` class
#' @export
summarize_bp <- function(incidences) {
  stopifnot(inherits(incidences, "bp-incidence"))

  cumulative <- apply(incidences, 2, cumsum)
  class(cumulative) <- append(class(incidences), "bp-cumulative")

  retval <- list(
    incidence = incidences,
    cumulative = cumulative,
    p_no_spread = sum(colSums(incidences) == 1) / ncol(incidences)
  )

  class(retval) <- append(class(retval), "bp-summary")

  return(retval)
}
