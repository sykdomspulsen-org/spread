#'
#'
#' Samples from a branching process model with a constant reproduction number
#'
#' R0 - basic reproduction number
#' serial_interval a distcrete object
#'
#' @export
branching_process <- function(initial_cases=1,
                              R0=3,
                              dispersion=1,
                              serial_interval=0,
                              end_time=10,
                              N=1000){
  incidences <- matrix(nrow=end_time, ncol=N)
  incidences[1,] <- initial_cases
  for(day in 2:end_time){
    if(day == 2){
      FI <- incidences[1:(day-1),] * rev(serial_interval$d(1:(day-1)))*R0
    }else{
      FI <- colSums(incidences[1:(day-1),] * rev(serial_interval$d(1:(day-1))))*R0
    }
    new <- rnbinom(N, mu=FI, size=dispersion)
    incidences[day, ] <- new
  }
  return(incidences)

}


#' Fit parameters
#'
#' @export
fit_params_bp <- function(cases_min, cases_max, param_list, N=100){

  run_bp <- function(param){
    incidences <- branching_process(param$initial_cases,
                                    param$R0,
                                    param$dispersion,
                                    param$serial_interval,
                                    param$end_time, N=N)
    cum <- colSums(incidences)
    df <- data.table(cummulative=cum,
                     R0=param$R0,
                     dispersion=param$dispersion,
                     initial_cases = param$initial_cases,
                     time= param$end_time,
                     serial_interval=param$serial_interval$parameters$shape
                     )


    return(df)
  }
#  results_list <- lapply(param_list,run_bp)
  results_list <- parallel::mclapply(param_list,run_bp, mc.cores=4)
  results <- rbindlist(results_list)


  fitting_results <- results[cummulative >= cases_min & cummulative < cases_max]

  return(fitting_results)

}





#' plot_quantiles
#'
#' @import ggplot2
#'
#' @export
plot_quantiles <- function(da, x=NULL, max_v=NULL, min=0){
  med <- matrixStats::rowQuantiles(da, p=c(0.05, 0.5, 0.95))
  if(is.null(x)){
    x <- 1:nrow(med)
    }
    if(!is.null(max_v)){
        med[med[, 3] > max_v, 3] <- max_v

    }else{
        max_v <- max(med[, 3])
    }
    ggplot() + geom_line(aes(x=x, y=med[,2], group=1)) + geom_ribbon(aes(x=x, ymin=med[, 1], ymax=med[,3], group=1), alpha=0.3) + ylim(min, max_v) + theme_minimal()
}


#'summarize_bp
#'
#' @export
summarize_bp <- function(incidences){
  return(list(
    incidence = incidences,
    cummulative = apply(incidences, 2, cumsum),
    p_no_spread = sum(colSums(incidences)==1)/ ncol(incidences)
    ))
}


