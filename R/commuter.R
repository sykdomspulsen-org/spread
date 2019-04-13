commuter_calculate_beta_from_r0 <- function(
  r0,
  gamma,
  symptomatic_prob,
  asymptomatic_prob,
  asymptomatic_relative_infectiousness){

  beta <- r0 / gamma / (symptomatic_prob + asymptomatic_prob * asymptomatic_relative_infectiousness)
  beta <- round(beta, 3)

  return(beta)
}

#' commuter
#' @param pop_wo_com Data frame
#' @param di_edge_list Data frame
#' @param start_points Data frame
#' @param r0 Float, basic reproduction number
#' @param beta Float, infection parameter, 0.6
#' @param a Float, 1/latent period, 1/1.9
#' @param gamma Float, 1/infectious period, 1/3
#' @param asymptomatic_prob Float, Proportion/probability of asymptomatic given infectious
#' @param asymptomatic_relative_infectiousness Float, Relative infectiousness of asymptomatic infectious
#' @param N Int = 1 int, Number of repetitions
#' @param M Int, Number of days
#' @import data.table
#' @export
commuter <- function(
  pop_wo_com=spread::norway_pop_wo_com_2017,
  di_edge_list=spread::norway_di_edge_list_2017,
  start_points=spread::start_points_oslo,
  r0=NULL,
  beta=NULL,
  a=1/1.9,
  gamma=1/3,
  asymptomatic_prob=0,
  asymptomatic_relative_infectiousness=0,
  N=1,
  M=7*8
  ) {

  if(!is.null(r0) & !is.null(beta)){
    stop("You cannot specify both 'r0' and 'beta' simultaneously")
  } else if(is.null(r0) & is.null(beta)){
    stop("You need to specify either 'r0' or 'beta'")
  } else if(!is.null(r0) & is.null(beta)){
    beta <- commuter_calculate_beta_from_r0(
      r0=r0,
      gamma=gamma,
      symptomatic_prob = 1-asymptomatic_prob,
      asymptomatic_prob = asymptomatic_prob,
      asymptomatic_relative_infectiousness = asymptomatic_relative_infectiousness
    )
  }

  d <- commuter_cpp(
    pop_wo_com=pop_wo_com,
    di_edge_list=di_edge_list,
    start_points=start_points$I,
    beta=beta,
    a=a,
    gamma=gamma,
    asymptomaticProb = asymptomatic_prob,
    asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness,
    N=N,
    M=M)

  return(d)

}



